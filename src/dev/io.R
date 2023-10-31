#! /usr/bin/env Rscript
# Input/Output of h5 format st data
library(hdf5r)

# Function 1: Write st data to h5 file
h5data.write <- function(data_filename,    ## String: file name of h5 data 
                         sample_size,      ## Integer: sample size
                         data_list,        ## List: st data (count and coord)
                         platform
){
  
  data_h5 <- H5File$new(data_filename, mode = "w")
  data_grp <- data_h5$create_group("count")
  count_all <- plyr::alply(c(1: sample_size), 1, function(a){
    grp_name <- paste0("count.s", a)
    data_grp[[grp_name]] <- data_list[[a]]$count
  })
  data_grp <- data_h5$create_group("cellID")
  count_all <- plyr::alply(c(1: sample_size), 1, function(a){
    grp_name <- paste0("cellID.s", a)
    data_grp[[grp_name]] <- colnames(data_list[[a]]$count)
  })
  data_grp <- data_h5$create_group("featureID")
  count_all <- plyr::alply(c(1: sample_size), 1, function(a){
    grp_name <- paste0("featureID.s", a)
    data_grp[[grp_name]] <- rownames(data_list[[a]]$count)
  })
  data_grp <- data_h5$create_group("meta.data")
  meta_data_all <- plyr::alply(c(1: sample_size), 1, function(a){
    grp_name <- paste0("meta.data.s", a)
    data_grp[[grp_name]] <- data_list[[a]]$meta.data
  })
  data_grp <- data_h5$create_group("image.data")
  image_data_all <- plyr::alply(c(1: sample_size), 1, function(a){
    grp_name <- paste0("image.data.s", a)
    data_grp[[grp_name]] <- data_list[[a]]$image.data
  })
  data_grp <- data_h5$create_group("image.meta")
  image_meta_all <- plyr::alply(c(1: sample_size), 1, function(a){
    grp_name <- paste0("image.meta.s", a)
    data_grp[[grp_name]] <- data_list[[a]]$image.meta
  })
  data_grp <- data_h5$create_group("platform")
  data_grp[["platform"]] <- platform
  data_h5$close_all()
  
  return(0)
}

# Function 2: Load h5 data
h5data.load <- function(data_filename,          ## String: file name of h5 data
                        sample_size = 1,        ## Integer: sample size
                        load_count = TRUE,      ## Boolean: load count matrix 
                        normalization = TRUE,   ## Boolean: normalize count matrix
                        load_coord = FALSE,     ## Boolean: load coordinate data
                        coordinate = TRUE,      ## Boolean: coordinate gene for each sample
                        image = F               ## Boolean: load image data
){
  
  data_h5 <- H5File$new(data_filename, mode = "r")
  ## load count matrix
  if (load_count == TRUE) {
    
    if (!"count" %in% list.groups(data_h5)){
      
      stop("h5 of QC data should include \"count\" group!")
    } else {
      
      count_list <- plyr::alply(c(1: sample_size), 1, function(a){
        
        count_grp_name <- paste0("count/count.s", a)
        count_mat <- data_h5[[count_grp_name]][, ] 
        cellID_grp_name <- paste0("cellID/cellID.s", a)
        colnames(count_mat) <- data_h5[[cellID_grp_name]][] 
        featureID_grp_name <- paste0("featureID/featureID.s", a)
        rownames(count_mat) <- data_h5[[featureID_grp_name]][] 
        
        if(normalization == TRUE){
          
          # log-normalization
          scale_factor <- 10000
          total_counts <- colSums(count_mat)
          norm_mat <- apply(count_mat, 1, function(x){
            log(x / total_counts*scale_factor + 1)
          }) %>% t()
          dimnames(norm_mat) <- dimnames(count_mat)
          return(norm_mat)
        } else {
          
          return(count_mat)
        }
      })
      cat("Count data loaded!\n")
      
      ### coordinate count list
      if (coordinate == TRUE){
        gene_list <- lapply(count_list, function(a){
          rownames(a)
        })
        gene_inter <- Reduce("intersect", gene_list)
        count_list <- lapply(seq_along(count_list), function(a){
          count_s <- count_list[[a]][rownames(count_list[[a]]) %in% gene_inter, ]
          return(count_s)
        })
      }

    }
   
  } else {
    
    count_list <- NULL
  }
  
  ## load coord list
  if (load_coord == TRUE) {
    
    if (!"meta.data" %in% list.groups(data_h5)){
      
      stop("h5 of QC data should include \"meta.data\" group!")
    } else {
      
      coord_list <-  plyr::alply(c(1: sample_size), 1, function(a){
        coord_grp_name <- paste0("meta.data/meta.data.s", a)
        coord_mat <- data_h5[[coord_grp_name]][, ] 
        cellID_grp_name <- paste0("cellID/cellID.s", a)
        dimnames(coord_mat) <- list(data_h5[[cellID_grp_name]][],
                                    c("x", "y")) 
        return(coord_mat)
      })
      cat("Coord data loaded!\n")
    }
  } else {
    
    coord_list <- NULL
  }
  
  ## load platform
  platform <- data_h5[["platform/platform"]][]
  
  ## load coord list
  if (platform == "Visium" & image == TRUE) {
    
    if (!"image.data" %in% list.groups(data_h5)){
      
      stop("h5 of QC data should include \"image.data\" group!")
    } else {
      
      image_list <-  plyr::alply(c(1: sample_size), 1, function(a){
        image_grp_name <- paste0("image.data/image.data.s", a)
        image_data <- data_h5[[image_grp_name]][, ,] 
        image_meta_grp_name <- paste0("image.meta/image.meta.s", a)
        image_meta <- data_h5[[image_meta_grp_name]][]
        scalef <- as.numeric(image_meta[-1])
        names(scalef) <- c("fiducial_diameter_fullres",
                           "spot_diameter_fullres",
                           "tissue_hires_scalef",
                           "tissue_lowres_scalef")
        return(list("image_data" = image_data,
                    "image_path" = image_meta[1],
                    "scale_factors" = scalef))
      })
      cat("Image data loaded!\n")
    }
  } else {
    
    image_list <- NULL
  }
  
  ## final return count list and coord list
  return(list("count_list" = count_list,
              "coord_list" = coord_list, 
              "platform" = platform,
              "image_list" = image_list))
}

