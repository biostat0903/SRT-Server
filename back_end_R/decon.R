#! /usr/bin/env Rscript
# deconvolution spot referring to scRNA-seq or marker genes
# up-stream procedure: qc.R
# down-stream procedure: dc_plt.R

# load packages
library(CARD)
library(hdf5r)
library(plyr)
library(gtools)
library(bigreadr)
library(dplyr)
library(SingleCellExperiment)

# method_path setting
method_path <- "/root/st_server/"
refer_path <- "/root/st_server/"

# Fix parematers
MINCountGene = 0
MINCountSpot = 0
MAXCellNum = 10000
N_CORE = 10

# Function 1: check function
decon.check <- function(data_path = NULL,               ## String: output path for qc procedure
                     species = "Human",              ## String: species: {Human} and {Mouse}                 
                     status = NULL,                  ## String: combine tissue and disease status
                     marker_list = NULL,             ## String: path of a list format RData with cell type specific gene markers in CARDfree version
                     num_grid = 2000,                ## Integer: initial number of grids to construact a refined spatial map
                     scMapping = FALSE,              ## Boolean: whether to perform single cell resolution mapping
                     out_path = NULL                 ## String: output path for dc procedure
                     
){
  ## check inputs
  if (is.null(species) & is.null(marker_list)){
    
    stop("No scRNA-seq reference or marker_list!")
  }
  check_file <- paste0(data_path, "/qc_call_file.txt")
  if(!file.exists(check_file)){
    
    stop("No QC file! Please run QC module!")
  } else {

    qc_param <- read.table(check_file)[, 1]
    spatial_data_filename <- qc_param[1]
    sample_size <- qc_param[2]
  }
  
  ## create CARD object
  if (!is.null(species) & !is.null(status)){
    
    ref_type <- 'sc'
    sc_dir <- paste0(refer_path, "/scRNA_ref/", species, "/", status)
  } else{

    ref_type <- 'free'
    sc_dir <- NULL
  }
  # output
  
  write.table(c(spatial_data_filename, sample_size, sc_dir, 
                ref_type, num_grid, scMapping),
              file = paste0(out_path, "/dc_check_file.txt"), 
              row.names = F, col.names = F, quote = F)
  return(0)
}


# Function 2: call function
decon.call <- function(out_path                         ## String: output path for dc procedure
){
  
  ## load io code
  source(paste0(method_path, "/io.R"))
  
  ## load inputs
  check_file <- fread2(paste0(out_path, "/dc_check_file.txt"), header = F)[, 1]
  spatial_data_filename <- check_file[1]
  sample_size <- check_file[2] %>% as.numeric()
  sc_dir <- check_file[3]
  ref_type <- check_file[4] 
  num_grid <- check_file[5] %>% as.numeric()
  scMapping <- ifelse(check_file[6] == "FALSE", F, T)

  ## build CARD object
  st_list <- h5data.load(spatial_data_filename,         
                         sample_size = sample_size,       
                         load_count = TRUE,     
                         normalization = FALSE,  
                         load_coord = TRUE,    
                         coordinate = TRUE )
  platform <- st_list[["platform"]]
  if (ref_type == "sc"){
    load(paste0(sc_dir, "/sc_count.RData"))
    load(paste0(sc_dir, "/sc_meta.RData"))
    if (nrow(sc_meta) > MAXCellNum) {
      
      warning("The cell number of reference sc data is large. We randomly select 10,000 cells!")
      ref_cell_sub <- sample(sc_meta$cellID, MAXCellNum, replace = F)
      sc_meta <- sc_meta[match(ref_cell_sub, sc_meta$cellID),]
      sc_count <- sc_count[, ref_cell_sub]
    }
    rownames(sc_count) <- toupper(rownames(sc_count))
    sc_meta$cellType <- gsub(" ", "_", sc_meta$cellType)
    CARD_obj_list <- lapply(1:sample_size, function(a){
      createCARDObject(sc_count = sc_count,
                       sc_meta = sc_meta,
                       spatial_count = st_list[["count_list"]][[a]],
                       spatial_location = as.data.frame(st_list[["coord_list"]][[a]]),
                       ct.varname = "cellType",
                       ct.select = unique(sc_meta$cellType),
                       sample.varname = "sampleInfo",
                       minCountGene = MINCountGene,
                       minCountSpot = MINCountSpot)
    })
  } else {
    
    load(markerList)
    CARD_obj_list <- plyr::llply(1:sample_size, function(a){
      createCARDfreeObject(markerList = markerList,
                           spatial_count = count_list[[a]],
                           spatial_location = coord_list[[a]],
                           minCountGene = MINCountGene,
                           minCountSpot = MINCountSpot)
    })
  }
  
  ## deconvolution
  CARD_obj_list <- plyr::llply(seq_along(CARD_obj_list), function(a){
    
    CARD_obj <- CARD_obj_list[[a]]
    if (ref_type == 'sc'){
      CARD_obj <- CARD_deconvolution(CARD_obj)
    } else{
      
      CARD_obj <- CARD_refFree(CARD_obj)
    }
    CARD_obj <- CARD.imputation(CARD_obj,
                                NumGrids = num_grid,
                                ineibor = 10,
                                exclude = NULL)
    
    return(CARD_obj)
  })
  
  ## format results
  decon_result_list <- plyr::llply(seq_along(CARD_obj_list), function(a){
    
    decon_result <- list()
    decon_result[["Proportion_CARD"]] <- CARD_obj_list[[a]]@Proportion_CARD
    decon_result[["location"]] <- CARD_obj_list[[a]]@spatial_location
    decon_result[["sample_info"]] <- 
      data.frame(sample = a,
                 cell = rownames(CARD_obj_list[[a]]@spatial_location))
    return(decon_result)
  })
  decon_refine_result_list <- plyr::llply(seq_along(CARD_obj_list), function(a){
    
    decon_refine_result <- list()
    decon_refine_result[["refined_expression"]] <- CARD_obj_list[[a]]@refined_expression
    decon_refine_result[["refined_prop"]] <- CARD_obj_list[[a]]@refined_prop
    
    return(decon_refine_result)
  })
  print("CARD analysis completed!")
  
  ## single cell resolution mapping
  if (scMapping == TRUE) {
    
    print("Perform single cell resolution mapping!")
    num_cell <- ifelse(platform == "Visium", 7,
                       ifelse(platform == "Slideseq", 2, 20))
    
    scMapping_list <- plyr::llply(CARD_obj_list, function(CARD_obj){
      CARD_SCMapping(CARD_object = CARD_obj,
                     shapeSpot = "Square",
                     numCell = num_cell,
                     ncore = N_CORE)
    })
    
    ## save scMapping result
    sce_count_list <- list()
    sce_coord_meta_list <- list()
    for (a in seq_along(scMapping_list)) {
      scMapping_a <- scMapping_list[[a]]
      sce_count_list[[a]] <- assays(scMapping_a)$counts %>% 
        as.matrix() %>% as.data.frame()
      sce_coord_meta_list[[a]] <- colData(scMapping_a) %>% 
        as.data.frame()
    }
    scMapping_path <- c(paste0(out_path, '/sce_count_list.RData'),
                        paste0(out_path, '/sce_coord_meta_list.RData'))
    
    save(sce_count_list, file = scMapping_path[1])
    save(sce_coord_meta_list, file = scMapping_path[2])
  } else {
    
    scMapping_path <- NULL
  }
    
  ## output
  resule_path <- paste0(out_path, '/de_result_list.RData')
  save(decon_result_list, file = resule_path)
  refine_resule_path <- paste0(out_path, '/decon_refine_result_list.RData')
  save(decon_refine_result_list, file = refine_resule_path)
  write.table(c(resule_path, refine_resule_path, ref_type, scMapping_path),
              file = paste0(out_path, "/decon_call_file.txt"), 
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  return(0)
}

# Function 3: post function
decon.post <- function(out_path                         ## String: output path for decon procedure
){
  
  ## load data
  call_file <- fread2(paste0(out_path, "/decon_call_file.txt"), header = F)[, 1]
  result_path <- call_file[1]
  ref_type <- call_file[3]
  scMapping_path <- NULL
  if (length(call_file) > 3){
    scMapping_path <- c(call_file[4], call_file[5])
  } 
  
  ## format deconv result
  load(result_path) # CARD_result_list
  deconv <- lapply(decon_result_list, function(decon_result){
    sample_info <- decon_result$sample_info
    deconv <- cbind(sample_info, 
                    decon_result$Proportion_CARD[sample_info$cell,])
  }) %>% Reduce("rbind", .)

  result_dir <- paste0(out_path, "/decon_result/")
  if (!file.exists(result_dir)) {
    system(paste0("mkdir -p ", result_dir))
  }
  write.table(deconv, file = paste0(result_dir, "/deconv_", ref_type, ".txt"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
  message("Deconvolution result saved!")
  
  # format 
  if (!is.null(scMapping_path)) {
    
    load(scMapping_path[1]) # sce_count_list
    load(scMapping_path[2]) # sce_coord_meta_list
    
    scmap_path <- lapply(seq_along(sce_count_list), function(a){
      #
      sce_count <- sce_count_list[[a]]
      sce_coord_meta <- sce_coord_meta_list[[a]]
      sce_coord <- sce_coord_meta[,c("x", "y")]
      sce_meta <- sce_coord_meta[, c("Cell", "CT")]
      colnames(sce_meta) <- c("cell", "celltype")
      
      sce_result_dir <- paste0(out_path, "/decon_result/sample", a, "/")
      if (!file.exists(sce_result_dir)) {
        system(paste0("mkdir -p ", sce_result_dir))
      }
      #
      fwrite2(sce_count, 
              file = paste0(sce_result_dir, "/count.txt"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
      system(paste0("gzip -f ", sce_result_dir, "/count.txt"))
      
      fwrite2(sce_coord, 
              file = paste0(sce_result_dir, "/coord.txt"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
      system(paste0("gzip -f ", sce_result_dir, "/coord.txt"))
      
      fwrite2(sce_meta, 
              file = paste0(sce_result_dir, "/annotation.txt"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
      system(paste0("gzip -f ", sce_result_dir, "/annotation.txt"))
      
      print(paste0("Mapping result for sample", a, " saved!"))
      
      return(paste0(sce_result_dir, "/count.txt.gz,",
                    sce_result_dir, "/coord.txt.gz,",
                    sce_result_dir, "/annotation.txt.gz"))
    }) %>% unlist()

  } else {
    scmap_path <- NULL
  }
  
  write.table(c(paste0(result_dir, "/deconv_",ref_type,".txt"), 
                scmap_path), 
              file = paste0(out_path, "/decon_post_file.txt"), 
              sep = "\t", row.names = F, col.names = F, quote = F)
  return(0)
}

###################
### test code
#data_path <- "/net/mulan/disk2/yasheng/stwebProject/test/qc"
#output_path <- "/net/mulan/disk2/yasheng/stwebProject/test/decon"
#decon.check(data_path = data_path,
#         species = "Human",
#         status = "Breast_BRCA",
         # marker_list = NULL,
#         num_grid = 20000,
#         scMapping = T,
#         out_path = output_path)
#decon.call(out_path = output_path)
#decon.post(out_path = output_path)



