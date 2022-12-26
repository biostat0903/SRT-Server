#! /usr/bin/env Rscript
# quality control for st data
# up-stream procedures: none
# down-stream procedures: dr_cl.R, dr_cl_wr.R, dr_cl_sp.R, svg.R, dc.R

# Set method_path 
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# Load packages
library(Seurat)
library(dplyr)
library(bigreadr)

# Function 1: qc check 
qc.check <- function(data_path = NULL,        ## String: path string of samples *****(required)
                     platform = NULL,         ## String: platform for st data *****(required)
                     ## four platforms: {Visium}, {Merfish}, {Seqfish}, {Slideseq}
                     min_cells = 10,          ## Integer: features detected in at least this many cells
                     ## default setting: 10
                     min_features = 100,      ## Integer: cells where at least this many features are detected
                     ## default setting: 100
                     out_path = NULL          ## String: path for each sample
){
  
  ## split path string
  data_dir_s <- strsplit(data_path, "\\|")[[1]]
  write.table(c(data_dir_s, platform,
                min_features, min_cells), 
              file = paste0(out_path, "/qc_check_file.txt"), 
              row.names = F, quote = F, col.names = F)
  return(0)
}

# Function 2: quality control MERFISH data
LoadMerfish_QC <- function(data_path,         ## path for st data
                           min_features,      ## minimum feature number for each cell
                           min_cells          ## minimum cell number for each feature
){
  
  ## check files
  if(!"cell_by_gene.csv" %in% list.files(data_path) | !"cell_metadata.csv" %in% list.files(data_path)){
    
    stop("Count matrix or meta_data does not exist! Please check!")
  } 
  ## read data
  count_data <- fread2(paste0(data_path, "/cell_by_gene.csv"))
  ### delete "blank" gene
  count_data <- count_data[, !grepl("Blank", colnames(count_data))]
  message("Before QC, the data contains ", dim(count_data)[1], 
          " spots and ", dim(count_data)[2], " genes.\n")
  cell_ID <- count_data[, 1]
  count_data <- count_data[, -1] %>% t
  colnames(count_data) <- cell_ID
  meta_data <- fread2(paste0(data_path, "/cell_metadata.csv"))
  meta_data <- meta_data[, -1]
  rownames(meta_data) <- cell_ID
  
  ## pre-process data
  seurat_obj <- CreateSeuratObject(counts = count_data, 
                                   min.cells = min_cells, 
                                   min.features = min_features, 
                                   meta.data = meta_data)
  count_data_qc <- as.matrix(seurat_obj@assays$RNA@counts)
  message("After QC, the data contains ", dim(count_data_qc)[2], 
          " spots and ", dim(count_data_qc)[1], " genes.\n")
  meta_data_qc <- cbind.data.frame(x = seurat_obj@meta.data$center_x, 
                                   y = seurat_obj@meta.data$center_y)
  
  # if(ncol(count_data_qc) > 20000){
  #   
  #   seed <- as.numeric(as.POSIXct(Sys.time()))
  #   set.seed(seed)
  #   cell_ind <- sort(sample(c(1: ncol(count_data_qc)), 20000))
  #   count_data_qc <- count_data_qc[, cell_ind]
  #   meta_data_qc <- meta_data_qc[cell_ind, ]
  #   warning("The cell number is larger than 20 thousands. We randomly select 20 thousands cells to perform the following analysis.")
  # }
  
  spatial_result <- list(count = as.matrix(count_data_qc), 
                         meta.data = as.matrix(meta_data_qc))
  
  return(spatial_result)
}  

# Function 3: quality control Visium data 
LoadVisium_QC <- function(data_path,          ## path for st data
                          min_features,       ## minimum feature number for each cell
                          min_cells           ## minimum cell number for each feature
){
  
  ## read data
  spatial_obj <- Load10X_Spatial(data_path) 
  count_data <- as.matrix(spatial_obj@assays$Spatial@counts)
  meta_data <- spatial_obj@images$slice1@coordinates
  ### delete miRNA, lncRNA and AC gene
  gene_del_label <- grepl("^MT-", rownames(count_data)) |
    grepl("^MIR-", rownames(count_data)) |
    grepl("^LINC", rownames(count_data)) |
    grepl("^A[A-Z]", rownames(count_data)) |
    grepl("^B[A-Z]", rownames(count_data)) |
    grepl("^C[A-Z]", rownames(count_data)) |
    grepl("Rik", rownames(count_data)) |
    grepl("^Gm", rownames(count_data)) |
    grepl("^Vmn", rownames(count_data)) |
    grepl("^Olfr", rownames(count_data)) |
    grepl("Trav", rownames(count_data))
  count_data <- count_data[!gene_del_label, ]
  message("Before QC, the data contains ", dim(count_data)[2], 
          " spots and ", dim(count_data)[1], " genes.\n")
  
  ### delete spots with no coord
  meta_data <- subset(meta_data,
                      rowSums(is.na(meta_data)) == 0)
  count_data <- count_data[, which(colnames(count_data) %in% rownames(meta_data))]
  meta_data <- meta_data[colnames(count_data), ]
  
  ## pre-process data
  seurat_obj <- CreateSeuratObject(counts = count_data, 
                                   min.cells = min_cells, 
                                   min.features = min_features, 
                                   meta.data = meta_data)
  count_data_qc <- as.matrix(seurat_obj@assays$RNA@counts)
  message("After QC, the data contains ", dim(count_data_qc)[2], 
          " spots and ", dim(count_data_qc)[1], " genes.\n")
  meta_data_qc <- seurat_obj@meta.data[, c(8, 7)]
  meta_data_qc[, 2] <- -meta_data_qc[, 2]
  spatial_result <- list(count = as.matrix(count_data_qc), 
                         meta.data = as.matrix(meta_data_qc))
  
  return(spatial_result)
}

# Function 4: quality control SEQFISH data
LoadSeqfish_QC <- function(data_path = NULL,  ## path for st data
                           min_features,      ## minimum feature number for each cell
                           min_cells          ## minimum cell number for each feature
){
  
  ## check files list
  file_list_ref <- c("cell_locations/centroids_annot.txt", 
                     "cell_locations/centroids_coord.txt",
                     # "cell_locations/centroids_offset.txt",
                     "count_matrix/expression.txt", 
                     "raw_data/location_fields.png")
  file_list <- list.files(data_path, recursive = T)
  if(!all(file_list_ref %in% file_list)){
    
    miss_file_str = setdiff(file_list_ref, file_list) %>% paste(., collapse = ", ")
    stop(paste0(miss_file_str, " not found! Please check!"))
  } 
  
  ## read count data
  count_data <- fread2(paste0(data_path, "/count_matrix/expression.txt"))
  if (colnames(count_data)[1] == "V1") {
    count_data <- tibble::column_to_rownames(count_data, "V1")
  }
  ### delete "blank" gene
  count_data <- count_data[, !grepl("Blank", colnames(count_data))]
  message("Before QC, the data contains ", dim(count_data)[2], 
          " spots and ", dim(count_data)[1], " genes.\n")
  
  ## pre-process count data
  seurat_obj <- CreateSeuratObject(counts = count_data, 
                                   min.cells = min_cells, 
                                   min.features = min_features, 
                                   meta.data = NULL)
  count_data_qc <- as.matrix(seurat_obj@assays$RNA@counts)
  message("After QC, the data contains ", dim(count_data_qc)[2], 
          " spots and ", dim(count_data_qc)[1], " genes.\n")
  
  ## read coordinate data
  coord_data <- fread2(paste0(data_path, "/cell_locations/centroids_coord.txt"))
  if (ncol(coord_data) < 3) {
    stop("The first three columns of the centroids_coord file are required to be cell ID, x and y!")
  } else {
    coord_data <- coord_data[,1:3]
    colnames(coord_data) <- c("ID", "X", "Y")
  }
  
  ## read FOV data
  annot_data <- fread2(paste0(data_path, "/cell_locations/centroids_annot.txt"))
  if (ncol(coord_data) < 2) {
    stop("The first two columns of the centroids_annot file are required to be cell ID and FOV!")
  } else {
    annot_data <- annot_data[,1:2]
    colnames(annot_data) <- c("ID", "FOV")
  }
  
  ## pre-process meta data
  ### merge coord and FOV data
  if (!all(annot_data$ID == coord_data$ID)) {
    
    stop('Cell ID in centroids_coord and centroids_annot files should be in the same order!')
  } else {
    
    coord_annot <- coord_data
    coord_annot$FOV <- annot_data$FOV
  }
  
  ### check number and add cell_ID column
  if(nrow(coord_annot) != ncol(count_data)) {
    
    stop('Cell number in coordinate file do not equal to that of expression file!')
  } else {
    
    rownames(coord_annot) <- colnames(count_data)
  }
  
  ### stitch field coordinates from different FOV
  offset_file <- paste0(data_path, "/cell_locations/centroids_offset.txt")
  if (file.exists(offset_file)) {
    #
    offset_data <- fread2(offset_file)
    if (ncol(offset_data) < 3) {
      
      stop("The first Third columns of the centroids_offset file are required to be FOV, x_offset and y_offset!")
    } else {
      
      colnames(offset_data) <- c("FOV", "x_offset", "y_offset")
    }
    if (!all(unique(coord_annot$FOV %in% offset_data$FOV))) {
      
      stop("Not all FOV from centroids_annot file in centroids_offset file!")
    } else {
      
      offset_match <- offset_data[match(coord_annot$FOV, offset_data$FOV),]
      coord_annot$X <- coord_annot$X + offset_match$x_offset
      coord_annot$Y <- coord_annot$Y + offset_match$y_offset
    }
  }
  
  meta_data_qc <- coord_annot[colnames(count_data_qc), c("X", "Y")]  
  
  ## output
  spatial_result <- list(count = as.matrix(count_data_qc), 
                         meta.data = as.matrix(meta_data_qc))
  
  return(spatial_result)
}  

# Function 5: quality control Slide-seq data
LoadSlideseq_QC <- function(data_path,        ## path for st data
                            min_features,     ## minimum feature number for each cell
                            min_cells         ## minimum cell number for each feature
){
  
  # data_path <- "/net/mulan/disk2/yasheng/stwebProject/02_data/04_slide-seq/Mouse_Hippocampus2/"
  ## check files
  file_names <- list.files(data_path)
  if(sum(grepl("location", file_names)) != 1 | sum(grepl("expression", file_names)) != 1){
    
    stop("Gene Expression matrix or location matrix does not exist! Please check!")
  } 
  ## read data
  count_data <- fread2(paste0(data_path, "/", file_names[grepl("expression", file_names)]))
  geneID <- count_data[, 1]
  count_data <- count_data[, -1]
  rownames(count_data) <- geneID
  ### delete gene
  gene_del_label <- grepl("Rik", rownames(count_data)) |
    grepl("A[A-Z]", rownames(count_data)) |
    grepl("^RP2", rownames(count_data)) |
    grepl("^WI", rownames(count_data)) 
  count_data <- count_data[!gene_del_label, ]
  message("Before QC, the data contains ", dim(count_data)[2], 
          " spots and ", dim(count_data)[1], " genes.\n")
  meta_data <- fread2(paste0(data_path, "/", file_names[grepl("locations", file_names)]), 
                      skip = 1)
  meta_data[, 2] <- ceiling(meta_data[, 2])
  meta_data[, 3] <- ceiling(meta_data[, 3])
  colnames(meta_data) <- c("barcodes", "xcoord", "ycoord")
  rownames(meta_data) <- meta_data$barcodes
  
  ## pre-process data
  seurat_obj <- CreateSeuratObject(counts = count_data, 
                                   min.cells = min_cells, 
                                   min.features = min_features, 
                                   meta.data = meta_data)
  count_data_qc <- as.matrix(seurat_obj@assays$RNA@counts)
  message("After QC, the data contains ", dim(count_data_qc)[2], 
          " spots and ", dim(count_data_qc)[1], " genes.\n")
  meta_data_qc <- cbind.data.frame(x = seurat_obj@meta.data$center_x, 
                                   y = seurat_obj@meta.data$center_y)
  spatial_result <- list(count = as.matrix(count_data_qc), 
                         meta.data = as.matrix(meta_data_qc))
  
  return(spatial_result)
}  

# Function 6: quality control general format
LoadGeneal_QC <- function(data_path,         ## path for st data
                          min_features,      ## minimum feature number for each cell
                          min_cells          ## minimum cell number for each feature
){
  
  ## check files
  if(!"gene_expression.csv" %in% list.files(data_path) | !"cell_location.csv" %in% list.files(data_path)){
    
    stop("Gene expression or location does not exist! Please check!")
  } 
  ## read data
  count_data <- fread2(paste0(data_path, "/gene_expression.csv"))
  ### delete "blank" gene
  message("Before QC, the data contains ", dim(count_data)[1], 
          " spots and ", dim(count_data)[2], " genes.\n")
  cell_ID <- count_data[, 1]
  count_data <- count_data[, -1] %>% t
  colnames(count_data) <- cell_ID
  meta_data <- fread2(paste0(data_path, "/cell_location.csv"))
  meta_data <- meta_data[, c(2, 3)]
  rownames(meta_data) <- cell_ID
  
  ## pre-process data
  seurat_obj <- CreateSeuratObject(counts = count_data, 
                                   min.cells = min_cells, 
                                   min.features = min_features, 
                                   meta.data = meta_data)
  count_data_qc <- as.matrix(seurat_obj@assays$RNA@counts)
  message("After QC, the data contains ", dim(count_data_qc)[2], 
          " spots and ", dim(count_data_qc)[1], " genes.\n")
  meta_data_qc <- cbind.data.frame(x = seurat_obj@meta.data$x, 
                                   y = seurat_obj@meta.data$y)
  spatial_result <- list(count = as.matrix(count_data_qc), 
                         meta.data = as.matrix(meta_data_qc))
  
  return(spatial_result)
}  



# Function 7: qc call
qc.call <- function(data_path,                ## String: path for st data
                    out_path                  ## String: path for check file
                    
){
  
  ## load io code
  source(paste0(method_path, "/io.R"))
  
  ## load spatial data
  check_file <- read.table(paste0(data_path, "/qc_check_file.txt"))[, 1]
  sample_size <- length(check_file)-3
  platform <- check_file[sample_size + 1]
  min_features <- check_file[sample_size + 2] %>% as.numeric
  min_cells <- check_file[sample_size + 3] %>% as.numeric
  
  ### Merfish
  if (platform == "Merfish"){
    
    spatial_data_all <- plyr::alply(c(1: sample_size), 1, function(a){
      spatial_data <- try(LoadMerfish_QC(check_file[a], min_features, min_cells), silent = T)
      if (inherits(spatial_data, "try-error")){
        stop(paste0("File of sample ", a, " is wrong!"))
      }
      return(spatial_data)
    })
  }
  ### Visium  
  if (platform == "Visium"){
    
    spatial_data_all <- plyr::alply(c(1: sample_size), 1, function(a){
      spatial_data <- try(LoadVisium_QC(check_file[a], min_features, min_cells), silent = T)
      if (inherits(spatial_data, "try-error")){
        stop(paste0("File of sample ", a, " is wrong!"))
      }
      return(spatial_data)
    })
  }
  ### Seqfish  
  if (platform == "Seqfish"){
    
    spatial_data_all <- plyr::alply(c(1: sample_size), 1, function(a){
      spatial_data <- try(LoadSeqfish_QC(check_file[a], min_features, min_cells), silent = T)
      if (inherits(spatial_data, "try-error")){
        stop(paste0("File of sample ", a, " is wrong!"))
      }
      return(spatial_data)
    })
  }
  ### Slide-seq 
  if (platform == "Slideseq"){
    
    spatial_data_all <- plyr::alply(c(1: sample_size), 1, function(a){
      spatial_data <- try(LoadSlideseq_QC(check_file[a], min_features, min_cells), silent = T)
      if (inherits(spatial_data, "try-error")){
        
        stop(paste0("File of sample ", a, " is wrong!"))
      }
      return(spatial_data)
    })
  }
  
  if (sample_size != 1){
    
    ## get intersection of genes
    gene_all <- plyr::alply(c(1: sample_size), 1, function(a){
      rownames(spatial_data_all[[a]]$count)
    })
    gene_inter <- Reduce(intersect, gene_all)
    spatial_data_all_inter <- plyr::llply(c(1: sample_size), function(a){
      ### format cell ID
      count_a <- spatial_data_all[[a]]$count
      colnames(count_a) <- paste0("sample", a, "_", colnames(count_a))
      count_a <- count_a[match(gene_inter, rownames(count_a)), , drop = F]
      ### intersect gene for each sample
      spatial_data_a <- list(count = count_a, 
                             meta.data = spatial_data_all[[a]]$meta.data)
      return(spatial_data_a)
    })
  } else {
    
    spatial_data_all_inter <- spatial_data_all
  }
  
  ## save to h5 file
  result_dir <- paste0(data_path, "/qc_result")
  if (!file.exists(result_dir)){
    
    system(paste0("mkdir ", result_dir))
  }
  spatial_data_filename <- paste0(result_dir, "/qc_data.h5")
  if (file.exists(spatial_data_filename)){
    
    system(paste0("rm -rf ", spatial_data_filename))
    warning("The h5 exists. We delete it!")
  }
  bool_write <- h5data.write(spatial_data_filename, sample_size,  
                             spatial_data_all_inter, platform)
  
  ## save QC check file
  write.table(c(spatial_data_filename, sample_size, 
                min_features, min_cells), 
              file = paste0(out_path, "/qc_call_file.txt"), 
              row.names = F, col.names = F, quote = F)
  
  return(0)
}

# ###################
# ### test code
# input_path <- "/net/mulan/disk2/yasheng/stwebProject/02_data/02_Visium/Parent_Visium_Human_BreastCancer"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/test/qc"
# qc.check(data_path = input_path,
#          platform = "Visium",
#          min_cells = 100,
#          min_features = 100,
#          out_path = output_path)
# qc.call(data_path = output_path,
#         out_path = output_path)
