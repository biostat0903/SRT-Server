#! /usr/bin/env Rscript
# Identify differentially expressed gene in st data
# up-stream code: ct.R, sdd.R
# down-stream code: de_plt.R, ora.R

# Set method_path
method_path <- "/public/home/biostat03/project/stwebProject/01_code/srt_server/dev"

# Load packages
library(Seurat)
library(hdf5r)
library(dplyr)
library(bigreadr)

# Function 1: check function
de.check <- function(data_path1 = NULL,                 ## String: output path of ct or sdd procedure
                     data_path2 = NULL,                 ## String: output path of qc procedure
                     submodule = NULL,                  ## String: cell typing sub-modules: {CT_PCA}, {SDD_sPCA}, {CL_jo}
                     methods = NULL,                    ## String: cell typing methods: {Seurat}, {SpatialPCA}, {BASS}
                     Seurat_resolutions = NULL,         ## String: resolution
                     Seurat_dims = NULL,                ## String: PCA number (re)
                     Seurat_anchors = NULL,             ## String: number of neighbors to use when picking anchors
                     BASS_cellNum = NULL,               ## String: number of cell types. 
                     BASS_domainNum = NULL,             ## String: number of spatial domains.
                     SpatialPCA_domainNum = NULL,       ## String: desired cluster number
                     de_method = NULL,                  ## String: de methods: Wilcox 
                     out_path = NULL                    ## String: output path of de procedure
                     
){
  
  ## load inputs
  call_file <- paste0(data_path2, "/qc_call_file.txt")
  if(!file.exists(call_file)){
    
    stop("No QC file! Please run QC module!")
  } else {
    
    ### load count matrix
    qc_param <- read.table(call_file)[, 1]
    spatial_data_filename <- qc_param[1]
    sample_size <- qc_param[2]
  }
  ## Set cell typing module parameters
  if(grepl("CT", submodule) & submodule != "CT_Annot") {
    
    if (methods == "Seurat"){
      
      grid_pc <- ifelse(sample_size == 1,
                        paste0("cluster"),
                        paste0("cluster_anchor", Seurat_anchors,
                               "_dim", Seurat_dims))
      grid_use <- paste0(grid_pc, "_res", Seurat_resolutions)
    }
    post_file <- read.table(paste0(data_path1, "/ct_post_file.txt"))[,1]
    check_file <- c(spatial_data_filename, submodule, methods,
                    sample_size, grid_use, de_method, 
                    post_file[1], "ct")
  }
  ## Set clustering joint module parameters
  if(grepl("jo", submodule)){
    
    if(methods == "BASS"){
      
      grid_use <- paste0("R", BASS_domainNum, "_C", BASS_cellNum)
    }
    if(file.exists(paste0(data_path1, "/ct_post_file.txt"))){
      
      post_file <- read.table(paste0(data_path1, "/ct_post_file.txt"))[,1]
      check_file <- c(spatial_data_filename, submodule, methods,
                      sample_size, grid_use, de_method, 
                      post_file[1], "ct")
    } else {
      
      post_file <- read.table(paste0(data_path1, "/sdd_post_file.txt"))[,1]
      check_file <- c(spatial_data_filename, submodule, methods,
                      sample_size, grid_use, de_method, 
                      post_file[1], "sdd")
    }
    
  }
  
  ## Set spatial domain detection modules
  if(grepl("SDD", submodule)){
    
    if(methods == "SpatialPCA"){
      
      grid_use <- paste0("clust_", SpatialPCA_domainNum)
    }
    post_file <- read.table(paste0(data_path1, "/sdd_post_file.txt"))[,1]
    check_file <- c(spatial_data_filename, submodule, methods,
                    sample_size, grid_use, de_method, 
                    post_file[2], "sdd")
  }  
  
  ## Set annotation-based modules
  if(submodule == "CT_Annot"){
    
    grid_use <- "cluster"
    post_file <- read.table(paste0(data_path1, "/ct_post_file.txt"))[,1]
    check_file <- c(spatial_data_filename, submodule, methods,
                    sample_size, grid_use, de_method, 
                    post_file[1], "ct")
  }  
  
  ## output file
  write.table(check_file, 
              file = paste0(out_path, "/de_check_file.txt"), 
              row.names = F, quote = F, col.names = F)
  
  return(0)
}

# Function 2: de in Seurat
de.func <- function(st_list, 
                    cluster_df,
                    grid_use, 
                    de_method
){
  
  seurat_obj <- purrr::reduce(st_list[["count_list"]], function(x, y) {
    cbind(x = x, y = y)
  }) %>% CreateSeuratObject
  cluster_df <- cluster_df[, c("sample", "cell", grid_use)]
  colnames(cluster_df)[3] <- "cluster_label"
  ## add cluster information
  seurat_obj@meta.data$cluster_label <- cluster_df$cluster_label %>%
    factor(., levels = sort(unique(.)))
  Idents(seurat_obj) <- seurat_obj@meta.data$cluster_label
  ## identify markers
  markers <- FindAllMarkers(seurat_obj,
                            test.use = de_method,
                            only.pos = T,
                            min.pct = 0.1, 
                            return.thresh = 0.1,
                            logfc.threshold = 0)
  
  return(markers)
}

# Function 3: call function
de.call <- function(out_path                ## String: output path of de procedure
){
  
  source(paste0(method_path, "/io.R"))
  ## load inputs
  check_file <- read.table(paste0(out_path, "/de_check_file.txt"))[, 1]
  spatial_data_filename <- check_file[1]
  submodule <- check_file[2]
  methods <- check_file[3]
  sample_size <- check_file[4] %>% as.numeric
  grid_use <- check_file[5]
  de_method  <- check_file[6]
  ## load st and cluster data
  st_list <- h5data.load(spatial_data_filename,         
                         sample_size = sample_size,       
                         load_count = TRUE,     
                         normalization = FALSE,  
                         load_coord = FALSE,    
                         coordinate = TRUE)
  ## define markers in different modules
  result_file <- paste0(out_path, "/de_call_result.RData")
  clus_file <- paste0(out_path, "/de_call_cluster.RData")
  if(grepl("CT", submodule)) {
    
    ct_df <- fread2(check_file[7])
    ct_markers <- de.func(st_list, ct_df, grid_use, de_method)
    save(ct_markers, file = result_file)
    save(ct_df, file = clus_file)
  }
  if(grepl("jo", submodule)) {
    
    if(check_file[8] == "ct"){
      
      ct_df <- fread2(check_file[7])
      ct_markers <- de.func(st_list, ct_df, grid_use, de_method)
      save(ct_markers, file = result_file)
      save(ct_df, file = clus_file)
    } else {
      
      sdd_df <- fread2(check_file[7])
      sdd_markers <- de.func(st_list, sdd_df, grid_use, de_method)
      save(ct_markers, file = result_file)
      save(ct_df, file = clus_file)
    }
  }
  if(grepl("SDD", submodule)) {
    
    sdd_df <- fread2(check_file[7])
    sdd_markers <- de.func(st_list, sdd_df, grid_use, de_method)
    save(sdd_markers, file = result_file)
    save(sdd_df, file = clus_file)
  }
  ## output
  write.table(c(result_file, clus_file), 
              file = paste0(out_path, "/de_call_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}

# Function 4: post function
de.post <- function(out_path
){
  
  ## load data
  check_file <- read.table(paste0(out_path, "/de_check_file.txt"))[, 1]
  submodule <- check_file[2]
  methods <- check_file[3]
  call_file <- read.table(paste0(out_path, "/de_call_file.txt"))[, 1]
  load(call_file[1])
  clus_file <- call_file[2]
  
  ## process and output deg
  ct_de_mat_file <- sdd_de_mat_file <- "0"
  result_dir <- paste0(out_path, "/de_result/", submodule, "_", methods, "/")
  if (!file.exists(result_dir)) {
    
    system(paste0("mkdir -p ", result_dir))
  }
  if(check_file[8] == "ct") {
    
    ct_de_mat <- data.frame(gene = ct_markers$gene, 
                            cluster = ct_markers$cluster, 
                            log2FC = ct_markers$avg_log2FC, 
                            P = ct_markers$p_val,
                            p_adj = ct_markers$p_val_adj)
    ct_de_mat_file <- paste0(result_dir, "/ct_de_mat.txt")
    write.table(ct_de_mat, file = ct_de_mat_file, 
                sep = "\t", row.names = F, quote = F)
  }
  if(check_file[8] == "sdd") {
    
    sdd_de_mat <- data.frame(gene = sdd_markers$gene, 
                             cluster = sdd_markers$cluster, 
                             log2FC = sdd_markers$avg_log2FC, 
                             P = sdd_markers$p_val,
                             p_adj = sdd_markers$p_val_adj)
    sdd_de_mat_file <- paste0(result_dir, "/sdd_de_mat.txt")
    write.table(sdd_de_mat, file = sdd_de_mat_file, 
                sep = "\t", row.names = F, quote = F)
  }
  ## output
  write.table(c(ct_de_mat_file, sdd_de_mat_file, clus_file, submodule), 
              file = paste0(out_path, "/de_post_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}
