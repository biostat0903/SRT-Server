#! /usr/bin/env Rscript
# Identify differentially expressed gene in st data
# up-stream code: dr_cl.R, dr_cl_wr.R, dr_cl_sp.R
# down-stream code: de_plt.R

# Set method_path
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# Load packages
library(Seurat)
library(hdf5r)
library(dplyr)
library(bigreadr)

# Function 1: check function
de.check <- function(data_path1 = NULL,               ## String: output path of dr_cl procedure
                     data_path2 = NULL,               ## String: output path of qc procedure
                     clus_mode = NULL,                ## String: dr_cl, dr_cl_wr or dr_cl_sp
                     resolution = NULL,               ## Float: resolution of dr_cl
                     anchor_num = NULL,               ## Integer: number of neighbors to use when picking anchors in dr_cl
                     dim_num = NULL,                  ## Integer: PCA number (re) of dr_cl
                     cluster_num = NULL,              ## Integer: number of cell types in dr_cl_wr or dr_cl_sp
                     domain_num = NULL,               ## Integer: number of spatial domains in dr_cl_warpping
                     spc_num = NULL,                  ## Integer: number of spatial PC in dr_cl_sp
                     de_method = "wilcox",            ## de methods: Wilcox 
                     out_path = NULL                  ## String: output path of de procedure
                     
){
  
  ## load inputs
  post_file <- read.table(paste0(data_path1, "/", clus_mode, "_post_file.txt"))[,1]
  dr_cl_cluster <- post_file[1]
  dr_cl_pc_path <- post_file[length(post_file)]
  cluster_df <- fread2(dr_cl_cluster)
  check_file <- paste0(data_path2, "/qc_call_file.txt")
  if(!file.exists(check_file)){
    
    stop("No QC file! Please run QC module!")
  } else {
    
    ### load count matrix
    qc_param <- read.table(check_file)[, 1]
    spatial_data_filename <- qc_param[1]
    sample_size <- qc_param[2]
  }

  if (!file.exists(dr_cl_cluster)) {
    
    stop("Cluster file is not found!")
  } else {
    
    if (clus_mode == "dr_cl") {
      
      grid_pc <- ifelse(sample_size == 1,
                        paste0("cluster"),
                        paste0("cluster_anchor", anchor_num,
                               "_dim", dim_num))
      grid_use <- paste0(grid_pc, "_res", resolution)
    } else {
      
      grid_use <- ifelse(clus_mode == "dr_cl_wr", 
                         paste0("R", domain_num, "_C", cluster_num),
                         paste0("clust_", cluster_num))
    }
    #
    if (!all(c("sample", "cell", grid_use) %in% colnames(cluster_df)))  {
      
      stop("Please ensure the cluster file to contain sample, cells and cluster column!")
    } else {
      
      cluster_df <- cluster_df[,c("sample", "cell", grid_use)]
      colnames(cluster_df)[3] <- "cluster_label"
    }
  }  


  ### output file
  write.table(c(spatial_data_filename, dr_cl_cluster, 
                sample_size, clus_mode, grid_use, de_method), 
              file = paste0(out_path, "/de_check_file.txt"), 
              row.names = F, quote = F, col.names = F)
  
  return(0)
}
  
# Function 2: call function
de.call <- function(out_path                ## String: output path of de procedure
){
  
  source(paste0(method_path, "/io.R"))
  ## load inputs
  check_file <- read.table(paste0(out_path, "/de_check_file.txt"))[, 1]
  spatial_data_filename <- check_file[1]
  dr_cl_cluster_file <- check_file[2]
  sample_size <- check_file[3]
  clus_mode <- check_file[4]
  grid_use <- check_file[5] 
  de_method <- check_file[6]

  ## load st and cluster data
  st_list <- h5data.load(spatial_data_filename,         
                         sample_size = sample_size,       
                         load_count = TRUE,     
                         normalization = FALSE,  
                         load_coord = FALSE,    
                         coordinate = TRUE)
  seurat_obj <- purrr::reduce(st_list[["count_list"]], function(x, y) {
    cbind(x = x, y = y)
  }) %>% CreateSeuratObject
  cluster_df <- fread2(dr_cl_cluster_file)
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
  
  ## output
  result_path <- paste0(out_path, "/de_call_result.RData")
  save(markers, file = result_path)
  clus_file <- paste0(out_path, "/de_call_cluster.RData")
  save(cluster_df, file = clus_file)
  write.table(c(result_path, clus_file, clus_mode), 
              file = paste0(out_path, "/de_call_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}

# Function 3: post function
de.post <- function(out_path
){
  
  ## load data
  check_file <- read.table(paste0(out_path, "/de_check_file.txt"))[, 1]
  clus_mode <- check_file[4]
  call_file <- read.table(paste0(out_path, "/de_call_file.txt"))[, 1]
  load(call_file[1])
  clus_file <- call_file[2]

  ## process deg
  de_mat <- data.frame(gene = markers$gene, 
                       cluster = markers$cluster, 
                       log2FC = markers$avg_log2FC, 
                       P = markers$p_val,
                       p_adj = markers$p_val_adj)
  de_mat_sig <- subset(de_mat,
                       abs(de_mat$log2FC) >= 0.25 & de_mat$p_adj < 0.05)
  
  ## output 
  result_dir <- paste0(out_path, "/de_result/", clus_mode)
  if (!file.exists(result_dir)) {
    
    system(paste0("mkdir -p ", result_dir))
  }
  de_mat_file <- paste0(result_dir, "/de_mat.txt")
  de_mat_sig_file <- paste0(result_dir, "/de_mat_sig.txt")
  
  ## output
  write.table(de_mat, file = de_mat_file, 
              sep = "\t", row.names = F, quote = F)
  write.table(de_mat_sig, file = de_mat_sig_file, 
              sep = "\t", row.names = F, quote = F)
  write.table(c(de_mat_file, de_mat_sig_file, clus_file, clus_mode), 
              file = paste0(out_path, "/de_post_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}


# ###################
# ### test code
# data_path1 <- "/net/mulan/disk2/yasheng/stwebProject/test/dr_cl_wr"
# data_path2 <- "/net/mulan/disk2/yasheng/stwebProject/test/qc"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/test/de"
# de.check (data_path1 = data_path1,
#           data_path2 = data_path2,
#           clus_mode = "dr_cl_wr",
#           anchor_num = NULL,
#           dim_num = NULL,
#           resolution = NULL,
#           cluster_num = 10,
#           domain_num = 5,
#           spc_num = NULL,
#           de_method = "wilcox",
#           out_path = output_path)
# de.call(out_path = output_path)
# de.post(out_path = output_path)

