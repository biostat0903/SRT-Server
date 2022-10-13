#! /usr/bin/env Rscript
# Trajectory inference
# up-stream code: dr_cl.R, dr_cl_wrapping.R, dr_cl_sp.R
# down-stream procedure: None

# Set method_path
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# Load packages
library(slingshot)
library(SingleCellExperiment)
library(tradeSeq)
library(BiocParallel)
library(Seurat)
library(hdf5r)
library(stringr)
library(dplyr)
library(tidyr)
library(bigreadr)

# Fix parameter
NCores <- 10

# Function 1: check function
traj.check <- function(data_path1 = NULL,               ## String: output path of clustering procedure
                       data_path2 = NULL,               ## String: output path of qc procedure
                       clus_mode = NULL,                ## String: cluster mode, including dr_cl, dr_cl_wr or dr_cl_sp
                       resolution = NULL,               ## Float: resolution of dr_cl
                       anchor_num = NULL,               ## Integer: number of neighbors to use when picking anchors in dr_cl
                       dim_num = NULL,                  ## Integer: PCA number (re) of dr_cl
                       cluster_num = NULL,              ## Integer: number of cell types in dr_cl_wr or dr_cl_sp
                       domain_num = NULL,               ## Integer: number of spatial domains in dr_cl_warpping
                       spc_num = NULL,                  ## Integer: number of spatial PC in dr_cl_sp
                       start_clust = NULL,              ## Integer: pecific start cluster in traj analysis
                       out_path = NULL                  ## String: output path of traj procedure
){
  
  ## check file
  check_file <- paste0(data_path2, "/qc_call_file.txt")
  if(!file.exists(check_file)){
    
    stop("No QC file! Please run QC module!")
  }else {
    
    qc_param <- read.table(check_file)[, 1]
    spatial_data_filename <- qc_param[1]
    sample_size <- qc_param[2] %>% as.numeric
  }
  dr_cl_post <- read.table(paste0(data_path1, "/", clus_mode,
                                  "_post_file.txt"))[,1]
  dr_cl_cluster_file <- dr_cl_post[1]
  dr_cl_pc_file <- dr_cl_post[length(dr_cl_post)]
  if (!file.exists(dr_cl_cluster_file) | !file.exists(dr_cl_pc_file)) {
    
    stop("dr_cl file is not found!")
  } 
  if(clus_mode == "dr_cl"){
    
    if (sample_size == 1){
      
      anchor_num <- "NULL"
      dim_num <- "NULL"
    } 
    cluster_num <- domain_num <- spc_num <- "NULL"
  }
  if(clus_mode == "dr_cl_sp"){
    
    anchor_num <- resolution <- dim_num <- "NULL"
    cluster_num <- domain_num <- "NULL"
  }
  if(clus_mode == "dr_cl_wr"){
    
    anchor_num <- resolution <- dim_num <- "NULL"
    spc_num <- "NULL"
  }

  ## output file
  write.table(c(spatial_data_filename, sample_size, clus_mode,
                anchor_num, dim_num, resolution,
                cluster_num, domain_num, spc_num,
                start_clust, dr_cl_cluster_file), 
              file = paste0(out_path, "/traj_check_file.txt"), 
              row.names = F, quote = F, col.names = F)
  
  return(0)
}

# Function 2: call function
traj.call <- function(data_path,                        ## String: output path of clustering procedure
                      out_path = NULL                   ## String: output path of traj procedure
){
  
  ## Load io code
  source(paste0(method_path, "/io.R"))
  
  ## load data
  check_file <- read.table(paste0(out_path, "/traj_check_file.txt"))[, 1]
  spatial_data_filename <- check_file[1]
  sample_size <- check_file[2] %>% as.numeric()
  clus_mode <- check_file[3]
  anchor_num <- ifelse(check_file[4] == "NULL", 0, check_file[4])
  dim_num <- ifelse(check_file[5] == "NULL", 0, as.numeric(check_file[5]))
  resolution <- ifelse(check_file[6] == "NULL", 0, check_file[6])
  cluster_num <- ifelse(check_file[7] == "NULL", 0, as.numeric(check_file[7]))
  domain_num <- ifelse(check_file[8] == "NULL", 0, as.numeric(check_file[8]))
  spc_num <- ifelse(check_file[9] == "NULL", 0, as.numeric(check_file[9]))
  start_clust <- check_file[10] %>% 
    strsplit(",") %>% unlist %>% as.numeric
  dr_cl_cluster_file <- check_file[11]
  
  ## build sce object
  st_list <- h5data.load(spatial_data_filename,         
                         sample_size = sample_size,       
                         load_count = TRUE,     
                         normalization = FALSE,  
                         load_coord = TRUE,    
                         coordinate = TRUE)
  count_merged <- purrr::reduce(st_list[["count_list"]], function(x, y) {
    cbind(x = x, y = y)
  })
  coord_merged <- Reduce(rbind, st_list[["coord_list"]]) %>% 
    as.data.frame()
  ## create sce object
  sce_obj <- SingleCellExperiment(assays = List(counts = count_merged))
  colData(sce_obj)$x <- coord_merged[,1]
  colData(sce_obj)$y <- coord_merged[,2]

  ## obtain cluster information
  cluster_df <- fread2(dr_cl_cluster_file)
  sample_size <- length(unique(cluster_df$sample))
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
  cluster_df <- cluster_df[,c("sample", "cell", grid_use)]
  colnames(cluster_df)[3] <- "cluster_label"
  if (!all(colnames(sce_obj) %in% cluster_df$cell)) {
    
    stop("Number of cells in expression matrix and cluster label file do not match!")
  } else {
    
    cluster_df <- cluster_df[match(colnames(sce_obj),
                                   cluster_df$cell),]
    colData(sce_obj)$sample <- cluster_df$sample
    colData(sce_obj)$cell <- cluster_df$cell
    colData(sce_obj)$cluster_label <- cluster_df$cluster_label
  }
  
  pc_data_str <- load(paste0(data_path, "/", clus_mode, "_pc.RData"))
  if (clus_mode == "dr_cl") {
    
    pc_list <- eval(parse(text = pc_data_str))
    if (!grid_pc %in% names(pc_list)) {
      
      stop(paste0(grid_pc," not in PCs file!"))
    } else {
      
      pc_mtx <- pc_list[[grid_pc]]
    }
  } else {
    
    pc_mtx <- eval(parse(text = pc_data_str))
  }
  pc_mtx <- pc_mtx[,match(colnames(sce_obj),
                          colnames(pc_mtx))]
  reducedDims(sce_obj) <- SimpleList(DRM = t(pc_mtx))
  
  ## add UMAP for dr_cl
  if (clus_mode == "dr_cl") {
    
    dr_cl_umap_path <- dr_cl_post_file[2]
    umap_list_str <- load(dr_cl_umap_path)
    # check gridx in pc_list
    umap_list <- eval(parse(text = umap_list_str))
    if (!grid_pc %in% names(pc_list)) {
      
      stop(paste0("ERROR: ", grid_pc," not in UMAP file!"))
    } else {
      
      umap_df <- umap_list[[grid_pc]]
    }
    
    ## check match between pc and expr data
    if (!all(colnames(sce_obj) %in% rownames(umap_df))) {
      
      stop("Number of cells in expression matrix and UMAP file do not match!")
    } else {
      
      umap_df <- umap_df[colnames(sce_obj),]
      colData(sce_obj)$UMAP_1 <- umap_df[,1]
      colData(sce_obj)$UMAP_2 <- umap_df[,2]
    }
    print("UMAP information added!")
  }
  
  ## perform traj on each start point
  if (all(is.na(start_clust))) {
    
    start_clust_all <- unique(cluster_df$cluster_label)
  } else {
    
    start_clust_all <- start_clust
  }
  pseudotime_list <- list()
  ATres_list <- list()
  for (start_clustx in start_clust_all) {
    
    sce_obj  <- slingshot(sce_obj, 
                          clusterLabels = 'cluster_label', 
                          reducedDim = 'DRM',
                          start.clus = start_clustx)
    traj_ind <- grep("slingPseudotime", names(sce_obj@colData@listData))
    traj_mat <- Reduce("cbind", sce_obj@colData@listData[traj_ind]) %>% 
      as.data.frame()
    dimnames(traj_mat) <- list(colnames(sce_obj), 
                               names(sce_obj@colData@listData)[traj_ind])
    pseudotime_list[[paste0("start", start_clustx)]] <- traj_mat
    
    # fit negative binomial GAM
    multicoreParam <- MulticoreParam(workers = NCores)
    sce_obj_gam <- fitGAM(sce_obj, 
                          parallel = T, 
                          BPPARAM = multicoreParam)
    # test for dynamic expression
    ATres_list[[paste0("start", start_clustx)]] <- associationTest(sce_obj_gam)
    
    print(paste0("Trajectory inference start from ", start_clustx, " is ok!"))
  }
  
  traj_result <- list()
  traj_result[["pseudotime"]] <- pseudotime_list
  traj_result[["ATres"]] <- ATres_list
  traj_result[["coord_df"]] <- colData(sce_obj)[,c("x","y")] %>% 
    as.data.frame()
  traj_result[["cluster_label"]] <- colData(sce_obj)$cluster_label
  traj_result[["sample"]] <- colData(sce_obj)[,c("sample","cell")] %>% 
    as.data.frame()
  if (clus_mode == "dr_cl") {
    
    traj_result[["umap"]] <- colData(sce_obj)[,c("UMAP_1","UMAP_2")] %>% 
      as.data.frame()
  }
  
  # output
  result_path <- paste0(out_path, "/traj_call_", clus_mode, "_result.RData")
  save(traj_result, file = result_path)
  write.table(c(result_path, clus_mode), 
              file = paste0(out_path, "/traj_call_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}

# Function 3: post function
traj.post <- function(out_path = NULL                   ## String: output path of traj procedure
){
  
  # load
  call_file <- read.table(paste0(out_path, "/traj_call_file.txt"))[, 1]
  traj_result_file <- call_file[1]
  clus_mode <- call_file[2]
  load(traj_result_file) # traj_result
  
  # output
  result_dir <- paste0(out_path, "/traj_result/", clus_mode, "/")
  if (!file.exists(result_dir)) {
    system(paste0("mkdir -p ", result_dir))
  }
  
  result_out <- plyr::laply(seq_along(traj_result[["pseudotime"]]), function(a){
    
    traj_a <- cbind(traj_result[["sample"]],
                    traj_result[["pseudotime"]][[a]])
    write.csv(traj_a, 
              file = paste0(result_dir, "traj_start", a,".csv"), 
              row.names = F, quote = F)
    return(paste0(result_dir, "traj_start", a,".csv"))
  })
  
  return(0)
}


# ###################
# ### test code
# data_path1 <- "/net/mulan/disk2/yasheng/stwebProject/test/dr_cl_wr"
# data_path2 <- "/net/mulan/disk2/yasheng/stwebProject/test/qc"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/test/ccc"
# traj.check(data_path1 = data_path1,
#            data_path2 = data_path2,
#            clus_mode = "dr_cl_wr",
#            anchor_num = NULL,
#            dim_num = NULL,
#            resolution = NULL,
#            cluster_num = 10,
#            domain_num = 5,
#            start_clust = "1,2",
#            out_path = output_path)
# traj.call(data_path = data_path1, out_path = output_path)
# traj.post(out_path = output_path)
