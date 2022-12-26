#! /usr/bin/env Rscript
# Trajectory inference
# up-stream code: ct.R, sdd.R
# down-stream procedure: None

# Set method_path
method_path <- "/applicatoins/docker-mnt/scripts/ST/dev"

# Load packages
library(slingshot)
library(SingleCellExperiment)
library(tradeSeq)
library(BiocParallel)
library(Seurat)
library(hdf5r)
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
library(bigreadr)

# Fix parameter
N_CORES <- 4

# Function 1: check function
traj.check <- function(data_path1 = NULL,               ## String: output path of clustering procedure
                       data_path2 = NULL,               ## String: output path of qc procedure
                       submodule = NULL,                ## String: cell typing sub-modules: {CT_PCA}, {SDD_sPCA}, {CL_jo}
                       methods = NULL,                  ## String: cell typing methods: {Seurat}, {SpatialPCA}, {BASS}
                       Seurat_resolutions = NULL,       ## String: resolution
                       Seurat_dims = NULL,              ## String: PCA number (re)
                       Seurat_anchors = NULL,           ## String: number of neighbors to use when picking anchors
                       BASS_cellNum = NULL,             ## String: number of cell types. 
                       BASS_domainNum = NULL,           ## String: number of spatial domains.
                       SpatialPCA_domainNum = NULL,     ## String: desired cluster number
                       start_ct = NULL,                 ## Integer: specific start cluster(s) in traj analysis
                       start_sdd = NULL,                ## Integer: specific start domain(s) in traj analysis
                       out_path = NULL                  ## String: output path of traj procedure
){
  
  ## check file
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
  
  #
  ## Set cell typing module parameters
  if(grepl("CT", submodule)) {
    
    if (methods == "Seurat"){
      post_file <- read.table(paste0(data_path1, "/ct_post_file.txt"))[,1]
      grid_pc <- ifelse(sample_size == 1,
                        paste0("cluster"),
                        paste0("cluster_anchor", Seurat_anchors,
                               "_dim", Seurat_dims))
      grid_use <- paste0(grid_pc, "_res", Seurat_resolutions)
      umap_file <- post_file[2]
      pc_file <- post_file[3]
    }
    check_file <- c(spatial_data_filename, submodule, methods,
                    sample_size, grid_use, 
                    start_ct, NA,
                    post_file[1], NA,
                    pc_file, umap_file)
  }
  ## Set clustering joint module parameters
  if(grepl("jo", submodule)){
    
    if(methods == "BASS"){
      
      grid_use <- paste0("R", BASS_domainNum, "_C", BASS_cellNum)
    }
    if(file.exists(paste0(data_path1, "/ct_post_file.txt"))){
      
      post_file <- read.table(paste0(data_path1, "/ct_post_file.txt"))[,1]
    } else {
      
      post_file <- read.table(paste0(data_path1, "/sdd_post_file.txt"))[,1]
    }
    check_file <- c(spatial_data_filename, submodule, methods,
                    sample_size, grid_use, 
                    start_ct, start_sdd,
                    post_file[1], post_file[2],
                    post_file[3], NA)
  }
  
  ## Set spatial domain detection modules
  if(grepl("SDD", submodule)){
    
    if(methods == "SpatialPCA"){
      
      grid_use <- paste0("clust_", SpatialPCA_domainNum)
    }
    post_file <- read.table(paste0(data_path1, "/sdd_post_file.txt"))[,1]
    check_file <- c(spatial_data_filename, submodule, methods,
                    sample_size, grid_use, 
                    NA, start_sdd, 
                    NA, post_file[2],
                    post_file[3], NA)
  }  
  
  ## output file
  write.table(check_file, file = paste0(out_path, "/traj_check_file.txt"), 
              row.names = F, quote = F, col.names = F)
  
  return(0)
}

# Function 2: traj test use slingshot
traj.func <- function(st_list = NULL,
                      cluster_df = NULL,
                      grid_use = NULL,
                      methods = NULL,
                      pc_file = NA,
                      umap_file = NA,
                      start_clust = NULL,
                      n_cores = N_CORES){
  cat(methods, "\n")
  ## 2.1 create sce object from st_list
  count_merged <- purrr::reduce(st_list[["count_list"]], function(x, y) {
    cbind(x = x, y = y)
  })
  coord_merged <- Reduce(rbind, st_list[["coord_list"]]) %>% 
    as.data.frame()
  
  sce_obj <- SingleCellExperiment(assays = List(counts = count_merged))
  colData(sce_obj)$x <- coord_merged[,1]
  colData(sce_obj)$y <- coord_merged[,2]
  
  ## 2.2 obtain cluster information
  if (!all(colnames(sce_obj) %in% cluster_df$cell)) {
    
    stop("Number of cells in expression matrix and cluster label file do not match!")
  } else {
    
    cluster_df <- cluster_df[match(colnames(sce_obj),
                                   cluster_df$cell),]
    colData(sce_obj)$sample <- cluster_df$sample
    colData(sce_obj)$cell <- cluster_df$cell
    colData(sce_obj)$cluster_label <- cluster_df[[grid_use]]
  }
  
  ## 2.3 obtain pc information
  # load pc matrix
  pc_data_str <- load(pc_file)
  if (methods == "Seurat") {
    
    pc_list <- eval(parse(text = pc_data_str))
    grid_pc <- strsplit(grid_use, "_res")[[1]][1]
    if (!grid_pc %in% names(pc_list)) {
      
      stop(paste0(grid_pc," not in PCs file!"))
    } else {
      
      pc_mtx <- pc_list[[grid_pc]]
    }
  } 
  
  if (methods %in% c("BASS", "SpatialPCA")){
    
    pc_mtx <- eval(parse(text = pc_data_str))
  }
  # add to sce object
  pc_mtx <- pc_mtx[, match(colnames(sce_obj), colnames(pc_mtx))]
  reducedDims(sce_obj) <- SimpleList(DRM = t(pc_mtx))
  
  ## 2.4 add UMAP for dr_cl
  if (methods == "Seurat" & !is.na(umap_file)) {
    
    umap_list_str <- load(umap_file)
    
    # check gridx in pc_list
    umap_list <- eval(parse(text = umap_list_str))
    grid_pc <- strsplit(grid_use, "_res")[[1]][1]
    
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
    
    start_clust_all <- unique(colData(sce_obj)$cluster_label)
  } else {
    
    start_clust_all <- start_clust
  }
  
  pseudotime_list <- list()
  ATres_list <- list()
  for (start_clustx in start_clust_all) {
    
    sce_objx  <- slingshot(sce_obj, 
                           clusterLabels = 'cluster_label', 
                           reducedDim = 'DRM',
                           start.clus = start_clustx)
    traj_ind <- grep("slingPseudotime", names(sce_objx@colData@listData))
    traj_mat <- Reduce("cbind", sce_objx@colData@listData[traj_ind]) %>% 
      as.data.frame()
    dimnames(traj_mat) <- list(colnames(sce_objx), 
                               names(sce_objx@colData@listData)[traj_ind])
    pseudotime_list[[paste0("start", start_clustx)]] <- traj_mat
    
    # fit negative binomial GAM
    multicoreParam <- MulticoreParam(workers = n_cores)
    sce_obj_gam <- fitGAM(sce_objx, 
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
  if (methods == "Seurat") {
    
    traj_result[["umap"]] <- colData(sce_obj)[,c("UMAP_1","UMAP_2")] %>% 
      as.data.frame()
  }
  
  return(traj_result)
}


# Function 2: call function
traj.call <- function(out_path                           ## String: output path of traj procedure
){
  
  ## Load io code
  source(paste0(method_path, "/io.R"))
  
  ## load inputs
  check_file <- read.table(paste0(out_path, "/traj_check_file.txt"))[, 1]
  spatial_data_filename <- check_file[1]
  submodule <- check_file[2]
  methods <- check_file[3]
  sample_size <- check_file[4] %>% as.numeric
  grid_use <- check_file[5]
  pc_file <- check_file[10]
  umap_file <- check_file[11]
  
  ## build sce object
  st_list <- h5data.load(spatial_data_filename,         
                         sample_size = sample_size,       
                         load_count = TRUE, 
                         normalization = FALSE,  
                         load_coord = TRUE,    
                         coordinate = TRUE)
  
  ## define ccc in different modules
  result_file <- paste0(out_path, "/traj_call_result.RData")
  clus_file <- paste0(out_path, "/traj_call_cluster.RData")
  if(grepl("CT", submodule)) {
    #
    start_ct <- check_file[6] %>% 
      strsplit(",") %>% unlist %>% as.numeric
    #
    ct_df <- fread2(check_file[8])
    ct_traj_result <- traj.func(st_list = st_list,
                                cluster_df = ct_df,
                                grid_use = grid_use,
                                methods = methods,
                                pc_file = pc_file,
                                umap_file = umap_file,
                                start_clust = start_ct)
    
    save(ct_traj_result, file = result_file)
    save(ct_df, file = clus_file)
  }
  if(grepl("jo", submodule)) {
    #
    start_ct <- check_file[6] %>% 
      strsplit(",") %>% unlist %>% as.numeric
    ct_df <- fread2(check_file[8])
    ct_traj_result <- traj.func(st_list = st_list,
                                cluster_df = ct_df,
                                grid_use = grid_use,
                                methods = methods,
                                pc_file = pc_file,
                                umap_file = umap_file,
                                start_clust = start_ct)
    #
    start_sdd <- check_file[7] %>% 
      strsplit(",") %>% unlist %>% as.numeric
    sdd_df <- fread2(check_file[9])
    sdd_traj_result <- traj.func(st_list = st_list,
                                 cluster_df = sdd_df,
                                 grid_use = grid_use,
                                 methods = methods,
                                 pc_file = pc_file,
                                 umap_file = umap_file,
                                 start_clust = start_sdd)
    save(ct_traj_result, sdd_traj_result, file = result_file)
    save(ct_df, sdd_df, file = clus_file)
  }
  if(grepl("SDD", submodule)) {
    
    #
    start_sdd <- check_file[7] %>% 
      strsplit(",") %>% unlist %>% as.numeric
    sdd_df <- fread2(check_file[9])
    sdd_traj_result <- traj.func(st_list = st_list,
                                 cluster_df = sdd_df,
                                 grid_use = grid_use,
                                 methods = methods, 
                                 pc_file = pc_file,
                                 start_clust = start_sdd)
    save(sdd_traj_result, file = result_file)
    save(sdd_df, file = clus_file)
  }
  
  ## output
  write.table(c(result_file, clus_file), 
              file = paste0(out_path, "/traj_call_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
  
}

# Function 3: post function
traj.post <- function(out_path = NULL                    ## String: output path of traj procedure
){
  
  ## 1. load data
  check_file <- read.table(paste0(out_path, "/traj_check_file.txt"))[, 1]
  submodule <- check_file[2]
  methods <- check_file[3]
  start_ct <- check_file[6]
  start_sdd <- check_file[7]
  
  call_file <- read.table(paste0(out_path, "/traj_call_file.txt"),
                          header = F, sep = "\t")[, 1]
  load(call_file[1])
  clus_file <- call_file[2]
  
  ## 2. output 
  ct_pseudo_result_file <- sdd_pseudo_result_file <- 
    ct_ATres_result_file <- sdd_ATres_result_file <- NA
  result_dir <- paste0(out_path, "/traj_result/", submodule, "/")
  if (!file.exists(result_dir)) {
    system(paste0("mkdir -p ", result_dir))
  }
  
  
  if(grepl("CT", submodule) | grepl("jo", submodule)) {
    
    # ct pseudo reslut
    ct_pseudo_result_file <- plyr::laply(names(ct_traj_result[["pseudotime"]]), function(a){
      
      ct_pseudo_result_file_a <- paste0(result_dir, "ct_pseudo_", a,".txt")
      pseudo_a <- cbind(ct_traj_result[["sample"]],
                        ct_traj_result[["coord_df"]],
                        ct_traj_result[["cluster_label"]],
                        ct_traj_result[["pseudotime"]][[a]])
      colnames(pseudo_a) <- c("sample", "cell", "x", "y", "cluster_label",
                              names(ct_traj_result[["pseudotime"]][[a]]))
      
      # add umap infor
      if (methods == "Seurat") {
        pseudo_a$UMAP_1 <- ct_traj_result[["umap"]][,1]
        pseudo_a$UMAP_2 <- ct_traj_result[["umap"]][,2]
      }
      
      write.table(pseudo_a, file = ct_pseudo_result_file_a, 
                  sep = "\t", row.names = F, quote = F)
      
      return(ct_pseudo_result_file_a)
    }) %>% paste(., collapse = ",")
    
    # ct ATres result
    ct_ATres_result_file <- plyr::laply(names(ct_traj_result[["ATres"]]), function(a){
      
      ct_ATres_result_file_a <- paste0(result_dir, "ct_ATres_", a,".txt")
      ATres_a <- rownames_to_column(ct_traj_result[["ATres"]][[a]], var = "gene")
      ATres_a <- ATres_a[, c("gene", "meanLogFC", "waldStat", "pvalue")]
      write.table(ATres_a, file = ct_ATres_result_file_a, 
                  sep = "\t", row.names = F, quote = F)
      
      return(ct_ATres_result_file_a)
    }) %>% paste(., collapse = ",")
    
  }
  if(grepl("SDD", submodule) | grepl("jo", submodule)) {
    
    sdd_pseudo_result_file <- plyr::laply(names(sdd_traj_result[["pseudotime"]]), function(a){
      
      sdd_pseudo_result_file_a <- paste0(result_dir, "sdd_pseudo_", a,".txt")
      pseudo_a <- cbind(sdd_traj_result[["sample"]],
                        sdd_traj_result[["coord_df"]],
                        sdd_traj_result[["cluster_label"]],
                        sdd_traj_result[["pseudotime"]][[a]])
      colnames(pseudo_a) <- c("sample", "cell", "x", "y", "cluster_label",
                              names(sdd_traj_result[["pseudotime"]][[a]]))
      write.table(pseudo_a, file = sdd_pseudo_result_file_a, 
                  sep = "\t", row.names = F, quote = F)
      
      return(sdd_pseudo_result_file_a)
    }) %>% paste(., collapse = ",")
    
    # ct ATres result
    sdd_ATres_result_file <- plyr::laply(names(sdd_traj_result[["ATres"]]), function(a){
      
      sdd_ATres_result_file_a <- paste0(result_dir, "sdd_ATres_", a,".txt")
      ATres_a <- rownames_to_column(sdd_traj_result[["ATres"]][[a]], var = "gene")
      ATres_a <- ATres_a[, c("gene", "meanLogFC", "waldStat", "pvalue")]
      write.table(ATres_a, file = sdd_ATres_result_file_a, 
                  sep = "\t", row.names = F, quote = F)
      return(sdd_ATres_result_file_a)
    }) %>% paste(., collapse = ",")
    
  }
  
  write.table(c(ct_pseudo_result_file, sdd_pseudo_result_file, 
                ct_ATres_result_file, sdd_ATres_result_file,
                start_ct, start_sdd,
                methods, submodule), 
              file = paste0(out_path, "/traj_post_file.txt"), 
              sep = "\t", col.names = F, row.names = F, quote = F)
  return(0)
}


# ###################
# ### test code
# data_path1 <- "/pt_data/494433291@qq.com/c8ab46c7beb047628f0926ea5f14bfbe/tools-output/wf-483597487144174144/job-SDD-483597743558754880"
# data_path2 <- "/pt_data/494433291@qq.com/c8ab46c7beb047628f0926ea5f14bfbe/tools-output/wf-483597487144174144/job-QC-483597512712651328"
# output_path <- "/pt_data/494433291@qq.com/c8ab46c7beb047628f0926ea5f14bfbe/tools-output/wf-483597487144174144/job-TRAJ-483643761960682048"
# traj.check(data_path1 = data_path1,
#           data_path2 = data_path2,
#           submodule = "SDD_sPCA",
#           methods = "SpatialPCA",
#           SpatialPCA_domainNum = 3,
#           start_ct = "1",
#           start_sdd = "3",
#           out_path = output_path)
# traj.call(out_path = output_path)
# traj.post(out_path = output_path)
