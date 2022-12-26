#! /usr/bin/env Rscript
# performing dimensional reduction and clustering for spatial transcriptomics data
# up-stream code: qc.R
# down-stream plot: ct_plt.R
# down-stream code: de.R, ccc.R, traj.R

# Set method_path
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# Load packages
library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(dplyr)
library(tidyr)
library(stringr)
library(glmGamPoi)
library(BASS)
library(SpatialPCA)
library(umap)
library(Rtsne)
library(ggplot2)
library(bigreadr)

# Fix parameters
DIMS = 20              ## BASS: number of dimension for PCA
GENENUM = 2000         ## BASS: number of gene for PCA 
BURNIN = 2000          ## BASS: number of burn-in
SAMPLES = 10000        ## BASS: number of MCMC iteration 
NORM = TRUE       
SCALE = TRUE
BATCH = FALSE
N_SPC = 20             ## SpatialPCA: number of pcs used in single sample
                       ##             (same as fixed in multiple samples)

# Function 1: ct.check
sdd.check <- function(data_path,                               ## String: output path of qc procedure
                      sdd_submodule = NULL,                    ## String: spatial domain detection sub-modules: {SDD_sPCA}, {CL_jo}
                      methods = NULL,                          ## String: spatial domain detection methods: {SpatialPCA}, {BASS}
                      SpatialPCA_geneType = "spatial",         ## String: The type of genes to be used
                                                               ##         Three selections: {spatial}, {hvg}, {custom}
                      SpatialPCA_customGene = "NULL",          ## String: path for user specified genes. 
                                                               ##         required in gene_type=="custom".
                      SpatialPCA_sparkVer = "NULL",            ## String: selection spark version: {spark}, {sparkx}
                                                               ##         required in gene_type=="spatial".
                      SpatialPCA_clusMethod = "louvain",       ## String: clustering method
                                                               ##         Two selections: {walktrap}, {louvain}
                      SpatialPCA_domainNum = NULL,             ## String: desired cluster number
                      SpatialPCA_kernel = "gaussian",          ## String: kernel for SpatialPCA
                                                               ##         Four selections: {gaussian}, {cauchy}, {quadratic}, {delaunday}
                      SpatialPCA_coreNum = 1,                  ## Integer: The core numbers for SpatialPCA
                      BASS_cellNum = NULL,                     ## String: number of cell types. 
                      BASS_domainNum = NULL,                   ## String: number of spatial domains.
                      BASS_betaMethod = "SW",                  ## String: fix the cell-cell interaction parameter. 
                                                               ##         to be beta {fix} or estimate the parameter based on the 
                                                               ##         data at hand {SW}.
                      BASS_initMethod = "mclust",              ## String: Initialize the cell type clusters and spatial domains 
                                                               ##         Two selections: {kmeans}, {mclust}.
                      BASS_geneSelect = "sparkx",              ## String: feature selection method: spatial gene {sparkx} or
                                                               ##         high variable gene {hvgs}.
                      annot_method = "None",                   ## String: annotation methods: {None}, {scSorter}, {Garnett}
                      out_path                                 ## String: output path for ct procedure
){
  
  ## Set parameters and check QC module
  sdd_methods <- paste(sdd_submodule, methods, sep="-")
  qc_param <- read.table(paste0(data_path, "/qc_call_file.txt"))[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  if(!file.exists(spatial_data_filename)){
    
    stop("Please run the QC module!")
  }
  
  ## Output the first method: CT_PCA-Seurat
  if (sdd_methods == "SDD_sPCA-SpatialPCA"){

    write.table(c(sdd_methods, annot_method, 
                  SpatialPCA_geneType, SpatialPCA_customGene, SpatialPCA_sparkVer,
                  SpatialPCA_clusMethod, SpatialPCA_domainNum, SpatialPCA_kernel,
                  SpatialPCA_coreNum), 
                file = paste0(out_path, "/sdd_check_file.txt"), 
                row.names = F, quote = F, col.names = F)
  }
  
  ## Output the second method: CL_jo-BASS
  if (sdd_methods == "CL_jo-BASS"){

    write.table(c(sdd_methods, annot_method, 
                  BASS_cellNum, BASS_domainNum,
                  BASS_initMethod, BASS_betaMethod, BASS_geneSelect),
                file = paste0(out_path, "/sdd_check_file.txt"), 
                row.names = F, col.names = F, quote = F)
  }
  
  return(0)
}

# Function 2: SpatialPCA call function
SpatialPCA.call.func <- function(st_list,                      ## String: output path of qc procedure
                                 out_path                      ## String: output path of sdd procedure
){
  
  ## Load io code
  source(paste0(method_path, "/io.R"))
  ## set parameters
  check_file <- read.table(file = paste0(out_path, "/sdd_check_file.txt"))[, 1]
  gene_type <- check_file[3]
  if (gene_type == "custom"){
    
    custom_gene <- strsplit(check_file[4], ",")
    sparkversion <- NULL
  } else {
    
    custom_gene <- NULL
    sparkversion <- check_file[5]
  }
  clus_method <- check_file[6]
  clus_num <- strsplit(check_file[7], ",")[[1]] %>% as.numeric()
  kernel <- check_file[8]
  core_num <- check_file[9] %>% as.numeric()
  sample_size <- length(st_list$count_list)
  platform <- st_list[["platform"]]
  
  ## run SpatialPCA
  if (sample_size == 1) {
    
    ## create SpatialPCA Object
    pca_obj <- CreateSpatialPCAObject(counts = st_list[["count_list"]][[1]],
                                      location = st_list[["coord_list"]][[1]],
                                      project = "SpatialPCA", 
                                      gene.type = gene_type,
                                      numCores_spark = core_num,
                                      sparkversion = sparkversion, 
                                      gene.number = 3000,
                                      customGenelist = custom_gene, 
                                      min.loctions = 0, 
                                      min.features = 0) 
    ## estimate spatial PCs
    n_cells <- nrow(st_list[["coord_list"]][[1]])
    bandwidthtype_use <- ifelse(n_cells <= 3000, "SJ", "Silverman")
    sparse_Kernel <- ifelse(n_cells <= 3000, FALSE, TRUE)
    pca_obj <- SpatialPCA_buildKernel(pca_obj, 
                                      kerneltype = kernel,
                                      bandwidthtype = bandwidthtype_use,
                                      bandwidth.set.by.user = NULL,
                                      sparseKernel = sparse_Kernel,
                                      sparseKernel_tol = 1e-20,
                                      sparseKernel_ncore = core_num)
    ## set extra parameters 
    if (n_cells <= 3000) {
      
      est_fast <- pca_fast <- FALSE
    } else {
      
      est_fast <- pca_fast <- TRUE
    }
    pca_obj <- SpatialPCA_EstimateLoading(pca_obj, 
                                          maxiter = 300,
                                          initial_tau = 1,
                                          fast = est_fast,
                                          SpatialPCnum = N_SPC)
    pca_obj <- SpatialPCA_SpatialPCs(pca_obj, 
                                     fast = pca_fast)
    ## sample information
    sample_info <- data.frame(sample = "Sample1",
                              cell = rownames(pca_obj@location))
  } else {
    
    ## wrapped multiple samples
    pca_obj <- SpatialPCA_Multiple_Sample(count_list = st_list[[1]],
                                          location_list = st_list[[2]],
                                          gene.type = gene_type,
                                          sparkversion = sparkversion,
                                          numCores_spark = core_num,
                                          gene.number = 3000, 
                                          customGenelist = custom_gene,
                                          min.loctions = 0, 
                                          min.features = 0,
                                          bandwidth_common = 0.1)
    
    # sample information
    sample_info <- lapply(seq_along(st_list[[2]]), function(a){
      data.frame(sample = a,
                 cell = rownames(st_list[[2]][[a]]))
    }) %>% Reduce("rbind", .)
  }
  
  ## find clusters 
  refine_shape <- ifelse(platform == "Visium", "hexagon", "square")
  cl_list <- list()
  for (n_clust in clus_num) {
    
    # Detect spatial domains
    if(clus_method == "louvain"){
      cluster_label <- louvain_clustering(clusternum = n_clust, 
                                          latent_dat = as.matrix(pca_obj@SpatialPCs),
                                          knearest = round(sqrt(dim(pca_obj@SpatialPCs)[2])))
    } else if(clus_method == "walktrap"){
      cluster_label <- walktrap_clustering(clusternum = n_clust, 
                                           latent_dat = as.matrix(pca_obj@SpatialPCs),
                                           knearest = round(sqrt(dim(pca_obj@SpatialPCs)[2])))
    }
    clusterlabel_refine  <- refine_cluster_10x(clusterlabels = cluster_label, 
                                               location = pca_obj@location, 
                                               shape = refine_shape)
    
    # rename last cluster
    gridx <- paste0("clust_", length(unique(clusterlabel_refine)))
    cl_list[[gridx]] <- clusterlabel_refine
    
    ## check max clusters 
    if (length(unique(clusterlabel_refine)) < n_clust) {
      
      message(paste0("Max number of clusters reached: ", 
                     length(unique(clusterlabel_refine)),"!"))
      
      break
    }
  }
  
  ## output 
  spatialPCA_result <- list()
  spatialPCA_result$pcs <- pca_obj@SpatialPCs %>% 
    as.data.frame()
  dimnames(spatialPCA_result$pcs) <- list(paste0("pc", 1: N_SPC),
                                        sample_info$cell)
  spatialPCA_result$normalized_expr <- pca_obj@normalized_expr %>% 
    as.matrix()
  spatialPCA_result$location <- pca_obj@location %>% 
    as.data.frame()
  spatialPCA_result$sample <- sample_info
  spatialPCA_result$clus <- cl_list[!duplicated(names(cl_list))]  # remove items with duplicated gridx names
  
  save(spatialPCA_result,
       file = paste0(out_path, "/sdd_result.RData"))
  write.table(c(paste0(out_path, "/sdd_result.RData"), 
                length(unique(clusterlabel_refine))), 
              file = paste0(out_path, "/sdd_call_file.txt"), 
              row.names = F, quote = F, col.names = F)
  
  return(0)
}

# Function 3: BASS.call.func
BASS.call.func <- function(st_list,                            ## String: output path of qc procedure
                           out_path                            ## String: output path of sdd procedure
){
  
  sample_size <- length(st_list$count_list)
  ## load settings
  check_file <- read.table(paste0(out_path, "/sdd_check_file.txt"))[, 1]
  Cs <- check_file[3] %>% 
    strsplit(",") %>% unlist %>% as.numeric
  Rs <- check_file[4] %>% 
    strsplit(",") %>% unlist %>% as.numeric
  init_method <- check_file[5] 
  beta_method <- check_file[6]
  gene_select <- check_file[7]
  
  ## run BASS on R*C grid
  cl_dom_list <- list()
  ct_pi_list <- list()
  for (Rx in Rs) {
    
    for (Cx in Cs) {
      
      gridx <- paste0("R",Rx, "_C", Cx)
      ## build and preprocess BASS object
      bass_obj <- createBASSObject(X = st_list[[1]], 
                                   xy = st_list[[2]], 
                                   C = Cx, 
                                   R = Rx, 
                                   init_method = init_method,
                                   beta_method = beta_method,
                                   nsample = SAMPLES, 
                                   burnin = BURNIN)
      bass_obj <- BASS.preprocess(bass_obj, 
                                  doLogNormalize = NORM,
                                  geneSelect = gene_select, 
                                  nSE = GENENUM, 
                                  doPCA = TRUE, 
                                  nPC = DIMS,
                                  scaleFeature = SCALE, 
                                  doBatchCorrect = BATCH)
      
      ## fit model
      bass_obj <- BASS.run(bass_obj)
      bass_obj <- BASS.postprocess(bass_obj)
      
      ## format output
      cell_type_label <- bass_obj@results$c   # cell type clusters
      domain_label <- bass_obj@results$z      # spatial domain labels
      #
      c_d_list <- plyr::alply(seq_along(cell_type_label), 1, function(a){
        tmp <- data.frame(sample = a,
                          cell = rownames(bass_obj@xy[[a]]), 
                          cluster_label = cell_type_label[[a]], 
                          domain = domain_label[[a]])
        return(tmp)
      })
      
      cl_dom_list[[gridx]] <- Reduce("rbind", c_d_list)
      
      # cell type composition matrix
      cell_type_pi_mtx <- bass_obj@results$pi
      dimnames(cell_type_pi_mtx) <- list(paste0("ct", 1:nrow(cell_type_pi_mtx)),
                                         paste0("d", 1:ncol(cell_type_pi_mtx)))
      ct_pi_list[[gridx]] <- cell_type_pi_mtx
    }
  }
  
  ## output
  cl_dom_file <- paste0(out_path, "/cl_jo_BASS_cl_dom.RData")
  save(cl_dom_list, file = cl_dom_file)
  ct_pi_file <- paste0(out_path, "/cl_jo_BASS_ct_pi.RData")
  save(ct_pi_list, file = ct_pi_file)
  
  # pc matrix with batch corrected
  bass_pc_df <- bass_obj@X_run %>% as.matrix()
  dimnames(bass_pc_df) <- list(paste0("pc", 1:DIMS),
                               cl_dom_list[[gridx]]$cell)
  bass_pc_file <- paste0(out_path, "/cl_jo_BASS_pc.RData")
  save(bass_pc_df, file = bass_pc_file)
  write.table(c(cl_dom_file, bass_pc_file, ct_pi_file), 
              file = paste0(out_path, "/sdd_call_file.txt"), 
              row.names = F, quote = F, col.names = F)
  return(0)
}

# Function 4: ct.call
sdd.call <- function(data_path,                                ## String: output path of qc procedure
                     out_path                                  ## String: output path of sdd procedure
){
  
  ## Load io code
  source(paste0(method_path, "/io.R"))
  
  ## load st data
  call_file <- paste0(data_path, "/qc_call_file.txt")
  qc_param <- read.table(call_file)[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  st_list <- h5data.load(spatial_data_filename,      
                         sample_size,      
                         load_count = TRUE,    
                         normalization = FALSE, 
                         load_coord = TRUE) 
  
  ## load ct check file
  check_file <- paste0(out_path, "/sdd_check_file.txt")
  sdd_param <- read.table(check_file)[, 1]
  sdd_methods <- sdd_param[1]
  
  ## choose different ct methods
  if (sdd_methods == "SDD_sPCA-SpatialPCA"){
    
    sdd_result <- SpatialPCA.call.func(st_list, out_path)
  }
  if (sdd_methods == "CL_jo-BASS"){
    
    sdd_result <- BASS.call.func(st_list, out_path)
  }
  
  return(0)
}

# Function 5: SpatialPCA.post.func
SpatialPCA.post.func <- function(data_path = NULL,             ## String: output path of annot procedure
                                 out_path                      ## String: output path of sdd procedure
){
  
  ## load data
  load(paste0(out_path, "/sdd_result.RData"))
  
  ## cluster label
  cl_list <- spatialPCA_result$clus
  domain_label <- lapply(seq_along(cl_list), function(a){
    clust_df <- data.frame(cl_list[[a]])
    colnames(clust_df) <- names(cl_list)[a]
    return(clust_df)
  }) %>% Reduce("cbind", .) %>%
    cbind.data.frame(spatialPCA_result$sample, .)
  
  ## save spatial PCs
  spc_df <- spatialPCA_result$pcs
  save(spc_df, file = paste0(out_path, "/sdd_pc.RData"))
  
  ## output
  result_dir <- paste0(out_path, "/sdd_result")
  if (!file.exists(result_dir)){
    system(paste0("mkdir ", result_dir))
  }
  ### get annotation
  check_file <- read.table(paste0(out_path, "/sdd_check_file.txt"))[, 1]
  annot_method <- check_file[8]
  if(annot_method %in% c("scSorter", "Garnett")){
    
    annot_file <- bigreadr::fread2(paste0(data_path, "/annot_result/celltype_", 
                                          annot_method, ".txt"), 
                                   header = F)
    
    annot_out <- plyr::aaply(c(3: ncol(domain_label)), 1, function(s){
      
      ss <- domain_label[, s]
      freq_tab <- as.matrix(table(ss, annot_file[, 2]))
      prop_tab <- apply(freq_tab, 1, function(x) x/sum(x))
      write.csv(prop_tab, 
                file = paste0(result_dir, "/", annot_method, "_annot_", 
                              colnames(cluster_label)[s], ".csv"), 
                quote = F)
      return(s)
    })
  }
  write.table(domain_label, 
              file = paste0(result_dir, "/sdd_spca_SpatialPCA_domain_label.txt"), 
              sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(c("0", 
                paste0(result_dir, "/sdd_spca_SpatialPCA_domain_label.txt"),
                paste0(out_path, "/sdd_pc.RData")), 
              file = paste0(out_path, "/sdd_post_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}

# Function 6: BASS.post.func
BASS.post.func <- function(data_path,                          ## String: output path of annot procedure
                           out_path                            ## String: output path of sdd procedure
){
  
  ## load
  call_file <- read.table(paste0(out_path, "/sdd_call_file.txt"))[, 1]
  cl_dom_file <- call_file[1]
  bass_pc_file <- call_file[2]
  load(cl_dom_file)
  ## format output
  ### cluster label
  cluster_label <- lapply(seq_along(cl_dom_list), function(a){
    clust_df <- cl_dom_list[[a]][, "cluster_label", drop = F]
    colnames(clust_df) <- names(cl_dom_list)[a]
    return(clust_df)
  }) %>% Reduce("cbind", .)
  cluster_label <- cbind(cl_dom_list[[1]][,c("sample", "cell")],
                         cluster_label)
  ### domain label
  domain_label <- lapply(seq_along(cl_dom_list), function(a){
    domain_df <- cl_dom_list[[a]][, "domain", drop = F]
    colnames(domain_df) <- names(cl_dom_list)[a]
    return(domain_df)
  }) %>% Reduce("cbind", .)
  domain_label <- cbind(cl_dom_list[[1]][,c("sample", "cell")],
                        domain_label)
  
  ## output
  result_dir <- paste0(out_path, "/sdd_result")
  if (!file.exists(result_dir)){
    system(paste0("mkdir ", result_dir))
  }
  ## get annotation
  check_file <- read.table(paste0(out_path, "/sdd_check_file.txt"))[, 1]
  annot_method <- check_file[6]
  if(annot_method %in% c("scSorter", "Garnett")){
    
    annot_file <- bigreadr::fread2(paste0(data_path, "/annot_result/celltype_", 
                                          annot_method, ".txt"), 
                                   header = F)
    
    annot_out <- plyr::aaply(c(3: ncol(cluster_label)), 1, function(s){
      
      ss <- cluster_label[, s]
      freq_tab <- as.matrix(table(ss, annot_file[, 2]))
      prop_tab <- apply(freq_tab, 1, function(x) x/sum(x))
      write.csv(prop_tab, 
                file = paste0(result_dir, "/cl_jo_BASS_", annot_method, "_annot_", 
                              colnames(cluster_label)[s], ".csv"), 
                quote = F)
      return(s)
    })
  }
  cluster_label_file <- paste0(result_dir, "/cl_jo_BASS_cluster_label.txt")
  domain_label_file <- paste0(result_dir, "/cl_jo_BASS_domain_label.txt")
  write.table(cluster_label, file = cluster_label_file, 
              sep = "\t", row.names = F, quote = F)
  write.table(domain_label, file = domain_label_file, 
              sep = "\t", row.names = F, quote = F)
  write.table(c(cluster_label_file, domain_label_file, bass_pc_file),
              file = paste0(out_path, "/sdd_post_file.txt"),
              row.names = F, quote = F, col.names = F)
  return(0)
}

# Function 7: ct.post.func
sdd.post <- function(data_path = NULL,                          ## String: output path of annot procedure
                     out_path                                   ## String: output path of annot procedure
){
  
  ## load ct check file
  sdd_param <- read.table(paste0(out_path, "/sdd_check_file.txt"))[, 1]
  sdd_methods <- sdd_param[1]
  ## choose different ct methods
  if (sdd_methods == "SDD_sPCA-SpatialPCA"){
    
    sdd_result <- SpatialPCA.post.func(data_path, out_path)
  }
  if (sdd_methods == "CL_jo-BASS"){
    
    sdd_result <- BASS.post.func(data_path, out_path)
  }
  
  return(0)
}

# ###################
# ### test code
# data_path <- "/net/mulan/disk2/yasheng/stwebProject/03_result/02_Visium/V1_Mouse_Brain_Sagittal_Posterior"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/03_result/02_Visium/V1_Mouse_Brain_Sagittal_Posterior"
# sdd.check(data_path = data_path,
#           sdd_submodule = "SDD_sPCA",
#            methods = "SpatialPCA",
#            spatialPCA_geneType = "spatial",
#            SpatialPCA_customGene = "NULL",
#            SpatialPCA_sparkVer = "sparkx",
#            SpatialPCA_clusMethod = "louvain",
#            SpatialPCA_domainNum = "6,8",
#            SpatialPCA_kernel = "gaussian",
#            SpatialPCA_coreNum = 5,
#            annot = "Garnett",
#            out_path = output_path)
# sdd.call(data_path = data_path, out_path = output_path)
# sdd.post(out_path = output_path)
