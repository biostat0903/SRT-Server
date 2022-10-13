#! /usr/bin/env Rscript
# Infer a low dimensional representation of the gene expression data in spatial transcriptomics

# Set method_path
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# load packages 
library(SpatialPCA)
library(hdf5r)
library(dplyr)
library(umap)
library(Rtsne)
library(ggplot2)
library(bigreadr)

# Fix parameters
N_SPC = 20 # number of pcs used in single sample (same as fixed in multiple samples)

# Function 1: check function
clsp.check <- function(data_path = NULL,              ## String: output path of qc procedure
                       gene_type = "spatial",         ## String: The type of genes to be used
                                                      ## Three selections: spatial, hvg, custom
                       custom_gene = "NULL",          ## String: path for user specified genes. 
                                                      ## required in gene_type=="custom".
                       spark_ver = "NULL",            ## String: selection spark, sparkx
                                                      ## required in gene_type=="spatial".
                       clus_method = "louvain",       ## String: clustering method
                                                      ## Two selections: walktrap, louvain
                       clus_num = "5,10,15",          ## String: desired cluster number
                       kernel = "gaussian",           ## String: kernel for SpatialPCA
                                                      ## Four selections: "gaussian", "cauchy", "quadratic", "delaunday"
                       core_num = 1,                  ## The core numbers for SpatialPCA
                       annot_method = "None",         ## String: annotation methods: {None}, {scSorter}, {Garnett}
                       out_path = "./"                ## String: output path for dr_cl_sp procedure
){
  
  ## check file
  check_file <- paste0(data_path, "/qc_call_file.txt")
  if(!file.exists(check_file)){
    
    stop("No QC file! Please run QC module!")
  } 
  
  ## check inputs
  if (gene_type == "spatial"){
    
    message("We use spatial variable genes to fit SpatialPCA.")
    if(spark_ver == "sparkx"){
      
      message("We use sparkx to identify svg.")
    } else {
      
      message("We use spark to identify svg.")
    }
  }
  if (gene_type == "custom"){
    
    if(is.null(custom_gene)){
      
      stop("We need the list for custom gene!")
    }
  }
  if (gene_type == "hvg"){
    
    message("We use Seurat to select the hvg.")
  }
  ## output
  write.table(c(gene_type, custom_gene, spark_ver,
                clus_method, clus_num, kernel, core_num, annot_method), 
              file = paste0(out_path, "/dr_cl_sp_check_file.txt"), 
              row.names = F, quote = F, col.names = F)
  return(0)
}

# Function 2: call function
clsp.call <- function(data_path = NULL,               ## String: output path of qc procedure
                      out_path                        ## String: output path of dr_cl_sp procedure
){
  
  ## Load io code
  source(paste0(method_path, "/io.R"))
  
  ## load st data
  check_file <- paste0(data_path, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  st_list <- h5data.load(spatial_data_filename,      
                         sample_size,      
                         load_count = TRUE,    
                         normalization = FALSE, 
                         load_coord = TRUE) 
  
  ## set parameters
  check_file <- read.table(file = paste0(out_path, "/dr_cl_sp_check_file.txt"))[, 1]
  gene_type <- check_file[1]
  if (gene_type == "custom"){

    custom_gene <- strsplit(check_file[2], ",")
    sparkversion <- NULL
  } else {
    
    custom_gene <- NULL
    sparkversion <- check_file[3]
  }
  clus_method <- check_file[4]
  clus_num <- strsplit(check_file[5], ",")[[1]] %>% as.numeric()
  kernel <- check_file[6]
  core_num <- check_file[7] %>% as.numeric()
  sample_size <- length(st_list[["count_list"]])
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
  dr_cl_sp_result <- list()
  dr_cl_sp_result$pcs <- pca_obj@SpatialPCs %>% 
    as.data.frame()
  dimnames(dr_cl_sp_result$pcs) <- list(paste0("pc", 1: N_SPC),
                                          sample_info$cell)
  dr_cl_sp_result$normalized_expr <- pca_obj@normalized_expr %>% 
    as.matrix()
  dr_cl_sp_result$location <- pca_obj@location %>% 
    as.data.frame()
  dr_cl_sp_result$sample <- sample_info
  dr_cl_sp_result$clus <- cl_list[!duplicated(names(cl_list))]  # remove items with duplicated gridx names
  
  save(dr_cl_sp_result,
       file = paste0(out_path, "/dr_cl_sp_result.RData"))

  return(0)
}

# Function 3: post function
clsp.post <- function(data_path = NULL,               ## String: output path of annot procedure
                      out_path                        ## String: output path of dr_cl_sp procedure
){
  
  ## load data
  load(paste0(out_path, "/dr_cl_sp_result.RData"))
 
  ## cluster label
  cl_list <- dr_cl_sp_result$clus
  cluster_label <- lapply(seq_along(cl_list), function(a){
    clust_df <- data.frame(cl_list[[a]])
    colnames(clust_df) <- names(cl_list)[a]
    return(clust_df)
  }) %>% Reduce("cbind", .) %>%
    cbind.data.frame(dr_cl_sp_result$sample, .)

  ## save spatial PCs
  spc_df <- dr_cl_sp_result$pcs
  save(spc_df, file = paste0(out_path, "/dr_cl_sp_pc.RData"))
  
  ## output
  result_dir <- paste0(out_path, "/dr_cl_sp_result")
  if (!file.exists(result_dir)){
    system(paste0("mkdir ", result_dir))
  }
  ### get annotation
  check_file <- read.table(paste0(out_path, "/dr_cl_sp_check_file.txt"))[, 1]
  annot_method <- check_file[8]
  if(annot_method %in% c("scSorter", "Garnett")){
    
    annot_file <- bigreadr::fread2(paste0(data_path, "/annot_result/celltype_", 
                                          annot_method, ".txt"), 
                                   header = F)
    
    annot_out <- plyr::aaply(c(3: ncol(cluster_label)), 1, function(s){
      
      ss <- cluster_label[, s]
      freq_tab <- as.matrix(table(ss, annot_file[, 2]))
      prop_tab <- apply(freq_tab, 1, function(x) x/sum(x))
      write.csv(prop_tab, 
                file = paste0(result_dir, "/", annot_method, "_annot_", 
                              colnames(cluster_label)[s], ".csv"), 
                quote = F)
      return(s)
    })
  }
  write.table(cluster_label, 
              file = paste0(result_dir, "/cluster_label.txt"), 
              sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(c(paste0(result_dir, "/cluster_label.txt"),
                paste0(out_path, "/dr_cl_sp_pc.RData")), 
              file = paste0(out_path, "/dr_cl_sp_post_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}


# ###################
# ### test code
# data_path1 <- "/net/mulan/disk2/yasheng/stwebProject/test/qc"
# data_path2 <- "/net/mulan/disk2/yasheng/stwebProject/test/annot"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/test/dr_cl_sp"
# clsp.check(data_path = data_path1,
#            gene_type = "spatial",
#            # custom_gene = NULL,
#            spark_ver = "sparkx",
#            clus_method = "louvain",
#            clus_num = "5,10",
#            kernel = "gaussian",
#            core_num = 10,
#            annot = "Garnett",
#            out_path = output_path)
# clsp.call(data_path = data_path1, out_path = output_path)
# clsp.post(data_path = data_path2, out_path = output_path)
