#! /usr/bin/env Rscript
# wrap dimensional reduction and clustering for st data
# up-stream procedure: qc.R
# down-stream procedure: de.R, ccc.R

# Set method_path
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# Load packages
library(hdf5r)
library(dplyr)
library(tidyr)
library(ggplot2)
library(BASS)

# Fix parameters
DIMS = 20              ## number of dimension for PCA
GENENUM = 2000         ## number of gene for PCA 
BURNIN = 20          ## number of burn-in
SAMPLES = 100        ## number of MCMC iteration 

NORM = TRUE
SCALE = TRUE
BATCH = FALSE

# Function 1: check function
clwr.check <- function(data_path = NULL,                ## String: output path of qc procedure
                       C_list = NULL,                   ## String: number of cell types. 
                       R_list = NULL,                   ## String: number of spatial domains.
                       beta_method = "SW",              ## String: fix the cell-cell interaction parameter. 
                                                        ## to be beta {fix} or estimate the parameter based on the 
                                                        ## data at hand {SW}.
                       init_method = "mclust",          ## String: Initialize the cell type clusters and spatial 
                                                        ## domains with either {kmeans} or {mclust}.
                       gene_select = "sparkx",          ## String: feature selection method: spatial gene {sparkx} or
                                                        ## high variable gene {hvgs}.
                       annot_method = "None",           ## String: annotation methods: {None}, {scSorter}, {Garnett}
                       out_path = NULL                  ## String: path for BASS object
){
  
  ## check inputs
  if(is.null(beta_method) | is.null(init_method) | is.null(gene_select)){
    
    stop("Please specify the method for initialization, beta or feature selection!")
  } else {
    
    if (!init_method %in% c("kmeans", "mclust")) {
      
      stop("Please choose initialization method from kmeans or mclust!")
    } else {
      
      cat(paste0("We use ", init_method, " method to initialize the cell type clusters and spatial domains!\n"))
    }
    
    if (!beta_method %in% c("fix", "SW")) {
      
      stop("Please choose beta method from fix or SW!")
    } else {
      
      cat(paste0("We set ", beta_method, " for the cell-cell interaction parameter!\n"))
    }
    
    if (!gene_select %in% c("sparkx", "hvgs")) {
      
      stop("Please choose feature selection method from sparkx or hvgs!")
    } else {
      
      cat(paste0("We use ", gene_select, " for feature selection!\n"))
    }
    
  }
  
  ## output
  write.table(c(C_list, R_list,
                init_method, beta_method, gene_select, annot_method),
              file = paste0(out_path, "/dr_cl_wr_check_file.txt"), 
              row.names = F, col.names = F, quote = F)
  
  return(0)
}


# Function 2: BASS_call
clwr.call <- function(data_path,                        ## String: output path of qc procedure
                      out_path                          ## String: output path of dr_cl_wr procedure
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
  
  ## load settings
  check_file <- read.table(paste0(out_path, "/dr_cl_wr_check_file.txt"))[, 1]
  Cs <- check_file[1] %>% 
    strsplit(",") %>% unlist %>% as.numeric
  Rs <- check_file[2] %>% 
    strsplit(",") %>% unlist %>% as.numeric
  init_method <- check_file[3] 
  beta_method <- check_file[4]
  gene_select <- check_file[5]
  
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
  cl_dom_file <- paste0(out_path, "/dr_cl_wr_cl_dom.RData")
  save(cl_dom_list, file = cl_dom_file)
  ct_pi_file <- paste0(out_path, "/dr_cl_wr_ct_pi.RData")
  save(ct_pi_list, file = ct_pi_file)
  
  # pc matrix with batch corrected
  bass_pc_df <- bass_obj@X_run %>% as.matrix()
  dimnames(bass_pc_df) <- list(paste0("pc", 1:DIMS),
                                cl_dom_list[[gridx]]$cell)
  
  bass_pc_file <- paste0(out_path, "/dr_cl_wr_pc.RData")
  save(bass_pc_df, file = bass_pc_file)
  write.table(c(cl_dom_file, bass_pc_file, ct_pi_file), 
              file = paste0(out_path, "/dr_cl_wr_call_file.txt"), 
              row.names = F, quote = F, col.names = F)
  return(0)
}

# Function 3: post function
clwr.post <- function(data_path = NULL,                 ## String: output path of annot procedure                             
                      out_path                          ## String: output path of dr_cl_wr procedure 
){
  
  ## load
  call_file <- read.table(paste0(out_path, "/dr_cl_wr_call_file.txt"))[, 1]
  cl_dom_file <- call_file[1]
  bass_pc_file <- call_file[2]
  load(cl_dom_file)
  
  ## format output
  # cluster label
  cluster_label <- lapply(seq_along(cl_dom_list), function(a){
    clust_df <- cl_dom_list[[a]][, "cluster_label", drop = F]
    colnames(clust_df) <- names(cl_dom_list)[a]
    return(clust_df)
  }) %>% Reduce("cbind", .)
  #
  cluster_label <- cbind(cl_dom_list[[1]][,c("sample", "cell")],
                         cluster_label)
  
  # domain label
  domain_label <- lapply(seq_along(cl_dom_list), function(a){
    domain_df <- cl_dom_list[[a]][, "domain", drop = F]
    colnames(domain_df) <- names(cl_dom_list)[a]
    return(domain_df)
  }) %>% Reduce("cbind", .)
  #
  domain_label <- cbind(cl_dom_list[[1]][,c("sample", "cell")],
                        domain_label)
  
  ## output
  result_dir <- paste0(out_path, "/dr_cl_wr_result")
  if (!file.exists(result_dir)){
    system(paste0("mkdir ", result_dir))
  }
  ### get annotation
  check_file <- read.table(paste0(out_path, "/dr_cl_wr_check_file.txt"))[, 1]
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
                file = paste0(result_dir, "/", annot_method, "_annot_", 
                              colnames(cluster_label)[s], ".csv"), 
                quote = F)
      return(s)
    })
  }
  write.table(cluster_label, 
              file = paste0(result_dir, "/cluster_label.txt"), 
              sep = "\t", row.names = F, quote = F)
  write.table(domain_label, 
              file = paste0(result_dir, "/domain_label.txt"), 
              sep = "\t", row.names = F, quote = F)
  write.table(c(paste0(result_dir, "/cluster_label.txt"), 
                paste0(result_dir, "/domain_label.txt"),
                bass_pc_file),
              file = paste0(out_path, "/dr_cl_wr_post_file.txt"),
              row.names = F, quote = F, col.names = F)
  return(0)
}

# ###################
# ### test code
# data_path1 <- "/net/mulan/disk2/yasheng/stwebProject/test/qc"
# data_path2 <- "/net/mulan/disk2/yasheng/stwebProject/test/annot"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/test/dr_cl_wr"
# clwr.check(data_path = data_path1,
#            C_list = "10",
#            R_list = "4,5",
#            beta_method = "SW",
#            init_method = "mclust",
#            gene_select = "sparkx",
#            annot = "Garnett",
#            out_path = output_path)
# clwr.call(data_path = data_path1, out_path = output_path)
# clwr.post(data_path = data_path2, out_path = output_path)
