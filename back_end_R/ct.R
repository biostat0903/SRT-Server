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

# Fix parameters
DIMS = 20              ## number of dimension for PCA
GENENUM = 2000         ## number of gene for PCA 
BURNIN = 2        ## number of burn-in
SAMPLES = 10        ## number of MCMC iteration 

NORM = TRUE
SCALE = TRUE
BATCH = FALSE

# Function 1: ct.check
ct.check <- function(data_path,                         ## String: output path of qc procedure
                     ct_submodule = NULL,               ## String: cell typing sub-modules: {CT_PCA}, {CL_jo}
                     methods = NULL,                    ## String: cell typing methods: {Seurat}, {BASS}
                     Seurat_resolutions = NULL,         ## String: resolution
                     Seurat_integration = NULL,         ## String: integration method: {integration} and {SCTransform}
                     Seurat_dims = NULL,                ## String: PCA number (re)
                     Seurat_anchors = NULL,             ## String: number of neighbors to use when picking anchors
                     BASS_cellNum = NULL,               ## String: number of cell types. 
                     BASS_domainNum = NULL,             ## String: number of spatial domains.
                     BASS_betaMethod = "SW",            ## String: fix the cell-cell interaction parameter. 
                                                        ##         to be beta {fix} or estimate the parameter based on the 
                                                        ##         data at hand {SW}.
                     BASS_initMethod = "mclust",        ## String: Initialize the cell type clusters and spatial domains 
                                                        ##         Two selections: {kmeans}, {mclust}.
                     BASS_geneSelect = "sparkx",        ## String: feature selection method: spatial gene {sparkx} or
                                                        ##         high variable gene {hvgs}.
                     annot_method = "None",             ## String: annotation methods: {None}, {scSorter}, {Garnett}
                     out_path                           ## String: output path for dr_cl procedure
){
  
  ## Set parameters and check QC module
  ct_methods <- paste(ct_submodule, methods, sep="-")
  qc_param <- read.table(paste0(data_path, "/qc_call_file.txt"))[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  if(!file.exists(spatial_data_filename)){
    
    stop("Please run the QC module!")
  }

  ## Output the first method: CT_PCA-Seurat
  if (ct_methods == "CT_PCA-Seurat"){
    
    if(is.null(Seurat_integration)){
      
      Seurat_dims <- NULL
      Seurat_anchors <- NULL
      message("We regard st data as one sample setting!")
    } else {
      
      message("We regard st data as multiple samples setting!")
    }
    if (sample_size != 1){
      
      write.table(c(ct_methods, annot_method, 
                    Seurat_resolutions, Seurat_integration, Seurat_dims, Seurat_anchors), 
                  file = paste0(out_path, "/ct_check_file.txt"), 
                  row.names = F, col.names = F, quote = F)
    } else {
      
      write.table(c(ct_methods, annot_method, 
                    Seurat_resolutions), 
                  file = paste0(out_path, "/ct_check_file.txt"), 
                  row.names = F, col.names = F, quote = F)
    }
  }
    
  ## Output the second method: CL_jo-BASS
  if (ct_methods == "CL_jo-BASS"){

    ## output
    write.table(c(ct_methods, annot_method, 
                  BASS_cellNum, BASS_domainNum,
                  BASS_initMethod, BASS_betaMethod, BASS_geneSelect),
                file = paste0(out_path, "/ct_check_file.txt"), 
                row.names = F, col.names = F, quote = F)
  }
  
  return(0)
}

# Function 2: single
single.seurat <- function(seurat_obj                    ## seurat object with one sample
){
  
  ## scale data
  seurat_obj <- NormalizeData(seurat_obj, 
                              normalization.method = "LogNormalize", 
                              scalBe.factor = 10000,
                              verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, 
                                     selection.method = "vst", 
                                     nfeatures = 2000, 
                                     verbose = FALSE)
  featureID <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, 
                          features = featureID, 
                          verbose = FALSE)
  
  ## output
  return(seurat_obj)
}

# Function 3: integration
integration.seurat <- function(seurat_obj_list,         ## seurat object with one sample
                               anchors,                 ## number of neighbors to use when picking anchors
                               dims                     ## dimdension selection
){
  
  # select features that are repeatedly variable across datasets for integration
  featureID <- SelectIntegrationFeatures(object.list = seurat_obj_list,
                                         nfeatures = 2000,
                                         fvf.nfeatures = 2000)
  seurat_obj_list <- lapply(X = seurat_obj_list, FUN = function(x) {
    x <- ScaleData(x, features = featureID, 
                   verbose = FALSE)
    x <- RunPCA(x, features = featureID, verbose = FALSE)
  })
  anchorID <- FindIntegrationAnchors(object.list = seurat_obj_list, 
                                     reduction = "rpca",
                                     k.anchor = anchors,
                                     k.filter = 200,
                                     k.score = 30,
                                     anchor.features = featureID)
  seurat_combined <- IntegrateData(anchorset = anchorID,
                                   dims = 1: dims,
                                   new.assay.name = "integrated",
                                   normalization.method = "LogNormalize")
  DefaultAssay(seurat_combined) <- "integrated"
  
  ## scale data
  seurat_combined <- ScaleData(seurat_combined, 
                               verbose = FALSE)
  
  ## output
  return(seurat_combined)
}

# Function 4: sctransform
sctransform.seurat <- function(seurat_obj_list,         ## seurat object with one sample
                               anchors,                 ## number of neighbors to use when picking anchors
                               dims                     ## dimension selection
){
  
  seurat_list <- lapply(seurat_obj_list, function(x){
    SCTransform(x, variable.features.n = 2000, 
                method = "glmGamPoi")
  })
  featureID <- SelectIntegrationFeatures(object.list = seurat_list,
                                         nfeatures = 2000)
  seurat_list <- PrepSCTIntegration(object.list = seurat_list,
                                    anchor.features = featureID)
  anchors_list <- FindIntegrationAnchors(object.list = seurat_list,
                                         normalization.method = "SCT",
                                         k.anchor = anchors,
                                         k.filter = 200,
                                         k.score = 30,
                                         anchor.features = featureID)
  seurat_sct <- IntegrateData(anchorset = anchors_list, 
                              normalization.method = "SCT",
                              dims = 1: dims)
  
  return(seurat_sct)
}

# Function 5: Seurat.call.func
Seurat.call.func <- function(st_list,                   ## String: output path of qc procedure
                             out_path                   ## String: output path for ct procedure
){

  sample_size <- length(st_list$count_list)
  ## Transform to Seurat format
  if (sample_size == 1){
    
    seurat_obj <- CreateSeuratObject(st_list$count_list[[1]])
  } else {
    
    seurat_obj <- plyr::llply(st_list$count_list, function(s){
      seurat_obj_s <- CreateSeuratObject(s)
      return(seurat_obj_s)
    })
  }
  
  ## Set parameters for different sample settings
  check_file <- read.table(paste0(out_path, "/ct_check_file.txt"))[, 1]
  ### one sample
  resolutions <- check_file[3] %>% 
    strsplit(",") %>% unlist %>% as.numeric
  
  ### two and more samples
  if(sample_size != 1){
    
    integration <- check_file[4]
    dims <- check_file[5] %>% 
      strsplit(",") %>% unlist %>% as.numeric
    anchors <- check_file[6] %>% 
      strsplit(",") %>% unlist %>% as.numeric
  }
  
  pc_list <- list()
  umap_list <- list()
  ## one sample
  if(sample_size == 1){
    
    ## Dimensional reduction
    ### PCA
    seurat_obj <- single.seurat(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj,
                         npcs = 50,
                         verbose = FALSE)
    eig_val <- (seurat_obj@reductions$pca@stdev)^2
    var_explained <- eig_val / sum(eig_val)
    pc_clust <- min(which(cumsum(var_explained)>0.8))
    cat(paste0("We select ", pc_clust, " PCs to explain 80% variance.\n"))
    ### UMAP
    seurat_obj <- RunUMAP(seurat_obj,
                          reduction = "pca",
                          dims = 1:pc_clust, 
                          verbose = FALSE)
    ### Clustering
    seurat_obj <- FindNeighbors(seurat_obj,
                                reduction = "pca",
                                dims = 1: pc_clust, 
                                verbose = FALSE)
    clust <- plyr::alply(resolutions, 1, function(s){
      
      clust_s <- FindClusters(seurat_obj, 
                              resolution = s,
                              verbose = FALSE) %>% Idents
      message(paste0("When the res = ", s, ", the spots cluster into ", 
                     length(unique(clust_s)), " clusters."))
      return(clust_s)
    }) %>% Reduce("cbind", .)
    clust_name <- paste0("cluster_res", resolutions)
    
    ## Output
    ### Extract PCs and UMAPs used
    pc_list[["cluster"]] <- Embeddings(seurat_obj[["pca"]])[, 1: pc_clust] %>% t()
    dimnames(pc_list[["cluster"]]) <- list(paste0("pc", 1: pc_clust),
                                           colnames(seurat_obj))
    umap_list[["cluster"]]  <- Embeddings(seurat_obj[["umap"]]) %>% as.data.frame()
    
    # Save seurat object
    SaveH5Seurat(seurat_obj, 
                 filename = paste0(out_path, "/ct_pca_seurat_obj.h5Seurat"), 
                 overwrite = TRUE)
    
  } else {
    
    clust_list <- list()
    clust_name <- vector()
    count <- 1
    for (anchor in anchors){
      
      for (dim in dims){
        
        gridx <- paste0("cluster_anchor", anchor, "_dim", dim)
        ### integration method
        if (integration == "Integration"){
          
          seurat_obj_norm <- lapply(X = seurat_obj, FUN = function(x) {
            x <- NormalizeData(x, verbose = FALSE)
            x <- FindVariableFeatures(x, 
                                      selection.method = "vst", 
                                      nfeatures = 2000, 
                                      verbose = FALSE)
          })
          seurat_obj <- integration.seurat(seurat_obj_norm, anchor, dim)
        }
        ### SCTransform method
        if (integration == "SCTransform"){
          
          seurat_obj <- sctransform.seurat(seurat_obj, anchor, dim)
        }
        
        ## set default assay to 'integrated'
        DefaultAssay(seurat_obj) <- "integrated"
        
        ### PCA
        seurat_obj <- RunPCA(seurat_obj,
                             npcs = 50,
                             verbose = FALSE)
        eig_val <- (seurat_obj@reductions$pca@stdev)^2
        var_explained <- eig_val / sum(eig_val)
        pc_clust <- min(which(cumsum(var_explained)>0.8))
        
        ## UMAP
        seurat_obj <- RunUMAP(seurat_obj,
                              reduction = "pca",
                              dims = 1:pc_clust, 
                              verbose = FALSE)
        
        ### Clustering
        seurat_obj <- FindNeighbors(seurat_obj,
                                    reduction = "pca",
                                    dims = 1: pc_clust)
        
        # 
        for (rr in resolutions){
          
          cluster <- FindClusters(seurat_obj, resolution = rr) %>% Idents
          clust_list[[count]] <- cluster
          count <- count + 1
          clust_name <- c(clust_name, paste0(gridx, "_res", rr))
        }
        
        ## Output
        ### Extract UMAP used
        umap_list[[gridx]] <- Embeddings(seurat_obj[["umap"]]) %>% as.data.frame()
        pc_list[[gridx]] <- Embeddings(seurat_obj[["pca"]])[, 1:pc_clust] %>% t()
        dimnames(pc_list[[gridx]]) <- list(paste0("pc", 1:pc_clust),
                                           colnames(seurat_obj))
        
        ## Save seurat object
        SaveH5Seurat(seurat_obj, 
                     filename = paste0(out_path, "/ct_pca_seurat_obj_anchor", anchor, "_dim", dim, ".h5Seurat"), 
                     overwrite = TRUE)
      }
    }
    clust <- Reduce("cbind", clust_list)
    
  }
  
  ### Output
  clust_outpath <- paste0(out_path, "/ct_pca_seurat_cluster.txt")
  umap_outpath <- paste0(out_path, "/ct_pca_seurat_umap_list.RData")
  pc_outpath <- paste0(out_path, "/ct_pca_seurat_pc_list.RData")
  write.table(clust, 
              file = clust_outpath,
              quote = F, row.names = F, col.names = F)
  save(umap_list, file = umap_outpath)
  save(pc_list, file = pc_outpath)
  write.table(c(clust_outpath, umap_outpath, pc_outpath, clust_name),
              file = paste0(out_path, "/ct_call_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}

# Function 6: BASS.call.func
BASS.call.func <- function(st_list,                     ## String: output path of qc procedure
                           out_path                     ## String: output path of ct procedure
){
  
  sample_size <- length(st_list$count_list)
  ## load settings
  check_file <- read.table(paste0(out_path, "/ct_check_file.txt"))[, 1]
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
              file = paste0(out_path, "/ct_call_file.txt"), 
              row.names = F, quote = F, col.names = F)
  return(0)
}

# Function 7: ct.call
ct.call <- function(data_path,                          ## String: output path of qc procedure
                    out_path                            ## String: output path of ct procedure
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
  check_file <- paste0(out_path, "/ct_check_file.txt")
  ct_param <- read.table(check_file)[, 1]
  ct_methods <- ct_param[1]
  
  ## choose different ct methods
  if (ct_methods == "CT_PCA-Seurat"){
    
    ct_result <- Seurat.call.func(st_list, out_path)
  }
  if (ct_methods == "CL_jo-BASS"){
    
    ct_result <- BASS.call.func(st_list, out_path)
  }

  return(0)
}

# Function 8: Seurat.post.func
Seurat.post.func <- function(data_path1,                ## String: output path of qc procedure
                             data_path2,                ## String: output path of annot procedure
                             out_path                   ## String: output path of ct procedure
){
  
  ## get cluster
  ### get sample and cell ID
  qc_file <- paste0(data_path1, "/qc_call_file.txt")
  qc_param <- read.table(qc_file)[, 1]
  spatial_data_filename <- qc_param[1]
  
  spatial_data_h5 <- H5File$new(spatial_data_filename, mode = "r")
  sample_size <- qc_param[2] %>% as.numeric
  cellID_list <- plyr::alply(c(1: sample_size), 1, function(a){
    cellID_grp_name <- paste0("cellID/cellID.s", a)
    cellID <- cbind(a, spatial_data_h5[[cellID_grp_name]][])
    return(cellID)
  })
  
  ### get cluster info
  call_file <- read.table(paste0(out_path, "/ct_call_file.txt"))[, 1]
  cluster_label <- read.table(call_file[1])
  cluster_label <- cbind.data.frame(Reduce("rbind", cellID_list),
                                    cluster_label)
  colnames(cluster_label) <- c("sample", "cell", call_file[-c(1:3)])
  umap_outpath <- call_file[2]
  pc_outpath <- call_file[3]
  
  ## output
  result_dir <- paste0(out_path, "/ct_result")
  if (!file.exists(result_dir)){

    system(paste0("mkdir ", result_dir))
  }
  ### get annotation
  check_file <- read.table(paste0(out_path, "/ct_check_file.txt"))[, 1]
  annot_method <- check_file[3]
  if(annot_method %in% c("scSorter", "Garnett")){
    
    annot_file <- bigreadr::fread2(paste0(data_path2, "/annot_result/celltype_", 
                                          annot_method, ".txt"), 
                                   header = F)
    annot_out <- plyr::aaply(c(3: ncol(cluster_label)), 1, function(s){
      
      ss <- cluster_label[, s]
      freq_tab <- as.matrix(table(ss, annot_file[, 2]))
      prop_tab <- apply(freq_tab, 1, function(x) x/sum(x))
      write.csv(prop_tab, file = paste0(result_dir, "/ct_pca_seurat_", annot_method, "_annot_", 
                                        colnames(cluster_label)[s], ".csv"), 
                quote = F)
      return(s)
    })
  }
  cluster_label_file <- paste0(result_dir, "/ct_pca_seurat_cluster_label.txt")
  write.table(cluster_label, file = cluster_label_file, 
              sep = "\t", row.names = F, quote = F)
  write.table(c(cluster_label_file, umap_outpath, pc_outpath), 
              file = paste0(out_path, "/ct_post_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}

# Function 9: Seurat.post.func
BASS.post.func <- function(data_path,                   ## String: output path of annot procedure
                           out_path                     ## String: output path of ct procedure
){
  
  ## load
  call_file <- read.table(paste0(out_path, "/ct_call_file.txt"))[, 1]
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
  #
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
  result_dir <- paste0(out_path, "/ct_result")
  if (!file.exists(result_dir)){
    system(paste0("mkdir ", result_dir))
  }
  ## get annotation
  check_file <- read.table(paste0(out_path, "/ct_check_file.txt"))[, 1]
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
              file = paste0(out_path, "/ct_post_file.txt"),
              row.names = F, quote = F, col.names = F)
  return(0)
}

# Function 10: ct.post
ct.post <- function(data_path1,                         ## String: output path of qc procedure
                    data_path2 = NULL,                  ## String: output path of annot procedure
                    out_path     
){
  
  ## load ct check file
  ct_param <- read.table(paste0(out_path, "/ct_check_file.txt"))[, 1]
  ct_methods <- ct_param[1]
  
  ## choose different ct methods
  if (ct_methods == "CT_PCA-Seurat"){
    
    ct_result <- Seurat.post.func(data_path1, data_path2, out_path)
  }
  if (ct_methods == "CL_jo-BASS"){
    
    ct_result <- BASS.post.func(data_path2, out_path)
  }
  
  return(0)
}

# ###################
# ### test code
# data_path1 <- "/net/mulan/disk2/yasheng/stwebProject/03_result/02_Visium/V1_Mouse_Brain_Sagittal_Posterior"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/03_result/02_Visium/V1_Mouse_Brain_Sagittal_Posterior"
# ct.check(data_path = data_path1,
#          ct_submodule = "CT_PCA",
#          methods = "Seurat",
#          Seurat_resolutions = "0.2,0.5",
#          out_path = output_path)
# ct.call(data_path = data_path1, out_path = output_path)
# ct.post(data_path1 = data_path1, data_path2 = NULL, out_path = output_path)
# 
# ct.check(data_path = data_path1,
#          ct_submodule = "CL_jo",
#          methods = "BASS",
#          BASS_cellNum = "6,10",
#          BASS_domainNum = "3,5",
#          out_path = output_path)
# ct.call(data_path = data_path1, out_path = output_path)
# ct.post(data_path1 = data_path1, data_path2 = NULL, out_path = output_path)
