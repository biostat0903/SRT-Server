#! /usr/bin/env Rscript
# performing dimensional reduction and clustering for spatial transcriptomics data
# up-stream code: qc.R
# down-stream plot: dr_cl_plt.R
# down-stream code: de.R, ccc.R

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
library(ggplot2)

# Function 1: cl.check
cl.check <- function(data_path,                      ## String: output path of qc procedure
                     resolutions = NULL,             ## String: resolution
                     integration = NULL,             ## String: integration method: {integration} and {SCTransform}
                     dims = NULL,                    ## String: PCA number (re)
                     anchors = NULL,                 ## String: number of neighbors to use when picking anchors
                     annot_method = "None",          ## String: annotation methods: {None}, {scSorter}, {Garnett}
                     out_path                        ## String: output path for dr_cl procedure
){
  
  ## check input
  if(is.null(resolutions)){
    
    cat("We regard st data as multiple samples setting!\n")
    if(is.null(integration) | is.null(dims) | is.null (anchors)){
      
      stop("The setting of \"integration\", \"integration\" and \"integration\" should not be NULL!")
    }
  } else {
    
    cat("We regard st data as one sample setting!\n")
  }
  call_file <- read.table(paste0(data_path, "/qc_call_file.txt"))[, 1]
  sample_size <- call_file[2]
  
  ## output 
  if (sample_size != 1){
    
    write.table(c(sample_size, resolutions, annot_method, integration, dims, anchors), 
                file = paste0(out_path, "/dr_cl_check_file.txt"), 
                row.names = F, col.names = F, quote = F)
  } else {
    
    write.table(c(sample_size, resolutions, annot_method), 
                file = paste0(out_path, "/dr_cl_check_file.txt"), 
                row.names = F, col.names = F, quote = F)
  }

  return(0)
}

# Function 2: single
cl.single.seurat <- function(seurat_obj              ## seurat object with one sample
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
cl.integration.seurat <- function(seurat_obj_list,   ## seurat object with one sample
                                  anchors,           ## number of neighbors to use when picking anchors
                                  dims               ## dimdension selection
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
cl.sctransform.seurat <- function(seurat_obj_list,   ## seurat object with one sample
                                  anchors,           ## number of neighbors to use when picking anchors
                                  dims               ## dimension selection
){
  
  seurat_list <- lapply(seurat_obj_list, function(x){
    SCTransform(x, variable.features.n = 2000, 
                method = "glmGamPoi")
  })
  featureID <- SelectIntegrationFeatures(object.list = seurat_list,
                                         nfeatures = 2000)
  seurat_list <- PrepSCTIntegration(object.list = seurat_list,
                                    anchor.features = featureID)
  anchorID <- FindIntegrationAnchors(object.list = seurat_list,
                                     normalization.method = "SCT",
                                     k.anchor = anchors,
                                     k.filter = 200,
                                     k.score = 30,
                                     anchor.features = featureID)
  seurat_sct <- IntegrateData(anchorset = anchorID, 
                              normalization.method = "SCT",
                              dims = 1: dims)
  
  ## output
  return(seurat_sct)
}

# Function 5: cl.call
cl.call <- function(data_path,                       ## String: output path of qc procedure
                    out_path                         ## String: output path for dr_cl procedure
){
  
  ## Load io code
  source(paste0(method_path, "/io.R"))
  
  ## Load st data
  check_file <- paste0(data_path, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  st_list <- h5data.load(spatial_data_filename,      
                         sample_size,      
                         load_count = TRUE,    
                         normalization = FALSE, 
                         load_coord = FALSE) 
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
  check_file <- read.table(paste0(out_path, "/dr_cl_check_file.txt"))[, 1]
  ### one sample
  sample_size <- as.numeric(check_file[1])
  resolutions <- check_file[2] %>% 
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
    seurat_obj <- cl.single.seurat(seurat_obj)
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
                 filename = paste0(out_path, "/dr_cl_obj.h5Seurat"), 
                 overwrite = TRUE)
    
  } else {
    
    seurat_obj_norm <- lapply(X = seurat_obj, FUN = function(x) {
      x <- NormalizeData(x, verbose = FALSE)
      x <- FindVariableFeatures(x, 
                                selection.method = "vst", 
                                nfeatures = 2000, 
                                verbose = FALSE)
    })
    
    clust_list <- list()
    clust_name <- vector()
    count <- 1
    for (anchor in anchors){
      
      for (dim in dims){
        
        gridx <- paste0("cluster_anchor", anchor, "_dim", dim)
        ### integration method
        if (integration == "Integration"){
          seurat_obj <- cl.integration.seurat(seurat_obj_norm, 
                                              anchor, 
                                              dim)
          DefaultAssay(seurat_obj) <- "integrated"
        }
        ### SCTransform method
        if (integration == "SCTransform"){
          
          seurat_obj <- cl.sctransform.seurat(seurat_obj, 
                                              anchor, 
                                              dim)
          DefaultAssay(seurat_obj) <- "SCT"
        }
        
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
                     filename = paste0(out_path, "/dr_cl_obj_anchor", anchor, "_dim", dim, ".h5Seurat"), 
                     overwrite = TRUE)
      }
    }
    clust <- Reduce("cbind", clust_list)

  }
  
  ### Output
  clust_outpath <- paste0(out_path, "/dr_cl_cluster.txt")
  umap_outpath <- paste0(out_path, "/dr_cl_umap_list.RData")
  pc_outpath <- paste0(out_path, "/dr_cl_pc_list.RData")
  write.table(clust, 
              file = clust_outpath,
              quote = F, row.names = F, col.names = F)
  save(umap_list, file = umap_outpath)
  save(pc_list, file = pc_outpath)
  write.table(c(clust_outpath, umap_outpath, pc_outpath, clust_name),
              file = paste0(out_path, "/dr_cl_call_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}

# Function 6: cl.post
cl.post <- function(data_path1,                      ## String: output path of qc procedure
                    data_path2 = NULL,               ## String: output path of annot procedure
                    out_path                         ## String: output path of dr_cl procedure
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
  call_file <- read.table(paste0(out_path, "/dr_cl_call_file.txt"))[, 1]
  cluster_label <- read.table(call_file[1])
  cluster_label <- cbind.data.frame(Reduce("rbind", cellID_list),
                                    cluster_label)
  colnames(cluster_label) <- c("sample", "cell", call_file[-c(1:3)])
  umap_outpath <- call_file[2]
  pc_outpath <- call_file[3]
  

  ## output
  result_dir <- paste0(out_path, "/dr_cl_result")
  if (!file.exists(result_dir)){
    system(paste0("mkdir ", result_dir))
  }
  ### get annotation
  check_file <- read.table(paste0(out_path, "/dr_cl_check_file.txt"))[, 1]
  annot_method <- check_file[3]
  if(annot_method %in% c("scSorter", "Garnett")){
    
    annot_file <- bigreadr::fread2(paste0(data_path2, "/annot_result/celltype_", 
                                          annot_method, ".txt"), 
                                   header = F)
    
    annot_out <- plyr::aaply(c(3: ncol(cluster_label)), 1, function(s){
      
      ss <- cluster_label[, s]
      freq_tab <- as.matrix(table(ss, annot_file[, 2]))
      prop_tab <- apply(freq_tab, 1, function(x) x/sum(x))
      write.csv(prop_tab, file = paste0(result_dir, "/", annot_method, "_annot_", 
                                        colnames(cluster_label)[s], ".csv"), 
                quote = F)
      return(s)
    })
  }
  write.table(cluster_label, 
              file = paste0(result_dir, "/cluster_label.txt"), 
              sep = "\t", row.names = F, quote = F)
  write.table(c(paste0(result_dir, "/cluster_label.txt"),
                umap_outpath,
                pc_outpath), 
              file = paste0(out_path, "/dr_cl_post_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}


# ###################
# ### test code
# data_path1 <- "/net/mulan/disk2/yasheng/stwebProject/test/qc"
# data_path2 <- "/net/mulan/disk2/yasheng/stwebProject/test/annot"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/test/dr_cl"
# cl.check(data_path = data_path1,
#          resolutions = "0.2,0.5,0.7",
#          out_path = output_path,
#          annot = "Garnett")
# cl.call(data_path = data_path1, out_path = output_path)
# cl.post(data_path1 = data_path1, data_path2 = data_path2, out_path = output_path)
