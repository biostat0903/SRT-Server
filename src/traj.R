#! /usr/bin/env Rscript
# Trajectory inference
# up-stream code: ct.R, sdd.R
# down-stream procedure: None

# Set method_path
method_path <- "/public/home/biostat03/project/stwebProject/01_code/srt_server/dev"

# Load packages
library(slingshot)
library(SingleCellExperiment)
library(tradeSeq)
library(BiocParallel)
library(Seurat)
library(harmony)
library(hdf5r)
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
library(bigreadr)

# Fix parameter
N_CORES = 2
USE_HVG = T
# Function 1: check function
traj.check <- function(data_path1 = NULL,               ## String: output path of clustering procedure
                       data_path2 = NULL,               ## String: output path of qc procedure
                       data_path3 = NULL,               ## String: output path of select_ct_file.txt
                       data_path4 = NULL,               ## String: output path of start_clus_file.txt
                       data_path5 = NULL,               ## String: output path of select_clus_file.txt
                       start_clus = NULL,
                       select_clus = NULL,
                       select_cell_type = NULL,
                       submodule = NULL,                ## String: cell typing sub-modules: {CT_PCA}, {SDD_sPCA}, {CL_jo}, {DECON}
                       methods = NULL,                  ## String: cell typing/decon methods: {Seurat}, {SpatialPCA}, {BASS}, {CARD}, {cell2location}, {tangram}
                       Seurat_resolutions = NULL,       ## String: resolution
                       Seurat_dims = NULL,              ## String: PCA number (re)
                       Seurat_anchors = NULL,           ## String: number of neighbors to use when picking anchors
                       BASS_cellNum = NULL,             ## String: number of cell types. 
                       BASS_domainNum = NULL,           ## String: number of spatial domains.
                       SpatialPCA_domainNum = NULL,     ## String: desired cluster number
                       out_path = NULL                  ## String: output path of traj procedure
){
  
  ##
  if (submodule == "DECON") {
    
    sel_ct_file <- NA
    decon_file <- paste0(data_path1, "/decon_post_file.txt")
    ## load inputs
    if (!file.exists(data_path3)) {
      
      stop("Please provide started clusters file!")
    } else {
      
      start_ct_file <- data_path3
    }
    
    if(!file.exists(decon_file)){
      
      stop("No DECON file! Please run DECON module!")
    } else {
      
      post_file <- read.table(decon_file)[,1]
      spatial_data_filename <- post_file[-c(1:2)] %>% paste(., collapse = "|")
      post_methods <- post_file[1]
      if (post_methods != methods) {
        stop(paste0("Please run DECON module in ", methods, " first!"))
      }
      sample_size <- length(spatial_data_filename)
      check_file <- c(spatial_data_filename, submodule, methods,
                      sample_size, NA, 
                      start_ct_file, sel_ct_file,
                      NA, NA,
                      NA, NA)
    }
    
  } else {
    
    call_file <- paste0(data_path2, "/qc_call_file.txt")
    if(!file.exists(call_file)){
      
      stop("No QC file! Please run QC module!")
    } else {
      
      qc_param <- read.table(call_file)[, 1]
      spatial_data_filename <- qc_param[1]
      sample_size <- qc_param[2]
    }
    
    ##
    if (!file.exists(data_path4) | !file.exists(data_path5)) {
      
      stop("Please provide selected clusters and strated clusters file!")
    } else {
      
      start_clus_file <- data_path4
      sel_clus_file <- data_path5
      start_clus <- read.table(start_clus_file, header = F, sep = "\t")[, 1, drop = T]
      sel_clus <- read.table(sel_clus_file, header = F, sep = "\t")[, 1, drop = T]
      # check lines
      if(length(start_clus) != length(sel_clus)){
        
        stop("Please ensure a consistant number of lines in the start clusters file and the selected clusters file!")
      } else {
        
        # check match
        check_match <- lapply(seq_along(start_clus), function(nn){
          start_clusx <- start_clus[nn] %>% as.character()
          sel_clusx <- sel_clus[nn] %>% 
            strsplit(., ",") %>% 
            unlist() %>% 
            as.character()
          return((any(toupper(sel_clusx) == "ALL") | start_clusx %in% sel_clusx))
        }) %>% unlist
        
        if (!all(check_match)) {
          
          stop("Please ensure each selected cluster group to containe the corresponding start cluster!")
        }
      }
    }
    
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
      if (methods %in% c("Garnett", "scSorter")){
        post_file <- read.table(paste0(data_path1, "/ct_post_file.txt"))[,1]
        grid_pc <- NA
        grid_use <- "cluster"
        umap_file <- NA
        pc_file <- NA
        # stop("slingshot should use the PCs!")
      }
      check_file <- c(spatial_data_filename, submodule, methods,
                      sample_size, grid_use, 
                      start_clus_file, sel_clus_file,
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
        check_file <- c(spatial_data_filename, submodule, methods,
                        sample_size, grid_use, 
                        start_clus_file, sel_clus_file,
                        post_file[1], NA,
                        post_file[3], NA)
      } else {
        
        post_file <- read.table(paste0(data_path1, "/sdd_post_file.txt"))[,1]
        check_file <- c(spatial_data_filename, submodule, methods,
                        sample_size, grid_use, 
                        start_clus_file, sel_clus_file,
                        NA, post_file[2],
                        post_file[3], NA)
      }
    }
    ## Set spatial domain detection modules
    if(grepl("SDD", submodule)){
      
      if(methods == "SpatialPCA"){
        
        grid_use <- paste0("clust_", SpatialPCA_domainNum)
      }
      post_file <- read.table(paste0(data_path1, "/sdd_post_file.txt"))[,1]
      check_file <- c(spatial_data_filename, submodule, methods,
                      sample_size, grid_use, 
                      start_clus_file, sel_clus_file,
                      NA, post_file[2],
                      post_file[3], NA)
    }  
  }
  ## output file
  write.table(check_file, file = paste0(out_path, "/traj_check_file.txt"), 
              row.names = F, quote = F, col.names = F)
  
  return(0)
}

# Function 2: pca using Seurat
pca.seurat <- function(counts = NULL,
                       anno_df = NULL,
                       sample_idx = "sample"){
  #
  rownames(anno_df) <- colnames(counts)
  seurat_obj <- CreateSeuratObject(counts = counts, 
                                   meta.data = anno_df, 
                                   min.cells = 10) %>%
    NormalizeData(., verbose = FALSE) %>%
    FindVariableFeatures() %>%
    ScaleData() %>% 
    RunPCA(verbose = FALSE)
  if (length(unique(anno_df[[sample_idx]])) == 1) {
    
    return(Embeddings(seurat_obj[["pca"]]))
  } else {
    
    seurat_obj <- RunHarmony(seurat_obj, group.by.vars = sample_idx)
    return(Embeddings(seurat_obj[["harmony"]]))
  }
}

# Function 3: traj test use slingshot
traj.func <- function(st_list = NULL,
                      cluster_df = NULL,
                      grid_use = NULL,
                      methods = NULL,
                      pc_file = NA,
                      umap_file = NA,
                      start_clust = NULL,
                      n_cores = N_CORES
){
  
  ## 3.1 create sce object from st_list
  count_merged <- purrr::reduce(st_list[["count_list"]], function(x, y) {
    cbind(x = x, y = y)
  })
  coord_merged <- lapply(seq_along(st_list[["coord_list"]]), function(samplex){
    coords <- st_list[["coord_list"]][[samplex]] %>% as.data.frame()
    coords$sample <- as.character(samplex)
    return(coords)
  }) %>% Reduce("rbind", .)
  rownames(coord_merged) <- colnames(count_merged)
  sce_obj <- SingleCellExperiment(assays = List(counts = count_merged))
  colData(sce_obj)$x <- coord_merged[,1]
  colData(sce_obj)$y <- coord_merged[,2]
  ## 3.2 obtain cluster information
  if (!all(colnames(sce_obj) %in% cluster_df$cell)) {
    
    stop("Number of cells in expression matrix and cluster label file do not match!")
  } else {
    
    cluster_df <- cluster_df[match(colnames(sce_obj),
                                   cluster_df$cell),]
    colData(sce_obj)$sample <- cluster_df$sample
    colData(sce_obj)$cell <- cluster_df$cell
    colData(sce_obj)$cluster_label <- cluster_df[[grid_use]]
  }
  ## 3.3 obtain pc information
  # load pc matrix
  if (!is.na(pc_file)) {
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
    pc_mtx <- pc_mtx[, match(colnames(sce_obj), colnames(pc_mtx))]
    reducedDims(sce_obj) <- SimpleList(DRM = t(pc_mtx))
    
    
  } else {
    # gene Filtering
    print("Perform PCA as no PCA file provided!")
    
    # if (nrow(sce_obj) < 20000 & ncol(sce_obj) < 5000 &
    #     length(unique(coord_merged$sample)) == 1) {
    #
    # geneFilter <- apply(assays(sce_obj)$counts,1,function(x){
    #   sum(x > 0) >= 10
    # })
    # sce_obj <- sce_obj[geneFilter, ]
    # #
    # FQnorm <- function(counts){
    #   rk <- apply(counts,2,rank,ties.method='min')
    #   counts.sort <- apply(counts,2,sort)
    #   refdist <- apply(counts.sort,1,median)
    #   norm <- apply(rk,2,function(r){ refdist[r] })
    #   rownames(norm) <- rownames(counts)
    #   return(norm)
    # }
    # assays(sce_obj)$norm <- FQnorm(assays(sce_obj)$counts)
    # sce_obj_pca <- prcomp(t(log1p(assays(sce_obj)$norm)), scale. = FALSE)$x
    # } else {
    
    sce_obj_pca <- pca.seurat(counts = count_merged,
                              anno_df = coord_merged,
                              sample_idx = "sample")
    reducedDims(sce_obj) <- SimpleList(DRM = sce_obj_pca)
    # }
  }
  ## 3.4 add UMAP for ct
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
    
    print("Running slingshot...")
    sce_objx  <- slingshot(sce_obj, 
                           clusterLabels = 'cluster_label', 
                           reducedDim = 'DRM',
                           start.clus = start_clustx)
    traj_ind <- grep("slingPseudotime", names(sce_objx@colData@listData))
    traj_mat <- Reduce("cbind", sce_objx@colData@listData[traj_ind]) %>% 
      as.data.frame()
    dimnames(traj_mat) <- list(colnames(sce_objx), 
                               names(sce_objx@colData@listData)[traj_ind])
    pseudotime_list[[paste0("start_", start_clustx)]] <- traj_mat
    
    # fit negative binomial GAM
    if (USE_HVG) {
      
      print("Use HVG for GAM fitting!")
      gam_genes <- CreateAssayObject(counts = count_merged, 
                                     min.cells = 10) %>%
        NormalizeData(., verbose = FALSE) %>%
        FindVariableFeatures() %>%
        VariableFeatures()
    } else {
      
      print("Use all genes for GAM fitting!")
      gam_genes <- rownames(count_merged)
    }
    print("Running GAM...")
    multicoreParam <- MulticoreParam(workers = n_cores)
    sce_obj_gam <- try(fitGAM(sce_objx, parallel = T, 
                              genes = gam_genes,
                              BPPARAM = multicoreParam), 
                       silent = T)
    # test for dynamic expression
    if(inherits(sce_obj_gam, "try-error")){
      
      ATres_list <- "NA"
      warning("Trajectory gene starting from ", start_clustx, " fails!")
    } else {
      
      ATres_list[[paste0("start", start_clustx)]] <- associationTest(sce_obj_gam)
      message(paste0("Trajectory inference starting from ", start_clustx, " is ok!"))
    }
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

# Function 4: call function
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
  start_clus_file <- check_file[6]
  sel_clus_file <- check_file[7]
  pc_file <- check_file[10]
  umap_file <- check_file[11]
  jo_model <- "NA"
  
  result_file <- paste0(out_path, "/traj_call_result.RData")
  clus_file <- paste0(out_path, "/traj_call_cluster.RData")
  
  # set start clusters
  start_clus <- read.table(start_clus_file, header = F, sep = "\t")[, 1, drop = T]
  start_clus <- gsub(" ", "_", start_clus)
  scene_idx <- seq_along(start_clus)
  
  if (submodule == "DECON") {
    
    # load st and annotate data
    spatial_data_list <- str_split(spatial_data_filename, "\\|") %>% unlist() %>%
      lapply(., function(x){
        str_split(x, ",") %>% unlist
      })
    st_list <- lapply(seq_along(spatial_data_list), function(samplex){
      counts <- fread2(spatial_data_list[[samplex]][1])
      coords <- fread2(spatial_data_list[[samplex]][2])
      annos <- fread2(spatial_data_list[[samplex]][3])
      #
      annos$sample <- samplex
      annos$cell <- gsub(" ", "_", annos$cell)
      annos$celltype <- gsub(" ", "_", annos$celltype)
      rownames(counts) <- annos$cell
      return(list("counts" = counts %>% t,
                  "coords" = coords,
                  "anno" = annos))
    })
    
    # 
    anno_df <- lapply(st_list, function(st_lists){
      st_lists[["anno"]]
    }) %>% Reduce("rbind", .)
    
    st_list <- list(
      "count_list" = lapply(st_list, function(st_lists){
        st_lists[["counts"]]
      }),
      "coord_list" = lapply(st_list, function(st_lists){
        st_lists[["coords"]]
      })
    )
    #
    if (!all(start_clus %in% unique(anno_df$celltype))) {
      stop("Not all start clusters in annotation file!")
    } else {
      sel_clus <- start_clus
    }
    #
    grid_use <- "celltype"
    
  } else {
    #
    if (submodule == "CL_jo"){
      
      jo_model <- ifelse(is.na(check_file[8]), "sdd", "ct")
    }
    ## build sce object
    st_list <- h5data.load(spatial_data_filename,         
                           sample_size = sample_size,       
                           load_count = TRUE, 
                           normalization = FALSE,  
                           load_coord = TRUE,    
                           coordinate = TRUE)
    ## define ccc in different modules
    if(grepl("CT", submodule) | (grepl("jo", submodule)&jo_model == "ct")) {
      
      anno_df <- fread2(check_file[8])
    } else if(grepl("SDD", submodule) | (grepl("jo", submodule)&jo_model == "sdd")) {
      
      anno_df <- fread2(check_file[9])
    }
    #
    anno_df[[grid_use]] <- gsub(" ", "_", anno_df[[grid_use]])
    
    sel_clus <- read.table(sel_clus_file, header = F, sep = "\t")[, 1, drop = T]
  }
  
  ##
  traj_result_list <- lapply(scene_idx, function(scenex){
    
    # set cluster
    start_clusx <- start_clus[scenex] %>% 
      as.character()
    sel_clusx <- sel_clus[scenex]
    if (toupper(sel_clusx) == "ALL"){
      sel_clusx <- unique(anno_df[[grid_use]]) %>% as.character() %>% sort
    } else {
      sel_clusx  <- sel_clusx %>% 
        strsplit(",") %>% 
        unlist %>% 
        gsub(" ", "_", .)
    }
    
    # subset on start and select clusters
    anno_dfx <- subset(anno_df, as.character(anno_df[[grid_use]]) %in% sel_clusx)
    st_listx <- lapply(seq_along(st_list[["count_list"]]), function(samplex){
      cellx <- anno_dfx$cell[anno_dfx$sample == samplex]
      idx <- match(cellx, colnames(st_list[["count_list"]][[samplex]]))
      return(list("counts" = st_list[["count_list"]][[samplex]][,idx],
                  "coords" = st_list[["coord_list"]][[samplex]][idx,]))
    })
    st_listx <- list("count_list" = lapply(st_listx, function(st_lists){st_lists[["counts"]]}),
                     "coord_list" = lapply(st_listx, function(st_lists){st_lists[["coords"]]}))
    # run traj analysis
    traj_resultx <- traj.func(st_list = st_listx,
                              cluster_df = anno_dfx,
                              grid_use = grid_use,
                              methods = methods,
                              pc_file = pc_file,
                              umap_file = umap_file,
                              start_clust = start_clusx)
    return(list("traj_result" = traj_resultx,
                "anno_df" = anno_dfx))
  })
  
  # format output
  traj_list <- lapply(traj_result_list, function(traj_result_listx){
    traj_result_listx[["traj_result"]]
  })
  cl_list <- lapply(traj_result_list, function(traj_result_listx){
    traj_result_listx[["anno_df"]]
  })
  names(traj_list) <- names(cl_list) <- paste0("scene", scene_idx)
  #
  save(traj_list, file = result_file)
  save(cl_list, file = clus_file)
  
  ## output
  write.table(c(result_file, clus_file), 
              file = paste0(out_path, "/traj_call_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
  
}

# Function 5: post function
traj.post <- function(out_path = NULL                    ## String: output path of traj procedure
){
  
  ## 1. load data
  check_file <- read.table(paste0(out_path, "/traj_check_file.txt"))[, 1]
  submodule <- check_file[2]
  methods <- check_file[3]
  start_clus_file <- check_file[6]
  sel_clus_file <- check_file[7]
  jo_model <- "NA"
  if (submodule == "CL_jo"){
    
    jo_model <- ifelse(is.na(check_file[8]), "sdd", "ct")  
  }
  call_file <- read.table(paste0(out_path, "/traj_call_file.txt"),
                          header = F, sep = "\t")[, 1]
  load(call_file[1])
  clus_file <- call_file[2]
  
  ## 2. output 
  pseudo_result_file <- ATres_result_file <- NA
  result_dir <- paste0(out_path, "/traj_result/", submodule, "_", methods, "/")
  if (!file.exists(result_dir)) {
    system(paste0("mkdir -p ", result_dir))
  }
  
  ##
  pseudo_result_file <- plyr::laply(names(traj_list), function(a){
    
    pseudo_result_file_a <- paste0(result_dir, "pseudo_", a,".txt")
    pseudo_a <- cbind(traj_list[[a]][["sample"]],
                      traj_list[[a]][["coord_df"]],
                      traj_list[[a]][["cluster_label"]],
                      traj_list[[a]][["pseudotime"]][[1]])
    colnames(pseudo_a) <- c("sample", "cell", "x", "y", "cluster_label",
                            names(traj_list[[a]][["pseudotime"]][[1]]))
    if (methods == "Seurat") {
      pseudo_a$UMAP_1 <- traj_list[[a]][["umap"]][,1]
      pseudo_a$UMAP_2 <- traj_list[[a]][["umap"]][,2]
    }
    write.table(pseudo_a, file = pseudo_result_file_a, 
                sep = "\t", row.names = F, quote = F)
    return(pseudo_result_file_a)
  }) %>% paste(., collapse = ",")
  
  ##
  ATres_result_file <- plyr::laply(names(traj_list), function(a){
    #
    if (traj_list[[a]][["ATres"]] == "NA"){
      
      ATres_result_file_a <- NA
    } else {
      ATres_result_file_a <- paste0(result_dir, "ATres_", a,".txt")
      ATres_a <- rownames_to_column(traj_list[[a]][["ATres"]][[1]], 
                                    var = "gene")
      ATres_a <- ATres_a[, c("gene", "meanLogFC", "waldStat", "pvalue")]
      write.table(ATres_a, file = ATres_result_file_a, 
                  sep = "\t", row.names = F, quote = F)
      
    }
    return(ATres_result_file_a)
  })
  if (all(is.na(unlist(ATres_result_file)))) {
    ATres_result_file <- "NA"
  } else {
    ATres_result_file <- paste(ATres_result_file, collapse = ",")
  }
  
  #  
  write.table(c(pseudo_result_file,
                ATres_result_file,
                start_clus_file, sel_clus_file,
                methods, submodule), 
              file = paste0(out_path, "/traj_post_file.txt"), 
              sep = "\t", col.names = F, row.names = F, quote = F)
  return(0)
}
