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
                       data_path3 = NULL,
                       submodule = NULL,                ## String: cell typing sub-modules: {CT_PCA}, {SDD_sPCA}, {CL_jo}, {DECON}
                       methods = NULL,                  ## String: cell typing/decon methods: {Seurat}, {SpatialPCA}, {BASS}, {CARD}, {cell2location}, {tangram}
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
  
  ## load inputs
  start_clus_file <- paste0(data_path3, "/start_clus_file.txt")
  sel_clus_file <- paste0(data_path3, "/sel_clus_file.txt")
  if (!file.exists(start_clus_file) | !file.exists(sel_clus_file)) {
    
    stop("No cluster file! Please provide start clusters file and selected clusters file!")
  } else {
    start_clus <- read.table(start_clus_file, header = F, sep = "\t")[, 1, drop = T]
    sel_clus <- read.table(sel_clus_file, header = F, sep = "\t")[, 1, drop = T]
    
    # check lines
    if(length(start_clus != sel_clus)){
      
      stop("Please ensure a consistant number of lines in the start clusters file and the selected clusters file!")
    } else {
      
      # check match
      check_match <- lapply(seq_along(start_clus), function(nn){
        start_clusx <- start_clus[nn]
        sel_clusx <- sel_clus[nn] %>% strsplit(., ",") %>% unlist()
        return(start_clusx %in% sel_clusx)
      })
      
      if (!all(check_match)) {
        
        stop("Please ensure each selected cluster group to containe the corresponding start cluster!")
      } else {
        cluster_file <- paste0(start_clus_file, ",", sel_clus_file)
      }
    }
  }
  ##
  if (submodule == "DECON") {
    
    decon_file <- paste0(data_path1, "/decon_post_file.txt")
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
      start_ct <- gsub(" ", "_", start_ct)
      check_file <- c(spatial_data_filename, submodule, methods,
                      sample_size, NA, 
                      start_ct, NA,
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
      start_ct <- gsub(" ", "_", start_ct)
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
        check_file <- c(spatial_data_filename, submodule, methods,
                        sample_size, grid_use, 
                        start_ct, NA,
                        post_file[1], NA,
                        post_file[3], NA)
      } else {
        
        post_file <- read.table(paste0(data_path1, "/sdd_post_file.txt"))[,1]
        check_file <- c(spatial_data_filename, submodule, methods,
                        sample_size, grid_use, 
                        NA, start_sdd,
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
                      NA, start_sdd, 
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
    coords$sample <- samplex
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
    seurat_pca <- pca.seurat(counts = count_merged,
                             anno_df = coord_merged,
                             sample_idx = "sample")
    reducedDims(sce_obj) <- SimpleList(DRM = seurat_pca)
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
  pc_file <- check_file[10]
  umap_file <- check_file[11]
  jo_model <- "NA"
  result_file <- paste0(out_path, "/traj_call_result.RData")
  clus_file <- paste0(out_path, "/traj_call_cluster.RData")
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
      return(list("counts" = counts%>% t,
                  "coords" = coords,
                  "anno" = annos))
    })
    # set start clusters
    anno_df <- lapply(st_list, function(st_lists){
      st_lists[["anno"]]
    }) %>% Reduce("rbind", .)
    
    start_ct <- check_file[6]
    if (toupper(start_ct) == "ALL"){
      start_ct <- unique(anno_df$celltype) %>% as.character() %>% sort
    } else {
      start_ct  <- start_ct %>% 
        strsplit(",") %>% unlist
    }
    # run traj
    decon_traj_result_list <- lapply(start_ct, function(start_ctx){
      
      st_listx <- list(
        "count_list" = lapply(st_list, function(st_lists){
          idx_ct <- st_lists[["anno"]][["celltype"]] %in% start_ctx
          st_lists[["counts"]][, idx_ct, drop = F]
        }),
        "coord_list" = lapply(st_list, function(st_lists){
          idx_ct <- st_lists[["anno"]][["celltype"]] %in% start_ctx
          st_lists[["coords"]][idx_ct, , drop = F]
        })
      )
      anno_dfx <- anno_df[anno_df$celltype %in% start_ct, ]
      decon_traj_resultx <- traj.func(st_list = st_listx,
                                      cluster_df = anno_dfx,
                                      grid_use = "celltype",
                                      methods = methods,
                                      pc_file = pc_file,
                                      umap_file = umap_file,
                                      start_clust = start_ctx)
      return(list("decon_traj_result" = decon_traj_resultx,
                  "anno_df" = anno_df))
    })
    # format output
    decon_traj_list <- lapply(decon_traj_result_list, function(decon_traj_result_listx){
      decon_traj_result_listx[["decon_traj_result"]]
    })
    decon_ct_list <- lapply(decon_traj_result_list, function(decon_traj_result_listx){
      decon_traj_result_listx[["anno_df"]]
    })
    names(decon_traj_list) <- names(decon_ct_list) <- start_ct
    #
    save(decon_traj_list, file = result_file)
    save(decon_ct_list, file = clus_file)
    
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
      # set start clusters
      ct_df <- fread2(check_file[8])
      start_ct <- check_file[6]
      if (toupper(start_ct) == "ALL"){
        
        start_ct <- unique(ct_df[[grid_use]]) %>% 
          as.character() %>% sort
      } else {
        
        start_ct  <- start_ct %>% 
          strsplit(",") %>% unlist
        #
        ct_df <- subset(ct_df, as.character(ct_df[[grid_use]]) %in% start_ct)
        cellx <- ct_df$cell
        st_list <- lapply(seq_along(st_list[["count_list"]]), function(samplex){
          idx <- match(cellx, colnames(st_list[["count_list"]][[samplex]]))
          return(list("counts" = st_list[["count_list"]][[samplex]][,idx],
                      "coords" = st_list[["coord_list"]][[samplex]][idx,]))
        })
        st_list <- list("count_list" = lapply(st_list, function(st_lists){st_lists[["counts"]]}),
                        "coord_list" = lapply(st_list, function(st_lists){st_lists[["coords"]]}))
      }
      #
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
    
    if(grepl("SDD", submodule) | (grepl("jo", submodule)&jo_model == "sdd")) {
      
      # set start clusters
      sdd_df <- fread2(check_file[9])
      start_sdd <- check_file[7]
      if (toupper(start_sdd) == "ALL"){
        
        start_sdd <- unique(sdd_df[[grid_use]]) %>% 
          as.character() %>% sort
      } else {
        
        start_sdd  <- start_sdd %>% 
          strsplit(",") %>% unlist
        #
        sdd_df <- subset(sdd_df, as.character(sdd_df[[grid_use]]) %in% start_sdd)
        cellx <- sdd_df$cell
        st_list <- lapply(seq_along(st_list[["count_list"]]), function(samplex){
          idx <- match(cellx, colnames(st_list[["count_list"]][[samplex]]))
          return(list("counts" = st_list[["count_list"]][[samplex]][,idx],
                      "coords" = st_list[["coord_list"]][[samplex]][idx,]))
        })
        st_list <- list("count_list" = lapply(st_list, function(st_lists){st_lists[["counts"]]}),
                        "coord_list" = lapply(st_list, function(st_lists){st_lists[["coords"]]}))
      }
      #
      sdd_traj_result <- traj.func(st_list = st_list,
                                   cluster_df = sdd_df,
                                   grid_use = grid_use,
                                   methods = methods, 
                                   pc_file = pc_file,
                                   start_clust = start_sdd)
      save(sdd_traj_result, file = result_file)
      save(sdd_df, file = clus_file)
    }
  }
  
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
  start_ct <- check_file[6]
  start_sdd <- check_file[7]
  jo_model <- "NA"
  if (submodule == "CL_jo"){
    
    jo_model <- ifelse(is.na(check_file[8]), "sdd", "ct")  
  }
  call_file <- read.table(paste0(out_path, "/traj_call_file.txt"),
                          header = F, sep = "\t")[, 1]
  load(call_file[1])
  clus_file <- call_file[2]
  
  ## 2. output 
  decon_pseudo_result_file <- ct_pseudo_result_file <- sdd_pseudo_result_file <- 
    decon_ATres_result_file <- ct_ATres_result_file <- sdd_ATres_result_file <- NA
  result_dir <- paste0(out_path, "/traj_result/", submodule, "_", methods, "/")
  if (!file.exists(result_dir)) {
    system(paste0("mkdir -p ", result_dir))
  }
  ##
  if (submodule == "DECON") {
    decon_pseudo_result_file <- plyr::laply(names(decon_traj_list), function(a){
      
      decon_pseudo_result_file_a <- paste0(result_dir, "decon_pseudo_", a,".txt")
      pseudo_a <- cbind(decon_traj_list[[a]][["sample"]],
                        decon_traj_list[[a]][["coord_df"]],
                        decon_traj_list[[a]][["cluster_label"]],
                        decon_traj_list[[a]][["pseudotime"]][[1]])
      colnames(pseudo_a) <- c("sample", "cell", "x", "y", "cluster_label",
                              names(decon_traj_list[[a]][["pseudotime"]][[1]]))
      write.table(pseudo_a, file = decon_pseudo_result_file_a, 
                  sep = "\t", row.names = F, quote = F)
      return(decon_pseudo_result_file_a)
    }) %>% paste(., collapse = ",")
    #
    decon_ATres_result_file <- plyr::laply(names(decon_traj_list), function(a){
      #
      if (decon_traj_list[[a]][["ATres"]] == "NA"){
        
        decon_ATres_result_file_a <- NA
      } else {
        decon_ATres_result_file_a <- paste0(result_dir, "decon_ATres_", a,".txt")
        ATres_a <- rownames_to_column(decon_traj_list[[a]][["ATres"]][[1]], 
                                      var = "gene")
        ATres_a <- ATres_a[, c("gene", "meanLogFC", "waldStat", "pvalue")]
        write.table(ATres_a, file = decon_ATres_result_file_a, 
                    sep = "\t", row.names = F, quote = F)
        
      }
      return(decon_ATres_result_file_a)
    })
    if (all(is.na(unlist(decon_ATres_result_file)))) {
      decon_ATres_result_file <- "NA"
    } else {
      decon_ATres_result_file <- paste(decon_ATres_result_file, collapse = ",")
    }
  }
  
  # ct pseudo reslut
  if(grepl("CT", submodule) | (grepl("jo", submodule)&jo_model == "ct")) {
    
    ct_pseudo_result_file <- plyr::laply(names(ct_traj_result[["pseudotime"]]), function(a){
      
      ct_pseudo_result_file_a <- paste0(result_dir, "ct_pseudo_", a,".txt")
      pseudo_a <- cbind(ct_traj_result[["sample"]],
                        ct_traj_result[["coord_df"]],
                        ct_traj_result[["cluster_label"]],
                        ct_traj_result[["pseudotime"]][[a]])
      colnames(pseudo_a) <- c("sample", "cell", "x", "y", "cluster_label",
                              names(ct_traj_result[["pseudotime"]][[a]]))
      if (methods == "Seurat") {
        pseudo_a$UMAP_1 <- ct_traj_result[["umap"]][,1]
        pseudo_a$UMAP_2 <- ct_traj_result[["umap"]][,2]
      }
      write.table(pseudo_a, file = ct_pseudo_result_file_a, 
                  sep = "\t", row.names = F, quote = F)
      return(ct_pseudo_result_file_a)
    }) %>% paste(., collapse = ",")
    if (all(ct_traj_result[["ATres"]] == "NA")){
      
      ct_traj_result <- "NA"
    } else {
      
      ct_ATres_result_file <- plyr::laply(names(ct_traj_result[["ATres"]]), function(a){
        
        ct_ATres_result_file_a <- paste0(result_dir, "ct_ATres_", a,".txt")
        ATres_a <- rownames_to_column(ct_traj_result[["ATres"]][[a]], var = "gene")
        ATres_a <- ATres_a[, c("gene", "meanLogFC", "waldStat", "pvalue")]
        write.table(ATres_a, file = ct_ATres_result_file_a, 
                    sep = "\t", row.names = F, quote = F)
        
        return(ct_ATres_result_file_a)
      }) %>% paste(., collapse = ",")
    }
  }
  if(grepl("SDD", submodule) | (grepl("jo", submodule)&jo_model == "sdd")) {
    
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
    if (all(sdd_traj_result[["ATres"]] == "NA")){
      
      sdd_ATres_result_file <- "NA"
    } else {
      
      sdd_ATres_result_file <- plyr::laply(names(sdd_traj_result[["ATres"]]), function(a){
        
        sdd_ATres_result_file_a <- paste0(result_dir, "sdd_ATres_", a,".txt")
        ATres_a <- rownames_to_column(sdd_traj_result[["ATres"]][[a]], var = "gene")
        ATres_a <- ATres_a[, c("gene", "meanLogFC", "waldStat", "pvalue")]
        write.table(ATres_a, file = sdd_ATres_result_file_a, 
                    sep = "\t", row.names = F, quote = F)
        return(sdd_ATres_result_file_a)
      }) %>% paste(., collapse = ",")
    }
  }
  
  write.table(c(decon_pseudo_result_file, ct_pseudo_result_file, sdd_pseudo_result_file, 
                decon_ATres_result_file, ct_ATres_result_file, sdd_ATres_result_file,
                start_ct, start_sdd,
                methods, submodule), 
              file = paste0(out_path, "/traj_post_file.txt"), 
              sep = "\t", col.names = F, row.names = F, quote = F)
  return(0)
}
