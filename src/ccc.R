#! /usr/bin/env Rscript
# Identify cell-cell communication considering the spatial location
# up-stream procedure: sdd.R, ct.R
# down-stream procedure: ccc_plt.R

# Set method_path
method_path <- "/public/home/biostat03/project/stwebProject/01_code/srt_server"
refer_path <- "/public/home/biostat03/project/stwebProject/02_data/reference_data/SRT-Server/"

# Load packages 
library(SpaTalk)
library(CellChat)
library(liana)
library(dplyr)
library(bigreadr)
library(tibble)
library(SingleCellExperiment)
library(Giotto)
library(future)

# Set parameters
N_CORES = 8
CELL_RES_PLATFORM <- c("MERFISH", "Seqfish", "Generic_Cell")
SPOT_RES_PLATFORM <- c("Visium", "SlideseqV2", "Generic_Spot")

# Function 1: ccc.check
ccc.check <- function(data_path1 = NULL,                 ## String: output path of ct or sdd procedure
                      data_path2 = NULL,                 ## String: output path of qc procedure
                      data_path3 = NULL,                 ## String: output path of decon procedure
                      submodule = NULL,                  ## String: cell typing sub-modules: {CT_PCA}, {SDD_sPCA}, {CL_jo}
                      methods = NULL,                    ## String: cell typing methods: {Seurat}, {SpatialPCA}, {BASS}
                      Seurat_resolutions = NULL,         ## String: resolution
                      Seurat_dims = NULL,                ## String: PCA number (re)
                      Seurat_anchors = NULL,             ## String: number of neighbors to use when picking anchors
                      BASS_cellNum = NULL,               ## String: number of cell types. 
                      BASS_domainNum = NULL,             ## String: number of spatial domains.
                      SpatialPCA_domainNum = NULL,       ## String: desired cluster number
                      ccc_species = "Human",             ## String: species to analysis on, could be "Human" or "Mouse" (required)******
                      ccc_method = "spatalk",            ## String: ccc methods
                      ccc_db = "CellPhoneDB",            ## String: LR database to be used for ccc analysis, default is "CellPhoneDB" 
                      out_path = "./"                    ## String: output path of ccc procedure              
){
  
  ## Load general parameters
  check_file <- paste0(data_path2, "/qc_check_file.txt")
  call_file <- paste0(data_path2, "/qc_call_file.txt")
  if(!file.exists(call_file) | !file.exists(check_file)){
    
    stop("No QC file! Please run QC module!")
  } else{
    
    qc_param <- read.table(call_file)[, 1]
    spatial_data_filename <- qc_param[1]
    sample_size <- qc_param[2]
    qc_param <- read.table(check_file)[, 1]
    platform <- qc_param[length(qc_param)-2]
  }
  ## Set parameters for different platforms
  if(platform %in% SPOT_RES_PLATFORM){
    
    post_file <- read.table(paste0(data_path1, "/sdd_post_file.txt"))[,1]
    ### Set clustering joint module parameters
    if(grepl("jo", submodule)){
      
      if(methods == "BASS"){
        
        grid_use <- paste0("R", BASS_domainNum, "_C", BASS_cellNum)
      }
      check_file <- c(spatial_data_filename, submodule, methods,
                      sample_size, grid_use, 
                      ccc_method, ccc_species, ccc_db, NA,
                      post_file[2])
    }
    ### Set spatial domain detection modules
    if(grepl("SDD", submodule)){
      
      if(methods == "SpatialPCA"){
        
        call_file <- read.table(paste0(data_path1, "/sdd_call_file.txt"))[,1]
        clut_num <- as.numeric(call_file[2])
        grid_use <- paste0("clust_", min(as.numeric(SpatialPCA_domainNum), clut_num))
      }
      
      check_file <- c(spatial_data_filename, submodule, methods,
                      sample_size, grid_use, 
                      ccc_method, ccc_species, ccc_db, NA,
                      post_file[2])
    } 
  }
  if(platform %in% CELL_RES_PLATFORM){
    
    post_file <- read.table(paste0(data_path1, "/ct_post_file.txt"))[, 1]
    ### Set clustering joint module parameters
    if(grepl("jo", submodule)){
      
      if(methods == "BASS"){
        
        grid_use <- paste0("R", BASS_domainNum, "_C", BASS_cellNum)
      }
      check_file <- c(spatial_data_filename, submodule, methods,
                      sample_size, grid_use, 
                      ccc_method, ccc_species, ccc_db, NA,
                      post_file[2], post_file[1])
    }
    ### Set PCA modules
    if(grepl("CT", submodule)){
      
      if (methods == "Seurat"){
        
        grid_pc <- ifelse(sample_size == 1,
                          paste0("cluster"),
                          paste0("cluster_anchor", Seurat_anchors,
                                 "_dim", Seurat_dims))
        grid_use <- paste0(grid_pc, "_res", Seurat_resolutions)
      }
      if (methods == "Garnett" | methods == "scSorter"){
        
        grid_use <- "cluster"
      }
      check_file <- c(spatial_data_filename, submodule, methods,
                      sample_size, grid_use, 
                      ccc_method, ccc_species, ccc_db, NA,
                      post_file[1])
    }
  }
  ## Set parameters for spatalk
  if (ccc_method == "spatalk"){
    
    if(platform %in% SPOT_RES_PLATFORM){
      
      if(is.null(data_path3) ) {
        
        stop("spatalk needs the result of DECON!")  
      } else {
        
        check_file[is.na(check_file)] <- data_path3
      }
    }
  }
  ## output file
  write.table(check_file, 
              file = paste0(out_path, "/ccc_check_file.txt"), 
              row.names = F, quote = F, col.names = F)
  
  return(0)
} 

# format_db <- function(ccc_db,
#                       method){
#   
# }

# Function 2: ccc test
ccc.func <- function(st_list = NULL,
                     cluster_df = NULL,
                     grid_use = NULL,
                     decon_path = NULL,
                     platform = NULL, 
                     sample_size = NULL,
                     ccc_type = "ccc",
                     ccc_method,
                     ccc_species,
                     ccc_db
){
  
  
  ## 1. SpaTalk
  if (ccc_method == c("spatalk")) {
    
    ## 0. format expression data and annotation data
    count_merged <- purrr::reduce(st_list[["count_list"]], function(x, y) {
      cbind(x = x, y = y)
    })
    coord_merged <- as.data.frame(Reduce(rbind, st_list[["coord_list"]]))
    colnames(count_merged) <- rownames(coord_merged) <- 
      make.names(colnames(count_merged))
    annot_dat <- data.frame(cell = rownames(coord_merged),
                            x = coord_merged[,1],
                            y = coord_merged[,2])
    annot_dat$celltype <- cluster_df[[grid_use]]
    
    # 1.1 format LR database for SpaTalk
    if (ccc_db == "giotto.mouse") {
      
      LR_data <- fread2(system.file("extdata", 
                                   "mouse_ligand_receptors.txt", 
                                   package = 'Giotto'))
    } else if (ccc_db %in% liana::show_resources()) {
      
      omni_resources <- readRDS(system.file(package = 'liana', "omni_resources.rds"))
      LR_data <- omni_resources[[ccc_db]][,1:2]
    } else if (ccc_db == "SpaTalk.DB"){
      
      LR_data <- SpaTalk::lrpairs
      LR_data <- subset(LR_data, LR_data$species == ccc_species,
                        select = c("ligand", "receptor"))
    }
    
    # 1.2 fit SpaTalk model
    message(paste0("We use ", ccc_db, " to fit ", ccc_method, "!"))
    colnames(LR_data) <- c("ligand", "receptor")
    LR_data$species <- ccc_species
    ccc_list <- list()
    ccc_list[[1]] <- SpaTalk.call(expr_dat = count_merged,
                                  annot_dat = annot_dat,
                                  ccc_species = ccc_species,
                                  decon_path = decon_path,
                                  platform = platform, 
                                  p_thresh = 1,
                                  ccc_method = ccc_method,
                                  ccc_db = LR_data,
                                  ccc_pathway = NULL)
  }
  ## CellChat ("For a SINGLE spatial imaging dataset using CellChat")
  if (ccc_method == "cellchat") {
    # format expression data and annotation data
    ccc_list <- lapply(1:sample_size, function(sample_idx){
      #
      expr_dat <- st_list[["count_list"]][[sample_idx]]
      coord_dat <- st_list[["coord_list"]][[sample_idx]]
      cluster_dat <- subset(cluster_df, cluster_df$sample == sample_idx)
      annot_dat <- data.frame(cell = rownames(coord_dat),
                              x = coord_dat[,1],
                              y = coord_dat[,2])
      annot_dat$celltype <- cluster_df[[grid_use]][match(annot_dat$cell, cluster_df$cell)]
      scale_factor <- st_list[["image_list"]][[sample_idx]][["scale_factors"]] %>% as.list()
      
      # 2.1 set LR database for CellChat
      if (ccc_species == "Human"){
        #
        LR_data <- CellChat::CellChatDB.human
      } else {
        #
        LR_data <- CellChat::CellChatDB.mouse
      }
      ## 2.2 fit CellChat model
      message(paste0("We use ", ccc_db, " to fit ", ccc_method, "!"))
      ccc_list_sample <- CellChat.call(expr_dat = expr_dat,
                                       annot_dat = annot_dat,
                                       ccc_species = ccc_species,
                                       scale_factor = scale_factor, 
                                       p_thresh = 1,
                                       ccc_method = ccc_method,
                                       ccc_db = LR_data)
      return(ccc_list_sample)
    })
  }
  ## output ccc_list
  return(ccc_list)
}

# Function 3: SpaTalk call function
SpaTalk.call <- function(expr_dat = NULL,
                         annot_dat = NULL,
                         p_thresh = 1,
                         ccc_species = NULL,
                         decon_path = NULL,
                         platform = NULL, 
                         ccc_method = NULL,
                         ccc_db = NULL,
                         ccc_pathway = NULL
){
  
  ## create SpaTalk data
  if (is.null(ccc_pathway)) {
    
    ccc_pathway <- SpaTalk::pathways
  }
  if (platform %in% SPOT_RES_PLATFORM){
    
    ## build SpaTalk object
    decon_post <- read.table(paste0(decon_path, "/decon_post_file.txt"))[, 1]
    decon_method <- decon_post[1]
    message(paste0("Using DECON results from ", decon_method, "!"))
    decon_prop <- fread2(decon_post[2])[,-1]
    decon_prop <- column_to_rownames(decon_prop, var = "cell")
    spot_label <- gsub("-", ".", rownames(decon_prop))
    rownames(decon_prop) <- spot_label
    colnames(annot_dat)[1] <- "spot"
    expr_dat_sub <- expr_dat[, colnames(expr_dat) %in% spot_label]
    annot_dat_sub <- annot_dat[annot_dat[, 1] %in% spot_label, ]
    
    ## perform deconvolution by DECON result
    ref_dat <- read.table(paste0(decon_path, "/decon_check_file.txt"))[3, 1]
    ref_dat_path <- sub(refer_path, refer_path, ref_dat)
    load(paste0(ref_dat_path, "/sc_count.RData"))
    load(paste0(ref_dat_path,"/sc_meta.RData"))
    inter_gene <- intersect(rownames(expr_dat_sub), rownames(sc_count))
    sc_count <- sc_count[match(inter_gene, rownames(sc_count)), ]
    expr_dat_sub <- expr_dat_sub[match(inter_gene, rownames(expr_dat_sub)), ]
    celltype <- gsub(" ", "_", sc_meta$cellType)
    obj <- createSpaTalk(st_data = as.matrix(expr_dat_sub),
                         st_meta = annot_dat_sub[, c(1, 2, 3)],
                         species = ccc_species,
                         if_st_is_sc = F,
                         spot_max_cell = 4)
    obj <- dec_celltype(obj, 
                        sc_data = as.matrix(sc_count),
                        sc_celltype = celltype, 
                        if_doParallel = T, 
                        use_n_cores = N_CORES, 
                        dec_result = as.matrix(decon_prop))

  }
  if (platform %in% CELL_RES_PLATFORM){
    
    obj <- createSpaTalk(st_data = as.matrix(expr_dat),
                         st_meta = annot_dat[, -4],
                         species = ccc_species,
                         if_st_is_sc = T,
                         spot_max_cell = 1,
                         celltype = annot_dat$celltype %>% as.character())
  }
  ## Filter LRIs with downstream targets
  obj <- find_lr_path(object = obj, 
                      lrpairs = ccc_db, 
                      pathways = ccc_pathway,
                      if_doParallel = T,
                      use_n_cores = N_CORES)
  
  ## Infer all cell-cell communications 
  if (is.character(obj)) {
    
    stop("ERROR: No ligand-recepotor pairs found in st_data!")
  } else {
    
    obj <- dec_cci_all(object = obj,
                       n_neighbor = 10,
                       min_pairs = 5,
                       min_pairs_ratio = 0,
                       per_num = 1000,
                       pvalue = p_thresh,
                       co_exp_ratio = 0.1,
                       if_doParallel = T,
                       use_n_cores = N_CORES)
    ccc_mat_lr <- obj@lrpair
    ccc_mat_lr <- ccc_mat_lr[,c("celltype_sender", "celltype_receiver",
                                "ligand", "receptor",
                                "lr_co_ratio", "lr_co_ratio_pvalue")]
    colnames(ccc_mat_lr) <- c("source", "target",
                              "ligand", "receptor",
                              "magnitude", "significance")
    obj_pathway <- obj@lr_path$pathways
    ccc_mat_lr$pathway <- obj_pathway$pathway[match(paste0(ccc_mat_lr$ligand, ":", ccc_mat_lr$receptor),
                                                    paste0(obj_pathway$src, ":", obj_pathway$dest))]
    ##
    obj_lr_path <- mclapply(1:nrow(ccc_mat_lr), function(x){
      obj_lr_pathx <- get_lr_path(object = obj,
                                  celltype_sender = ccc_mat_lr$source[x], 
                                  celltype_receiver = ccc_mat_lr$target[x],
                                  ligand = ccc_mat_lr$ligand[x],
                                  receptor = ccc_mat_lr$receptor[x])
      # print(x)
      path_pvaluex <- obj_lr_pathx$path_pvalue
      if (nrow(path_pvaluex) > 0) {
        path_pvaluex <- path_pvaluex[order(path_pvaluex$gene_count, decreasing = T),]
        path_pvaluex <- path_pvaluex[!duplicated(paste0(path_pvaluex$celltype_sender, 
                                                        ":",
                                                        path_pvaluex$celltype_receiver,
                                                        ":",
                                                        path_pvaluex$receptor_pathways)),]
        return(path_pvaluex[,c("celltype_sender", "celltype_receiver", 
                               "receptor_pathways", "pvalue", "gene_count")])
      } else {
        return(NULL)
      }
    }, mc.cores = N_CORES)
    
    #
    ccc_mat_path <- Reduce("rbind", obj_lr_path)
    ccc_mat_path <- ccc_mat_path[order(ccc_mat_path$gene_count, decreasing = T),]
    ccc_mat_path <- ccc_mat_path[!duplicated(paste0(ccc_mat_path$celltype_sender, 
                                                    ":",
                                                    ccc_mat_path$celltype_receiver,
                                                    ":",
                                                    ccc_mat_path$receptor_pathways)),]
    colnames(ccc_mat_path) <- c("source", "target", "pathway", 
                                "significance", "magnitude")
    return(list("ccc_mat_lr" = ccc_mat_lr,
                "ccc_mat_path" = ccc_mat_path))
  }
}

# Function 4: CellChat call function
CellChat.call <- function(expr_dat = NULL,
                          annot_dat = NULL,
                          ccc_species = NULL,
                          scale_factor = NULL, 
                          p_thresh = 1,
                          ccc_method = NULL,
                          ccc_db = NULL
){

  ## Create a CellChat object
  cellchat_obj <- createCellChat(object = expr_dat,
                                 meta = annot_dat[,c("cell", "celltype")],
                                 group.by = "celltype",
                                 datatype = "spatial",
                                 coordinates = annot_dat[,c("x", "y")],
                                 scale.factors = list(spot.diameter = 65, 
                                                      spot = scale_factor$spot_diameter_fullres,
                                                      fiducial = scale_factor$fiducial_diameter_fullres, 
                                                      hires = scale_factor$tissue_hires_scalef, 
                                                      lowres = scale_factor$tissue_lowres_scalef),
                                 assay = NULL,
                                 do.sparse = T)

  ## Set the ligand-receptor interaction database
  cellchat_obj@DB <- ccc_db
  ## Preprocessing the expression data for cell-cell communication analysis
  # subset the expression data of signaling genes for saving computation cost
  cellchat_obj <- subsetData(cellchat_obj) 

  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj,
                                             only.pos = TRUE,
                                             return.object = TRUE,
                                             thresh.pc = 0,
                                             thresh.fc = 0,
                                             thresh.p = p_thresh)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj,
                                                    return.object = TRUE)
  ## Compute the communication probability and infer cellular communication network
  cellchat_obj <- computeCommunProb(cellchat_obj, 
                                    type = "truncatedMean", 
                                    trim = 0.1,
                                    distance.use = TRUE,
                                    interaction.length = 200,
                                    scale.distance = 0.01)
  cellchat_obj <- filterCommunication(cellchat_obj, 
                                      min.cells = 10)
  ## Infer the cell-cell communication at a signaling pathway level
  cellchat_obj <- computeCommunProbPathway(cellchat_obj,
                                           net = NULL,
                                           pairLR.use = NULL,
                                           thresh = p_thresh)
  
  ## extract ccc info
  ccc_mat_lr <- subsetCommunication(cellchat_obj, 
                                    slot.name = "net",
                                    sources.use = NULL,
                                    targets.use = NULL,
                                    signaling = NULL,
                                    pairLR.use = NULL,
                                    thresh = p_thresh)
  ccc_mat_lr <- ccc_mat_lr[,c(1:6, 9)]
  colnames(ccc_mat_lr) <- c("source", "target",
                            "ligand", "receptor",
                            "magnitude", "significance",
                            "pathway")
  #
  ccc_mat_path <- subsetCommunication(cellchat_obj, 
                                      slot.name = "netP",
                                      sources.use = NULL,
                                      targets.use = NULL,
                                      signaling = NULL,
                                      pairLR.use = NULL,
                                      thresh = p_thresh)
  colnames(ccc_mat_path) <- c("source", "target", "pathway", 
                              "magnitude", "significance")
  return(list("ccc_mat_lr" = ccc_mat_lr,
              "ccc_mat_path" = ccc_mat_path))
  
}

# Function 5: ccc.call
ccc.call <- function(out_path            ## String: output path of ccc procedure      
){
  
  ## Load io code
  source(paste0(method_path, "/dev/io.R"))
  
  ## Load inputs
  check_file <- read.table(paste0(out_path, "/ccc_check_file.txt"))[, 1]
  spatial_data_filename <- check_file[1]
  submodule <- check_file[2]
  methods <- check_file[3]
  sample_size <- check_file[4] %>% as.numeric
  grid_use <- check_file[5]
  ccc_method <- check_file[6]
  ccc_species <- check_file[7]
  ccc_db <- check_file[8]
  decon_path <- check_file[9]
  
  ## Load st data
  image_set <- norm_set <- ifelse(ccc_method %in% c("cellchat"),
                                  TRUE, FALSE)
  st_list <- h5data.load(spatial_data_filename,         
                         sample_size = sample_size,       
                         load_count = TRUE,     
                         normalization = norm_set,  
                         load_coord = TRUE,    
                         coordinate = TRUE, 
                         image = image_set)
  platform <- st_list$platform
  ## define ccc in different modules
  result_file <- paste0(out_path, "/ccc_call_result.RData")
  clus_file <- paste0(out_path, "/ccc_call_cluster.RData")
  if(!grepl("CT", submodule)) {
    #
    sdd_df <- fread2(check_file[10])
    sdd_ccc_list <- ccc.func(st_list = st_list,
                             cluster_df = sdd_df,
                             grid_use = grid_use,      
                             decon_path = decon_path,
                             platform = st_list$platform, 
                             sample_size = sample_size,
                             ccc_type = "ddc",
                             ccc_method = ccc_method,
                             ccc_species = ccc_species,
                             ccc_db = ccc_db)
    #
    ct_df <- ct_ccc_list <- NULL
    if(grepl("jo", submodule) & length(check_file) > 10){
      
      ct_df <- fread2(check_file[length(check_file)])
      ct_ccc_list <- ccc.func(st_list = st_list,
                              cluster_df = ct_df,
                              grid_use = grid_use,
                              decon_path = decon_path,
                              platform = st_list$platform,
                              sample_size = sample_size,
                              ccc_type = "ccc",
                              ccc_method = ccc_method,
                              ccc_species = ccc_species,
                              ccc_db = ccc_db)
    }
    save(ct_ccc_list, sdd_ccc_list, file = result_file)
    save(ct_df, sdd_df, file = clus_file)
  } 
  if(grepl("CT", submodule)){
    #
    ct_df <- fread2(check_file[10])
    ct_ccc_list <- ccc.func(st_list = st_list,
                            cluster_df = ct_df,
                            grid_use = grid_use,
                            decon_path = decon_path,
                            platform = st_list$platform,
                            sample_size = sample_size,
                            ccc_type = "ccc",
                            ccc_method = ccc_method,
                            ccc_species = ccc_species,
                            ccc_db = ccc_db)
    save(ct_ccc_list, file = result_file)
    save(ct_df, file = clus_file)
  }
  
  ## output
  write.table(c(result_file, clus_file), 
              file = paste0(out_path, "/ccc_call_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}

# Function 6: ccc post
ccc.post <- function(out_path            ## String: output path of ccc procedure      
){
  
  ## load data
  check_file <- read.table(paste0(out_path, "/ccc_check_file.txt"))[, 1]
  submodule <- check_file[2]
  ccc_method <- check_file[6]
  ccc_db <- check_file[8]
  call_file <- read.table(paste0(out_path, "/ccc_call_file.txt"),
                          header = F, sep = "\t")[, 1]
  load(call_file[1])
  clus_file <- call_file[2]
  
  ## output 
  ct_ccc_file <- sdd_ccc_file <- c(NA, NA)
  result_dir <- paste0(out_path, "/ccc_result/", submodule, "/")
  if (!file.exists(result_dir)) {
    system(paste0("mkdir -p ", result_dir))
  }
  
  if(!is.null(ct_ccc_list)) {

    ## for SpaTalk or CellChat with only one sample
    if(length(ct_ccc_list) == 1){
      
      ## set ccc out file
      ct_ccc_lr_file <- paste0("/ct_ccc_", ccc_method,  "_", ccc_db, "_lr_mat.txt")
      ct_ccc_path_file <- paste0("/ct_ccc_", ccc_method,  "_", ccc_db, "_path_mat.txt")
      
      # output for lr and pathway results
      write.table(ct_ccc_list[[1]][["ccc_mat_lr"]], 
                  file = paste0(result_dir, ct_ccc_lr_file), 
                  sep = "\t", row.names = F, quote = F)
      write.table(ct_ccc_list[[1]][["ccc_mat_path"]], 
                  file = paste0(result_dir, ct_ccc_path_file), 
                  sep = "\t", row.names = F, quote = F)
      pathway_active <- intersect(ct_ccc_list[[1]][["ccc_mat_lr"]]$pathway, 
                                  ct_ccc_list[[1]][["ccc_mat_path"]]$pathway)
      write.table(pathway_active, 
                  file = paste0(result_dir, "/pathway_active.txt"), 
                  sep = "\t", row.names = F, col.names = F, quote = F)
      
      # set ct_ccc_file
      ct_ccc_file <- paste0(result_dir, c(ct_ccc_lr_file, ct_ccc_path_file))
      
    } else {
      
      ## for CellChat with multiple samples
      result_dir_all <- lapply(1:length(ct_ccc_list), function(sample_idx){
        
        ## set ccc out file
        ct_ccc_lr_file <- paste0("/sample", sample_idx, "_ct_ccc_", ccc_method,  "_", 
                                  ccc_db, "_lr_mat.txt")
        ct_ccc_path_file <- paste0("/sample", sample_idx, "_ct_ccc_", ccc_method,  "_", 
                                    ccc_db, "_path_mat.txt")
        
        # output for lr and pathway results
        write.table(ct_ccc_list[[sample_idx]][["ccc_mat_lr"]], 
                    file = paste0(result_dir, ct_ccc_lr_file), 
                    sep = "\t", row.names = F, quote = F)
        write.table(ct_ccc_list[[sample_idx]][["ccc_mat_path"]], 
                    file = paste0(result_dir, ct_ccc_path_file), 
                    sep = "\t", row.names = F, quote = F)
        pathway_active <- intersect(ct_ccc_list[[sample_idx]][["ccc_mat_lr"]]$pathway, 
                                    ct_ccc_list[[sample_idx]][["ccc_mat_path"]]$pathway)
        write.table(pathway_active, 
                    file = paste0(result_dir, "/sample", sample_idx, "_pathway_active.txt"), 
                    sep = "\t", row.names = F, col.names = F, quote = F)
        # return ct_ccc_file
        return(result_dir)
        
      }) %>% unlist
      
      # re-format ct_ccc_file for CellChat with multiple samples
      ct_ccc_file <- c(paste0(result_dir, ct_ccc_lr_file) %>% 
                         paste(., collapse = ","),
                       paste0(result_dir, ct_ccc_path_file) %>% 
                         paste(., collapse = ","))
    }
    
  }
  ##
  if(grepl("SDD", submodule) | grepl("jo", submodule)) {
    
    ## for SpaTalk or CellChat with only one sample
    if(length(sdd_ccc_list) == 1){
      
      ## set ccc out file
      sdd_ccc_lr_file <- paste0("/sdd_ccc_", ccc_method,  "_", ccc_db, "_lr_mat.txt")
      sdd_ccc_path_file <- paste0("/sdd_ccc_", ccc_method,  "_", ccc_db, "_path_mat.txt")
      
      # output for lr and pathway results
      write.table(sdd_ccc_list[[1]][["ccc_mat_lr"]], 
                file = paste0(result_dir, sdd_ccc_lr_file), 
                sep = "\t", row.names = F, quote = F)
      write.table(sdd_ccc_list[[1]][["ccc_mat_path"]], 
                file = paste0(result_dir, sdd_ccc_path_file), 
                sep = "\t", row.names = F, quote = F)
      pathway_active <- intersect(sdd_ccc_list[[1]][["ccc_mat_lr"]]$pathway, 
                                  sdd_ccc_list[[1]][["ccc_mat_path"]]$pathway)
      write.table(pathway_active, 
                  file = paste0(result_dir, "/pathway_active.txt"), 
                  sep = "\t", row.names = F, col.names = F, quote = F)
      # set ct_ccc_file
      sdd_ccc_file <- paste0(result_dir, c(sdd_ccc_lr_file, sdd_ccc_path_file))
      
    } else {
      
      ## for CellChat with multiple samples
      result_dir_all <- lapply(1:length(sdd_ccc_list), function(sample_idx){
        
        ## set ccc out file
        sdd_ccc_lr_file <- paste0("/sample", sample_idx, "_sdd_ccc_", ccc_method,  "_", 
                                  ccc_db, "_lr_mat.txt")
        sdd_ccc_path_file <- paste0("/sample", sample_idx, "_sdd_ccc_", ccc_method,  "_", 
                                    ccc_db, "_path_mat.txt")

        # output for lr and pathway results
        write.table(sdd_ccc_list[[sample_idx]][["ccc_mat_lr"]], 
                    file = paste0(result_dir, sdd_ccc_lr_file), 
                    sep = "\t", row.names = F, quote = F)
        write.table(sdd_ccc_list[[sample_idx]][["ccc_mat_path"]], 
                    file = paste0(result_dir, sdd_ccc_path_file), 
                    sep = "\t", row.names = F, quote = F)
        pathway_active <- intersect(sdd_ccc_list[[sample_idx]][["ccc_mat_lr"]]$pathway, 
                                    sdd_ccc_list[[sample_idx]][["ccc_mat_path"]]$pathway)
        write.table(pathway_active, 
                    file = paste0(result_dir, "/sample", sample_idx, "_pathway_active.txt"), 
                    sep = "\t", row.names = F, col.names = F, quote = F)
        
        # return ct_ccc_file
        return(result_dir)
      }) %>% unlist
      
      # re-format ct_ccc_file for CellChat with multiple samples
      sdd_ccc_file <- c(paste0(result_dir, sdd_ccc_lr_file) %>% 
                          paste(., collapse = ","),
                        paste0(result_dir, sdd_ccc_path_file) %>% 
                          paste(., collapse = ","))
    }
  }
  #
  write.table(c(ct_ccc_file, sdd_ccc_file, 
                ccc_method, ccc_db, submodule), 
              file = paste0(out_path, "/ccc_post_file.txt"), 
              sep = "\t", col.names = F, row.names = F, quote = F)
  return(0)
  
}
