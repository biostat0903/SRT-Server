#! /usr/bin/env Rscript
# Identify cell-cell communication considering the spatial location
# up-stream procedure: sdd.R, ct.R
# down-stream procedure: ccc_plt.R

# Set method_path
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# Load packages 
library(SpaTalk)
library(liana)
library(dplyr)
library(bigreadr)
library(SingleCellExperiment)
library(Giotto)

# Set parameters
N_CORES = 8
CELL_RES_PLATFORM <- c("Merfish", "Seqfish", "Generic_Cell")
SPOT_RES_PLATFORM <- c("Visium", "Slideseq", "Generic_Spot")

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

# Function 2: ccc test
ccc.func <- function(st_list = NULL,
                     cluster_df = NULL,
                     grid_use = NULL,
                     decon_path = NULL,
                     platform = NULL, 
                     ccc_type = "ccc",
                     ccc_method,
                     ccc_species,
                     ccc_db
){
  
  ## format LR database
  if (ccc_method %in% c("spatalk")) {
    
    if (ccc_db == "giotto.mouse") {
      
      LR_data <- fread(system.file("extdata", 
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
  }
  message(paste0("We use ", ccc_db, " to fit ", ccc_method, "!"))
  ## format expression data and annotation data
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
  ## fit SpaTalk model
  if (ccc_method == "spatalk") {
    
    colnames(LR_data) <- c("ligand", "receptor")
    LR_data$species <- ccc_species
    ccc_mat <- SpaTalk.test(expr_dat = count_merged,
                            annot_dat = annot_dat,
                            ccc_species = ccc_species,
                            decon_path = decon_path,
                            platform = platform, 
                            p_thresh = 1,
                            ccc_method = ccc_method,
                            ccc_db = LR_data,
                            ccc_pathway = NULL)
  }
  ## output ccc_mat_list and write call_file
  return(ccc_mat)
}

# Function 2.1: SpaTalk test
SpaTalk.test <- function(expr_dat = NULL,
                         annot_dat = NULL,
                         p_thresh = 1,
                         ccc_species = NULL,
                         decon_path = NULL,
                         platform = NULL, 
                         ccc_method = NULL,
                         ccc_db = NULL,
                         ccc_pathway = NULL,
                         n_cores = N_CORES
){
  
  ## create SpaTalk data
  if (is.null(ccc_pathway)) {
    
    ccc_pathway <- SpaTalk::pathways
  }
  if (platform %in% SPOT_RES_PLATFORM){
    
    ## build SpaTalk object
    load(paste0(decon_path, "/decon_result_list.RData"))
    decon_prop <- plyr::alply(1: length(decon_result_list), 1, function(pp)
      decon_result_list[[pp]]$Proportion_CARD) %>%
      Reduce(rbind, .)
    spot_label <- gsub("-", ".", rownames(decon_prop))
    rownames(decon_prop) <- spot_label
    colnames(annot_dat)[1] <- "spot"
    expr_dat_sub <- expr_dat[, colnames(expr_dat) %in% spot_label]
    annot_dat_sub <- annot_dat[annot_dat[, 1] %in% spot_label, ]

    ## perform deconvolution by DECON result
    ref_dat <- read.table(paste0(decon_path, "/decon_check_file.txt"))[3, 1]
    load(paste0(ref_dat, "/sc_count.RData"))
    load(paste0(ref_dat,"/sc_meta.RData"))
    inter_gene <- intersect(rownames(expr_dat_sub), rownames(sc_count))
    sc_count <- sc_count[match(inter_gene, rownames(sc_count)), ]
    expr_dat_sub <- expr_dat_sub[match(inter_gene, rownames(expr_dat_sub)), ]
    celltype <- gsub(" ", "_", sc_meta$cellType)
    obj <- createSpaTalk(st_data = as.matrix(expr_dat_sub),
                         st_meta = annot_dat_sub[, c(1, 2, 3)],
                         species = ccc_species,
                         if_st_is_sc = F,
                         spot_max_cell = 4)
    # start_time <- Sys.time()
    obj <- dec_celltype(obj, 
                        sc_data = as.matrix(sc_count),
                        sc_celltype = celltype, 
                        if_doParallel = T, 
                        use_n_cores = N_CORES, 
                        dec_result = decon_prop)
    # end_time <- Sys.time()
    
    # save(obj, sc_count, celltype,decon_prop, file = "/net/mulan/disk2/yasheng/stwebProject/01_code/02_run/02_Visium/test.RData")
    # save(decon_prop, file = "/net/mulan/disk2/yasheng/stwebProject/01_code/02_run/02_Visium/test_decon.RData")
    
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
    ccc_mat <- obj@lrpair
    colnames(ccc_mat)[4:5] <- c("source", "target")
    return(ccc_mat)
  }
}

# Function 3: ccc.call
ccc.call <- function(out_path            ## String: output path of ccc procedure      
){
  
  ## Load io code
  source(paste0(method_path, "/io.R"))
  
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
  st_list <- h5data.load(spatial_data_filename,         
                         sample_size = sample_size,       
                         load_count = TRUE,     
                         normalization = FALSE,  
                         load_coord = TRUE,    
                         coordinate = TRUE)
  platform <- st_list$platform
  ## define ccc in different modules
  result_file <- paste0(out_path, "/ccc_call_result.RData")
  clus_file <- paste0(out_path, "/ccc_call_cluster.RData")
  if(!grepl("SDD", submodule)) {

    sdd_df <- sdd_ccc_mat <- NULL
    if(grepl("jo", submodule) & length(check_file) > 10){

      sdd_df <- fread2(check_file[length(check_file)])
      sdd_ccc_mat <- ccc.func(st_list = st_list,
                              cluster_df = sdd_df,
                              grid_use = grid_use,      
                              decon_path = decon_path,
                              platform = st_list$platform, 
                              ccc_type = "ddc",
                              ccc_method = ccc_method,
                              ccc_species = ccc_species,
                              ccc_db = ccc_db)
    }
    ct_df <- fread2(check_file[10])
    ct_ccc_mat <- ccc.func(st_list = st_list,
                           cluster_df = ct_df,
                           grid_use = grid_use,
                           decon_path = decon_path,
                           platform = st_list$platform,
                           ccc_type = "ccc",
                           ccc_method = ccc_method,
                           ccc_species = ccc_species,
                           ccc_db = ccc_db)
    save(ct_ccc_mat, sdd_ccc_mat, file = result_file)
    save(ct_df, sdd_df, file = clus_file)
  } 
  if(grepl("SDD", submodule)){

    sdd_df <- fread2(check_file[10])
    sdd_ccc_mat <- ccc.func(st_list = st_list,
                            cluster_df = sdd_df,
                            grid_use = grid_use,      
                            decon_path = decon_path,
                            platform = st_list$platform, 
                            ccc_type = "ddc",
                            ccc_method = ccc_method,
                            ccc_species = ccc_species,
                            ccc_db = ccc_db)
    save(sdd_ccc_mat, file = result_file)
    save(sdd_df, file = clus_file)
  }
  
  ## output
  write.table(c(result_file, clus_file), 
              file = paste0(out_path, "/ccc_call_file.txt"), 
              col.names = F, row.names = F, quote = F)
  return(0)
}

# Function 4: ccc post
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
  ct_ccc_mat_file <- sdd_ccc_mat_file <- NA
  result_dir <- paste0(out_path, "/ccc_result/", submodule, "/")
  if (!file.exists(result_dir)) {
    system(paste0("mkdir -p ", result_dir))
  }
  
  if(grepl("CT", submodule)) {
    
    ct_ccc_mat_file <- paste0(result_dir, 
                              "/ct_ccc_", ccc_method,  "_", ccc_db, "_mat.txt")

    write.table(ct_ccc_mat, file = ct_ccc_mat_file, 
                sep = "\t", row.names = F, quote = F)
  }
  if(grepl("SDD", submodule) | grepl("jo", submodule)) {
    
    sdd_ccc_mat_file <- paste0(result_dir, 
                              "/sdd_ccc_", ccc_method,  "_", ccc_db, "_mat.txt")
    write.table(sdd_ccc_mat, file = sdd_ccc_mat_file, 
                sep = "\t", row.names = F, quote = F)
  }
  
  write.table(c(ct_ccc_mat_file, sdd_ccc_mat_file, 
                ccc_method, ccc_db, submodule), 
              file = paste0(out_path, "/ccc_post_file.txt"), 
              sep = "\t", col.names = F, row.names = F, quote = F)
  return(0)
  
}
