#! /usr/bin/env Rscript
# Identify cell-cell communication considering the spatial location
# up-stream procedure: dr_cl.R, dr_cl_wr.R, dr_cl_sp.R
# down-stream procedure: ccc_plt.R

# Set method_path
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# Load packages 
library(hdf5r)
library(SpaTalk)
library(liana)
library(bigreadr)
library(dplyr)
library(tibble)
library(Matrix)
library(SingleCellExperiment)
library(scuttle)
library(CellChat)
library(data.table)
library(Giotto)

N_CORES = 4

# Function 1: ccc.check
ccc.check <- function(data_path1 = NULL,                 ## String: output path of ct or sdd procedure
                      data_path2 = NULL,                 ## String: output path of qc procedure
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
      
      grid_pc <- ifelse(sample_size == 1,
                        paste0("cluster"),
                        paste0("cluster_anchor", Seurat_anchors,
                               "_dim", Seurat_dims))
      grid_use <- paste0(grid_pc, "_res", Seurat_resolutions)
    }
    post_file <- read.table(paste0(data_path1, "/ct_post_file.txt"))[,1]
    check_file <- c(spatial_data_filename, submodule, methods,
                    sample_size, grid_use, 
                    ccc_method, ccc_species, ccc_db, 
                    post_file[1])
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
                    ccc_method, ccc_species, ccc_db, 
                    post_file[1], post_file[2])
  }
  
  ## Set spatial domain detection modules
  if(grepl("SDD", submodule)){
    
    if(methods == "SpatialPCA"){
      
      call_file <- read.table(paste0(data_path1, "/sdd_call_file.txt"))[,1]
	  cat(call_file[1], "\n")
      clut_num <- as.numeric(call_file[2])
      grid_use <- paste0("clust_", min(as.numeric(SpatialPCA_domainNum), clut_num))
    }
    post_file <- read.table(paste0(data_path1, "/sdd_post_file.txt"))[,1]
    check_file <- c(spatial_data_filename, submodule, methods,
                    sample_size, grid_use, 
                    ccc_method, ccc_species, ccc_db,  
                    post_file[2])
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
                     ccc_method,
                     ccc_species,
                     ccc_db){
  
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
  # expression data
  count_merged <- purrr::reduce(st_list[["count_list"]], function(x, y) {
    cbind(x = x, y = y)
  })
  coord_merged <- as.data.frame(Reduce(rbind, st_list[["coord_list"]]))
  
  # make names
  colnames(count_merged) <- rownames(coord_merged) <- 
    make.names(colnames(count_merged))
  
  # annotation data
  annot_dat <- data.frame(cell = rownames(coord_merged),
                          x = coord_merged[,1],
                          y = coord_merged[,2])
  annot_dat$celltype <- cluster_df[[grid_use]]#[match(annot_dat$cell, 
                                                     # make.names(cluster_df$cell))]
  
  ## fit SpaTalk model
  if (ccc_method == "spatalk") {
    colnames(LR_data) <- c("ligand", "receptor")
    LR_data$species <- ccc_species
    ccc_mat <- SpaTalk.test(expr_dat = count_merged,
                            annot_dat = annot_dat,
                            ccc_species = ccc_species,
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
                         ccc_method = NULL,
                         ccc_db = NULL,
                         ccc_pathway = NULL,
                         n_cores = N_CORES
){
  
  ## create SpaTalk data
  if (is.null(ccc_pathway)) {
    
    ccc_pathway <- SpaTalk::pathways
  }
  obj <- createSpaTalk(st_data = as.matrix(expr_dat),
                       st_meta = annot_dat[, -4],
                       species = ccc_species,
                       if_st_is_sc = T,
                       spot_max_cell = 1,
                       celltype = annot_dat$celltype %>% as.character())
  
  ## Filter LRIs with downstream targets
  obj <- tryCatch({
    
    find_lr_path(object = obj, 
                 lrpairs = ccc_db, 
                 pathways = ccc_pathway,
                 if_doParallel = T,
                 use_n_cores = n_cores)
  }, error = function(e){
    
    print("ERROR")
  })
  
  
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
                       use_n_cores = n_cores)
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
  
  ## load inputs
  check_file <- read.table(paste0(out_path, "/ccc_check_file.txt"))[, 1]
  spatial_data_filename <- check_file[1]
  submodule <- check_file[2]
  methods <- check_file[3]
  sample_size <- check_file[4] %>% as.numeric
  grid_use <- check_file[5]
  ccc_method <- check_file[6]
  ccc_species <- check_file[7]
  ccc_db <- check_file[8]
    
  ## load st data
  st_list <- h5data.load(spatial_data_filename,         
                         sample_size = sample_size,       
                         load_count = TRUE,     
                         normalization = FALSE,  
                         load_coord = TRUE,    
                         coordinate = TRUE)
    
  ## define ccc in different modules
  result_file <- paste0(out_path, "/ccc_call_result.RData")
  clus_file <- paste0(out_path, "/ccc_call_cluster.RData")
  if(grepl("CT", submodule)) {
    
    ct_df <- fread2(check_file[9])
    ct_ccc_mat <- ccc.func(st_list = st_list,
                       cluster_df = ct_df,
                       grid_use = grid_use,
                       ccc_method = ccc_method,
                       ccc_species = ccc_species,
                       ccc_db = ccc_db)
    
    save(ct_ccc_mat, file = result_file)
    save(ct_df, file = clus_file)
  }
  if(grepl("jo", submodule)) {
    #
    ct_df <- fread2(check_file[9])
    ct_ccc_mat <- ccc.func(st_list = st_list,
                           cluster_df = ct_df,
                           grid_use = grid_use,
                           ccc_method = ccc_method,
                           ccc_species = ccc_species,
                           ccc_db = ccc_db)
    #
    sdd_df <- fread2(check_file[10])
    sdd_ccc_mat <- ccc.func(st_list = st_list,
                            cluster_df = sdd_df,
                            grid_use = grid_use,
                            ccc_method = ccc_method,
                            ccc_species = ccc_species,
                            ccc_db = ccc_db)
    save(ct_ccc_mat, sdd_ccc_mat, file = result_file)
    save(ct_df,sdd_df, file = clus_file)
  }
  if(grepl("SDD", submodule)) {
    
    sdd_df <- fread2(check_file[9])
    sdd_ccc_mat <- ccc.func(st_list = st_list,
                            cluster_df = sdd_df,
                            grid_use = grid_use,
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
  
  ## 1. load data
  check_file <- read.table(paste0(out_path, "/ccc_check_file.txt"))[, 1]
  submodule <- check_file[2]
  ccc_method <- check_file[6]
  ccc_db <- check_file[8]
  call_file <- read.table(paste0(out_path, "/ccc_call_file.txt"),
                              header = F, sep = "\t")[, 1]
  load(call_file[1])
  clus_file <- call_file[2]
  
  ## 2. output 
  ct_ccc_mat_file <- sdd_ccc_mat_file <- NA
  result_dir <- paste0(out_path, "/ccc_result/", submodule, "/")
  if (!file.exists(result_dir)) {
    system(paste0("mkdir -p ", result_dir))
  }
  
  if(grepl("CT", submodule) | grepl("jo", submodule)) {
    
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


# ###################
# ### test code
# data_path1 <- "/pt_data/494433291@qq.com/0a4d7f5a254e45028b4d74a998958526/tools-output/wf-483227122920325696/job-SDD-483233367291068992"
# data_path2 <- "/pt_data/494433291@qq.com/0a4d7f5a254e45028b4d74a998958526/tools-output/wf-483227122920325696/job-QC-483227148551717440"
# output_path <- "/pt_data/494433291@qq.com/0a4d7f5a254e45028b4d74a998958526/tools-output/wf-483227122920325696/job-CCC-483276849707745856"
# ccc.check(data_path1 = data_path1,
          # data_path2 = data_path2,
          # submodule = "SDD_sPCA",
          # methods = "SpatialPCA",
          # SpatialPCA_domainNum = "10",
          # ccc_method = "spatalk",
          # ccc_species = "Mouse",
          # ccc_db = "SpaTalk.DB",
          # out_path = output_path)
# ccc.call(out_path = output_path)
# ccc.post(out_path = output_path)
