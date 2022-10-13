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

N_CORES = 10

# Function 1: ccc.check
ccc.check <- function(data_path1 = NULL,               ## String: output path of dr_cl procedure
                      data_path2 = NULL,               ## String: output path of qc procedure
                      clus_mode = NULL,                ## String: cluster mode, including dr_cl, dr_cl_wr or dr_cl_sp
                      resolution = NULL,               ## Float: resolution of dr_cl
                      anchor_num = NULL,               ## Integer: number of neighbors to use when picking anchors in dr_cl
                      dim_num = NULL,                  ## Integer: PCA number (re) of dr_cl
                      cluster_num = NULL,              ## Integer: number of cell types in dr_cl_wr or dr_cl_sp
                      domain_num = NULL,               ## Integer: number of spatial domains in dr_cl_warpping
                      spc_num = NULL,                  ## Integer: number of spatial PC in dr_cl_sp
                      ccc_species = "Human",           ## String: species to analysis on, could be "Human" or "Mouse" (required)******
                      ccc_method = "spatalk",          ## String: ccc methods
                      ccc_db = "CellPhoneDB",          ## String: LR database to be used for ccc analysis, default is "CellPhoneDB" 
                      out_path = "./"                  ## String: output path of ccc procedure              
){
  
  ## check file
  check_file <- paste0(data_path2, "/qc_call_file.txt")
  if(!file.exists(check_file)){
    
    stop("No QC file! Please run QC module!")
  }else {
    
    qc_param <- read.table(check_file)[, 1]
    spatial_data_filename <- qc_param[1]
    sample_size <- qc_param[2] %>% as.numeric
  }
  dr_cl_post <- read.table(paste0(data_path1, "/", clus_mode,
                                  "_post_file.txt"))[,1]
  dr_cl_cluster_file <- dr_cl_post[1]
  dr_cl_pc_file <- dr_cl_post[length(dr_cl_post)]
  if (!file.exists(dr_cl_cluster_file) | !file.exists(dr_cl_pc_file)) {
    
    stop("dr_cl file is not found!")
  } 
  if(clus_mode == "dr_cl"){
    
    if (sample_size == 1){
      
      anchor_num <- "NULL"
      dim_num <- "NULL"
    } 
    cluster_num <- domain_num <- spc_num <- "NULL"
  }
  if(clus_mode == "dr_cl_sp"){
    
    anchor_num <- resolution <- dim_num <- "NULL"
    cluster_num <- domain_num <- "NULL"
  }
  if(clus_mode == "dr_cl_wr"){
    
    anchor_num <- resolution <- dim_num <- "NULL"
    spc_num <- "NULL"
  }
  
  ## output file
  write.table(c(spatial_data_filename, sample_size, clus_mode,
                anchor_num, dim_num, resolution,
                cluster_num, domain_num, spc_num,
                ccc_method, ccc_species, ccc_db, 
                dr_cl_cluster_file), 
              file = paste0(out_path, "/ccc_check_file.txt"), 
              row.names = F, quote = F, col.names = F)
  return(0)
  
} 

######### SpaTalk #########
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
                 pathways = ccc_pathway)
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

# Function 2: ccc.call
ccc.call <- function(out_path                          ## String: output path of ccc procedure      
){
  
  ## Load io code
  source(paste0(method_path, "/io.R"))
  
  ## load data
  check_file <- read.table(paste0(out_path, "/ccc_check_file.txt"))[, 1]
  spatial_data_filename <- check_file[1]
  sample_size <- check_file[2] %>% as.numeric()
  clus_mode <- check_file[3]
  anchor_num <- ifelse(check_file[4] == "NULL", 0, check_file[4])
  dim_num <- ifelse(check_file[5] == "NULL", 0, as.numeric(check_file[5]))
  resolution <- ifelse(check_file[6] == "NULL", 0, check_file[6])
  cluster_num <- ifelse(check_file[7] == "NULL", 0, as.numeric(check_file[7]))
  domain_num <- ifelse(check_file[8] == "NULL", 0, as.numeric(check_file[8]))
  spc_num <- ifelse(check_file[9] == "NULL", 0, as.numeric(check_file[9]))
  ccc_method <- check_file[10]
  ccc_species <- check_file[11]
  ccc_db <- check_file[12]
  dr_cl_cluster_file <- check_file[13]
  
  ## check LR database
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
  
  ## load st data
  st_list <- h5data.load(spatial_data_filename,         
                         sample_size = sample_size,       
                         load_count = TRUE,     
                         normalization = FALSE,  
                         load_coord = TRUE,    
                         coordinate = TRUE)
  count_merged <- purrr::reduce(st_list[["count_list"]], function(x, y) {
    cbind(x = x, y = y)
  })
  coord_merged <- Reduce(rbind, st_list[["coord_list"]]) %>% 
    as.data.frame()
  
  # make names
  colnames(count_merged) <- rownames(coord_merged) <- 
    make.names(colnames(count_merged))
  
  ## obtain cluster information
  cluster_df <- fread2(dr_cl_cluster_file)
  if (clus_mode == "dr_cl") {
    
    grid_pc <- ifelse(sample_size == 1,
                      paste0("cluster"),
                      paste0("cluster_anchor", anchor_num,
                             "_dim", dim_num))
    grid_use <- paste0(grid_pc, "_res", resolution)
  } else {
    
    grid_use <- ifelse(clus_mode == "dr_cl_wr", 
                       paste0("R", domain_num, "_C", cluster_num),
                       paste0("clust_", cluster_num))
  }
  cluster_df <- cluster_df[,c("sample", "cell", grid_use)]
  colnames(cluster_df)[3] <- "cluster_label"

  ## check meta column
  annot_dat <- data.frame(cell = rownames(coord_merged),
                          x = coord_merged[,1],
                          y = coord_merged[,2],
                          celltype = cluster_df$cluster_label)

  ## fit SpaTalk model
  res_file <- paste0(out_path, "/ccc_call_", 
                     ccc_method, "_", ccc_db, ".RData")
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
    ## output ccc_mat_list and write call_file
    save(ccc_mat, file = res_file)
  }
  write.table(c(res_file, ccc_method, ccc_db),
              file = paste0(out_path, "/ccc_call_file.txt"), 
              sep = "\t", col.names = F, row.names = F, quote = F)
}

# Function 3: ccc post
ccc.post <- function(out_path                          ## String: output path of ccc procedure      
){
  
  ## 1. load
  call_file <- read.table(paste0(out_path, "/ccc_call_file.txt"),
                              header = F, sep = "\t")[, 1]
  ccc_result <- call_file[1]
  ccc_method <- call_file[2]
  ccc_db <- call_file[3]
  check_file <- read.table(file = paste0(out_path, "/ccc_check_file.txt"),
                            header = F, sep = "\t")[, 1]
  clus_mode <- check_file[3]
  
  ## 2. output
  load(ccc_result)
  result_dir <- paste0(out_path, "/ccc_result/", clus_mode, "/")
  if (!file.exists(result_dir)) {
    system(paste0("mkdir -p ", result_dir))
  }
  
  result_path <- paste0(result_dir, 
                        "/ccc_", ccc_method,  "_", ccc_db, "_result.csv")
  write.csv(ccc_mat, 
            file = result_path, 
            row.names = F, quote = F)
  write.table(c(result_path, ccc_method, ccc_db, clus_mode), 
              file = paste0(out_path, "/ccc_post_file.txt"), 
              sep = "\t", col.names = F, row.names = F, quote = F)
  return(0)
}


# ###################
# ### test code
# data_path1 <- "/net/mulan/disk2/yasheng/stwebProject/test/dr_cl_wr"
# data_path2 <- "/net/mulan/disk2/yasheng/stwebProject/test/qc"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/test/ccc"
# ccc.check(data_path1 = data_path1,
#           data_path2 = data_path2,
#           clus_mode = "dr_cl_wr",
#           anchor_num = NULL,
#           dim_num = NULL,
#           resolution = NULL,
#           cluster_num = 10,
#           domain_num = 5,
#           ccc_method = "spatalk",
#           ccc_species = "Human",
#           ccc_db = "SpaTalk.DB",
#           out_path = output_path)
# ccc.call(out_path = output_path)
# ccc.post(out_path = output_path)
