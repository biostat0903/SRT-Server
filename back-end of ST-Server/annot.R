#! /usr/bin/env Rscript
# unsupervised annotation using the marker database
# up-stream procedure: qc.R

# Load packages
library(Seurat)
library(scSorter)
library(monocle)
library(garnett)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(dplyr)

# Set refer_path and method_path
refer_path <- "/net/mulan/disk2/yasheng/stwebProject/02_data/00_ref"
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# Function 1: annot.check
annot.check <- function(data_path,                         ## String: output path of qc procedure
                        method_use = "scSorter",           ## String: annotation methods: {scSorter}, {Garnett}
                        species = "Human",                 ## String: species: {Human} and {Mouse}                 
                        status = NULL,                     ## String: combine tissue and disease status
                        out_path                           ## String: output path for annot procedure
){
  
  ## Check data
  check_file <- paste0(data_path, "/qc_call_file.txt")
  if(!file.exists(check_file)){
    
    stop("Can't find \"qc_call_file.txt\" file! Please run QC module!")
  } 
  ## Check parameters
  if(method_use == "scSorter"){
    
    marker_file <- paste0(refer_path, "/scSorter/", species, "/", status, "_Marker.csv")
    if(!file.exists(marker_file)){
      
      stop("Marker file does not exist! Please check!")
    } 
    ## Output
    write.table(c(method_use, species, status, 
                  marker_file, "NULL", "NULL"), 
                file = paste0(out_path, "/annot_check_file.txt"), 
                row.names = F, col.names = F, quote = F)
  }
  if(method_use == "Garnett"){
    
    marker_file <- paste0(refer_path, "/scRNA_ref/", species, "/", status, "/marker.txt")
    classifer_file <- paste0(refer_path, "/scRNA_ref/", species, "/", status, "/classifer.RData")
    if(!file.exists(marker_file) & !file.exists(classifer_file)){
      
      stop("We do not have marker or classifer file for your setting!")
    }
    if(!file.exists(classifer_file)){
      
      message("We train the classifer using marker file!")
      ## Output
      write.table(c(method_use, species, status, 
                    "NULL", marker_file, "TRUE"), 
                  file = paste0(out_path, "/annot_check_file.txt"), 
                  row.names = F, col.names = F, quote = F)
    } else {
      
      message("We use the trained classifer!")
      ## Output
      write.table(c(method_use, species, status, 
                    "NULL", classifer_file, "FALSE"), 
                  file = paste0(out_path, "/annot_check_file.txt"), 
                  row.names = F, col.names = F, quote = F)
    }
  }
  
  return(0)
}

# Function 2: scSorter
scSorter.fit <- function(markers, 
                         species, 
                         status,
                         count_dat
){
  
  ## process markers
  markers_num <- plyr::count(markers, "Type")
  message(paste0("From marker dataset, we collect ", 
                 nrow(markers_num), " cell types in ", 
                 status, " of ", species, "."))
  ## fit model
  seurat_obj <- CreateSeuratObject(count_dat)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", 
                              scale.factor = 10000, verbose = F)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                     nfeatures = 2000, verbose = F)
  hvgs <- head(VariableFeatures(seurat_obj), 2000)
  expr <- GetAssayData(seurat_obj)
  topgene_filter <- rowSums(as.matrix(expr)[hvgs, ]!=0) > ncol(expr)*.1
  hvgs <- hvgs[topgene_filter]
  picked_genes <- intersect(markers$Marker, hvgs)
  markers_sel <- markers[markers$Marker %in% picked_genes, ]
  message("After intersecting hvgs and marker gene, we obtain ", 
          nrow(markers_sel), " genes.")
  ##
  expr_marker_genes <- setdiff(rownames(expr), markers_sel$Marker)
  if (length(expr_marker_genes) == 0){
    
    stop("Gene number of st data should be larger than that of marker data.")
  } 
  rts <- scSorter::scSorter(expr, markers_sel)
  cell_type <- rts$Pred_Type
    
  return(cell_type)
}

# Function 3: Garnett
Garnett.fit <- function(train = FALSE, 
                        garnett_file,
                        species, 
                        count_dat
){
  
  ## Process data
  cell_id <- data.frame(cellID = colnames(count_dat),
                        row.names = colnames(count_dat))
  gene_name <- data.frame(gene_short_name = rownames(count_dat), 
                          row.names = rownames(count_dat))
  pdata <- new("AnnotatedDataFrame",data = cell_id)
  fdata <- new("AnnotatedDataFrame",data = gene_name)
  sc_cds <- newCellDataSet(as(count_dat, "dgCMatrix"), 
                           phenoData = pdata, 
                           featureData = fdata)
  sc_cds <- estimateSizeFactors(sc_cds)
  
  ## Two different selections for Garnett model
  if (train == TRUE){
    
    if (species == "Human"){
      
      sc_classifier <- train_cell_classifier(cds = sc_cds,
                                             marker_file = garnett_file,
                                             db=org.Hs.eg.db,
                                             cds_gene_id_type = "SYMBOL",
                                             marker_file_gene_id_type = "SYMBOL")
    } else {
      
      sc_classifier <- train_cell_classifier(cds = sc_cds,
                                             marker_file = garnett_file,
                                             db=org.Mm.eg.db,
                                             cds_gene_id_type = "SYMBOL",
                                             marker_file_gene_id_type = "SYMBOL")
    }
    
  } else {
    
    load(garnett_file)
  }
  
  ## Fit classifier 
  sc_cds <- classify_cells(sc_cds, 
                           sc_classifier,
                           db = org.Hs.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "SYMBOL")
  
  return(pData(sc_cds)$cluster_ext_type)
}

# Function 4: annot.call
annot.call <- function(data_path,                          ## String: output path of qc procedure 
                       out_path                            ## String: output path of annot procedure 
){
  
  ## load st data
  source(paste0(method_path, "/io.R"))
  check_file <- paste0(data_path, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  st_list <- h5data.load(data_filename = spatial_data_filename,
                         sample_size = sample_size,
                         load_count = TRUE,
                         normalization = TRUE, 
                         load_coord = FALSE, 
                         coordinate = TRUE)
  # count_dat <- Reduce("rbind", st_list[["count_list"]]) %>% t
  count_dat <- Reduce("cbind", st_list[["count_list"]]) # transpose is canceled in io
  
  ## load markers
  check_file <- read.table(paste0(out_path, "/annot_check_file.txt"))[, 1]
  method_use <- check_file[1]
  if(method_use == "scSorter"){
    
    markers <- bigreadr::fread2(check_file[5])
    markers <- markers[order(markers$Type), ]
    markers <- markers[markers$Type %in% c("Cancer cell", "B cell",
                                           "Endothelial cell", "Epithelial cell"), ]
    cell_type <- scSorter.fit(markers, count_dat)
  } else {
   
    species <- check_file[2]
    garnett_file <- check_file[5]
    train <- as.logical(check_file[6])
    cell_type <- Garnett.fit(train, garnett_file, species, count_dat)
  }
  
  ## Output
  if (!file.exists(paste0(out_path, "/annot_result"))) {
    system(paste0("mkdir '", out_path, "/annot_result'"))
  }
  cell_type_dat <- cbind(colnames(count_dat), cell_type)
  write.table(cell_type_dat, 
              file = paste0(out_path, "/annot_result/celltype_", 
                            method_use, ".txt"), sep = "\t",
              row.names = F, col.names = F, quote = F)
  
  return(0)
}


# ###################
# ### test code
# data_path <- "/net/mulan/disk2/yasheng/stwebProject/test/qc"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/test/annot"
# annot.check(data_path = data_path,
#             method_use = "Garnett",
#             species = "Human",
#             status = "Breast_BRCA",
#             out_path = output_path)
# annot.call(data_path = data_path, out_path = output_path)
