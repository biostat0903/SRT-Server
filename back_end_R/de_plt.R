#! /usr/bin/env Rscript
# Visulize DEGs
# up-stream code: de.R
# down-stream procedure: none

# Set method_path
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# load packages
library(dplyr)
library(bigreadr)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(reshape2)

COLOR_USE = c("#98C1D9", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51", 
              "#E0FBFC", "#A8DADC", "#3D5A80", "#81B29A", "#E07A5F", 
              "#DBC9D8", "#b388eb", "#A4277C", "#BC93B2", "#0077b6", 
              "#BB3E03", "#FFDDD2", "#F19C79", "#006D77", "#6A3569",
              "#E6194B", "#4363D8", "#FFE119", "#3CB44B", "#F58231", 
              "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", 
              "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000", 
              "#AAFFC3", "#808000", "#FFD8B1", "#000075", "#808080")

# Function 0: extract top genes from DE result
get.TopDEgene <- function(data_path,
                          n_top = 10){
  ## load DE result
  de_post <- paste0(data_path, "/de_post_file.txt")
  if (!file.exists(de_post)) {
    stop("Please run de procedure first!")
  }
  de_post_file <- read.table(de_post)[, 1]
  if(de_post_file[1] == "0"){
  
    de_mat <- fread2(de_post_file[2])  
  } else {
    de_mat <- fread2(de_post_file[1])
  }
  ## top genes
  top_gene <- de_mat %>%
    group_by(cluster) %>%
    top_n(n = n_top, wt = -p_adj) %>%
    top_n(n = n_top, wt = log2FC)
  marker_gene <- top_gene[!duplicated(top_gene$gene), ]$gene
  
  return(marker_gene)
}

# Function 1: Visualize heatmap plot
heatmap.visualize <- function(datt,
                              cluster_info,
                              type, 
                              show_genes = NULL,
                              cols
){
  
  datt <- as.matrix(datt)
  ## row annotation
  if (is.null(show_genes)) {
    show_genes <- rownames(datt)
  }
  gene_pos <- match(show_genes,rownames(datt))
  row_anno <- rowAnnotation(gene=anno_mark(at = gene_pos,
                                           labels = show_genes))
  ## top annotation
  cluster_info <- as.factor(cluster_info)
  names(cols) <- levels(cluster_info)
  cols <- cols[levels(cluster_info)]
  if (type == "ct"){
    
    top_anno <- HeatmapAnnotation(Celltype = cluster_info,
                                  col = list(Celltype = cols),
                                  annotation_legend_param = list(Celltype = list(title = "Celltype"),
                                                                 just = c("right", "bottom")),
                                  show_legend = T,
                                  show_annotation_name = F,
                                  annotation_height = 0.1)
  } else {
    
    top_anno <- HeatmapAnnotation(Domain = cluster_info,
                                  col = list(Domain = cols),
                                  annotation_legend_param = list(Domain = list(title = "Domain"),
                                                                 just = c("right", "bottom")),
                                  show_legend = T,
                                  show_annotation_name = F,
                                  annotation_height = 0.1)
  }
  ## plot heatmap
  col_fun = colorRamp2(c(-4,-2, 0, 1, 2), viridis(5))
  plt <- Heatmap(datt,
                 na_col = "red",
                 cluster_rows = F,
                 cluster_columns = F,
                 show_column_names = F,
                 show_row_names = F,
                 column_split = cluster_info,
                 top_annotation = top_anno, 
                 column_title = NULL,
                 right_annotation = row_anno,
                 # raster_by_magick = F,
                 heatmap_legend_param = list(
                   title = "Expression", 
                   at = c(-2, 0, 2)
                 ),
                 col = col_fun,
                 use_raster = TRUE,
                 raster_resize = TRUE,
                 raster_device = "png")
  return(plt)  
}

# Function 2: Process heatmap data
heatmap.process <- function(data_path1,                       ## String: output path of de procedure
                            data_path2,                       ## String: output path of qc procedure
                            marker_gene
){
  
  ## Load io code
  source(paste0(method_path, "/io.R"))
  
  ## load cluster label
  call_file <- read.table(paste0(data_path1, "/de_call_file.txt"))[, 1]
  load(call_file[2])
  post_file <- read.table(paste0(data_path1, "/de_post_file.txt"))[, 1]
  submodule <- post_file[length(post_file)]
  check_file <- read.table(paste0(data_path1, "/de_check_file.txt"))[, 1]
  jo_model <- "NA"
  if (length(check_file) == 8){
    
    jo_model <- check_file[8]
  }
  ## form count
  call_file <- paste0(data_path2, "/qc_call_file.txt")
  qc_param <- read.table(call_file)[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  ## norm data with transpose
  norm_list <- h5data.load(data_filename = spatial_data_filename, 
                           sample_size = sample_size,
                           load_count = T,
                           normalization = T, 
                           load_coord = F)[["count_list"]]
  # combine norm data and scale
  norm_mat_f <- lapply(norm_list, function(a){
    a <- t(a) %>% as.data.frame() # new added as transpose is canceled in io
    a[, marker_gene]
  }) %>% Reduce("rbind", .)
  scale_df <- scale(norm_mat_f, center = T, scale = T) %>%
    t() %>% as.data.frame()
  heatmap_list <- list(ct_datt = "NULL", ct_annot_df = "NULL", 
                       sdd_datt = "NULL", sdd_annot_df = "NULL")
  if (grepl("CT", submodule) | (grepl("jo", submodule)&jo_model=="ct") ){
    
    heatmap_list[[1]] <- scale_df[, ct_df$cell]
    heatmap_list[[2]] <- ct_df
    
  }
  if (grepl("SDD", submodule) | (grepl("jo", submodule)&jo_model == "sdd") ){
    
    heatmap_list[[3]] <- scale_df[, sdd_df$cell]
    heatmap_list[[4]] <- sdd_df
  }
  
  return(heatmap_list)
}

# Function 3: heatmap in general format
heatmap.plot.func <- function(datt, 
                              annot_df, 
                              grid_use, 
                              type = "ct",
                              marker_gene,
                              out_path, 
                              out_figure, 
                              zip_figure
){

  ## build heatmap object
  annot_df <- annot_df[, c("sample", "cell", grid_use)]
  colnames(annot_df)[3] <- "cluster_label"
  n_ct <- length(unique(annot_df$cluster_label))
  # typex <- ifelse(type == "ct", "Cell Type", "Spatial Domain")
  heatmap_plt <- heatmap.visualize(datt = datt,
                                   cluster_info = annot_df$cluster_label,
                                   type = type, 
                                   show_genes = NULL,
                                   cols = COLOR_USE)
  save(heatmap_plt, file = paste0(out_path, "/de_result/plot.RData"))
  if (out_figure == TRUE){
    
    hmp_ht <- 1.5 + round(length(marker_gene) * 0.18)
    hmp_wt <- 4 + round(n_ct * 1)
    file_name <- paste0(out_path, "/de_result/Heatmap_", type, ".tiff")
    tiff(file_name, height = hmp_ht, width = hmp_wt, units = "in", 
         res = 300, compression = "lzw")
    draw(heatmap_plt)
    dev.off()
    if(zip_figure == TRUE){
      
      system(paste0("gzip -f ", file_name))
    }
  } 
  
  return(0)
}

# Function 4: Visualize heatmap after data process
de_plt.plot <- function(data_path1,                       ## String: output path of de procedure
                        data_path2,                       ## String: output path of qc procedure
                        out_path,                         ## String: output path of de_plt procedure
                        n_top = 5,
                        out_figure = FALSE,
                        zip_figure = FALSE
){
  
  result_dir <- paste0(out_path, "/de_result")
  if (!file.exists(result_dir)){
    
    system(paste0("mkdir ", result_dir))
  }
  ## inputs
  check_file <- read.table(paste0(data_path1, "/de_check_file.txt"))[, 1]
  submodule <- check_file[2]
  grid_use <- check_file[5]
  jo_model <- "NA"
  if (length(check_file) == 8){
    
    jo_model <- check_file[8]
  }
  
  if(grepl("CT", submodule) | (grepl("jo", submodule)&jo_model == "ct")){
    
    marker_gene <- get.TopDEgene(data_path1, n_top = n_top)
    heatmap_res <- heatmap.process(data_path1 = data_path1, 
                                   data_path2 = data_path2,
                                   marker_gene = marker_gene)
    ct_datt <- heatmap_res[[1]]
    ct_df <- heatmap_res[[2]]
    ct_heatmap <- heatmap.plot.func(datt = ct_datt, 
                                    annot_df = ct_df, 
                                    grid_use, 
                                    type = "ct",
                                    marker_gene,
                                    out_path, 
                                    out_figure,
                                    zip_figure)
  }
  if(grepl("SDD", submodule) | (grepl("jo", submodule)&jo_model == "sdd")){
    
    marker_gene <- get.TopDEgene(data_path1, n_top = n_top)
    heatmap_res <- heatmap.process(data_path1 = data_path1, 
                                   data_path2 = data_path2,
                                   marker_gene = marker_gene)
    sdd_datt <- heatmap_res[[3]]
    sdd_df <- heatmap_res[[4]]
    sdd_heatmap <- heatmap.plot.func(datt = sdd_datt, 
                                     annot_df = sdd_df, 
                                     grid_use, 
                                     type = "sdd",
                                     marker_gene,
                                     out_path, 
                                     out_figure,
                                     zip_figure)
  }

  return(0)
}
