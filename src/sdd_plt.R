#! /usr/bin/env Rscript
# visualize UMAP and location plot after using sdd analysis
# up-stream code: sdd.R

# Load packages
library(dplyr)
library(bigreadr)
library(ggplot2)
library(Seurat)

# Set method_path
method_path <- "/public/home/biostat03/project/stwebProject/01_code/srt_server/dev"

# Set colors
CL_COLS <- c("#FD7446", "#709AE1", "#31A354", "#9EDAE5", "#DE9ED6",
             "#BCBD22", "#CE6DBD", "#DADAEB", "#FFFF00", "#FF9896",
             "#91D1C2", "#C7E9C0" , "#6B6ECF", "#7B4173", "#FABEBE", 
             "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000", 
             "#AAFFC3", "#808000", "#FFD8B1", "#000075", "#808080")

# Function 1: de in Seurat
de.func <- function(st_list, 
                    cluster_df,
                    grid_use, 
                    de_method
){
  
  seurat_obj <- purrr::reduce(st_list[["count_list"]], function(x, y) {
    cbind(x = x, y = y)
  }) %>% CreateSeuratObject
  cluster_df <- cluster_df[, c("sample", "cell", grid_use)]
  colnames(cluster_df)[3] <- "cluster_label"
  ## add cluster information
  seurat_obj@meta.data$cluster_label <- cluster_df$cluster_label %>%
    factor(., levels = sort(unique(.)))
  Idents(seurat_obj) <- seurat_obj@meta.data$cluster_label
  ## identify markers
  markers <- FindAllMarkers(seurat_obj,
                            test.use = de_method,
                            only.pos = T,
                            min.pct = 0.1, 
                            return.thresh = 0.1,
                            logfc.threshold = 0)
  
  return(markers)
}

# Function 2: define de genes
define.de <- function(data_path1,                                   ## String: output path for sdd procedure
                      data_path2,                                   ## String: output path for qc procedure
                      result_dir                                    ## String: output path for sdd_plt procedure
){
 
  cl_clus <- fread2(paste0(data_path1, "/domain_label.txt"))
  # cl_dir <- ifelse(file.exists(paste0(data_path1, "/sdd_result")), 
  #                  paste0(data_path1, "/sdd_result"), 
  #                  paste0(data_path1, "/ct_result"))
  # cl_clus <- list.files(cl_dir)[grep("domain_label", list.files(cl_dir))] %>%
  #   paste0(cl_dir, "/", .) %>%
  #   fread2()
  cl_type <- ifelse(grepl("sdd_result", data_path1), "sd", "ct")
  ## define deg 
  source(paste0(method_path, "/io.R"))
  call_file <- read.table(paste0(data_path2, "/qc_call_file.txt"))[, 1]
  spatial_data_filename <- call_file[1]
  sample_size <- call_file[2]
  st_list <- h5data.load(spatial_data_filename,         
                         sample_size = sample_size,       
                         load_count = TRUE,     
                         normalization = FALSE,  
                         load_coord = FALSE,    
                         coordinate = TRUE)
  ## output
  for (grid in colnames(cl_clus)[3: ncol(cl_clus)]){

    cl_markers <- de.func(st_list = st_list, 
                          cluster_df = cl_clus,
                          grid_use = grid, 
                          de_method = "wilcox")
    cl_de_mat <- data.frame(gene = cl_markers$gene,
                            cluster = cl_markers$cluster,
                            log2FC = cl_markers$avg_log2FC,
                            P = cl_markers$p_val,
                            p_adj = cl_markers$p_val_adj)
    cl_de_mat <- cl_de_mat[cl_de_mat$p_adj < 0.01, ]
    write.csv(cl_de_mat, row.names = F, quote = F,
              file = paste0(result_dir, "/de_", grid, ".csv"))
    ### color
    cl_num <- length(unique(cl_markers$cluster))
    cl_col <- data.frame(cluster = paste0("cluster", c(1: cl_num)),
                         color = CL_COLS[c(1: cl_num)])
    write.csv(cl_col, row.names = F, quote = F,
              file = paste0(result_dir, "/clus_col_match_", grid, ".csv"))
  }
  cl_grid <- colnames(cl_clus)[3: ncol(cl_clus)]
  de_loc <- data.frame(de = paste0("de_", cl_grid, ".csv"), 
                       fig = paste0("Location_", cl_type, "_", cl_grid, ".tiff"))
  write.csv(de_loc, file = paste0(result_dir, "/de_loc_match.csv"), 
            row.names = F, quote = F)
  return(0)
}

# Function 3: sdd plot
sdd_plt.plot <- function(data_path1,                                ## String: output path for sdd procedure
                         data_path2,                                ## String: output path for qc procedure
                         out_path,                                  ## String: output path for sdd_plt procedure
                         de_active = FALSE,                         ## Boolean: 
                         ft_marker_gene_list = NULL,
                         ft_marker_num = NULL,
                         bb_marker_gene_list = NULL,
                         bb_marker_num = NULL,
                         out_figures, 
                         zip_figures = FALSE
){
  
  source(paste0(method_path, "/plt_utils.R"))
  ## load check file
  sdd_param <- read.table(paste0(data_path1, "/sdd_check_file.txt"))[, 1]
  sdd_methods <- sdd_param[1]
  result_dir <- paste0(out_path, "/sdd_result/", gsub("-", "_", sdd_methods), "/")
  if (!file.exists(result_dir)){
    
    system(paste0("mkdir ", result_dir))
  }
  ## location plot
  loc_plt2 <- loc.plot(data_path1 = data_path1, 
                       data_path2 = data_path2, 
                       mode_usage = "sdd", 
                       out_path = result_dir, 
                       vis_type = "spatial_domain",
                       out_figure = out_figures, 
                       zip_figure = zip_figures)
  if (de_active == TRUE){
    
    result_dir_sdd <- paste0(out_path, "/sdd_result/", 
                             gsub("-", "_", sdd_methods), "/")
    define.de(data_path1 = result_dir_sdd, 
              data_path2 = data_path2, 
              result_dir = result_dir)
  }
  ## feature plot
  if (!is.null(ft_marker_gene_list) & !is.null(ft_marker_num)){
    
    ft_marker_num <- NULL
  }
  ft_plt <- feature.plot(data_path1 = data_path1,
                         data_path2 = data_path2,
                         mode_usage = "sdd",
                         marker_gene_list = ft_marker_gene_list,
                         marker_num = ft_marker_num,
                         out_path = result_dir,
                         out_figure = out_figures, 
                         zip_figure = zip_figures)
  
  ## 
  if (!is.null(bb_marker_gene_list) & !is.null(bb_marker_num)){
    
    bb_marker_num <- NULL
  }
  bb_plt <- bubble.plot(data_path1 = data_path1,
                        data_path2 = data_path2,
                        mode_usage = "sdd",
                        marker_gene_list = bb_marker_gene_list,
                        marker_num = bb_marker_num,
                        out_path = result_dir,
                        out_figure = out_figures, 
                        zip_figure = zip_figures)
  ## choose different ct methods
  save(loc_plt2, ft_plt, bb_plt, 
       file = paste0(out_path, "/sdd_result/plot.RData"))
  
  return(0)
}
