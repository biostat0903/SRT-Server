#! /usr/bin/env Rscript
# visualize UMAP and location plot after using sdd analysis
# up-stream code: sdd.R

# Load packages
library(dplyr)
library(bigreadr)
library(ggplot2)

# Set method_path
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

CT_COLS <- c("#98C1D9", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51", 
             "#E0FBFC", "#A8DADC", "#3D5A80", "#81B29A", "#E07A5F", 
             "#DBC9D8", "#b388eb", "#A4277C", "#BC93B2", "#0077b6", 
             "#BB3E03", "#FFDDD2", "#F19C79", "#006D77", "#6A3569",
             "#E6194B", "#4363D8", "#FFE119", "#3CB44B", "#F58231", 
             "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", 
             "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000", 
             "#AAFFC3", "#808000", "#FFD8B1", "#000075", "#808080")

# Function 1
sdd_plt.plot <- function(data_path1,                                ## String: output path for sdd procedure
                         data_path2,                                ## String: output path for qc procedure
                         out_path,                                  ## String: output path for sdd_plt procedure
                         ft_marker_gene_list = NULL,
                         ft_marker_num,
                         bb_marker_gene_list = NULL,
                         bb_marker_num,
                         out_figures, 
                         zip_figures = FALSE
){
  
  source(paste0(method_path, "/plt_utils.R"))
  ## load ct check file
  sdd_param <- read.table(paste0(data_path1, "/sdd_check_file.txt"))[, 1]
  sdd_methods <- sdd_param[1]
  result_dir <- paste0(out_path, "/sdd_result")
  if (!file.exists(result_dir)){
    
    system(paste0("mkdir ", result_dir))
  }
  loc_plt2 <- loc.plot(data_path1 = data_path1, 
                       data_path2 = data_path2, 
                       mode_usage = "sdd", 
                       out_path = out_path, 
                       vis_type = "spatial_domain",
                       out_figure = out_figures, 
                       zip_figure = zip_figures)
  if (!is.null(ft_marker_gene_list) & !is.null(ft_marker_num)){
    
    ft_marker_num <- NULL
  }
  ft_plt <- feature.plot(data_path1 = data_path1,
                         data_path2 = data_path2,
                         mode_usage = "sdd",
                         marker_gene_list = ft_marker_gene_list,
                         marker_num = ft_marker_num,
                         out_path = out_path,
                         out_figure = out_figures, 
                         zip_figure = zip_figures)
  if (!is.null(bb_marker_gene_list) & !is.null(bb_marker_num)){
    
    bb_marker_num <- NULL
  }
  bb_plt <- bubble.plot(data_path1 = data_path1,
                        data_path2 = data_path2,
                        mode_usage = "sdd",
                        marker_gene_list = bb_marker_gene_list,
                        marker_num = bb_marker_num,
                        out_path = out_path,
                        out_figure = out_figures, 
                        zip_figure = zip_figures)
  ## choose different ct methods
  save(loc_plt2, ft_plt, bb_plt, 
       file = paste0(out_path, "/sdd_result/plot.RData"))

  
  return(0)
}
