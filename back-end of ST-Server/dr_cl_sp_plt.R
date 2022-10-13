#! /usr/bin/env Rscript
# visualize UMAP and location plot after using dr_cl_sp analysis
# up-stream code: dr_cl_sp.R

library(dplyr)
library(bigreadr)
library(ggplot2)
library(hdf5r)
library(glue)

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
dr_cl_sp_plt.plot <- function(data_path1,                              ## String: output path of dr_cl procedure
                              data_path2,                              ## String: output path of qc procedure
                              out_path,                                ## String: output path for dr_cl_sp procedure
                              ft_marker_gene_list = NULL,
                              ft_marker_num,
                              bb_marker_gene_list = NULL,
                              bb_marker_num,
                              out_figures
){
  
  source(paste0(method_path, "/plt_utils.R"))
  result_dir <- paste0(out_path, "/dr_cl_sp_result")
  if (!file.exists(result_dir)){
    
    system(paste0("mkdir ", result_dir))
  }
  loc_plt <- loc.plot(data_path1 = data_path1,
                      data_path2 = data_path2,
                      mode_usage = "dr_cl_sp", 
                      out_path = out_path, 
                      vis_type = "cell_type",
                      out_figure = TRUE)
  if (!is.null(ft_marker_gene_list) & !is.null(ft_marker_num)){
    
    ft_marker_num <- NULL
  }
  ft_plt <- feature.plot(data_path1 = data_path1,
                         data_path2 = data_path2,
                         mode_usage = "dr_cl_sp",
                         marker_gene_list = ft_marker_gene_list,
                         marker_num = ft_marker_num,
                         out_path = out_path,
                         out_figure = T)
  if (!is.null(bb_marker_gene_list) & !is.null(bb_marker_num)){
    
    bb_marker_num <- NULL
  }
  bb_plt <- bubble.plot(data_path1 = data_path1,
                        data_path2 = data_path2,
                        mode_usage = "dr_cl_sp",
                        marker_gene_list = bb_marker_gene_list,
                        marker_num = bb_marker_num,
                        out_path = out_path,
                        out_figure = T)

  save(loc_plt, ft_plt, bb_plt, 
       file = paste0(out_path, "/dr_cl_sp_result/plot.RData"))
  
  return(0)
}


# ###################
# ### test code
# data_path <- "/net/mulan/disk2/yasheng/stwebProject/test/qc"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/test/dr_cl_sp"
# dr_cl_sp_plt.plot(data_path1 = output_path,
#                   data_path2 = data_path,
#                   out_path = output_path,
#                   ft_marker_gene_list = NULL,
#                   ft_marker_num = 10,
#                   bb_marker_gene_list = NULL,
#                   bb_marker_num = 50,
#                   out_figures = TRUE)
