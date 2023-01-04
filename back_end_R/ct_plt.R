#! /usr/bin/env Rscript
# visualize UMAP and location plot after using seurat analysis
# up-stream code: dr_cl.R

# Load pacakges
library(ggplot2)
library(dplyr)
library(plyr)
library(bigreadr)

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

# Function 1: Visualize location plot
umap.visualize <- function(datt, 
                           pointsize, 
                           color_in
){
  
  plt <- ggplot(datt, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
    geom_point(alpha = 1, size = pointsize) + 
    scale_color_manual("Cluster", values = color_in) + 
    guides(color = guide_legend(byrow = T, nrow = 9,
                                override.aes = list(size = 2)))+
    
    theme_bw() +
    theme(legend.position = "right",
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 6),
          panel.grid = element_blank())
  return(plt)
}

# Function 2: Process UMAP plot data
umap.process <- function(data_path 
){
  
  ## load umap file
  cl_post_file <- read.table(paste0(data_path, "/ct_post_file.txt"))[,1]
  umap_outpath <- cl_post_file[2]
  load(umap_outpath)
  
  ## load labels
  cluster_label <- read.table(paste0(data_path, "/ct_result/ct_pca_seurat_cluster_label.txt"), 
                              head = T) 
  umap_cluster_list <- lapply(names(umap_list), function(a){
    umap_df_a <- umap_list[[a]][cluster_label$cell, ]
    # corresponding cluster colnames
    cluster_name_a <- colnames(cluster_label)[grep(a, colnames(cluster_label))]
    # combine data
    umap_cluster <- cbind(umap_df_a,
                          cluster_label[,c("cell", "sample", cluster_name_a)])
    colnames(umap_cluster) <- c("UMAP_1", "UMAP_2",
                                "cell", "sample", cluster_name_a)
    return(umap_cluster)
  })
  names(umap_cluster_list) <- names(umap_list)
  return(umap_cluster_list)
}

# Function 3: Visualize UMAP plot after data process
umap.plot <- function(data_path, 
                      out_path, 
                      out_figure = FALSE, 
                      zip_figure = FALSE
){
  
  ## process data
  datt_list <- umap.process(data_path)
  
  ## UMAP plot
  ct_umap_list <- lapply(names(datt_list), function(a){
    
    # plot on each cluster label
    cluster_name_a <- colnames(datt_list[[a]])[-c(1:4)]
    lapply(cluster_name_a, function(cl){
      
      datt <- data.frame(datt_list[[a]][, 1:4],
                         cluster = datt_list[[a]][[cl]] %>% 
                           as.factor())
      ## umap: classify umap by cell type
      ct_umap_plt <- umap.visualize(datt = datt, 
                                    pointsize = 0.3, 
                                    color_in = CT_COLS)
      if(out_figure == TRUE){
        
        plt_ht <- 3
        plt_wt <- 4
        ggsave(filename = paste0(out_path, "/ct_result/UMAP_ct_", cl, ".tiff"),
               plot = ct_umap_plt,
               height = plt_ht, width = plt_wt, units = "in", dpi = 300)
        if(zip_figure == TRUE){
          
          system(paste0("gzip -f ", out_path, "/ct_result/UMAP_ct_",
                        cl, ".tiff"))
        }
      }
      return(ct_umap_plt)
    })
  })
  
  return(ct_umap_list)
}

# Function 4
Seurat.plot.func <- function(data_path1,                              ## String: output path of ct procedure
                             data_path2,                              ## String: output path of qc procedure
                             out_path,                                ## String: output path of ct_plt procedure
                             ft_marker_gene_list = NULL,              ## String: gene list for feature plot
                             ft_marker_num,                           ## String: gene number for feature plot
                             bb_marker_gene_list = NULL,              ## String: gene list for bubble plot
                             bb_marker_num,                           ## String: gene number for bubble plot
                             out_figures,                             ## Boolean: output figure
                             zip_figures
){
  
  ## Load plt_utils
  source(paste0(method_path, "/plt_utils.R"))
  result_dir <- paste0(out_path, "/ct_result")
  if (!file.exists(result_dir)){
    
    system(paste0("mkdir ", result_dir))
  }
  umap_plt <- umap.plot(data_path = data_path1,
                        out_path = out_path,
                        out_figure = out_figures, 
                        zip_figure = zip_figures)
  loc_cl_plt <- loc.plot(data_path1 = data_path1,
                         data_path2 = data_path2,
                         mode_usage = "ct",
                         vis_type = "cell_type",
                         out_path = out_path,
                         out_figure = out_figures, 
                         zip_figure = zip_figures)
  if (!is.null(ft_marker_gene_list) & !is.null(ft_marker_num)){
    
    ft_marker_num <- NULL
  }
  ft_plt <- feature.plot(data_path1 = data_path1,
                         data_path2 = data_path2,
                         mode_usage = "ct",
                         marker_gene_list = NULL,
                         marker_num = ft_marker_num,
                         out_path = out_path,
                         out_figure = out_figures, 
                         zip_figure = zip_figures)
  if (!is.null(bb_marker_gene_list) & !is.null(bb_marker_num)){
    
    bb_marker_num <- NULL
  }
  bb_plt <- bubble.plot(data_path1 = data_path1,
                        data_path2 = data_path2,
                        mode_usage = "ct",
                        marker_gene_list = NULL,
                        marker_num = bb_marker_num,
                        out_path = out_path,
                        out_figure = out_figures, 
                        zip_figure = zip_figures)
  save(umap_plt, loc_cl_plt, ft_plt, bb_plt,
       file = paste0(out_path, "/ct_result/plot.RData"))
  
  return(0)
}

# Function 5
BASS.plot.func <- function(data_path1,                                ## String: output path for ct procedure
                           data_path2,                                ## String: output path for qc procedure
                           out_path,                                  ## String: output path for ct procedure
                           ft_marker_gene_list = NULL,
                           ft_marker_num,
                           bb_marker_gene_list = NULL,
                           bb_marker_num,
                           out_figures,                             ## Boolean: output figure
                           zip_figures
){
  
  source(paste0(method_path, "/plt_utils.R"))
  result_dir <- paste0(out_path, "/ct_result")
  if (!file.exists(result_dir)){
    
    system(paste0("mkdir ", result_dir))
  }
  loc_plt1 <- loc.plot(data_path1 = data_path1, 
                       data_path2 = data_path2, 
                       mode_usage = "ct", 
                       out_path = out_path, 
                       vis_type = "cell_type",
                       out_figure = out_figures, 
                       zip_figure = zip_figures)
  loc_plt2 <- loc.plot(data_path1 = data_path1, 
                       data_path2 = data_path2,  
                       mode_usage = "ct", 
                       out_path = out_path, 
                       vis_type = "spatial_domain",
                       out_figure = out_figures, 
                       zip_figure = zip_figures)
  if (!is.null(ft_marker_gene_list) & !is.null(ft_marker_num)){
    
    ft_marker_num <- NULL
  }
  ft_plt <- feature.plot(data_path1 = data_path1,
                         data_path2 = data_path2,
                         mode_usage = "ct",
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
                        mode_usage = "ct",
                        marker_gene_list = bb_marker_gene_list,
                        marker_num = bb_marker_num,
                        out_path = out_path,
                        out_figure = out_figures, 
                        zip_figure = zip_figures)
  save(loc_plt1, loc_plt2, ft_plt, bb_plt, 
       file = paste0(out_path, "/ct_result/plot.RData"))
  
  return(0)
}

# Function 5
Annot.plot.func <- function(data_path1,                                ## String: output path for ct procedure
                            data_path2,                                ## String: output path for qc procedure
                            out_path,                                  ## String: output path for ct procedure
                            ft_marker_gene_list = NULL,
                            ft_marker_num,
                            bb_marker_gene_list = NULL,
                            bb_marker_num,
                            out_figures,                               ## Boolean: output figure
                            zip_figures
){
  
  source(paste0(method_path, "/plt_utils.R"))
  result_dir <- paste0(out_path, "/ct_result")
  if (!file.exists(result_dir)){
    
    system(paste0("mkdir ", result_dir))
  }
  loc_plt1 <- loc.plot(data_path1 = data_path1, 
                       data_path2 = data_path2, 
                       mode_usage = "ct", 
                       out_path = out_path, 
                       vis_type = "cell_type",
                       out_figure = out_figures, 
                       zip_figure = zip_figures)
  if (!is.null(ft_marker_gene_list) & !is.null(ft_marker_num)){
    
    ft_marker_num <- NULL
  }
  ft_plt <- feature.plot(data_path1 = data_path1,
                         data_path2 = data_path2,
                         mode_usage = "ct",
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
                        mode_usage = "ct",
                        marker_gene_list = bb_marker_gene_list,
                        marker_num = bb_marker_num,
                        out_path = out_path,
                        out_figure = out_figures, 
                        zip_figure = zip_figures)
  save(loc_plt1, ft_plt, bb_plt, 
       file = paste0(out_path, "/ct_result/plot.RData"))
  
  return(0)
}

# Function 6
ct_plt.plot <- function(data_path1,                                ## String: output path for ct procedure
                        data_path2,                                ## String: output path for qc procedure
                        out_path,                                  ## String: output path for ct procedure
                        ft_marker_gene_list = NULL,
                        ft_marker_num,
                        bb_marker_gene_list = NULL,
                        bb_marker_num,
                        out_figures, 
                        zip_figures
){
  
  ## load ct check file
  ct_param <- read.table(paste0(data_path1, "/ct_check_file.txt"))[, 1]
  ct_methods <- ct_param[1]
  
  ## choose different ct methods
  if (ct_methods == "CT_PCA-Seurat"){
    
    ct_result <- Seurat.plot.func(data_path1,                          
                                  data_path2,                          
                                  out_path,                            
                                  ft_marker_gene_list,
                                  ft_marker_num,
                                  bb_marker_gene_list,
                                  bb_marker_num,
                                  out_figures, 
                                  zip_figures)
  }
  if (ct_methods == "CL_jo-BASS"){
    
    ct_result <- BASS.plot.func(data_path1,          
                                data_path2,          
                                out_path,            
                                ft_marker_gene_list,
                                ft_marker_num,
                                bb_marker_gene_list,
                                bb_marker_num,
                                out_figures, 
                                zip_figures)
  }
  
  if (ct_methods %in% c("CT_Annot-Garnett", "CT_Annot-scSorter")){
    
    ct_result <- Annot.plot.func(data_path1,          
                                 data_path2,          
                                 out_path,            
                                 ft_marker_gene_list,
                                 ft_marker_num,
                                 bb_marker_gene_list,
                                 bb_marker_num,
                                 out_figures, 
                                 zip_figures)
  }
  
  return(0)
}
