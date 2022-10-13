#! /usr/bin/env Rscript
# Visualize spatially variable gene (SVG) for ST data
# We use pattern plot, qq plot and histogram plot to visualize the result of svg

# Set method_path
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# Load packages 
library(hdf5r)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(tidyr)
library(reshape2)
library(scales) 
library(patchwork)

COLOR_USE = c("#98C1D9", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51", 
              "#E0FBFC", "#A8DADC", "#3D5A80", "#81B29A", "#E07A5F", 
              "#DBC9D8", "#b388eb", "#A4277C", "#BC93B2", "#0077b6", 
              "#BB3E03", "#FFDDD2", "#F19C79", "#006D77", "#6A3569",
              "#E6194B", "#4363D8", "#FFE119", "#3CB44B", "#F58231", 
              "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", 
              "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000", 
              "#AAFFC3", "#808000", "#FFD8B1", "#000075", "#808080")

# Function 1: Process trajectory base data
trajBasedata.process <- function(data_path,
                                 location_type = NULL,
                                 start_clust = NULL
){
  ##
  call_file <- read.table(file = paste0(data_path, "/traj_call_file.txt"))[, 1]
  traj_result_file <- call_file[1]
  clust_mode <- call_file[2]
  load(traj_result_file) # traj_result
  pseudo_name <- paste0("start", start_clust)
  
  ##
  if (!pseudo_name %in% names(traj_result$pseudotime)) {
    stop("Target trajectory not in the result!")
  }
  
  # combind data
  clusterlabels <- traj_result$cluster_label %>% 
    as.factor() %>% as.data.frame()
  sample_info <- traj_result$sample %>% 
    as.data.frame()
  pseudotime <- traj_result$pseudotime[[pseudo_name]][,1, drop = F] %>% 
    as.data.frame()
  traj_df <- cbind(sample_info, clusterlabels, pseudotime)
  colnames(traj_df) <- c("sample", "cell",
                         "cluster_label", 
                         "pseudotime")
  traj_df$sample <- paste0("Sample", traj_df$sample)
  
  if (!is.null(location_type)) {
    
    if (location_type == "coord") {
      
      location <- traj_result$coord_df %>% 
        as.data.frame()
      colnames(location) <- c("x", "y")
      
    } else if (location_type == "umap") {
      
      if (clust_mode != "dr_cl") {
        stop("UMAP can only be plottd for dr_cl!")
      } else {
        location <- traj_result$umap %>% 
          as.data.frame()
        colnames(location) <- c("UMAP_1", "UMAP_2")
      }
      
    } else {
      
      stop("Unrecognized location index!")
      
    }
    
    traj_df <- cbind(traj_df, location)
  }
  
  return(list("traj_df" = traj_df,
              "clust_mode" = clust_mode))
}


# Function 2: Visualize trajectory plot
traj.visualize <- function(datt,
                      arrow_df,
                      color_in,
                      pointsize = 1,
                      arrowlength = 0.2,
                      arrowsize = 1){
  # plot  
  datt$pseudotime <- datt$pseudotime + 0.01
  plt1 <- ggplot(datt, aes(x = x, y = y, color = pseudotime)) +
    geom_point(alpha = 1,size = pointsize) +
    scale_color_viridis_c("Pseudotime") +
    theme_void()+
    theme(plot.title = element_text(size = 15,  face = "bold"),
          text = element_text(size = 15),
          legend.position = "bottom") +
    facet_wrap(~sample, ncol = 1)
  
  plt2 <- ggplot(datt, aes(x = x, y = y)) +
    geom_point(alpha =1,size = pointsize,
               aes(color = cluster_label)) +
    theme_void()+
    scale_colour_manual("Cluster", values = color_in)+
    theme(plot.title = element_text(size = 15,  face = "bold"),
          text = element_text(size = 15),
          legend.position = "bottom")+
    geom_segment(aes(x = start_x,
                     y = start_y,
                     xend = end_x,
                     yend = end_y,
                     colour = "segment"),
                 arrow = arrow(length = unit(arrowlength, "cm")),
                 size = arrowsize,
                 color = "grey20",
                 data = arrow_df)+
    facet_wrap(~sample, ncol = 1)+
    guides(color = guide_legend(ncol = 5, byrow = T, 
                                override.aes = list(size = 2)))
  # output
  return(plt1+plt2)
}

# Function 3: Process trajectory plot data
traj.process <- function(data_path,
                         gridnum = 10,
                         start_clust = NULL
){
  ##
  base_res <- trajBasedata.process(data_path = data_path,
                                   location_type = "coord",
                                   start_clust = start_clust)
  traj_df <- base_res[["traj_df"]]
  clust_mode <- base_res[["clust_mode"]]
  
  ##
  traj_list <- split(traj_df, traj_df$sample)
  arrow_list <- lapply(seq_along(traj_list), function(a){
    # define data
    traj_list_a <- traj_list[[a]]
    location_a <- traj_list_a[,c("x", "y")] %>% as.data.frame()
    
    # get anchors
    space_x <- diff(range(location_a$x))/gridnum
    space_y <- diff(range(location_a$x))/gridnum
    x_anchor <- min(location_a$x) + c(0:gridnum) * space_x
    y_anchor <- min(location_a$y) + c(0:gridnum) * space_y
    
    # label square by num_x, num_y
    pseudotime_use <- traj_list_a$pseudotime
    ar_list <- list()
    for(num_x in 1:gridnum){
      for(num_y in 1:gridnum){
        
        # find points in each grid
        points_in_grid <- which(location_a$x >= x_anchor[num_x] & 
                                  location_a$x <= x_anchor[num_x+1] &
                                  location_a$y >= y_anchor[num_y] & 
                                  location_a$y <= y_anchor[num_y+1])
        pseudotime_in_grid <- pseudotime_use[points_in_grid]
        
        # find min pseudotime and max pseudotime in each grid
        if(length(points_in_grid)>1 & 
           sum(!is.na(pseudotime_in_grid))>0){
          
          min_pseudotime_ind <- points_in_grid[which.min(pseudotime_in_grid)]
          max_pseudotime_ind <- points_in_grid[which.max(pseudotime_in_grid)]
          
          start_x <- location_a[min_pseudotime_ind, "x"]
          start_y <- location_a[min_pseudotime_ind, "y"]
          end_x <- location_a[max_pseudotime_ind, "x"]
          end_y <- location_a[max_pseudotime_ind, "y"]
          ar_list <- c(ar_list,
                       list(c(start_x, start_y, end_x, end_y)))
        }
      }
    }
    ar_df <- Reduce("rbind", ar_list) %>% as.data.frame()
    colnames(ar_df) <- c("start_x", "start_y", "end_x", "end_y")
    ar_df$sample <- paste0("Sample", a)
    return(ar_df)
  })
  
  arrow_df <- Reduce("rbind", arrow_list)
  return(list("traj_df" = traj_df,
              "arrow_df" = arrow_df,
              "clust_mode" = clust_mode))
}

# Function 4: Visualize trajectory plot after data process
traj.plot <- function(data_path, 
                      out_path, 
                      start_clust = NULL,
                      gridnum = 10,
                      out_figure = FALSE
){
  
  ## inputs
  traj_res <- traj.process(data_path = data_path, 
                           gridnum = gridnum,
                           start_clust = start_clust)
  traj_df <- traj_res[[1]]
  arrow_df <- traj_res[[2]]
  clust_mode <- traj_res[[3]]
  sample_size <- length(unique(traj_df$sample))
  n_cl <- length(unique(traj_df$cluster_label))
  
  ## traj plot 
  traj_plt <- traj.visualize(datt = traj_df,
                             arrow_df = arrow_df,
                             color_in = COLOR_USE,
                             pointsize = 1,
                             arrowlength = 0.3,
                             arrowsize = 0.8)
  
  ## output figure
  if (out_figure == TRUE){
    
    traj_ht <- sample_size*4 + ceiling(n_cl/5)*0.2
    traj_wt <- 10
    plt_name <- paste0(out_path, "/traj_result/",clust_mode,"/TrajPlot_", 
                       clust_mode, "_start", start_clust,".tiff")
    ggsave(filename = plt_name, 
           plot = traj_plt, 
           height = traj_ht, width = traj_wt, 
           units = "in", dpi = 300)
    system(paste0("gzip -f ", plt_name))
  } 
  return(traj_plt)
}

# Function 5: plot UMAP plot
trajUMAP.visualize <- function(datt,
                               pointsize = 1
){
  # plot  
  datt$pseudotime <- datt$pseudotime + 0.01
  plt <- ggplot(datt, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
    geom_point(alpha = 1,size = pointsize) +
    scale_color_viridis_c("Pseudotime") +
    theme_bw() +
    theme(legend.position = "right",
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13),
          panel.grid = element_blank())
  
  # output
  return(plt)
}

# Function 6: Visualize UMAP plot after data process
trajUMAP.plot <- function(data_path, 
                          out_path, 
                          start_clust = NULL,
                          out_figure = FALSE
){
  
  ## data format
  base_res <- trajBasedata.process(data_path = data_path,
                                   location_type = "umap",
                                   start_clust = start_clust)
  traj_df <- base_res[["traj_df"]]
  clust_mode <- base_res[["clust_mode"]]
  
  ## traj plot 
  traj_plt <- trajUMAP.plot(datt = traj_df,
                            pointsize = 1)
  
  ## output figure
  if (out_figure == TRUE){
    
    traj_ht <- 7
    traj_wt <- 8
    plt_name <- paste0(out_path, "/traj_result/",clust_mode,"/TrajUMAP_", 
                       clust_mode, "_start", start_clust,".tiff")
    ggsave(filename = plt_name, 
           plot = traj_plt, 
           height = traj_ht, width = traj_wt, units = "in", dpi = 30)
    system(paste0("gzip -f ", plt_name))
  }
  return(traj_plt)
}

# Function 7: Visualize heatmap plot
heatmap.plot <- function(datt,
                         cluster_info,
                         cell_reorder = FALSE,
                         show_genes = NULL,
                         cols){
  
  datt <- as.matrix(datt)
  
  # row annotation
  if (is.null(show_genes)) {
    show_genes <- rownames(datt)
  }
  
  gene_pos <- match(show_genes,rownames(datt))
  row_anno <- rowAnnotation(gene=anno_mark(at = gene_pos,
                                           labels = show_genes))
  # top annotation
  cluster_info <- as.factor(cluster_info)
  names(cols) <- levels(cluster_info)
  cols <- cols[levels(cluster_info)]

  top_anno <- HeatmapAnnotation(Celltype = cluster_info,
                                col = list(Celltype = cols),
                                annotation_legend_param = list(
                                  Celltype = list(
                                    title = "Celltype"),
                                  just = c("right", "bottom")
                                ),
                                show_legend = T,
                                show_annotation_name = F,
                                annotation_height = 0.1)


  # plot heatmap
  col_fun = colorRamp2(c(-4,-2, 0, 1, 2), viridis(5))
  
  if (cell_reorder == TRUE) {
    column_split <- NULL
  } else {
    column_split <- cluster_info
  }
  
  plt <- Heatmap(datt,
                 na_col = "red",
                 cluster_rows = T,
                 cluster_columns = F,
                 show_column_names = F,
                 show_row_names = F,
                 column_split = column_split,
                 top_annotation = top_anno, 
                 column_title = NULL,
                 right_annotation = row_anno,
                 heatmap_legend_param = list(
                   title = "Expression", 
                   at = c(-2, 0, 2)
                 ),
                 col = col_fun,
                 use_raster = T)
  return(plt)  
}

# Function 8: extract top genes from DE result
get.TopTrajgene <- function(data_path,
                            start_clust,
                            n_top = 50){
  ##
  call_file <- read.table(file = paste0(data_path, "/traj_call_file.txt"))[, 1]
  traj_result_file <- call_file[1]
  load(traj_result_file) # traj_result
  pseudo_name <- paste0("start", start_clust)
  
  if (!pseudo_name %in% names(traj_result$pseudotime)) {
    stop("Target trajectory not in the result!")
  }
  
  
  ## top genes
  top_gene <- traj_result$ATres[[pseudo_name]] %>%
    top_n(n = n_top, wt = -pvalue) %>%
    top_n(n = n_top, wt = abs(meanLogFC))
  
  traj_gene <- rownames(top_gene) 
  message(paste0("We use top ", length(traj_gene), " trajectory genes!"))
  return(traj_gene)
}

# Function 9: Process heatmap data
heatmap.process <- function(data_path, 
                            retain_cell,
                            marker_gene
){
  
  ## form count
  check_file <- paste0(data_path, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  
  # norm data with transpose
  norm_list <- h5data.load(data_filename = spatial_data_filename, 
                           sample_size = sample_size,
                           load_count = T,
                           normalization = T, 
                           load_coord = F)[["count_list"]]
  
  # combine norm data and scale
  norm_mat_f <- lapply(norm_list, function(a){
    a <- t(a) %>% as.data.frame() # new added as transpose is canceled in io
    a[,marker_gene]
  }) %>% Reduce("rbind", .)
  
  scale_df <-scale(norm_mat_f, center = T, scale = T) %>%
    t() %>% as.data.frame()
  
  if (!is.null(retain_cell)) {
    scale_df <- scale_df[,retain_cell]
  }
  message(paste0("Scaled mat on ", ncol(scale_df), " cells and ", nrow(scale_df), " genes retained!"))
  return(list("datt" = scale_df))
}

# Function 10: Visualize heatmap after data process
trajHeatmap.plot <- function(data_path1, 
                             data_path2,
                             out_path, 
                             start_clust = NULL,
                             n_top = 50,
                             out_figure = FALSE
){
  
  ## inputs
  base_res <- trajBasedata.process(data_path = data_path1,
                                   location_type = NULL,
                                   start_clust = start_clust)
  traj_df <- base_res[["traj_df"]]
  clust_mode <- base_res[["clust_mode"]]
  traj_gene <- get.TopTrajgene(data_path1,
                               start_clust = start_clust,
                               n_top = n_top)
  traj_df_order <- traj_df[order(traj_df$pseudotime, na.last = NA), ]
  color_in <- COLOR_USE[sort(unique(traj_df_order$cluster_label))]
  traj_df_order$cluster_label <- factor(traj_df_order$cluster_label,
                                        levels = sort(unique(traj_df_order$cluster_label)))
  #
  heatmap_res <- heatmap.process(data_path = data_path2, 
                                 retain_cell = traj_df_order$cell,
                                 marker_gene = traj_gene)
  datt <- heatmap_res[[1]] %>% as.data.frame()
  sample_size <- length(unique(traj_df$sample))
  n_ct <- length(unique(traj_df$cluster_label))
  
  ## traj plot 
  heatmap_plt <- heatmap.plot(datt = datt,
                              cluster_info = traj_df_order$cluster_label,
                              cell_reorder = TRUE,
                              show_genes = NULL,
                              cols = color_in)
  
  ## output figure
  if (out_figure == TRUE){
    hmp_ht <- 1.5 + round(length(traj_gene) * 0.18)
    hmp_wt <- 4 + round(n_ct * 1)
    tiff(paste0(out_path, "/traj_result/", clust_mode, "/Heatmap_", 
                clust_mode, "_start", start_clust, ".tiff"),
         height = hmp_ht, width = hmp_wt, units = "in", res = 300, compression = "lzw")
    print(heatmap_plt)
    dev.off()
    system(paste0("gzip -f ", out_path, "/traj_result/",clust_mode,"/Heatmap_", 
                  clust_mode, "_start", start_clust,".tiff"))
  } 
  
  return(heatmap_plt)
}

# Function 11: Process scatter plot data
trajScatter.process <- function(data_path1,
                                data_path2,
                                n_top,
                                start_clust
){
  
  ## inputs
  base_res <- trajBasedata.process(data_path = data_path1,
                                   location_type = NULL,
                                   start_clust = start_clust)
  traj_df <- base_res[["traj_df"]]
  clust_mode <- base_res[["clust_mode"]]
  
  traj_gene <- get.TopTrajgene(data_path1,
                               start_clust = start_clust,
                               n_top = n_top)
  
  traj_dff <- traj_df[!is.na(traj_df$pseudotime), ]
  color_in <- COLOR_USE[sort(unique(traj_dff$cluster_label))]
  traj_dff$cluster_label <- factor(traj_dff$cluster_label,
                                   levels = sort(unique(traj_dff$cluster_label)))
  
  ## form count
  check_file <- paste0(data_path2, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  
  # norm data with transpose
  norm_list <- h5data.load(data_filename = spatial_data_filename, 
                           sample_size = sample_size,
                           load_count = T,
                           normalization = T, 
                           load_coord = F)[["count_list"]]
  
  # combine norm data and scale
  norm_mat_f <- lapply(norm_list, function(a){
    a <- t(a) %>% as.data.frame() # new added as transpose is canceled in io
    a[,traj_gene]
  }) %>% Reduce("rbind", .)
  
  norm_mat_f <- norm_mat_f[traj_dff$cell,]
  
  message(paste0("We normalized ", nrow(norm_mat_f), " cells and ", 
                 ncol(norm_mat_f), " genes retained!"))
  
  # format input
  scatter_mat <- cbind(traj_dff, norm_mat_f)
  scatter_df <- reshape2::melt(scatter_mat, 
                               id.vars = c("sample", "cell", "cluster_label", "pseudotime"), 
                               variable.name = "Gene",
                               value.name = "Expression")
  scatter_df$Gene <- as.character(scatter_df$Gene )
  
  return(list("scatter_df" = scatter_df,
              "color_in" = color_in,
              "clust_mode" = clust_mode))
}

# Function 12: plot scatter plot
trajScatter.visulize <- function(datt,
                         color_in,
                         pointsize = 1
){
  
  # plot  
  datt$pseudotime <- datt$pseudotime + 0.01
  plt <- ggplot(datt, aes(x = pseudotime, y = Expression, color = cluster_label)) +
    geom_point(alpha = 1,size = pointsize) +
    scale_colour_manual("Cluster", values = color_in)+
    theme_bw()+
    theme(legend.position = "right",
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13),
          strip.placement = "outside", 
          strip.background = element_blank(),
          strip.text.x = element_text(size = 12,color = "black"),
          panel.grid = element_blank())+
    facet_wrap(~datt$Gene, ncol = 2)
  
  # output
  return(plt)
}

# Function 13: Visualize scatter after data process
trajScatter.plot <- function(data_path1,
                             data_path2,
                             out_path, 
                             start_clust = NULL,
                             n_top = 4,
                             out_figure = FALSE
){
  
  ## inputs
  scatter_res <- trajScatter.process(data_path1 = data_path1, 
                                     data_path2 = data_path2, 
                                     n_top = n_top,
                                     start_clust = start_clust)
  scatter_df <- scatter_res[["scatter_df"]]
  color_in <- scatter_res[["color_in"]]
  clust_mode <- scatter_res[["clust_mode"]]
  n_gene <- length(unique(scatter_df$Gene))
  
  ## traj plot 
  scatter_plt <- trajScatter.visulize(datt = scatter_df,
                                      pointsize = 1,
                                      color_in = color_in)
  
  ## output figure
  if (out_figure == TRUE){
    sct_ht <- ceiling(n_gene/2) * 3
    sct_wt <- min(n_gene, 2) * 4
    tiff(paste0(out_path, "/traj_result/",clust_mode,
                "/Scatter_", clust_mode, "_start", start_clust,".tiff"),
         height = sct_ht, width = sct_wt, units = "in", res = 300, compression = "lzw")
    print(scatter_plt)
    dev.off()
    system(paste0("gzip -f ", out_path, "/traj_result/",clust_mode,
                  "/Scatter_", clust_mode, "_start", start_clust,".tiff"))
  } 
  
  return(scatter_plt)
}

# Function 14
traj_plt.plot <- function(data_path1,                ## String: output path of traj procedure
                          data_path2,                ## String: output path of qc procedure
                          out_path,                  ## String: output path of qc procedure
                          start_clust,
                          hm_gene_num = 10,
                          sc_gene_num = 4, 
                          clus_mode,
                          out_figures
                          
){
  
  source(paste0(method_path, "/io.R"))
  traj_plt <- traj.plot(data_path = data_path1,
                        out_path = out_path,
                        start_clust = start_clust,
                        gridnum = 10,
                        out_figure = out_figures)
  trajHeatmap <- trajHeatmap.plot(data_path1 = data_path1,
                                  data_path2 = data_path2,
                                  out_path = out_path,
                                  start_clust = start_clust,
                                  n_top = hm_gene_num,
                                  out_figure = out_figures)
  trajScatter <- trajScatter.plot(data_path1 = data_path1,
                                  data_path2 = data_path2,
                                  out_path = out_path,
                                  start_clust = start_clust,
                                  n_top = sc_gene_num,
                                  out_figure = out_figures)
  if(clus_mode == "dr_cl"){
    
    trajUMAP_plt <- trajUMAP.plot(data_path = data_path1,
                                  out_path = out_path,
                                  start_clust = start_clust,
                                  out_figure = out_figures)
    save(traj_plt, trajUMAP_plt, trajHeatmap, trajScatter,
         file = paste0(out_path, "/traj_result/plot.RData"))
  } else {
    
    save(traj_plt, trajHeatmap, trajScatter,
         file = paste0(out_path, "/traj_result/plot.RData"))
  }
  
  return(0)
}


# ###################
# ### test code
# data_path1 <- "/net/mulan/disk2/yasheng/stwebProject/test/traj"
# data_path2 <- "/net/mulan/disk2/yasheng/stwebProject/test/qc"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/test/traj_plt"
# traj_plt.plot(data_path1 = data_path1,
#               data_path2 = data_path2,
#               out_path = output_path,
#               start_clust = 1,
#               hm_gene_num = 10,
#               sc_gene_num = 4,
#               clus_mode = "dr_cl_wr",
#               out_figures = T)

