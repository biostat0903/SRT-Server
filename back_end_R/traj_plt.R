#! /usr/bin/env Rscript
# Visualize trajectory
# up-stream code: traj.R
# down-stream procedure: none

# Set method_path
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# Load packages 
library(hdf5r)
library(ggplot2)
library(bigreadr)
library(dplyr)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(tidyr)
library(reshape2)
library(scales) 
library(ggtree)
library(aplot)
library(patchwork)

# Set colors
CL_COLS <- c("#FD7446", "#709AE1", "#31A354", "#9EDAE5", "#DE9ED6",
             "#BCBD22", "#CE6DBD", "#DADAEB", "yellow", "#FF9896",
             "#91D1C2", "#C7E9C0" , "#6B6ECF", "#7B4173", "#FABEBE", 
             "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000", 
             "#AAFFC3", "#808000", "#FFD8B1", "#000075", "#808080")
CT_COLS <- c("#98C1D9", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51", 
             "#E0FBFC", "#A8DADC", "#3D5A80", "#81B29A", "#E07A5F", 
             "#DBC9D8", "#b388eb", "#A4277C", "#BC93B2", "#0077b6", 
             "#BB3E03", "#FFDDD2", "#F19C79", "#006D77", "#6A3569",
             "#E6194B", "#4363D8", "#FFE119", "#3CB44B", "#F58231", 
             "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", 
             "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000", 
             "#AAFFC3", "#808000", "#FFD8B1", "#000075", "#808080")

# Function 1: Process trajectory base data
trajBasedata.process <- function(traj_result_file,
                                 location_type = NULL
){
  
  ## load traj_result
  traj_result <- fread2(traj_result_file)
  base_df <- traj_result[,c("sample", "cell", "cluster_label")]
  base_df$sample <- paste0("Sample", base_df$sample)
  ## set base column
  pseudo_col <- paste0("slingPseudotime_1")
  if (!pseudo_col %in% colnames(traj_result)) {
    
    stop("Target trajectory not in the result!")
  } else {
    base_df$pseudotime <- traj_result[[pseudo_col]]
  }
  
  loc_col <- NULL
  ## set location column
  if (!is.null(location_type)) {
    # set
    if (location_type == "coord") {
      
      loc_col <- c("x", "y")
      
    } else if (location_type == "umap") {
      
      # if (clust_mode != "Seurat") {
      #   
      #   stop("UMAP can only be plottd for dr_cl!")
      # } else {
      
      loc_col <- c("UMAP_1", "UMAP_2")
      # }
      
    } else {
      
      stop("Unrecognized location index!")
      
    }
    # check
    if (!all(loc_col %in% colnames(traj_result))) {
      
      stop("Target location columns not in the result!")
    } 
  }
  
  traj_df <- cbind(base_df,
                   traj_result[,loc_col])
  
  return(traj_df)
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
  datt$cluster_label <- as.factor(datt$cluster_label)
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
traj.process <- function(traj_result_file,
                         gridnum = 10
){
  
  traj_df <- trajBasedata.process(traj_result_file = traj_result_file,
                                  location_type = "coord")
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
              "arrow_df" = arrow_df))
}

# Function 4: Visualize trajectory plot after data process
traj.plot <- function(traj_file, 
                      out_path, 
                      start_clust = NULL,
                      submodule,
                      jo_model,
                      gridnum = 10,
                      out_figure = FALSE, 
                      zip_figure = FALSE
){
  
  traj_plt <- list()
  for (i in seq_along(start_clust)) {
    
    traj_file_i <- traj_file[i]
    start_clust_i <- start_clust[i]
    namex <- paste0("start", start_clust_i)
    
    traj_res_i <- traj.process(traj_result_file = traj_file_i, 
                               gridnum = gridnum)
    traj_df_i <- traj_res_i[[1]]
    arrow_df_i <- traj_res_i[[2]]
    
    ## traj plot 
    if(grepl("SDD", submodule) | (grepl("jo", submodule)&jo_model == "sdd")){
      
      color_in <- CL_COLS
    }
    if(grepl("CT", submodule) | (grepl("jo", submodule)&jo_model == "ct")){
      
      color_in <- CT_COLS
    }
    traj_plt[[namex]] <- traj.visualize(datt = traj_df_i,
                                        arrow_df = arrow_df_i,
                                        color_in = color_in,
                                        pointsize = 1,
                                        arrowlength = 0.3,
                                        arrowsize = 0.8)
    
    ## output figure
    if (out_figure == TRUE){
      
      # set plot parameters
      sample_size_i <- length(unique(traj_df_i$sample))
      n_cl_i <- length(unique(traj_df_i$cluster_label))
      traj_ht <- sample_size_i*4 + ceiling(n_cl_i/5)*0.2
      traj_wt <- 10
      plt_name <- paste0(out_path, "/TrajPlot_", submodule, "_", namex,".tiff")
      ggsave(filename = plt_name, 
             plot = traj_plt[[namex]], 
             height = traj_ht, width = traj_wt, 
             units = "in", dpi = 300)
      if(zip_figure == TRUE){
        
        system(paste0("gzip -f ", plt_name))
      }
    } 
    
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
trajUMAP.plot <- function(traj_file, 
                          out_path, 
                          submodule,
                          start_clust = NULL,
                          out_figure = FALSE, 
                          zip_figure = FALSE
){
  
  trajUMAP_plt <- list()
  for (i in seq_along(start_clust)) {
    
    traj_file_i <- traj_file[i]
    start_clust_i <- start_clust[i]
    namex <- paste0("start", start_clust_i)
    ## data format
    traj_df_i <- trajBasedata.process(traj_result_file = traj_file_i,
                                      location_type = "umap")
    
    ## traj plot 
    trajUMAP_plt[[namex]] <- trajUMAP.visualize(datt = traj_df_i,
                                                pointsize = 1)
    
    ## output figure
    if (out_figure == TRUE){
      traj_ht <- 7
      traj_wt <- 8
      plt_name <- paste0(out_path, "/TrajUMAP_", submodule, "_", namex,".tiff")
      ggsave(filename = plt_name, 
             plot =   trajUMAP_plt[[namex]], 
             height = traj_ht, width = traj_wt, units = "in", dpi = 300)
      if(zip_figure == TRUE){
        
        system(paste0("gzip -f ", plt_name))
      } 
    }
  }
  return(trajUMAP_plt)
}

# Function 7: Visualize heatmap plot
heatmap.plot <- function(datt,
                         cluster_info,
                         cols
){
  
  # format plot mat
  datt <- as.matrix(datt)
  plot_df <- t(datt) %>% 
    as.data.frame() %>% 
    as.data.frame() %>%
    mutate(cell = row.names(.)) %>% 
    melt(., id.vars = "cell", 
         value.name = "Expression", 
         variable.name = "Gene")
  
  # make cluster on genes
  gene_clust <- hclust(dist(datt))
  phr <- ggtree(gene_clust, 
                layout="rectangular",
                branch.length="none")
  
  # make top annotation
  group_df <- data.frame(cell = colnames(datt),
                         group = cluster_info,
                         p = "group")
  top_anno <-  ggplot(group_df, aes(x = cell, y = p, fill = group))+
    geom_tile() + 
    scale_x_discrete(limits = colnames(datt))+
    scale_y_discrete(limits = "group")+
    scale_fill_manual(values = cols)+
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())+
    labs(fill = "Cluster")
  
  # base plot
  fill_limit <- range(plot_df$Expression) %>% 
    min() %>% abs() %>% max(.,2)
  plot_df$Expression[plot_df$Expression < -fill_limit] <- -fill_limit
  plot_df$Expression[plot_df$Expression > fill_limit] <- fill_limit
  plt_base <- ggplot(plot_df, aes(x = cell, y = Gene, fill = Expression))+
    geom_raster()+
    scale_fill_viridis_c()+
    scale_y_discrete(limits = gene_clust$labels[gene_clust$order],
                     position="right")+
    scale_x_discrete(limits = colnames(datt))+
    theme(axis.title = element_blank(),
          axis.text.x = element_blank()) 
  
  # combine plot to aplot object (could be convert to ggplot via ggplotify::ggplot2)
  plt <- plt_base %>%
    insert_top(top_anno, height = 0.05) %>%
    insert_left(phr, width = 0.1)
  
  return(plt)  
}

# Function 8: extract top genes from DE result
get.TopTrajgene <- function(ATres_result_file,
                            n_top = 50
){
  
  ATres_result <- fread2(ATres_result_file)
  top_gene <- ATres_result %>%
    filter(!is.na(waldStat)) %>%
    top_n(n = n_top, wt = -pvalue) %>%
    top_n(n = n_top, wt = abs(meanLogFC))
  traj_gene <- top_gene$gene 
  message(paste0("We use top ", length(traj_gene), " trajectory genes!"))
  return(traj_gene)
}

# Function 9: Process heatmap data
heatmap.process <- function(norm_list, 
                            retain_cell,
                            marker_gene
){
  
  # combine norm data and scale
  norm_mat_f <- lapply(norm_list, function(a){
    a <- t(a) %>% as.data.frame() # new added as transpose is canceled in io
    a[,marker_gene]
  }) %>% Reduce("rbind", .)
  scale_df <- scale(norm_mat_f, center = T, scale = T) %>%
    t() %>% as.data.frame()
  if (!is.null(retain_cell)) {
    scale_df <- scale_df[, retain_cell]
  }
  message(paste0("Scaled mat on ", ncol(scale_df), " cells and ", nrow(scale_df), " genes retained!"))
  return(scale_df)
}

# Function 10: Visualize heatmap after data process
trajHeatmap.plot <- function(traj_file, 
                             ATres_file,
                             start_clust,
                             norm_list,
                             out_path, 
                             n_top = 50,
                             submodule,
                             jo_model,
                             out_figure = FALSE, 
                             zip_figure = FALSE
){
  
  ## inputs
  heatmap_plt <- list()
  for (i in seq_along(start_clust)) {
    
    traj_file_i <- traj_file[i]
    ATres_file_i <- ATres_file[i]
    start_clust_i <- start_clust[i]
    namex <- paste0("start", start_clust_i)
    
    ## format plot inputs
    traj_df_i <- trajBasedata.process(traj_result_file = traj_file_i,
                                      location_type = NULL)
    
    traj_gene_i <- get.TopTrajgene(ATres_result_file = ATres_file_i,
                                   n_top = n_top)
    
    traj_df_order_i <- traj_df_i[order(traj_df_i$pseudotime, na.last = NA), ]
    traj_df_order_i$cluster_label <- as.factor(traj_df_order_i$cluster_label)
    if(grepl("SDD", submodule) | (grepl("jo", submodule)&jo_model == "sdd")){
      
      color_in_i <- CL_COLS[1:nlevels(traj_df_order_i$cluster_label)]
    }
    if(grepl("CT", submodule) | (grepl("jo", submodule)&jo_model == "ct")){
      
      color_in_i <- CT_COLS[1:nlevels(traj_df_order_i$cluster_label)]
    }
    datt <- heatmap.process(norm_list = norm_list, 
                            retain_cell = traj_df_order_i$cell,
                            marker_gene = traj_gene_i)%>% as.data.frame()
    heatmap_plt[[namex]] <- heatmap.plot(datt = datt,
                                         cluster_info = traj_df_order_i$cluster_label,
                                         cols = color_in_i)
    
    ## output figure
    if (out_figure == TRUE){
      
      n_ct <- length(unique(traj_df_i$cluster_label))
      hmp_ht <- 1.5 + round(length(traj_gene_i) * 0.18)
      hmp_wt <- 4 + round(n_ct * 1)
      tiff(paste0(out_path, "/Heatmap_", submodule, "_", namex, ".tiff"),
           height = hmp_ht, width = hmp_wt, units = "in", res = 300, compression = "lzw")
      print(heatmap_plt[[namex]])
      dev.off()
      if (zip_figure == TRUE){
        
        system(paste0("gzip -f ", out_path, "/Heatmap_", submodule, "_", namex, ".tiff"))
      }
    } 
  }
  return(heatmap_plt)
}

# Function 11: Process scatter plot data
trajScatter.process <- function(traj_result_file,
                                ATres_result_file,
                                norm_list,
                                submodule,
                                jo_model,
                                n_top
){
  
  ## inputs
  traj_df <- trajBasedata.process(traj_result_file = traj_result_file,
                                  location_type = NULL)
  traj_dff <- traj_df[!is.na(traj_df$pseudotime), ]
  traj_dff$cluster_label <- as.factor(traj_dff$cluster_label)
  if(grepl("SDD", submodule) | (grepl("jo", submodule)&jo_model == "sdd")){
    
    color_in <- CL_COLS[1:nlevels(traj_dff$cluster_label)]
  }
  if(grepl("CT", submodule) | (grepl("jo", submodule)&jo_model == "ct")){
    
    color_in <- CT_COLS[1:nlevels(traj_dff$cluster_label)]
  }
  traj_gene <- get.TopTrajgene(ATres_result_file = ATres_result_file,
                               n_top = n_top)
  ## combine norm data and scale
  norm_mat_f <- lapply(norm_list, function(a){
    a <- t(a) %>% as.data.frame() # new added as transpose is canceled in io
    a[, traj_gene]
  }) %>% Reduce("rbind", .)
  norm_mat_f <- norm_mat_f[traj_dff$cell, ]
  message(paste0("We normalized ", nrow(norm_mat_f), " cells and ", 
                 ncol(norm_mat_f), " genes retained!"))
  ## format input
  scatter_mat <- cbind(traj_dff, norm_mat_f)
  scatter_df <- reshape2::melt(scatter_mat, 
                               id.vars = c("sample", "cell", "cluster_label", "pseudotime"), 
                               variable.name = "Gene",
                               value.name = "Expression")
  scatter_df$Gene <- as.character(scatter_df$Gene)
  
  return(list("scatter_df" = scatter_df,
              "color_in" = color_in))
}

# Function 12: plot scatter plot
trajScatter.visulize <- function(datt,
                                 color_in,
                                 submodule, 
                                 jo_model,
                                 pointsize = 1
){
  
  # plot  
  datt$pseudotime <- datt$pseudotime + 0.01
  datt$cluster_label <- as.factor(datt$cluster_label)
  plt <- ggplot(datt, aes(x = pseudotime, y = Expression, color = cluster_label)) +
    geom_point(alpha = 1,size = pointsize) +
    theme_bw()+
    theme(legend.position = "right",
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13),
          strip.placement = "outside", 
          strip.background = element_blank(),
          strip.text.x = element_text(size = 12,color = "black"),
          panel.grid = element_blank())+
    facet_wrap(~datt$Gene, ncol = 2)
  
  if(grep("SDD", submodule) | (grepl("jo", submodule)&jo_model == "ct")) {
    
    plt <- plt + scale_colour_manual("Cell Types", values = color_in)
  } 
  if(grepl("CT", submodule) | (grepl("jo", submodule)&jo_model == "sdd")){
    
    plt <- plt + scale_colour_manual("Domains", values = color_in)
  }
  
  return(plt)
}

# Function 13: Visualize scatter after data process
trajScatter.plot <- function(traj_file,
                             ATres_file,
                             start_clust = NULL,
                             norm_list,
                             out_path, 
                             submodule,
                             jo_model, 
                             n_top = 4,
                             out_figure = FALSE, 
                             zip_figure = FALSE
){
  
  ## inputs
  scatter_plt <- list()
  for (i in seq_along(start_clust)) {
    
    traj_file_i <- traj_file[i]
    ATres_file_i <- ATres_file[i]
    start_clust_i <- start_clust[i]
    namex <- paste0("start", start_clust_i)
    scatter_res_i <- trajScatter.process(traj_result_file = traj_file_i, 
                                         ATres_result_file = ATres_file_i, 
                                         norm_list = norm_list,
                                         submodule = submodule,
                                         jo_model = jo_model,
                                         n_top = n_top)
    
    if (is.null(scatter_res_i) == FALSE){
      
      scatter_df_i <- scatter_res_i[["scatter_df"]]
      color_in_i <- scatter_res_i[["color_in"]]
      ## traj plot 
      scatter_plt[[namex]] <- trajScatter.visulize(datt = scatter_df_i,
                                                   submodule = submodule,
                                                   jo_model = jo_model,
                                                   pointsize = 1,
                                                   color_in = color_in_i)
      ## output figure
      if (out_figure == TRUE){
        
        n_gene_i <- length(unique(scatter_df_i$Gene))
        sct_ht <- ceiling(n_gene_i/2) * 3
        sct_wt <- min(n_gene_i, 2) * 4
        plt_name <- paste0(out_path, "/Scatter_", submodule, "_", namex,".tiff")
        ggsave(filename = plt_name, 
               plot = scatter_plt[[namex]], 
               height = sct_ht, width = sct_wt, 
               units = "in", dpi = 300)
        if(zip_figure == TRUE){
          
          system(paste0("gzip -f ", plt_name))
        }
      }
    }
  }
  return(scatter_plt)
}

# Function 14
traj_plt.plot <- function(data_path1,                ## String: output path of traj procedure
                          data_path2,                ## String: output path of qc procedure
                          out_path,                  ## String: output path of qc procedure
                          start_ct_plot = NULL,
                          start_sdd_plot = NULL,
                          hm_gene_num = 10,
                          sc_gene_num = 4, 
                          out_figures,
                          zip_figures
){
  
  ## Load normalized count list
  source(paste0(method_path, "/io.R"))
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
  ## 2. load and choose traj file
  check_file <- read.table(paste0(data_path1, "/traj_check_file.txt"))[, 1]
  submodule <- check_file[2]
  methods <- check_file[3]
  jo_model <- "NA"
  if (submodule == "CL_jo"){
    
    jo_model <- ifelse(is.na(check_file[8]), "sdd", "ct")  
  }
  post_file <- read.table(paste0(data_path1, "/traj_post_file.txt"))[, 1]
  if(grepl("CT", submodule) | (grepl("jo", submodule)&jo_model == "ct")) {
    
    ct_pseudo_file <- post_file[1] %>% 
      strsplit(",") %>% unlist
    if (is.na(post_file[4])){
      
      ct_ATres_file <- "NA"
    } else {
      
      ct_ATres_file <- post_file[4] %>% strsplit(",") %>% unlist
    }
    start_ct <- post_file[5] %>% 
      strsplit(",") %>% unlist %>% as.numeric()
    start_ct_use <- strsplit(start_ct_plot, ",") %>% 
      unlist %>% as.numeric() %>% intersect(., start_ct)
    
    ## Output
    result_dir <- paste0(out_path, "/traj_result/", submodule, "/")
    if (!file.exists(result_dir)){
      
      system(paste0("mkdir ", out_path, "/traj_result"))
      system(paste0("mkdir ", result_dir))
    }
    if (length(start_ct_use) == 0) {
      
      stop("Start CT cluster specificed for plot not in trajectory result!")
    } else {
      
      ct_pseudo_file_use <- ct_pseudo_file[match(start_ct_use, start_ct)]
      if(all(ct_ATres_file == "NA")){
        
        ct_ATres_file_use <- "NA"
      } else {
        
        ct_ATres_file_use <- ct_ATres_file[match(start_ct_use, start_ct)]
      }
    }
    traj_plt <- traj.plot(traj_file = ct_pseudo_file_use, 
                          out_path = result_dir, 
                          start_clust = start_ct_use,
                          submodule = submodule,
                          jo_model = jo_model,
                          gridnum = 10,
                          out_figure = out_figures, 
                          zip_figure = zip_figures)
    if(all(ct_ATres_file_use =="NA")){
      
      trajScatter <- trajHeatmap <- NULL
      
    } else {
      
      trajScatter <- trajScatter.plot(traj_file = ct_pseudo_file_use,
                                      ATres_file = ct_ATres_file_use,
                                      start_clust = start_ct_use,
                                      norm_list = norm_list,
                                      submodule = submodule,
                                      n_top = sc_gene_num,
                                      out_path = result_dir,
                                      jo_model = jo_model, 
                                      out_figure = out_figures, 
                                      zip_figure = zip_figures)
      trajHeatmap <- trajHeatmap.plot(traj_file = ct_pseudo_file_use,
                                      ATres_file = ct_ATres_file_use,
                                      start_clust = start_ct_use,
                                      norm_list = norm_list,
                                      out_path = result_dir,
                                      submodule = submodule,
                                      jo_model = jo_model, 
                                      n_top = hm_gene_num,
                                      out_figure = out_figures, 
                                      zip_figure = zip_figures)
    }
    if(methods == "Seurat"){
      
      trajUMAP_plt <- trajUMAP.plot(traj_file = ct_pseudo_file_use, 
                                    out_path = result_dir, 
                                    start_clust = start_ct_use,
                                    submodule = submodule,
                                    out_figure = out_figures, 
                                    zip_figure = zip_figures)
    } else {
      
      trajUMAP_plt <- NULL
    }
    save(traj_plt, trajHeatmap, trajScatter, trajUMAP_plt,
         file = paste0(result_dir, "/", submodule,"_plot.RData"))
    
  }
  
  if(grepl("SDD", submodule) | (grepl("jo", submodule)&jo_model == "sdd")) {
    
    sdd_pseudo_file <- post_file[2] %>% 
      strsplit(",") %>% unlist
    if (is.na(post_file[4])){
      
      sdd_ATres_file <- "NA"
    } else {
      
      sdd_ATres_file <- post_file[4] %>% strsplit(",") %>% unlist
    }
    start_sdd <- post_file[6] %>% 
      strsplit(",") %>% unlist %>% as.numeric()
    start_sdd_use <- strsplit(start_sdd_plot, ",") %>% 
      unlist %>% as.numeric() %>% intersect(., start_sdd)
    result_dir <- paste0(out_path, "/traj_result/", submodule, "/")
    if (!file.exists(result_dir)){
      
      system(paste0("mkdir ", out_path, "/traj_result"))
      system(paste0("mkdir ", result_dir))
    }
    if (length(start_sdd_use) == 0) {
      
      stop("Start SDD cluster specificed for plot not in trajectory result!")
    } else {
      
      sdd_pseudo_file_use <- sdd_pseudo_file[match(start_sdd_use, start_sdd)]
      if(all(sdd_ATres_file == "NA")){
        
        sdd_ATres_file_use <- "NA"
      } else {
        
        sdd_ATres_file_use <- sdd_ATres_file[match(start_sdd_use, start_sdd)]
      }
    }
    traj_plt <- traj.plot(traj_file = sdd_pseudo_file_use, 
                          out_path = result_dir, 
                          start_clust = start_sdd_use,
                          submodule = submodule,
                          jo_model = jo_model, 
                          gridnum = 10,
                          out_figure = out_figures, 
                          zip_figure = zip_figures)
    if (all(sdd_ATres_file_use =="NA")){
      
      trajScatter <- trajHeatmap <- NULL
    } else {
      
      trajScatter <- trajScatter.plot(traj_file = sdd_pseudo_file_use,
                                      ATres_file = sdd_ATres_file_use,
                                      start_clust = start_sdd_use,
                                      norm_list = norm_list,
                                      submodule = submodule,
                                      jo_model = jo_model, 
                                      n_top = sc_gene_num,
                                      out_path = result_dir,
                                      out_figure = out_figures, 
                                      zip_figure = zip_figures)
      trajHeatmap <- trajHeatmap.plot(traj_file = sdd_pseudo_file_use,
                                      ATres_file = sdd_ATres_file_use,
                                      start_clust = start_sdd_use,
                                      norm_list = norm_list,
                                      out_path = result_dir,
                                      submodule = submodule,
                                      jo_model = jo_model, 
                                      n_top = hm_gene_num,
                                      out_figure = out_figures, 
                                      zip_figure = zip_figures)
    }
    save(traj_plt, trajHeatmap, trajScatter,
         file = paste0(result_dir, "/", submodule,"_plot.RData"))
  }
  
  return(0)
}
