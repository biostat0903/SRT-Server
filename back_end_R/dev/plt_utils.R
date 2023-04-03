#! /usr/bin/env Rscript
# visualize feature plot and bubble plot after all dr_cl analysis
# up-stream code: cl.R, sdd.R

# Set method_path
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# Load packages
library(Seurat)
library(dplyr)
library(bigreadr)
library(ggplot2)
library(hdf5r)

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

# Function 1: Visualize feature plot
feature.visualize <- function(datt,
                              pointsize = 0.2
){
  
  sample_size <- datt$sample %>% unique %>% length
  plt <- ggplot(datt, aes(x = loc_x, y = loc_y, color = Expression)) + 
    geom_point(alpha = 1, size = pointsize) + 
    scale_colour_gradient(low = "lightgrey", high = "#FF4500")+
    theme_bw()
  
  if (sample_size == 1) {
    
    plt <- plt + facet_wrap(~Gene, ncol = 4)
  } else {
    
    plt <- plt + facet_wrap(~sample + Gene, ncol = 4, dir = "v")
  }
  
  return(plt)
}

# Function 2: Process feature plot data
feature.process <- function(data_path, 
                            marker_gene
){
  
  ## Load count and coordinate data 
  check_file <- paste0(data_path, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  ### count data
  h5_dat <- h5data.load(data_filename = spatial_data_filename,
                        sample_size = sample_size,
                        load_count = TRUE,
                        normalization = TRUE, 
                        load_coord = TRUE, 
                        coordinate = TRUE)
  norm_list <- h5_dat[["count_list"]]
  coord_list <- h5_dat[["coord_list"]]
  
  ## build list
  datt <- plyr::alply(c(1: sample_size), 1, function(a){
    
    norm_dat <- norm_list[[a]] %>%
      t() %>% as.data.frame() # new added as transpose is canceled in io
    # choose marker genes
    norm_dat <- norm_dat[, colnames(norm_dat) %in% marker_gene[[a]]]
    # add annotation
    norm_dat$cell <- rownames(norm_dat)
    norm_dat$loc_x <- coord_list[[a]][, 1]
    norm_dat$loc_y <- coord_list[[a]][, 2]
    norm_dat$sample <- paste0("Sample", a)
    plot_dat <- reshape2::melt(norm_dat, 
                               id.vars = c("cell", "loc_x", "loc_y", "sample"), 
                               variable.name = "Gene",
                               value.name = "Expression")
    return(plot_dat)
  }) %>% Reduce("rbind", .)
  return(datt)
}

# Function 3: Visualize feature plot after data process
feature.plot <- function(data_path1,                   ## String: output path of dr_cl procedure 
                         data_path2,                   ## String: output path of qc procedure
                         mode_usage = NULL,
                         marker_gene_list = NULL, 
                         marker_num = NULL,
                         out_path, 
                         out_figure = FALSE, 
                         zip_figure = FALSE
){
  
  ## load io code
  source(paste0(method_path, "/io.R"))
  
  ## Get featureID of st data
  check_file <- paste0(data_path2, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  featureID_data_h5 <- H5File$new(spatial_data_filename, mode = "r")
  sample_size <- qc_param[2] %>% as.numeric
  featureID <- plyr::alply(c(1: sample_size), 1, function(a){
    featureID_grp_name <- paste0("featureID/featureID.s", a)
    featureID <- featureID_data_h5[[featureID_grp_name]][] 
    return(featureID)
  }) %>% Reduce("rbind", .) %>% unique
  
  ## Get maker gene
  if (is.null(marker_gene_list)){
    
    if(!is.null(marker_num)){
      
      marker_num_in <- marker_num
    } else {
      
      marker_num_in <- 10
      warning("When \"marker_num\" is null, we select top ten HVGs.")
    }
    count_list <- h5data.load(data_filename = spatial_data_filename, 
                              sample_size = sample_size,
                              load_count = TRUE,
                              normalization = FALSE, 
                              load_coord = FALSE)[["count_list"]]

    marker_gene_sel <- plyr::llply(count_list, function(a){
      
      seurat_obj_s <- CreateSeuratObject(a)
      seurat_obj_s <- FindVariableFeatures(seurat_obj_s, 
                                           selection.method = "vst", 
                                           nfeatures = marker_num_in, 
                                           verbose = FALSE)
      hvg_s <- seurat_obj_s@assays$RNA@var.features[1: marker_num_in]
      return(hvg_s)
    })
  } else {
    
    marker_gene <- strsplit(marker_gene_list, ",")[[1]]
    marker_gene_inter <- marker_gene[marker_gene %in% featureID]
    if (length(marker_gene_inter) != 0){
      
      if (length(marker_gene_inter) != length(marker_gene)){
        
        marker_gene_without <- setdiff(marker_gene, marker_gene_inter)
        cat(marker_gene_without, " is not included in ST data!\n")
      } 
    } else {
      
      stop("All marker genes are not inlcuded in ST data!")
    }
    marker_gene_sel <- plyr::alply(c(1: sample_size), 1, function(a) marker_gene_inter)
  }
  
  ## process feature plot data
  ft_datt <- feature.process(data_path = data_path2,
                             marker_gene = marker_gene_sel)
  ft_plt <- feature.visualize(ft_datt)
  
  ## output figure
  if (out_figure == TRUE){
    
    result_dir <- paste0(out_path, "/", mode_usage, "_result")
    if (!file.exists(result_dir)){
      system(paste0("mkdir ", result_dir))
    }
    
    ft_plt <- ft_plt  +
      theme(legend.position = "bottom",
            legend.title = element_text(size = 6), 
            legend.text = element_text(size = 6),
            panel.background = element_rect(fill = "white", colour = "black", size = 1),
            panel.spacing.x = unit(0, "cm"), 
            panel.spacing.y = unit(0, "cm"), 
            strip.placement = "outside", 
            strip.background = element_blank(),
            strip.text.x = element_text(size = 6, color = "black", angle = 0),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank()) 
    
    # ft_ht <- ceiling(sample_size/8)*1.5 + 0.5
    ft_ht <- floor(length(marker_gene_sel[[1]]) / 4)*sample_size*1.5
    ft_wt <- min(4, length(marker_gene_sel[[1]]))
    ggsave(filename = paste0(result_dir, "/Feature_plt.tiff"), 
           plot = ft_plt,
           height = ft_ht, width = ft_wt, units = "in", dpi = 300)
    if(zip_figure == TRUE){
      
      system(paste0("gzip -f ", result_dir, "/Feature_plt.tiff"))
    }
  }
  
  return(ft_plt)
}

# Function 4: Visualize bubble plot
bubble.visualize <- function(datt
){
  
  sample_size <- unique(datt$sample) %>% length
  clus_num <- unique(datt$cluster_setting) %>% length
  plt <- ggplot(data = datt, 
                aes(x = cluster_label, y = gene)) + 
    geom_point(aes(size = pct_exp, 
                   color = avg_exp_scaled)) + 
    scale_colour_viridis_c() + 
    guides(size = guide_legend(title = "Percentage"),
           color = guide_colorbar(title = "Averge Exp."), 
           direction = "vertical") + 
    theme_bw()
  
  if (sample_size == 1) {
    
    plt <- plt + facet_wrap(~ cluster_setting, ncol = clus_num, 
                            dir = "v", labeller = label_parsed)
  } else {
    
    plt <- plt + facet_wrap(~ sample + cluster_setting, ncol = clus_num, 
                            dir = "v", labeller = label_parsed)
  }
  
  return(plt)
  
}

# Function 5: Process bubble plot data
bubble.process <- function(data_path, 
                           clus, 
                           marker_gene,
                           do_scale
){
  
  check_file <- paste0(data_path, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  spatial_data_h5 <- H5File$new(spatial_data_filename, mode = "r")
  sample_size <- qc_param[2] %>% as.numeric
  count_list <- h5data.load(data_filename = spatial_data_filename, 
                            sample_size = sample_size,
                            load_count = TRUE,
                            normalization = FALSE, 
                            load_coord = FALSE)[["count_list"]]
  
  clus_list <- split(clus, clus$sample)
  # format average expression and percent of expressed cell 
  datt <- lapply(c(1: sample_size), function(s){
    clus_s <- clus_list[[s]]
    count <- count_list[[s]][marker_gene[[s]],]
    datt_s <- lapply(unique(clus_s[, 3]), function(a) {
      
      count_s <- count[, clus_s[, 3] == a, drop = F]
      avg_exp <- rowMeans(count_s)
      if (do_scale == T) {
        
        avg_exp_scaled <- scale(avg_exp)
        avg_exp_scaled[avg_exp_scaled > 2.5] <- 2.5
        avg_exp_scaled[avg_exp_scaled < -2.5] <- -2.5
      } else {
        
        avg_exp_scaled <- log1p(avg_exp)
      }
      
      pct_exp <- apply(X = count_s, 1, function(x){
        return(sum(x > 0) / length(x))
      })
      gene_mat_s <- data.frame(gene = rownames(count_s),
                               avg_exp = avg_exp,
                               avg_exp_scaled = avg_exp_scaled,
                               pct_exp = pct_exp * 100,
                               cluster_label = a)
      return(gene_mat_s)
    }) %>% Reduce("rbind", .)
    datt_s$sample <- paste0("Sample", s)
    return(datt_s)
  }) %>% Reduce("rbind", .)
  
  datt$cluster_label <- datt$cluster_label %>%
    factor(., levels = sort(unique(.)))
  datt$cluster_setting <- sub("cluster_", "", colnames(clus)[3])
  
  return(datt)
}

# Function 6: Visualize bubble plot after data process
bubble.plot <- function(data_path1,                   ## String: output path of dr_cl procedure 
                        data_path2,                   ## String: output path of qc procedure
                        mode_usage = NULL,
                        marker_gene_list = NULL, 
                        marker_num = NULL, 
                        do_scale = TRUE,
                        out_path, 
                        out_figure = FALSE,
                        zip_figure = FALSE
){
  
  ## Load io code
  source(paste0(method_path, "/io.R"))
  
  ## Load cluster 
  cl_post_file <- paste0(data_path1, "/", mode_usage, "_post_file.txt")
  cl_clus_file <- read.table(cl_post_file)[1, 1]
  if (cl_clus_file == "0"){
    
    cl_clus_file <- read.table(cl_post_file)[2, 1]
  }
  cl_clus <- fread2(cl_clus_file)
  setting_num <- ncol(cl_clus) - 2
  setting_str <- sub("cluster_", "", colnames(cl_clus)[3: ncol(cl_clus)])
  
  ## Get featureID of st data
  check_file <- paste0(data_path2, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  featureID_data_h5 <- H5File$new(spatial_data_filename, mode = "r")
  sample_size <- qc_param[2] %>% as.numeric
  featureID <- plyr::alply(c(1: sample_size), 1, function(a){
    featureID_grp_name <- paste0("featureID/featureID.s", a)
    featureID <- featureID_data_h5[[featureID_grp_name]][] 
    return(featureID)
  }) %>% Reduce("rbind", .) %>% unique
  
  ## Get maker gene
  if (is.null(marker_gene_list)){
    
    if(!is.null(marker_num)){
      
      marker_num_in <- marker_num
    } else {
      
      marker_num_in <- 50
      warning("When {marker_num} is null, we select top 50 HVGs.")
    }
    count_list <- h5data.load(data_filename = spatial_data_filename, 
                              sample_size = sample_size,
                              load_count = TRUE,
                              normalization = FALSE, 
                              load_coord = FALSE)[["count_list"]]
    marker_gene_sel <- plyr::llply(count_list, function(a){
      seurat_obj_s <- CreateSeuratObject(a)
      seurat_obj_s <- FindVariableFeatures(seurat_obj_s, 
                                           selection.method = "vst", 
                                           nfeatures = marker_num_in, 
                                           verbose = FALSE)
      hvg_s <- seurat_obj_s@assays$RNA@var.features[1: marker_num_in]
      return(hvg_s)
    })
  } else {
    
    marker_gene <- strsplit(marker_gene_list, ",")[[1]]
    marker_gene_inter <- marker_gene[marker_gene %in% featureID]
    if (length(marker_gene_inter) != 0){
      
      if (length(marker_gene_inter) != length(marker_gene)){
        
        marker_gene_without <- setdiff(marker_gene, marker_gene_inter)
        marker_gene_sel <- plyr::alply(c(1: sample_size), 1, function(a){
          marker_gene_without
        })
        messae(paste0(marker_gene_without, " is not included in ST data!"))
      } else {
        
        marker_gene_sel <- plyr::alply(c(1: sample_size), 1, function(a){
          marker_gene_inter
        })
      } 
    } else {
      
      stop("Marker genes are not inlcuded in ST data!")
    }
  }
  
  ## process feature plot data
  bb_datt <- plyr::alply(c(1: setting_num), 1, function(a){
    clus_s <- cl_clus[, c(1, 2, 2+a)]
    bb_datt_s <- bubble.process(data_path = data_path2, 
                                clus = clus_s, 
                                marker_gene = marker_gene_sel,
                                do_scale = T)
  }) %>% Reduce("rbind", .)
  bb_plt <- bubble.visualize(bb_datt)
  
  ## output figure
  if (out_figure == TRUE){
    
    result_dir <- paste0(out_path, "/", mode_usage, "_result")
    if (!file.exists(result_dir)){
      
      system(paste0("mkdir ", result_dir))
    }
    
    gene_usage <- length(marker_gene_sel[[1]])
    bb_plt <- bb_plt  +
      scale_radius(range = c(0, 1)) +
      theme(axis.title = element_blank(),
            panel.grid = element_blank(),
            axis.text = element_text(size = 4),
            legend.position = "top",
            legend.title = element_text(size = 6), 
            legend.text = element_text(size = 4), 
            legend.direction = "horizontal", 
            legend.box = "vertical")
    
    bb_ht <- min(4, gene_usage/2) + 2
    bb_wt <- min(4, sample_size*4) + 0.5
    ggsave(filename = paste0(result_dir, "/Bubble_plt.tiff"), 
           plot = bb_plt,
           height = bb_ht, width = bb_wt, units = "in", dpi = 300)
    if(zip_figure == TRUE){
      
      system(paste0("gzip -f ", result_dir, "/Bubble_plt.tiff"))  
    }
    
  }
  
  return(bb_plt)
  
}

# Function 7: Visualize location plot
loc.visualize <- function(datt, 
                          pointsize, 
                          vis_type,
                          title_in, 
                          color_in
){
  
  sample_size <- length(unique(datt$sample))
  plt <- ggplot(datt, aes(x = loc_x, y = loc_y, color = cluster)) + 
    geom_point(alpha = 1, size = pointsize) + 
    ggtitle(title_in) + 
    theme_bw() + 
    theme(legend.position = "bottom",
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 5),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) 
  if (max(nchar(as.character(datt$cluster))) > 10){
    
    plt <- plt + guides(color = guide_legend(byrow = TRUE, 
                                             ncol = 2, 
                                             keyheight = 0.4))
  } else {
    
    plt <- plt + guides(color = guide_legend(byrow = TRUE, 
                                             ncol = 6, 
                                             keyheight = 0.4))
  }
  if (sample_size != 1){
    
    plt <- plt + facet_wrap(~ sample, ncol = 2)
  }
  if (vis_type == "cell_type"){
    
    plt <- plt + scale_color_manual("Types", values = color_in)
  } else {
    
    plt <- plt + scale_color_manual("Domains", values = color_in)
  }
  
  return(plt)
}

# Function 8: Process location plot data 
loc.process <- function(data_path 
){
  
  ## load file
  check_file <- paste0(data_path, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  spatial_data_h5 <- H5File$new(spatial_data_filename, mode = "r")
  sample_size <- qc_param[2] %>% as.numeric
  
  coord_df <- plyr::alply(c(1: sample_size), 1, function(a){
    coord_grp_name <- paste0("meta.data/meta.data.s", a)
    coord_mat_a <- spatial_data_h5[[coord_grp_name]][, ] %>% 
      as.data.frame()
    
    cellID_grp_name <- paste0("cellID/cellID.s", a)
    rownames(coord_mat_a) <- spatial_data_h5[[cellID_grp_name]][] 
    colnames(coord_mat_a) <- c("loc_x", "loc_y") 
    coord_mat_a$sample <- paste0("Sample", a)
    
    return(coord_mat_a)
  }) %>% Reduce("rbind", .)
  
  return(coord_df)
}

# Function 9: Visualize location plot data 
loc.plot <- function(data_path1,                   ## String: output path of ct procedure 
                     data_path2,                   ## String: output path of qc procedure
                     mode_usage = NULL,            ## String: "cl", "sdd"
                     vis_type = "cell_type",       ## String: "cell_type", "spatial_domain"
                     out_path,                     ## String: output path of cl/sdd procedure 
                     out_figure = FALSE,           ## Boolean: output figure
                     zip_figure = FALSE
){
  
  ## load io code
  source(paste0(method_path, "/io.R"))
  ## Load cluster 
  cl_post_file <- paste0(data_path1, "/", mode_usage, "_post_file.txt")
  cl_clus_file <- read.table(cl_post_file)[, 1]
  # if(mode_usage == "CT_PCA"){
  #   
  #   if (vis_type != "cell_type"){
  #     
  #     vis_type <- "cell_type"
  #     stop("The setting of \"vis_type\" should be \"cell_type\". We coerce to set it as \"cell_type\".\n")
  #   }
  # }
  if(vis_type == "cell_type"){
    
    cl_clus <- fread2(cl_clus_file[1])
  } else {
    
    cl_clus <- fread2(cl_clus_file[2])
  }
  clus_name <- colnames(cl_clus)[-c(1, 2)]
  setting_num <- ncol(cl_clus) - 2
  setting_str <- sub("cluster_", "", colnames(cl_clus)[3: ncol(cl_clus)])
  
  ## process data
  datt <- loc.process(data_path2)
  sample_size <- length(unique(cl_clus$sample))
  ## location plot
  p_size <- ifelse(nrow(datt)/sample_size > 10000, 0.03, 0.3)
  ct_loc_list <- lapply(clus_name, function(cl){
    datt$cluster <- cl_clus[, cl] %>% 
      factor(., levels = sort(unique(.)))
    clus_num <- length(unique(datt$cluster))
    if (vis_type == "cell_type"){
      
      ct_loc_plt <- loc.visualize(datt = datt, 
                                  title_in = "", 
                                  vis_type = vis_type,
                                  pointsize = p_size,
                                  color_in = CT_COLS) 
      loc_name <- paste0(out_path, "/", mode_usage, 
                         "_result/Location_ct_", cl, ".tiff")
    } else {
      
      ct_loc_plt <- loc.visualize(datt = datt, 
                                  title_in = "", 
                                  vis_type = vis_type,
                                  pointsize = p_size,
                                  color_in = CL_COLS)    
      loc_name <- paste0(out_path, "/", mode_usage, 
                         "_result/Location_sd_", cl, ".tiff")
    }
    if(out_figure == TRUE){
      
      plt_ht <- min(ceiling(sample_size/2)*2, 6)+1
      plt_wt <- max(plt_ht+1, sample_size*floor(sqrt(clus_num)/1.5))
      ggsave(filename = loc_name, plot = ct_loc_plt,
             height = plt_ht, width = plt_wt, units = "in",
             dpi = 300)
      if(zip_figure == TRUE){
        
        system(paste0("gzip -f ", loc_name))  
      }
    }
  })
  return(ct_loc_list)
}
