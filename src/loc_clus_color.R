#! /usr/bin/env Rscript
# visualize feature plot and bubble plot after all dr_cl analysis
# up-stream code: cl.R, sdd.R

# Set method_path
method_path <- "/public/home/biostat03/project/stwebProject/01_code/srt_server/dev"

# Load packages
library(Seurat)
library(dplyr)
library(bigreadr)
library(ggplot2)
library(hdf5r)

# Function 1: Visualize location plot
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

# Function 2: Process location plot data 
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
  
  return(list(coord_df, sample_size))
}

# Function 3: Visualize location plot data 
loc_clus_color.plot <- function(data_path1,                   ## String: output path of sdd/ct procedure 
                                data_path2,                   ## String: output path of qc procedure
                                out_path                      ## String: output path of location plot
){
  
  ## load io code
  source(paste0(method_path, "/io.R"))
  ## process data
  loc_datt <- loc.process(data_path2)
  datt <- loc_datt[[1]]
  sample_size <- loc_datt[[2]]

  ## load color
  cl_col <- read.csv(paste0(out_path, "/clus_col_match_updated.csv"), 
                     header = T)
  ## load cluster 
  cl_dir <- ifelse(file.exists(paste0(data_path1, "/sdd_result")), 
                   paste0(data_path1, "/sdd_result"), 
                   paste0(data_path1, "/ct_result"))

  if(file.exists(paste0(data_path1, "/sdd_result"))){
    
    cl_clus <- list.files(cl_dir, recursive = T)[grep("domain_label", list.files(cl_dir, recursive = T))] %>%
      paste0(cl_dir, "/", .) %>%
      fread2(.)
    cl_sel <- read.table(paste0(out_path, "/clus_plt_sel.txt"))[1, 1]  %>%
      gsub(".tiff", "", .) %>% 
      gsub("Location_sd_", "", .) 
    loc_name <- paste0(out_path, "/Location_sd_", cl_sel, "_updated_color.tiff")
  } else {
    
    cl_clus <- list.files(cl_dir, recursive = T)[grep("cluster_label", list.files(cl_dir, recursive = T))] %>%
      paste0(cl_dir, "/", .) %>%
      fread2()
    cl_sel <- read.table(paste0(out_path, "/clus_plt_sel.txt"))[1, 1] %>%
      gsub(".tiff", "", .) %>% 
      gsub("Location_ct_", "", .) 
    loc_name <-  paste0(out_path, "/Location_ct_", cl_sel, "_updated_color.tiff")
  }
  datt$cluster <- cl_clus[, colnames(cl_clus) %in% cl_sel] %>% 
    factor(., levels = sort(unique(.)))
  clus_num <- length(unique(datt$cluster))

  ## plot
  p_size <- ifelse(nrow(datt)/sample_size > 10000, 0.03, 0.3)
  vis_type <- ifelse(file.exists(paste0(data_path1, "/sdd_result")), 
                     "/spatial_domain", "cell_type")
  ct_loc_plt <- loc.visualize(datt = datt, 
                              title_in = "", 
                              vis_type = vis_type,
                              pointsize = p_size,
                              color_in = cl_col[, 2])
  ## output
  plt_ht <- min(ceiling(sample_size/2)*2, 6)+1
  plt_wt <- max(plt_ht+1, sample_size*floor(sqrt(clus_num)/1.5))
  ggsave(filename = loc_name, plot = ct_loc_plt,
         height = plt_ht, width = plt_wt, units = "in",
         dpi = 300)
  return(0)
}
