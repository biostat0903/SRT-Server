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

# Set colors
CL_COLS <- c("#FD7446", "#709AE1", "#31A354", "#9EDAE5", "#DE9ED6",
             "#BCBD22", "#CE6DBD", "#DADAEB", "#FFFF00", "#FF9896",
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


# Function 1: Process location plot data 
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

# Function 2: Visualize location plot data with selected cluster
loc_clus_sel.plot <- function(data_path1,                   ## String: output path of sdd/ct procedure 
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
  cl_selection <- read.csv(paste0(out_path, "/clus_cluster_selection.csv"), 
                           header = F)[, 1]
  cl_selection_label <- paste0(cl_selection, collapse = "&")
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
    loc_name <- paste0(out_path, "/Location_sd_", cl_sel, "_ct", cl_selection_label, ".tiff")
    cols <- CL_COLS
  } else {
    
    cl_clus <- list.files(cl_dir, recursive = T)[grep("cluster_label", list.files(cl_dir, recursive = T))] %>%
      paste0(cl_dir, "/", .) %>%
      fread2()
    cl_sel <- read.table(paste0(out_path, "/clus_plt_sel.txt"))[1, 1] %>%
      gsub(".tiff", "", .) %>% 
      gsub("Location_ct_", "", .) 
    loc_name <-  paste0(out_path, "/Location_ct_", cl_sel, "_ct", cl_selection_label, ".tiff")
    cols <- CT_COLS
  }
  datt$cluster <- cl_clus[, colnames(cl_clus) %in% cl_sel] %>% 
    factor(., levels = sort(unique(.)))
  clus_num <- length(unique(datt$cluster))
  datt$cluster <- as.character(datt$cluster)
  datt$cluster[!datt$cluster %in% cl_selection] <- NA
  
  # plot  
  datt_na <- subset(datt, is.na(datt$cluster))
  datt_a <- subset(datt, !is.na(datt$cluster))
  
  ## plot
  p_size <- ifelse(nrow(datt)/sample_size > 10000, 0.03, 0.3)
  vis_type <- ifelse(file.exists(paste0(data_path1, "/sdd_result")), 
                     "/spatial_domain", "cell_type")
  plt <- ggplot(datt, aes(x = loc_x, y = loc_y, color = cluster)) +
    geom_point(data = datt_na, color = "grey", size = p_size, show.legend = F)+
    geom_point(data = datt_a, alpha = 1,size = p_size, show.legend = T) +
    scale_color_manual(values = cols) +
    theme_void()+
    theme(plot.title = element_text(size = 15,  face = "bold"),
          text = element_text(size = 15),
          legend.position = "bottom",
          legend.text = element_text(size = 8)) 
  if (sample_size > 1){
    
    plt <-  plt + facet_wrap(~sample_size, ncol = 2)
  }

  ## output
  plt_ht <- min(ceiling(sample_size/2)*2, 6)+1
  plt_wt <- max(plt_ht+1, sample_size*floor(sqrt(clus_num)/1.5))
  ggsave(filename = loc_name, plot = plt,
         height = plt_ht, width = plt_wt, units = "in",
         dpi = 300)
  return(0)
}


loc_clus_sel.plot(data_path2,data_path2,data_path2)
