#! /usr/bin/env Rscript
# load packages 
# up-stream procedure: ccc.R

library(bigreadr)
library(dplyr)
library(CellChat)

#

# Function 1: Process plot data
cccdata.process <- function(ccc_mat = NULL, 
                            magin_thresh = 0,
                            target_file = NULL
){
  
  ##format input
  ccc_mat <- ccc_mat[,c("source", "target",
                        "ligand", "receptor",
                        "magnitude", "significance")]
  ccc_mat_sig <- subset(ccc_mat, ccc_mat$significance < 0.05 &
                          ccc_mat$magnitude > magin_thresh)
  ccc_mat_sig$significance <- with(ccc_mat_sig, 
                                   ifelse(significance > 0.01, 1,
                                          ifelse(significance > 0.001, 2,
                                                 ifelse(significance >0.0001, 3, 4)
                                          )))
  values <- c(1, 2, 3, 4)
  names(values) <- c("<0.05", "<0.01", "<0.001", "<0.0001")
  
  ##target pairs
  if (nrow(ccc_mat_sig) == 0) {
    
    stop("WARNING: Not significant communication iterms left!")
  } else {
    
    all_source <- unique(ccc_mat_sig$source) %>% sort()
    all_target <- unique(ccc_mat_sig$target) %>% sort()
    if (is.null(target_file)) {
      
      targ_pairs_f <- data.frame(source = rep(all_source, each = length(all_target)), 
                                 target = rep(all_target, times = length(all_source))) 
    } else {
      targ_pairs <- fread2(target_file)
      colnames(targ_pairs) <- c("source", "target")
      targ_pairs_f <- subset(targ_pairs, targ_pairs$source %in% all_source &
                               targ_pairs$target %in% all_target)
    }
  }
  return(list("ccc_mat_sig" = ccc_mat_sig,
              "values" = values,
              "targ_pairs" = targ_pairs_f))
}

# Function 2: circle plot
pathway.circle.visualize <- function(ccc_mat = NULL,
                                     pathway = NULL,
                                     out_path = NULL){
  
  all_ct <- unique(c(ccc_mat$source, ccc_mat$target))
  group_size <- length(all_ct)
  net_mat <- matrix(0, ncol = group_size, nrow = group_size,
                    dimnames = list(all_ct, all_ct))
  ## with rows correspond to sources, and columns correspond to targets
  for (i in seq_along(all_ct)) {
    for (j in seq_along(all_ct)) {
      net_mat[i,j] <- subset(ccc_mat, 
                             ccc_mat$source == all_ct[i] & 
                               ccc_mat$target == all_ct[j]) %>% nrow()
    }
  }
  ## output circle plot 1
  circle_plt_path <- paste0(out_path, "/circle_for_", pathway, ".tiff")
  tiff(circle_plt_path,
       height = 8, width = 8, units = "in",
       res = 300, compression = "lzw")
  netVisual_circle(as.matrix(net_mat), 
                   vertex.weight = group_size,
                   weight.scale = T, 
                   label.edge= T,
                   arrow.size = 0.5,
                   title.name = paste0("Number of interactions for ", pathway, " signaling pathway"))
  dev.off()
  return(circle_plt_path)
}


# Function 3: Visualize feature plot after data process
pathway.circle.plot <- function(ccc_mat, 
                                out_path,
                                pathway,
                                ccc_method
){
  ## inputs
  
  magin_thresh <- ifelse(ccc_method == "cellchat",
                         0.02, 0.2)
  circle_res <- cccdata.process(ccc_mat = ccc_mat, 
                                magin_thresh = magin_thresh)
  ccc_mat_sig <- circle_res[[1]]
  
  ## circle plot 
  circle_plt_path <- pathway.circle.visualize(ccc_mat = ccc_mat_sig,
                                              pathway = pathway,
                                              out_path = out_path)
  return(circle_plt_path)
}

# Function 4: Process feature plot data
pathway.cccdot.process <- function(ccc_mat, 
                                   n_show = 50
){
  
  ## Load
  dot_df <- with(ccc_mat,
                 data.frame(LR = paste0(ligand, ":", receptor),
                            CC = paste0(source, "->", target),
                            significance = significance,
                            magnitude = magnitude))
  
  # show top LR pairs
  if (length(unique(dot_df$LR)) > n_show) {
    top_LR <- table(dot_df$LR) %>% 
      sort(.,decreasing = T) %>% 
      head(.,n_show) %>%
      names()
    dot_df <- subset(dot_df, dot_df$LR %in% top_LR)
  }
  
  return(list("dot_df" = dot_df))
}

# Function 5: plot dotplot
pathway.cccdot.visualize <- function(dot_df = NULL,
                                     values = NULL
){
  # plot
  plt <- ggplot(dot_df, aes(CC,LR))+ 
    geom_point(aes(color=magnitude, size=significance)) +
    labs(x= "", y = "")+
    scale_radius(range = c(min(dot_df$significance), max(dot_df$significance)), 
                 breaks = sort(unique(dot_df$significance)),
                 labels = names(values)[values %in% sort(unique(dot_df$significance))], 
                 name = "Significance")+
    scale_colour_viridis_c(name = "magnitude")
  
  return(plt)
}

# Function 6: Visualize dotplot after data process
pathway.cccdot.plot <- function(ccc_mat, 
                                out_path, 
                                ccc_method,
                                pathway
){
  
  ## inputs
  magin_thresh <- ifelse(ccc_method == "cellchat",
                         0.02, 0.2)
  source(paste0(method_path, "/ccc_plt.R"))
  dot_res <- cccdata.process(ccc_mat = ccc_mat,
                             magin_thresh = magin_thresh)
  ccc_mat_sig <- dot_res[[1]]
  values <- dot_res[[2]]
  dot_df <- cccdot.process(ccc_mat = ccc_mat_sig, 
                           n_show = 50)[[1]]
  dot_plt <- cccdot.visualize(dot_df = dot_df,
                              values = values)
  # save path
  dot_ht <- length(unique(dot_df[,1]))*0.2 + 3
  dot_wt <- length(unique(dot_df[,2]))*0.5 + 3
  dot_plt <- dot_plt + 
    theme_bw()+
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.8),
          axis.title = element_text(size = 15, color = "black"),
          axis.text = element_text(size = 12, color = "black"))
  #
  dot_plt_path <- paste0(out_path, "/dot_for_", pathway, ".tiff")
  ggsave(filename = dot_plt_path, 
         plot = dot_plt,
         height = dot_ht, width = dot_wt, units = "in", dpi = 300)
  
  return(dot_plt_path)
}

# Function 7
pathway_plt.plot <- function(data_path,                         ## String: path of ccc procedure
                             out_path                           ## String: path of ccc_plt procedure
){
  
  ## with pathway_mat_file, ccc_method and pathway_id pre-wirtten in data_path
  post_file <- read.table(paste0(data_path, "/ccc_post_file.txt"),
                                 sep = "\t", header = F)[, 1]
  pathway_mat_file <- ifelse(is.na(post_file[2]), post_file[4], post_file[2])
  ccc_method <- read.table(paste0(out_path, "/ccc_method.txt"),
                           sep = "\t", header = F)[1, 1]
  ccc_mat_file <- gsub("_path_mat.txt$", "_lr_mat.txt",  pathway_mat_file)
  ccc_mat <- fread2(ccc_mat_file)
  #
  pathway_id <- read.table(paste0(out_path, "/sel_pathway.txt"),
                           sep = "\t", header = F)[1, 1]
  ccc_mat <- subset(ccc_mat, ccc_mat$pathway == pathway_id)
  
  ### circle plot
  circle_plt <- pathway.circle.plot(ccc_mat = ccc_mat, 
                                    out_path = out_path, 
                                    pathway = pathway_id,
                                    ccc_method = ccc_method)
  ### dot plot
  dot_plt <- pathway.cccdot.plot(ccc_mat = ccc_mat, 
                                 out_path = out_path, 
                                 pathway = pathway_id,
                                 ccc_method = ccc_method)
  
  
  return(0)  
}

# data_path="/public/home/biostat03/project/stwebProject/02_data/Case/case_study1_out/"
# pathway.ccc_plt.plot(data_path = data_path,
#                      out_path = data_path)