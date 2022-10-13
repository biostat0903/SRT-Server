#! /usr/bin/env Rscript
# load packages 
# up-stream procedure: ccc.R

library(bigreadr)
library(dplyr)
library(CellChat)

# Function 1: Process plot data
cccdata.process <- function(data_path, 
                           target_file = NULL
){
  
  ## 1. load
  post_file <- read.table(paste0(data_path, "/ccc_post_file.txt"), 
                          header = F, sep = "\t")[, 1]
  ccc_result <- post_file[1]
  clus_mode <- post_file[4]
  ccc_mat <- fread2(ccc_result)
  
  ## 2. format input
  ccc_mat <- ccc_mat[,c("source", "target",
                        "ligand", "receptor",
                        "lr_co_ratio", "lr_co_ratio_pvalue")]
  colnames(ccc_mat)[5:6] <- c("maginitude", "significance")
  ccc_mat_sig <- subset(ccc_mat, ccc_mat$significance < 0.05 &
                          ccc_mat$maginitude > 0.2)
  ccc_mat_sig$significance <- with(ccc_mat_sig, 
                                   ifelse(significance > 0.01, 1,
                                          ifelse(significance > 0.001, 2,
                                                 ifelse(significance >0.0001, 3, 4)
                                          )))
  values <- c(1, 2, 3, 4)
  names(values) <- c("<0.05", "<0.01", "<0.001", "<0.0001")
  
  ## 3. target pairs
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
              "targ_pairs" = targ_pairs_f,
              "clus_mode" = clus_mode))
}

# Function 2: circle plot
circle.visualize <- function(ccc_mat = NULL,
                             targ_pairs = NULL,
                             out_path = NULL){
  #
  all_ct <- unique(c(ccc_mat$source, ccc_mat$target))
  group_size <- length(all_ct)
  net_mat <- matrix(0, ncol = group_size, nrow = group_size,
                    dimnames = list(all_ct, all_ct))
  
  # with rows correspond to sources, and columns correspond to targets
  for (i in seq_along(all_ct)) {
    for (j in seq_along(all_ct)) {
      net_mat[i,j] <- subset(ccc_mat, 
                             ccc_mat$source == all_ct[i] & 
                               ccc_mat$target == all_ct[j]) %>% nrow()
    }
  }
  
  # output circle plot 1
  tiff(paste0(out_path, "/circle_all.tiff"),
       height = 8, width = 8, units = "in",
       res = 300, compression = "lzw")
  netVisual_circle(as.matrix(net_mat), 
                   vertex.weight = group_size,
                   weight.scale = T, 
                   label.edge= F,
                   arrow.size = 0.5,
                   title.name = "Number of interactions")
  dev.off()
  system(paste0("gzip -f ", out_path, "/circle_all.tiff"))
  
  # output circle plot 2
  n_row <- ceiling(nrow(net_mat)/5)
  circle_sub_wt <- min(nrow(net_mat), 5) * 3
  circle_sub_ht <- n_row * 3
  
  tiff(paste0(out_path, "/circle_by_celltypes.tiff"),
       height = circle_sub_ht, width = circle_sub_wt, units = "in",
       res = 300, compression = "lzw")
  par(mfrow = c(n_row, 5))
  for (i in sort(rownames(net_mat))) {
    targ_pairs_sub <- subset(targ_pairs, targ_pairs[,1] == i)
    netVisual_circle(net_mat, 
                     vertex.weight = group_size, 
                     targets.use = as.character(targ_pairs_sub$target),
                     sources.use = as.character(targ_pairs_sub$source),
                     weight.scale = T, 
                     label.edge= T,
                     arrow.size = 0.5,
                     edge.weight.max = max(net_mat), 
                     title.name = paste0("source ", i))
  }
  dev.off()
  system(paste0("gzip -f ", out_path, "/circle_by_celltypes.tiff"))
  return(c(paste0("Circle plot saved in: ", out_path, "circle_all.tiff"),
           paste0("Circle subplot saved in: ", out_path, "circle_by_celltypes.tiff")))
}


# Function 3: Visualize feature plot after data process
circle.plot <- function(data_path, 
                        out_path, 
                        target_file
){
  
  ## inputs
  circle_res <- cccdata.process(data_path = data_path, 
                                target_file = target_file)
  ccc_mat_sig <- circle_res[[1]]
  targ_pairs <- circle_res[[3]]
  clus_mode <- circle_res[[4]]
  
  ## circle plot 
  result_dir <- paste0(out_path, "/ccc_result/", clus_mode, "/")
  if (!file.exists(result_dir)) {
    system(paste0("mkdir -p ", result_dir))
  }
  circle_plt_path <- circle.visualize(ccc_mat = ccc_mat_sig,
                                      targ_pairs = targ_pairs,
                                      out_path = result_dir)
  return(circle_plt_path)
}


# Function 4: Process feature plot data
cccdot.process <- function(ccc_mat, 
                           n_show = 50
){
  
  ## Load
  dot_df <- with(ccc_mat,
                 data.frame(LR = paste0(ligand, ":", receptor),
                            CC = paste0(source, "->", target),
                            significance = significance,
                            maginitude = maginitude))
  
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
cccdot.visualize <- function(dot_df = NULL,
                             values = NULL
){
  # plot
  plt <- ggplot(dot_df, aes(CC,LR))+ 
    geom_point(aes(color=maginitude, size=significance)) +
    labs(x= "", y = "")+
    scale_radius(range = c(min(dot_df$significance), max(dot_df$significance)), 
                 breaks = sort(unique(dot_df$significance)),
                 labels = names(values)[values %in% sort(unique(dot_df$significance))], 
                 name = "Significance")+
    scale_colour_viridis_c(name = "Maginitude")
  
  return(plt)
}

# Function 6: Visualize dotplot after data process
cccdot.plot <- function(data_path, 
                        out_path, 
                        target_file,
                        out_figure = TRUE
){
  
  ## inputs
  circle_res <- cccdata.process(data_path = data_path, 
                                target_file = target_file)
  ccc_mat_sig <- circle_res[[1]]
  values <- circle_res[[2]]
  targ_pairs <- circle_res[[3]]
  clus_mode <- circle_res[[4]]
  
  ## circle plot 
  dot_plt <- list()
  for (i in unique(targ_pairs[,1])) {
    
    ccc_mat_sub <- subset(ccc_mat_sig, 
                          ccc_mat_sig$source == i | 
                          ccc_mat_sig$target == i)
    dot_df_i <- cccdot.process(ccc_mat = ccc_mat_sub, 
                               n_show = 50)[[1]]
    dot_plt_i <- cccdot.visualize(dot_df = dot_df_i,
                                  values = values)
    dot_plt[[i]] <- dot_plt_i
    if (out_figure == T) {
      
      # save path
      dot_ht <- length(unique(dot_df_i[,1]))*0.2 + 3
      dot_wt <- length(unique(dot_df_i[,2]))*0.5 + 3
      dot_plt_i <- dot_plt_i + 
                    theme_bw()+
                    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.8),
                          axis.title = element_text(size = 15, color = "black"),
                          axis.text = element_text(size = 12, color = "black"))
      ggsave(filename = paste0(out_path, "/ccc_result/dotplot_", i, ".tiff"), 
             plot = dot_plt_i,
             height = dot_ht, width = dot_wt, units = "in", dpi = 300)
      system(paste0("gzip -f ", out_path, "/ccc_result/dotplot_", i, ".tiff"))
    } 
    message(paste0("Dotplot saved for ", i))
  }
  
  return(dot_plt)
}

# Function 7
ccc_plt.plot <- function(data_path, 
                         out_path, 
                         target_file = NULL,
                         out_figures
){
  
  result_dir <- paste0(out_path, "/ccc_result")
  if (!file.exists(result_dir)){
    
    system(paste0("mkdir ", result_dir))
  }
  circle_plt <- circle.plot(data_path = data_path, 
                            out_path = out_path, 
                            target_file = target_file)
  dot_plt <- cccdot.plot(data_path = data_path, 
                         out_path = out_path, 
                         target_file = target_file,
                         out_figure = TRUE)
  save(dot_plt,
       file = paste0(out_path, "/ccc_result/plot.RData"))

  return(0)  
}


# ###################
# ### test code
# data_path <- "/net/mulan/disk2/yasheng/stwebProject/test/ccc"
# output_path <- "/net/mulan/disk2/yasheng/stwebProject/test/ccc_plt"
# ccc_plt.plot(data_path= data_path,
#              out_path= output_path,
#              target_file = NULL,
#              out_figures = T)

