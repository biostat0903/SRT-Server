# codePath setting
method_path <- "/public/home/biostat03/project/stwebProject/01_code/srt_server/dev"

# load packages 
library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyr)
library(reshape2)
library(scales) 

# Function 1: Visualize location plot with gene expression
locgene.visualize <- function(datt,
                              pointsize
){
  
  # plot
  sample_size <- unique(datt$sample) %>% length
  pal <- colorRampPalette(c("antiquewhite", viridis_pal()(10)[c(6,2)]))
  plt <- ggplot(datt, aes(x = x, y = y, color = Expression)) + 
    geom_point(alpha = 1, size = pointsize) + 
    scale_color_gradientn(colours = pal(5)) + 
    theme_bw() +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white",
                                          colour = "black", 
                                          size = 1),
          panel.spacing.x = unit(0, "cm"), 
          strip.placement = "outside", 
          strip.background = element_blank(),
          strip.text.x = element_text(size = 4, color = "black", angle = 0),
          strip.text.y = element_text(size = 4, color = "black"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
  if (sample_size != 1) {
    
    plt <- plt + facet_wrap(~sample, ncol = 2)
  }
  return(plt)
}

# Function 2: Process plot data
locgene.process <- function(data_path,                       ## String: output path for qc procedure
                            gene_name
){
  
  ## format norm count data
  source(paste0(method_path, "/io.R"))
  check_file <- paste0(data_path, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  spatial_data <- h5data.load(spatial_data_filename,
                              sample_size,
                              load_count = TRUE,
                              normalization = FALSE,
                              load_coord = TRUE, 
                              coordinate = TRUE)
  count_norm_list <- plyr::llply(spatial_data[["count_list"]], function(a){
    NormalizeData(t(a)) %>% as.matrix()
  })
  location_list <- spatial_data[["coord_list"]]
  
  # extract gene expression and location of SVGs
  dat_list <- lapply(1:sample_size, function(a){

    featureID <- colnames(count_norm_list[[a]])
    count_a <- count_norm_list[[a]][, featureID %in% gene_name] %>%
      as.data.frame()
    location_a <- location_list[[a]]
    count_a$cell <- rownames(count_a)
    count_a <- cbind(count_a, location_a)
    dat_a <- reshape2::melt(count_a, 
                            id.vars = c("cell", "x", "y"), 
                            variable.name = "Gene",
                            value.name = "Expression")
    dat_a$sample <- paste0("Sample", a)
    dat_a$Gene <- gene_name

    return(dat_a)
  })
  datt <- Reduce(rbind, dat_list)
  
  return(datt)
}

# Function 3: Plot
loc_gene.plot <- function(data_path,                          ## String: output path for qc procedure
                          out_path                            ## String: output path for loc_gene procedure
){
  
  ## inputs
  gene_id <- read.table(paste0(out_path, "/loc_gene.txt"))[1, 1]
  exp_datt <- locgene.process(data_path = data_path, gene_name = gene_id)
  sample_size <- length(unique(exp_datt$sample))
  ## plot 
  lgpt_plt <- locgene.visualize(datt = exp_datt, pointsize = 0.1)
  ## output figure
  plt_ht <- max(2*sample_size, 2)
  plt_wt <- max(4*sample_size, 4)
  ggsave(filename = paste0(out_path, "/loc_gene_plot.tiff"), 
         plot = lgpt_plt, height = plt_ht, width = plt_wt, units = "in", dpi = 300)
  
  return(0)
}
data_path="/public/home/biostat03/project/stwebProject/02_data/Case/case_study3_outx"
loc_gene.plot(data_path = data_path,
              out_path  = data_path)

