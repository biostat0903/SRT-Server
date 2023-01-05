#! /usr/bin/env Rscript
# visualize  and bubble plot after ora
# up-stream code: ora.R

# Load packages
library(ggplot2)
library(clusterProfiler)

# Function 1: 
bubble.process <- function(data_path
){
  
  ## load ora file
  post_file <- read.table(paste0(data_path, "/ora_post_file.txt"))[,1]
  ora_mat <- bigreadr::fread2(post_file[1])
  set.seed(20210826)
  id <- sample(c(1: nrow(ora_mat)))
  ora_datt <- data.frame(IDNum = id, 
                         Log_p = -log10(ora_mat$pvalue), 
                         Category = factor(ora_mat$PathwayInfo, 
                                           levels = unique(ora_mat$PathwayInfo)),
                         Count = ora_mat$Count)
  return(ora_datt)
}

# Function 2: 
bubble.visualize <- function(datt
){
  
  plt <- ggplot(datt, aes(x = IDNum, y = Log_p, color = Category))+
    geom_point(shape = 19,aes(fill = Category, size = Count, shape = 19),alpha=0.8)+
    scale_radius()+
    labs(x = "", y = expression(paste(bold(-log[10]),bold("("),bolditalic(p),bold("-value)"))))+
    theme(axis.text.x=element_blank(),
          plot.title = element_text(lineheight = .8, face="bold"),
          axis.text = element_text(size = 6),
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'grey80'),
          axis.title = element_text(size = 6, face = 'bold'),
          legend.title = element_text(size = 3,face = 'bold'),
          legend.text = element_text(size = 3),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = 'white'))+
    geom_hline(yintercept = 1.82, col = 'black',linetype = 2,size=2)+
    guides(color = guide_legend(order = 1, override.aes = list(alpha = 1,size=3)),
           size = guide_legend(order = 2, override.aes = list(alpha = 1,shape=21)),
           fill = FALSE)+ 
    labs(size = "Gene set size")+
    scale_color_manual(values=c("salmon","gold2","#42d4f4","#3cb44b",
                                "chocolate2","#4363d8","#bfef45","#911eb4","#f032e6","#a9a9a9"))+
    scale_fill_manual(values = c("salmon","gold2","#42d4f4","#3cb44b",
                                 "chocolate2","#4363d8","#bfef45","#911eb4","#f032e6","#a9a9a9"))+
    theme(legend.direction = "horizontal")+
    theme(legend.position = c(0.5, 0.85))+
    theme(legend.box = "horizontal")+
    theme(legend.title.align = 0)
  
  return(plt)
}

# Function 3: Visualize bubble plot after data process
bubble.plt <- function(data_path,                  ## String: output path of ora procedure
                       out_path,                   ## String: output path of ora_plt procedure
                       out_figure, 
                       zip_figure = FALSE
){
  
  bubble_datt <- bubble.process(data_path)
  bubble_plt <- bubble.visualize(bubble_datt)
  
  ## Output
  if(out_figure == TRUE){
    
    ggsave(filename = paste0(out_path, "/ora_result/Bubble_plot.tiff"), 
           plot = bubble_plt, height = 3, width = 6,
           units = "in", dpi = 300)
    if(zip_figure == TRUE){
      
      system(paste0("gzip -f ", out_path, 
                    "/ora_result/Bubble_plot.tiff"))
    }
  }
  return(0)
}

# Function 4: ora_plt.plot
ora_plt.plot <- function(data_path,                ## String: output path of ora procedure
                         out_path,                 ## String: output path of ora_plt procedure
                         out_figures, 
                         zip_figures
){
  
  ## Load file
  check_file <- read.table(paste0(data_path, "/ora_check_file.txt"))[, 1]
  pathway_db <- check_file[3]
  call_file <- read.table(paste0(out_path, "/ora_call_file.txt"))[, 1]
  
  ## plot
  result_dir <- paste0(out_path, "/ora_result")
  if (!file.exists(result_dir)){
    
    system(paste0("mkdir ", result_dir))
  }
  bubble_plt <- dot_plt <- NULL
  if (pathway_db == "ALL"){
    
    bubble_plt <- bubble.plt(data_path = data_path, 
                             out_path = out_path, 
                             out_figure = out_figures, 
                             zip_figure = zip_figures)
  } else {
    
    load(call_file[1])
    dot_plt <- dotplot(ora_res, showCategory = 10)
    cat("ok\n")
    if(out_figures == TRUE){
      
      tiff(filename = paste0(result_dir, "/dotplt_", pathway_db, ".tiff"), 
           width = 4, height = 4, units = "in", res = 300)
      dot_plt
      dev.off()
      # if(zip_figures == TRUE){
      #   
      #   system(paste0("gzip -f ", result_dir, "/dotplt_", pathway_db, ".tiff"))
      # }
    }
  }
  
  ## save plot data
  save(bubble_plt, dot_plt, 
       file = paste0(out_path, "/ora_result/ora_plot.RData"))
  
  return(0)
}
