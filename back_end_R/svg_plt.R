#! /usr/bin/env Rscript
# Visualize spatially variable gene (SVG) for ST data
# We use pattern plot, qq plot and histogram plot to visualize the result of svg
# up-stream procedure: svg.R

# codePath setting
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# load packages 
library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyr)
library(reshape2)
library(scales) 

# Function 1: Visualize pattern plot
svgpattern.visualize <- function(datt,
                                 pointsize
){
  
  # plot
  sample_size <- unique(datt$sample) %>% length
  datt$label_facet <- paste0(datt$Gene, "~(", round(datt$p_num, 2), 
                             "%*%", 10, "^", datt$p_exp, ")")
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
  if (sample_size == 1) {
    
    plt <- plt + facet_wrap(~label_facet, ncol = 4,
                            labeller = label_parsed)
  } else {
    
    plt <- plt + facet_wrap(~sample + label_facet, ncol = 4, 
                            labeller = label_parsed)
  }
  return(plt)
}

# Function 2: Process pattern plot data
svgpattern.process <- function(data_path1,                       ## String: output path for qc procedure
                               data_path2,                       ## String: output path for svg procedure
                               n_svg = NULL, 
                               gene_name = NULL
){
  
  ## load file
  call_file <- read.table(file = paste0(data_path2, "/svg_call_file.txt"))[, 1]
  svg_result_file <- call_file[1]
  svg_method <- call_file[2]
  load(svg_result_file)
  sample_size <- length(tot_svg_result)
  
  ## format norm count data
  check_file <- paste0(data_path1, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  spatial_data <- h5data.load(spatial_data_filename,
                              sample_size ,
                              load_count = TRUE,
                              normalization = FALSE,
                              load_coord = TRUE, 
                              coordinate = TRUE)
  count_norm_list <- plyr::llply(spatial_data[["count_list"]], function(a){
    NormalizeData(t(a)) %>% as.matrix()
  })
  location_list <- spatial_data[["coord_list"]]
  
  ## extract selected svgs and their p values
  if(!is.null(n_svg)){
    
    svg_list <- lapply(1: sample_size, function(a){
      svg_result_s <- tot_svg_result[[a]] %>% 
        top_n(n = n_svg, wt = -adjPval)
      return(svg_result_s$Gene[1: n_svg])
    })
  } else {
    
    svg_list <- lapply(1: sample_size, function(a){
      svg_result_s <- tot_svg_result[[a]]$Gene
      return(svg_result_s[svg_result_s %in% gene_name])
    })
  }
  svg_pval_list <- lapply(1: sample_size, function(a){
    svg_result_s <- tot_svg_result[[a]]
    return(svg_result_s$Pval[svg_result_s$Gene %in% svg_list[[a]]])
  })
  
  # extract gene expression and location of SVGs
  dat_list <- lapply(1:sample_size, function(a){
    
    show_genes_a <- svg_list[[a]]
    featureID <- colnames(count_norm_list[[a]])
    count_a <- count_norm_list[[a]][, featureID %in% show_genes_a] %>% 
      as.data.frame()
    location_a <- location_list[[a]]
    count_a$cell <- rownames(count_a)
    count_a <- cbind(count_a, location_a)
    
    dat_a <- reshape2::melt(count_a, 
                            id.vars = c("cell", "x", "y"), 
                            variable.name = "Gene",
                            value.name = "Expression")
    dat_a$sample <- paste0("Sample", a)
    
    p_val_a <- svg_pval_list[[a]] %>% 
      format(scientific = T) %>%
      as.character %>% 
      strsplit('e') %>% 
      unlist %>% 
      as.numeric
    p_val <- data.frame(Gene = show_genes_a,
                        p_num = p_val_a[seq(1, length(p_val_a), 2)],
                        p_exp = round(p_val_a[-seq(1, length(p_val_a), 2)]))
    rownames(p_val) <- p_val[,"Gene"]
    dat_a <- merge(dat_a, p_val, by="Gene", all.x = TRUE)
    
    return(dat_a)
  })
  
  datt <- Reduce(rbind, dat_list)
  return(datt)
}

# Function 3: Visualize pattern plot after data process
svgpattern.plot <- function(data_path1,                          ## String: output path for qc procedure
                            data_path2,                          ## String: output path for svg procedure
                            out_path, 
                            svg_num = NULL,
                            gene_name = NULL, 
                            out_figure = FALSE,
                            zip_figure = FALSE
){
  
  ## inputs
  if (is.null(svg_num)){
    
    if(is.null(gene_name)){
      
      stop("You should input eithter \"svg_num\" or \"gene_name\"!")
    } else {
      
      message("We show the genes: ", gene_name)
      gene_str <- strsplit(gene_name, ",")[[1]]
      pattern_datt <- svgpattern.process(data_path1 = data_path1, 
                                         data_path2 = data_path2,
                                         gene_name = gene_str)
    }
  } else {
    
    if(!is.null(gene_name)){
      
      stop("You should input eithter \"svg_num\" or \"gene_name\"!")
    } else {
      pattern_datt <- svgpattern.process(data_path1 = data_path1, 
                                         data_path2 = data_path2, 
                                         n_svg = svg_num)
    }
  }
  sample_size <- length(unique(pattern_datt$sample))
  
  ## pattern plot 
  svgpt_plt <- svgpattern.visualize(datt = pattern_datt, 
                                    pointsize = 0.1)
  
  ## output figure
  if (out_figure == TRUE){
    
    svgpt_ht <- max(2*sample_size, 2)
    svgpt_wt <- max(4*sample_size, 4)
    ggsave(filename = paste0(out_path, "/svg_result/Svg_pattern_plot.tiff"), 
           plot = svgpt_plt, 
           height = svgpt_ht, width = svgpt_wt, units = "in", dpi = 300)
    if(zip_figure == TRUE){
      
      system(paste0("gzip -f ", out_path, 
                    "/svg_result/Svg_pattern_plot.tiff"))
    }
  }
  
  return(svgpt_plt)
}

# Function 4: Visualize pattern plot
pattern.plot <- function(data_path, 
                         out_path, 
                         pointsize = 0.01, 
                         out_figure = FALSE,
                         zip_figure = FALSE
){
  
  ## inputs
  call_file <- read.table(file = paste0(data_path, "/svg_call_file.txt"))[, 1]
  svg_result_file <- call_file[1]
  svg_method <- call_file[2]
  load(svg_result_file)
  
  # plot
  sample_size <- length(tot_pattern_result)
  pattern_num <- length(unique(tot_pattern_result[[1]]$patternType))
  tot_pattern_datt <- plyr::alply(1: sample_size, 1, function(a){
    tot_pattern_result_s <- tot_pattern_result[[a]]
    tot_pattern_result_s$sample <- paste0("Sample", a)
    return(tot_pattern_result_s)
  }) %>% Reduce("rbind", .)
  pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
  pt_plt <- ggplot(tot_pattern_datt, 
                   aes(x = x, y = y, color = patternValue)) + 
    geom_point(alpha = 1, size = pointsize) + 
    scale_color_gradientn(colours = pal(5)) + 
    theme_bw() +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white",
                                          colour = "black", 
                                          size = 1),
          panel.spacing.x = unit(0, "cm"), 
          # strip.placement = "outside", 
          # strip.background = element_blank(),
          # strip.text.x = element_text(size = 4, color = "black", angle = 0),
          # strip.text.y = element_text(size = 4, color = "black"),
          strip.text.y = element_blank(),
          strip.text.x = element_blank(),
          # axis.title = element_blank(),
          axis.title.x = element_text(size = 8),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
  if (sample_size == 1) {
    
    pt_plt <- pt_plt + facet_wrap(~patternType, ncol = 3) + 
      theme(axis.title.y = element_blank())+
      xlab("Patterns")
  } else {
    
    pt_plt <- pt_plt + facet_wrap(~sample + patternType, ncol = 3)+ 
      theme(axis.title.y = element_text(size = 6))+
      xlab("Patterns") +
      ylab("Samples")
  }
  ## output figure
  if (out_figure == TRUE){
    
    pt_ht <- max(2*pattern_num*sample_size, 2)
    pt_wt <- max(4*pattern_num*sample_size, 4)
    ggsave(filename = paste0(out_path, "/svg_result/Pattern_plot.tiff"), 
           plot = pt_plt, limitsize = FALSE,
           height = pt_ht, width = pt_wt, units = "in", dpi = 300)
    if(zip_figure == TRUE){
      
      system(paste0("gzip -f ", out_path, "/svg_result/Pattern_plot.tiff"))
    }
  }
  
  return(pt_plt)
}

# Function 5: Visualize qq plot
## qplot_gg (copy from https://github.com/xzhoulab/SPARK-X-Analysis/blob/main/funcs/utilities.R)
qq.visualize <- function(datt
){
  
  ## parameters
  col.base <- c("hotpink", "mediumorchid2", "lightskyblue")
  ylab <- expression(paste("Observed ", -log[10], 
                           "(", italic(p), "-value)"))
  xlab <- expression(paste("Expected ", -log[10], 
                           "(", italic(p), "-value)"))
  
  ## format plot data
  n <- nrow(datt)
  datt <- datt[order(datt$P),]
  datt_s <- data.frame(
    expected = -log10(ppoints(n)),
    observed = -log10(datt$P),
    clower = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:n, shape2 = n:1)),
    sample = datt$sample
  )
  sample_size <- unique(datt$sample) %>% length
  
  ## plot
  if (sample_size == 1){
    
    plt <- ggplot(datt_s) + 
      geom_ribbon(aes(x = expected, ymax = cupper, ymin = clower),
                  fill = "darkgrey", alpha = 0.5)+
      geom_abline(intercept = 0, slope = 1, 
                  size = 0.8, linetype = 2, col = "blue")+
      geom_point(aes(expected, observed), 
                 size = 0.2, col = "#E9C46A", alpha = 1)+
      scale_x_continuous(xlab) +
      scale_y_continuous(ylab) +
      theme_bw()+
      theme(axis.text = element_text(size = 6), 
            axis.title = element_text(size = 6, face = "bold"), 
            legend.text = element_text(size=rel(1)),
            legend.title = element_text(size=rel(1)),
            panel.grid = element_blank(),
            strip.placement = "outside", 
            strip.background = element_blank(),
            strip.text.x = element_text(size = 4, face = "bold",
                                        color = "black", angle = 0),
            axis.line = element_line(colour = "black"), 
            legend.position = "none")
  }else{
    plt_list <- lapply(1: sample_size, function(a){
      
      plt_a <- ggplot(datt_s[datt_s$sample == paste0("Sample", a),]) + 
        geom_ribbon(aes(x = expected, ymax = cupper, ymin = clower),
                    fill = "darkgrey", alpha = 0.5)+
        geom_abline(intercept = 0, slope = 1, 
                    size = 0.8, linetype = 2, col = "blue")+
        geom_point(aes(expected, observed), 
                   size = 0.2, col = "#E9C46A", alpha = 1)+
        scale_x_continuous(xlab) +
        scale_y_continuous(ylab) +
        theme_bw()+
        theme(axis.text = element_text(size = 6), 
              axis.title = element_text(size = 6, face = "bold"), 
              plot.title = element_text(size = 8, face = "bold"),
              legend.text = element_text(size=rel(1)),
              legend.title = element_text(size=rel(1)),
              panel.grid = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(size = 4, face = "bold",
                                          color = "black", angle = 0),
              axis.line = element_line(colour = "black"), 
              legend.position = "none")+
        labs(title = paste0("Sample", a))
      return(plt_a)
    })
    plt <- patchwork::wrap_plots(plt_list, 
                                 nrow = floor(sqrt(sample_size)))
  }
  
  return(plt)
}

# Function 6: Process qq plot data
qq.process <- function(data_path
){
  
  # load file
  call_file <- read.table(file = paste0(data_path, "/svg_call_file.txt"))[, 1]
  svg_result_file <- call_file[1]
  load(svg_result_file)
  sample_size <- length(tot_svg_result_perm)
  
  # extract p values
  p_val_list <- lapply(1: sample_size, function(a){
    p_a <- tot_svg_result_perm[[a]] %>% 
      as.data.frame()
    dat_a <- p_a[, c(1, 2)]
    dat_a$sample <- paste0("Sample", a)
    dat_a[, 2] <- as.numeric(dat_a[, 2])
    return(dat_a)
  })
  
  datt <- Reduce(rbind, p_val_list) %>% data.frame
  names(datt)[2] <- 'P'
  
  return(datt)
}

# Function 7: Visualize qq plot after data process
qq.plot <- function(data_path, 
                    out_path,
                    out_figure = FALSE,
                    zip_figure = FALSE
){
  
  ## inputs
  p_val_data <- qq.process(data_path = data_path)
  sample_size <- length(unique(p_val_data$sample))
  
  ## pattern plot on top svg
  qq_plt <- qq.visualize(datt = p_val_data)
  if(out_figure== TRUE){
    
    qq_ht <- ifelse(sample_size == 1, 2, 1+2*floor(sqrt(sample_size)))
    qq_wt <- max(2.5, floor(sample_size/floor(sqrt(sample_size))) * 2.5)
    ggsave(filename = paste0(out_path, "/svg_result/qq_plot.tiff"),
           plot = qq_plt,
           height = qq_ht, width = qq_wt, units = "in", dpi = 300)
    if(zip_figure == TRUE){
      
      system(paste0("gzip -f ", out_path, "/svg_result/qq_plot.tiff"))
    }
  }
  
  return(qq_plt)
}

# Function 8: Visualize histogram
hist.visualize <- function (datt,
                            cutoff
){
  
  sample_size <- length(unique(datt$sample))
  datt_s <- datt[datt$P < cutoff,]
  adj_cutoff <- 0.01 / nrow(datt)
  x_lab <- expression(paste(italic(p), "-value"))
  
  if (sample_size == 1) {
    
    # Histogram
    plt <- ggplot(datt_s, aes(x = P, y=..scaled..)) +
      geom_histogram(aes(y = (..density..)/sum(..density..)), 
                     color="darkblue", fill="lightblue", bins = 20) +
      geom_density(color = "#000000", fill = "white", alpha = 0.6) +
      xlab(x_lab) + 
      ylab("Percentage") + 
      geom_vline(xintercept = adj_cutoff, linetype = "dashed", color = "red") + 
      scale_color_brewer(palette = "Paired") + 
      theme_classic() + 
      theme(legend.position = "top",
            axis.text = element_text(size = 6),
            axis.title = element_text(size = 6))
    n_s <- nrow(datt_s)
    p_s <- paste0('(', round(n_s/nrow(datt)*100, 2), '%)')
    n_as <- nrow(datt[datt$P < adj_cutoff,])
    p_as <- paste0('(', round(n_as/nrow(datt)*100, 2), '%)')
    plt <- plt + 
      annotate("text", x = 0.65*cutoff, y = 0.8, size = 1.5,
               label = bquote(Number~of~Genes[P < .(cutoff)]:.(n_s)~.(p_s))) +
      annotate("text", x = 0.65*cutoff, y = 0.7, size = 1.5,
               label = bquote(Number~of~Genes[P-adj < .(cutoff)]:.(n_as)~.(p_as)))
  } else {
    
    plt_list <- lapply(1: sample_size, function(a){
      # Histogram
      plt_a <- ggplot(datt_s[datt_s$sample == paste0("Sample", a),], aes(x = P, y=..scaled..)) +
        geom_histogram(aes(y = (..density..)/sum(..density..)), 
                       color="darkblue", fill="lightblue", bins = 20) +
        geom_density(color = "#000000", fill = "white", alpha = 0.6) +
        xlab(x_lab) + 
        ylab("Percentage") + 
        geom_vline(xintercept = adj_cutoff, linetype = "dashed", color = "red") + 
        scale_color_brewer(palette = "Paired") + 
        theme_classic() + 
        theme(legend.position = "top",
              axis.text = element_text(size = 6),
              plot.title = element_text(size = 6, face = "bold"))+
        ggtitle(paste0("Sample", a))
      n_a <- nrow(datt[datt$sample == paste0("Sample", a),])
      n_s <- nrow(datt_s[datt_s$sample == paste0("Sample", a),])
      p_s <- paste0('(', round(n_s/n_a*100, 2), '%)')
      n_as <- nrow(datt[datt$sample == paste0("Sample", a) & datt$P < adj_cutoff,])
      p_as <- paste0('(', round(n_as/n_a*100, 2), '%)')
      plt_a <- plt_a + 
        annotate("text", x = 0.65*cutoff, y = 0.8, size = 1.5,
                 label = bquote(Number~of~Genes[P < .(cutoff)]:.(n_s)~.(p_s))) +
        annotate("text", x = 0.65*cutoff, y = 0.7, size = 1.5,
                 label = bquote(Number~of~Genes[P-adj < .(cutoff)]:.(n_as)~.(p_as)))
      
      return(plt_a)
    })
    plt <- patchwork::wrap_plots(plt_list, 
                                 nrow = floor(sqrt(sample_size)))
  }
  
  return(plt)
}

# Function 9: Process histogram data
hist.process <- function(data_path
){
  
  # load file
  call_file <- read.table(file = paste0(data_path, "/svg_call_file.txt"))[, 1]
  svg_result_file <- call_file[1]
  load(svg_result_file)
  sample_size <- length(tot_svg_result)
  
  # extract p values
  p_val_list <- lapply(1: sample_size, function(a){
    dat_a <- tot_svg_result[[a]] %>% 
      as.data.frame()
    dat_a$sample <- paste0("Sample", a)
    dat_a[, 2] <- as.numeric(dat_a[, 2])
    return(dat_a[, c(2, 4)])
  })
  
  datt <- Reduce(rbind, p_val_list) %>% data.frame
  names(datt)[1] <- 'P'
  
  return(datt)
}

# Function 10: Visualize histogram after data process
hist.plot <- function(cutoff = 0.01,
                      data_path, 
                      out_path, 
                      out_figure = FALSE,
                      zip_figure = FALSE
){
  
  ## inputs
  p_val_dat <- hist.process(data_path = data_path)
  sample_size <- length(unique(p_val_dat$sample))
  
  ## histogram
  hist_plt <- hist.visualize(datt = p_val_dat,
                             cutoff = cutoff)
  if(out_figure == TRUE){
    
    hist_ht <- ifelse(sample_size == 1, 2, 1+2*floor(sqrt(sample_size)))
    hist_wt <- max(3, floor(sample_size/floor(sqrt(sample_size))) * 3)
    ggsave(filename = paste0(out_path, "/svg_result/Histogram.tiff"),
           plot = hist_plt,
           height = hist_ht, width = hist_wt, units = "in", dpi = 300)
    if(zip_figure == TRUE){
      
      system(paste0("gzip -f ", out_path, "/svg_result/Histogram.tiff"))
    }
  }
  return(hist_plt)
}

# Function 11: total plot
svg_plt.plot <- function(data_path1,                       ## String: output path for qc procedure
                         data_path2,                       ## String: output path for svg procedure
                         out_path, 
                         svg_num, 
                         qqplot = FALSE,                   ## Boolean: qqplot
                         patternplot = FALSE,              ## Boolean: pattern plot
                         out_figures,
                         zip_figures
){
  
  ## load io code
  source(paste0(method_path, "/io.R"))
  check_file <- read.table(paste0(data_path2, "/svg_check_file.txt"))[, 1]
  result_dir <- paste0(out_path, "/svg_result")
  if (!file.exists(result_dir)){
    
    system(paste0("mkdir ", result_dir))
  }
  
  ## qq plot 
  qq_plt <- NULL
  if(qqplot == TRUE){
    
    if (length(check_file) == 2 | length(check_file) == 5){
      
      qq_plt <- qq.plot(data_path = data_path2,
                        out_path = out_path,
                        out_figure = out_figures, 
                        zip_figure = zip_figures)
      
    }else{
      
      stop("Please perform permuation test first!")
    }
  }
  
  ## pattern plot 
  pattern_plt <- NULL
  if(patternplot == TRUE){
    
    if (length(check_file) == 4 | length(check_file) == 5){
      
      pattern_plt <- pattern.plot(data_path = data_path2,
                                  out_path = out_path,
                                  pointsize = 0.2,
                                  out_figure = out_figures, 
                                  zip_figure = zip_figures)
    }else{
      
      stop("Please perform pattern estimation first!")
    }
  }
  
  ## plot
  svgpattern_plt <- svgpattern.plot(data_path1 = data_path1,
                                    data_path2 = data_path2,
                                    out_path = out_path,
                                    svg_num = svg_num,
                                    out_figure = out_figures, 
                                    zip_figure = zip_figures)
  hist_plt <- hist.plot(data_path = data_path2,
                        out_path = out_path,
                        out_figure = out_figures, 
                        zip_figure = zip_figures)
  
  ## save plot data
  save(svgpattern_plt, pattern_plt, qq_plt, hist_plt, 
       file = paste0(out_path, "/svg_result/svg_plot.RData"))
  
  return(0)
}
