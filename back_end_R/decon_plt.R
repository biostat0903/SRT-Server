#! /usr/bin/env Rscript
# Visualize the results of cell type deconvolution
# up-stream procedure: decon.R

# method_path setting
method_path <- "/net/mulan/disk2/yasheng/stwebProject/01_code/01_method"

# load packages
library(hdf5r)
library(dplyr)
library(bigreadr)
library(ggplot2)
library(scatterpie)
library(SingleCellExperiment)

# Set colors
COLOR_USE <- c("#98C1D9", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51", 
               "#E0FBFC", "#A8DADC", "#3D5A80", "#81B29A", "#E07A5F", 
               "#DBC9D8", "#b388eb", "#A4277C", "#BC93B2", "#0077b6", 
               "#BB3E03", "#FFDDD2", "#F19C79", "#006D77", "#6A3569",
               "#E6194B", "#4363D8", "#FFE119", "#3CB44B", "#F58231", 
               "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", 
               "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000", 
               "#AAFFC3", "#808000", "#FFD8B1", "#000075", "#808080")

# Function 1: Process DECON data
base.process <- function(data_path){
  
  # 1. load data
  call_path <- paste0(data_path, "/decon_call_file.txt")
  if (!file.exists(call_path)) {
    
    stop("Please run CARD first!")
  } else {
    
    call_file <- fread2(call_path, header = F)[, 1]
  }
  
  resule_path <- call_file[1]
  # ref_type <- call_file[2]
  load(resule_path) # CARD_result_list
  
  base_df <- lapply(decon_result_list, function(decon_result){
    sample_info <- decon_result[["sample_info"]]
    sample_info$sample <- paste0("Sample", sample_info$sample)
    base_df_a <- cbind(sample_info,
                       decon_result[["location"]][sample_info$cell, , drop = F],
                       decon_result[["Proportion_CARD"]][sample_info$cell, , drop = F])
    colnames(base_df_a)[1:4] <- c("sample", "cell", "loc_x", "loc_y")
    
    return(base_df_a)
  }) %>% Reduce("rbind", .)
  
  print("Base DECON data formatted!")
  return(base_df)
  
}

# Function 2: Process refine DECON data
refine.process <- function(data_path,
                           expr_norm = FALSE
){
  
  ## load data
  call_path <- paste0(data_path, "/decon_call_file.txt")
  if (!file.exists(call_path)) {
    
    stop("Please run CARD first!")
  } else {
    
    call_file <- fread2(call_path, header = F)[, 1]
  }
  
  refine_resule_path <- call_file[2]
  load(refine_resule_path) # CARD_refine_result_list
  
  ## format base df
  base_df <- lapply(seq_along(decon_refine_result_list), function(a){
    
    decon_result_a <- decon_refine_result_list[[a]]
    refined_prop <- decon_result_a[["refined_prop"]]
    location <- stringr::str_split(rownames(refined_prop), "x", simplify = T) %>% 
      as.data.frame()
    sample_info <- data.frame(sample = paste0("Sample", a),
                              cell = rownames(refined_prop))
    base_df_a <- cbind(sample_info,
                       location,
                       refined_prop)
    colnames(base_df_a)[1:4] <- c("sample", "cell", "loc_x", "loc_y")
    
    return(base_df_a)
  }) %>% Reduce("rbind", .)
  
  ## format
  if (expr_norm == TRUE) {
    
    expr_refine <- lapply(seq_along(CARD_refine_result_list), function(a){
      CARD_refine_result_list[[a]][["refined_expression"]]
    }) %>% Reduce("cbind", .)
    expr_refine <- expr_refine[, base_df$cell]
    expr_refine <- t(expr_refine) %>% as.data.frame()
  } else {
    
    expr_refine <- NULL
  }
  
  print("Base redined DECON data formatted!")
  
  return(list("base_df" = base_df,
              "expr_refine" = expr_refine))
}

# Function 3: Visualize Pie plot
pie.visualize <- function(datt,
                          pie_scale,
                          color_in,
                          add_theme = FALSE
){
  
  ct.select <- colnames(datt)[-c(1:4)]
  sample_size <- length(unique(datt$sample))
  plt <- ggplot() + 
    geom_scatterpie(data = datt, aes(x = loc_x, y = loc_y),
                    pie_scale = 0.3,
                    cols = ct.select, color = NA) + 
    coord_equal() + 
    scale_fill_manual(values =  color_in) +
    theme_void() +
    guides(fill = guide_legend(title = "Cell Type")) 
  if (sample_size != 1){
    
    plt <- plt + facet_wrap(~sample, ncol = 2)
  }
  if (add_theme == TRUE) {
    
    plt <- plt  + 
      theme(legend.position="right",
            legend.title = element_text(size = 13, face="bold"),
            legend.text = element_text(size = 12),
            strip.placement = "outside", 
            strip.background = element_blank(),
            strip.text.x = element_text(size = 12,color = "black"),
            panel.border = element_rect(colour = "grey20", fill = NA),
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
  }
  
  return(plt)
}

# Function 4: Visualize Pie plot after data process
pie.plot <- function(data_path,
                     pie_scale,
                     out_path,
                     out_figure = FALSE, 
                     zip_figure = FALSE
){
  
  # process plot data
  pie_df <- base.process(data_path = data_path)
  sample_size <- length(unique(pie_df$sample))
  
  # plot
  decon_pie_plt <- pie.visualize(datt = pie_df,
                                pie_scale = pie_scale,
                                color_in = COLOR_USE,
                                add_theme = out_figure)
  if(out_figure == TRUE){
    
    pie_ht <- ceiling(sample_size / 2) * 5.5
    pie_wt <- min(2, sample_size) * 5 + 3
    ggsave(filename = paste0(out_path, "/decon_result/decon_pie.tiff"),
           plot = decon_pie_plt,
           height = pie_ht, width = pie_wt, units = "in", 
           dpi = 300, compression = "lzw")
    if (zip_figure == TRUE){
      
      system(paste0("gzip -f ", out_path, "/decon_result/decon_pie.tiff"))
    }
  }
  
  return(decon_pie_plt)
  
}

# Function 5: Process proportion pattern data
propPattern.process <- function(base_df,
                                ct_sel_str = NULL
){
  
  ##
  if (!is.null(ct_sel_str)) {
    ct_sel <-  strsplit(ct_sel_str, ",") %>% 
      unlist %>% as.numeric
  } else {
    ct_sel <- colnames(base_df)[-c(1:4)]
  }
  
  ct_scale <- sapply(ct_sel, function(a){
    prop_a <-  base_df[,a]
    (prop_a - min(prop_a)) / diff(range(prop_a))
  }) %>% as.data.frame()
  
  #
  base_scale_df <- cbind(base_df[, 1:4], ct_scale)
  m_scale_df <- reshape2::melt(base_scale_df,
                               id.vars = c("sample", "cell", "loc_x", "loc_y"),
                               variable.name = "celltype",
                               value.name = "value")
  
  print(paste0("DECON data formatted for ", paste(ct_sel, collapse = ", "), "!"))
  return(m_scale_df)
}

# Function 6: plot proportion pattern
propPattern.visualize <- function(datt,
                                  point_size, 
                                  add_theme = FALSE
){
  
  # 
  sample_size <- datt$sample %>% unique %>% length
  plt <- ggplot(datt, aes(x = loc_x, y = loc_y)) + 
    geom_point(aes(colour = value),size = point_size) +
    scale_color_gradientn(colours = c("lightblue", "lightyellow", "red")) + 
    coord_equal() +
    theme_void() +
    guides(fill = guide_legend(title="value"))
  
  #
  if (sample_size == 1) {
    plt <- plt + facet_wrap(~celltype, ncol = 5,
                            labeller = label_parsed)
  } else {
    plt <- plt + facet_wrap(~celltype * sample , ncol = 5, 
                            dir = "h", labeller = label_parsed)
  }
  
  #  
  if (add_theme == TRUE) {
    
    plt <- plt + 
      theme(legend.position = "right",
            legend.title = element_text(size = 13, face="bold"),
            legend.text = element_text(size = 12),
            strip.placement = "outside", 
            strip.background = element_blank(),
            strip.text.x = element_text(size = 12,color = "black"),
            panel.border = element_rect(colour = "grey20", fill = NA),
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
  }
  return(plt)
}

# Function 7: visualize proportion pattern after data process
propPattern.plot <- function(data_path,
                             out_path,
                             ct_sel,
                             point_size,
                             out_figure = FALSE, 
                             zip_figure = FALSE
){
  # process plot data
  base_df <- base.process(data_path = data_path)
  m_scale_df <- propPattern.process(base_df = base_df)
  n_ct <- length(unique(m_scale_df$celltype))
  
  # plot
  decon_propPattern <- propPattern.visualize(datt = m_scale_df,
                                             point_size = 0.5,
                                             add_theme = out_figure)
  if(out_figure == TRUE){
    
    pp_ht <- ceiling(n_ct / 5) * 5.5
    pp_wt <- min(5, n_ct) * 3 + 1
    ggsave(filename = paste0(out_path, "/decon_result/decon_propPattern.tiff"),
           plot = decon_propPattern,
           height = pp_ht, width = pp_wt, units = "in", 
           dpi = 300, compression = "lzw")
    if(zip_figure == TRUE){
      
      system(paste0("gzip -f ", out_path, "/decon_result/decon_propPattern.tiff"))
    }
  }
  
  return(decon_propPattern)
}

# Function 8: visualize proportion pattern after data process
refine.propPattern.plot <- function(data_path,
                                    out_path,
                                    ct_sel,
                                    point_size,
                                    out_figure = FALSE, 
                                    zip_figure = FALSE
){
  # process plot data
  base_df <- refine.process(data_path,
                            expr_norm = FALSE)[["base_df"]]
  m_scale_df <- propPattern.process(base_df = base_df)
  n_ct <- length(unique(m_scale_df$celltype))
  
  # plot
  decon_propPattern <- propPattern.visualize(datt = m_scale_df,
                                             point_size = point_size,
                                             add_theme = out_figure)
  if(out_figure == TRUE){
    
    pp_ht <- ceiling(n_ct / 5) * 5.5
    pp_wt <- min(5, n_ct) * 3 + 1
    ggsave(filename = paste0(out_path, "/decon_result/decon_refine_propPattern.tiff"),
           plot = decon_propPattern,
           height = pp_ht, width = pp_wt, units = "in", 
           dpi = 300, compression = "lzw")
    if (zip_figure == TRUE){
      
      system(paste0("gzip -f ", out_path, "/decon_result/decon_refine_propPattern.tiff"))
    }
  }
  
  return(decon_propPattern)
}

# Function 9: Visualize refine feature plot after data process 
refine.feature.plot <- function(data_path, 
                                out_path, 
                                point_size,
                                marker_gene = NULL,
                                out_figure = FALSE, 
                                zip_figure = FALSE
){
  
  ## inputs
  decon_refine_list <- refine.process(data_path,
                                   expr_norm = TRUE)
  base_df <- decon_refine_list[["base_df"]]
  expr_refine <- decon_refine_list[["expr_refine"]]
  norm_list <- split(as.data.frame(expr_refine), base_df$sample)
  coord_list <- split(base_df[,c("loc_x", "loc_y")], base_df$sample)
  
  marker_gene_f <- intersect(marker_gene, colnames(expr_refine))
  if (length(marker_gene_f) == 0) {
    stop("No target genes found in expression matrix!")
  } else {
    print(paste0("Plot feature plot on ", paste(marker_gene_f, collapse = ", ", "!")))
  }
  ## inputs
  mrefine_feature_df <- feature.process(norm_list = norm_list, 
                                        coord_list = coord_list,
                                        marker_gene = marker_gene_f)
  n_ct <- length(marker_gene_f)
  feature_plt <- feature.plot(datt = mrefine_feature_df,
                              pointsize = point_size,
                              add_theme = out_figure)
  
  ## output figure
  if (out_figure == TRUE){
    
    ft_ht <- ceiling(n_ct / 5) * 6
    ft_wt <- min(5, n_ct) * 3 + 1
    ggsave(filename = paste0(out_path, "/decon_result/decon_refine_featurePlot.tiff"), 
           plot = feature_plt,
           height = ft_ht, width = ft_wt, units = "in", dpi = 300, compression = "lzw")
    if(zip_figure == TRUE){
    
      system(paste0("gzip -f ", out_path, "/decon_result/decon_refine_featurePlot.tiff"))  
    }
    print("Feature plot saved!")
  }
  
  return(feature_plt)
}

# Function 10: Visualize mapped location plot after data process 
deconMapping.loc.plot <- function(data_path, 
                                  out_path, 
                                  point_size,
                                  out_figure = FALSE, 
                                  zip_figure = FALSE
){
  
  # 1. load data
  call_path <- paste0(data_path, "/decon_call_file.txt")
  if (!file.exists(call_path)) {
    stop("Please run CARD first!")
  } else {
    
    call_file <- fread2(call_path, header = F)[, 1]
  }
  
  map_meta_path <- call_file[5]
  load(map_meta_path) # sce_coord_meta_list
  
  # 2. format input
  sce_coord_meta <- lapply(seq_along(sce_coord_meta_list), function(a){
    
    sce_coord_meta_a <- sce_coord_meta_list[[a]][,c("x", "y", "CT", "Cell")]
    sce_coord_meta_a$sample <- paste0("Sample", a)
    return(sce_coord_meta_a)
  }) %>% Reduce("rbind", .)
  colnames(sce_coord_meta) <- c("loc_x", "loc_y", "celltype", "cell", "sample")
  sample_size <- length(unique(sce_coord_meta$sample))
  n_ct <- length(unique(sce_coord_meta$celltype))
  
  # plot
  ct_plt <- ggplot(sce_coord_meta, aes(x = loc_x, y = loc_y, colour = celltype)) +
    geom_point(size = point_size) +
    scale_colour_manual(values =  COLOR_USE) +
    guides(color = guide_legend(title="Cell Type", byrow = T, nrow = 9,
                                override.aes = list(size = 2)))+
    facet_wrap(~ sample, ncol = 2) +
    theme_bw()
  
  ## output figure
  if (out_figure == TRUE){
    
    ct_plt <- ct_plt +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "white", colour = "black", size = 1),
            panel.spacing.x = unit(0, "cm"), 
            strip.placement = "outside", 
            strip.background = element_blank(),
            strip.text.x = element_text(size = 12,color = "black"),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
    
    ct_ht <- ceiling(sample_size / 2) * 5
    ct_wt <- min(2, sample_size) * 5 + ceiling(n_ct / 9) * 2
    ggsave(filename =  paste0(out_path, "/decon_result/decon_mapping_locPlot.tiff"), 
           plot = ct_plt,
           height = ct_ht, width = ct_wt, units = "in", dpi = 300, compression = "lzw")
    if(zip_figure == TRUE){
    
      system(paste0("gzip -f ", out_path, "/decon_result/decon_mapping_locPlot.tiff"))  
    }
    print("Feature plot saved!")
  }
  
  return(ct_plt)
}

# Function 11
decon_plt.plot <- function(data_path1,                     ## String: output path for decon procedure 
                           data_path2,                     ## String: output path for qc procedure 
                           out_path,                       ## String: output path for decon_plt procedure
                           scMapping = FALSE,              ## Boolean: scMapping
                           out_figures, 
                           zip_figures
){
  
  source(paste0(method_path, "/plt_utils.R"))
  ## plot
  decon_pie_plt <- pie.plot(data_path = data_path1,
                            pie_scale = 0.3,
                            out_path = out_path,
                            out_figure = out_figures, 
                            zip_figure = zip_figures)
  prop_pattern_list <- propPattern.plot(data_path = data_path1,
                                        out_path = out_path,
                                        ct_sel = NULL,
                                        point_size = 0.2,
                                        out_figure = out_figures, 
                                        zip_figure = zip_figures)
  decon_feature_plt <- feature.plot(data_path1 = data_path1,
                                    data_path2 = data_path2,
                                    out_path = out_path,
                                    mode_usage = "decon",
                                    out_figure = out_figures, 
                                    zip_figure = zip_figures)
  ## save data
  save(prop_pattern_list, decon_pie_plt, decon_feature_plt, 
       file = paste0(out_path, "/decon_result/plot.RData"))
  
  if (scMapping == TRUE){
    
    refine_prop_pattern_list <- refine.propPattern.plot(data_path = data_path1,
                                                        out_path = out_path,
                                                        ct_sel = NULL,
                                                        point_size = 0.2,
                                                        out_figure = out_figures, 
                                                        zip_figure = zip_figures)
    save(prop_pattern_list, decon_pie_plt, decon_feature_plt, 
         refine_prop_pattern_list,
         file = paste0(out_path, "/decon_result/plot.RData"))
  } 
  
  return(0)
}
