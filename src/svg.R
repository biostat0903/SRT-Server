#! /usr/bin/env Rscript
# Identify spatially variable genes
# up-stream procedure: qc.R
# down-stream procedure: svg_plt.R

# method_path setting
method_path <- "/public/home/biostat03/project/stwebProject/01_code/srt_server/dev"

# load packages 
library(SPARK)
library(dplyr)
library(amap)

# Function 1: check function
svg.check <- function(data_path = NULL,              ## String: output path of qc procedure
                      svg_method = NULL,             ## String: method for svg: {SPARK}, {SPARK-X}, {SpatialDE}
                      pattern_est = FALSE,           ## Boolean: estimate svg pattern 
                      pval_thr = 0.01,               ## Float: P threshold for select significant svg
                      clus_num = 3,                  ## Integer: the number of cluster
                      perm = FALSE,                  ## Boolean: estimate permutation test
                      out_path = NULL                ## String: output path of svg procedure
){
  
  ## load st data
  check_file <- paste0(data_path, "/qc_call_file.txt")
  if(!file.exists(check_file)){
    
    stop("Can't find \"qc_call_file.txt\" file! Please run QC module!")
  } 
  
  ### output file
  if(perm == FALSE){
    
    if(pattern_est == FALSE){
      
      write.table(c(svg_method), 
                  file = paste0(out_path, "/svg_check_file.txt"), 
                  row.names = F, quote = F, col.names = F)
    } else {
      
      write.table(c(svg_method, pattern_est, pval_thr, clus_num), 
                  file = paste0(out_path, "/svg_check_file.txt"), 
                  row.names = F, quote = F, col.names = F)
    }
  } else {
    
    if(svg_method != "SPARK-X"){
      
      warning("The permutation procedure might spent some time! We recommend use \"SPARK-X\".")
    }
    
    if(pattern_est == FALSE){
      
      write.table(c(svg_method, perm), 
                  file = paste0(out_path, "/svg_check_file.txt"), 
                  row.names = F, quote = F, col.names = F)
    } else {
      
      write.table(c(svg_method, pattern_est, perm, pval_thr, clus_num), 
                  file = paste0(out_path, "/svg_check_file.txt"), 
                  row.names = F, quote = F, col.names = F)
    }
  }
  
  return(0)
}

# Function 2: summarize pattern
## Function 2.1: variance stabilizing transformation: NB
var.stabilize <- function(x, sv = 1) {
  
  varx <- apply(x, 1, var)
  meanx <- apply(x, 1, mean)
  phi <- coef(nls(varx ~ meanx + phi * meanx^2, 
                  start = list(phi = sv)))
  return(log(x + 1/(2 * phi)))
}
## Function 2.2: relative expression estimation
relExpr.estimate <- function(expr_mat){
  
  expr_range <- max(expr_mat) - min(expr_mat)
  rel_expr <- (expr_mat - min(expr_mat))/expr_range
  return(rel_expr)
}
## Function 2.3: pattern definition
pattern.define <- function(exp_mat,
                           loc_info,
                           svg_result,
                           pval_thr,
                           clus_num
){
  
  ## 
  LMReg <- function(ct, LB) {
    return(lm(ct ~ LB)$residuals)
  }
  vst_count <- var.stabilize(exp_mat)
  sig_vst_idx <- which(svg_result$adjPval < pval_thr)
  if (length(sig_vst_idx) > 10000){
    
    warning("We select 10,000 svgs to perform clustering!")
    sig_vst_count <- vst_count[order(svg_result$adjPval), ][1:10000, ]
  } else {
    
    sig_vst_count <- vst_count[sig_vst_idx, ]
  }
  lib_size <- apply(exp_mat, 2, sum)
  sig_vst_res <- t(apply(sig_vst_count, 1, LMReg, 
                         LB = log(lib_size)))
  ##
  hc <- hcluster(sig_vst_res, method = "euc",
                 link = "ward", nbproc = 5, 
                 doubleprecision = TRUE)
  memb <- cutree(hc, k = clus_num)
  cent <- NULL
  for (k in 1: clus_num) {
    cent <- cbind(cent, colMeans(sig_vst_res[memb == k, , drop = FALSE]))
  }
  rownames(loc_info) <- rownames(cent)
  rel_cent <- t(apply(cent, 1, relExpr.estimate))
  pattern_dat <- setNames(cbind.data.frame(rownames(loc_info), loc_info, rel_cent), 
                          c("cell", "x", "y", paste0("Pattern", c(1: clus_num))))
  
  return(pattern_dat)
}

# Function 3: SPARK
Spark.svg.test <- function(exp_mat, 
                           loc_info,
                           percentage = 0.1,       
                           min_counts = 100        
){
  
  ## test
  if (!identical(colnames(exp_mat), rownames(loc_info))){
    colnames(exp_mat) <- rownames(loc_info)
  }
  spark <- CreateSPARKObject(counts = exp_mat, 
                             location = loc_info[,1:2],
                             percentage = percentage, 
                             min_total_counts = min_counts)
  spark@lib_size <- apply(spark@counts, 2, sum)
  spark <- spark.vc(spark, 
                    covariates = NULL, 
                    lib_size = spark@lib_size, 
                    num_core = 1,
                    verbose = F)
  spark <- spark.test(spark, 
                      check_positive = T, 
                      verbose = F)
  res_spark <- spark@res_mtest[, c("combined_pvalue", "adjusted_pvalue")]
  res_spark <- cbind.data.frame(row.names(res_spark), res_spark)
  return(res_spark)
}

# Function 4: SpatialDE
SpatialDE.svg.test <- function(exp_mat, 
                               loc_info
){
  
  ## test  
  norm_expr <- stabilize(exp_mat)
  resid_expr <- regress_out(norm_expr, sample_info = loc_info)
  res_spatialDE <- spatialDE::run(resid_expr, 
                                  coordinates = loc_info[, c(1, 2)])
  return(res_spatialDE[, c(3, 17, 18)])
}

# Function 5: SPARK-X
SparkX.svg.test <- function(sp_exp_mat, 
                            location
){
  
  ## test
  sparkX <- sparkx(sp_exp_mat, 
                   location, 
                   numCores = 1,
                   option = "mixture")
  res_sparkx <- sparkX$res_mtest
  res_sparkx <- cbind.data.frame(row.names(res_sparkx),
                                 res_sparkx)
  return(res_sparkx)
}

# Function 6: call function
svg.call <- function(data_path = NULL,                ## String: output path of qc procedure
                     out_path = NULL                  ## String: output path of svg procedure
){
  
  ## load io code
  source(paste0(method_path, "/io.R"))
  
  ## load st data
  check_file <- paste0(data_path, "/qc_call_file.txt")
  qc_param <- read.table(check_file)[, 1]
  spatial_data_filename <- qc_param[1]
  sample_size <- qc_param[2] %>% as.numeric
  spatial_data <- h5data.load(spatial_data_filename,
                              sample_size ,
                              load_count = TRUE,
                              normalization = FALSE,
                              load_coord = TRUE,
                              coordinate = TRUE)
  count_list <- spatial_data[["count_list"]]
  coord_list <- spatial_data[["coord_list"]]
  
  ## load input
  check_file <- read.table(file = paste0(out_path, "/svg_check_file.txt"))[, 1]
  svg_method <- check_file[1]
  if (length(check_file) == 1){
    
    pattern_est <- FALSE
    perm <- FALSE
  }
  if(length(check_file) == 2){
    
    pattern_est <- FALSE
    perm <- TRUE
  } 
  if (length(check_file) == 4){
    
    pattern_est <- TRUE
    perm <- FALSE
  }
  if (length(check_file) == 5){
    
    pattern_est <- TRUE
    perm <- TRUE
  }
  if (pattern_est == TRUE){
    
    pval_thr <- check_file[length(check_file)-1] %>% as.numeric()
    clus_num <- check_file[length(check_file)] %>% as.numeric()
  }
  
  ## construct data for three methods
  ### SPARK and SpatialDE
  if (svg_method %in% c("SPARK", "SpatialDE")){
    
    #### total data
    tot_list <- plyr::alply(c(1: sample_size), 1, function(a){
      
      exp_mat_s <- as.matrix(count_list[[a]])
      loc_info_s <- data.frame(x = coord_list[[a]][, 1],
                               y = coord_list[[a]][, 2],
                               total_counts = colSums(exp_mat_s)) %>% as.matrix
      colnames(loc_info_s) <- c("x", "y", "total_counts")
      return(list(exp_mat_s, loc_info_s))
    })
  }
  ### SPARK-X
  if (svg_method == "SPARK-X"){
    
    #### total data
    tot_list <- plyr::alply(c(1: sample_size), 1, function(a){
      
      exp_mat_s <- as.matrix(t(count_list[[a]]))
      loc_info_s <- data.frame(x = coord_list[[a]][, 1],
                               y = coord_list[[a]][, 2]) %>% as.matrix
      colnames(loc_info_s) <- c("x", "y")
      rownames(loc_info_s) <- rownames(exp_mat_s)
      return(list(exp_mat_s, loc_info_s))
    })
  }
  
  ## define svg using different methods
  tot_svg_result <- plyr::llply(tot_list, function(a){
    
    if (svg_method == "SPARK-X"){
      
      exp_mat_s <- t(a[[1]])
    } else {
      
      exp_mat_s <- a[[1]]
    }
    loc_info_s <- a[[2]]
    if (svg_method == "SPARK-X"){
      
      svg_s <- SparkX.svg.test(exp_mat_s, loc_info_s)
    }
    if (svg_method == "SPARK"){
      
      svg_s <- Spark.svg.test(exp_mat_s, as.data.frame(loc_info_s))
    }
    if (svg_method == "SpatialDE"){
      
      svg_s <- SpatialDE.svg.test(exp_mat_s, loc_info_s)
    }
    colnames(svg_s) <- c("Gene", "Pval", "adjPval")
    return(svg_s) 
  })
  
  ## permuate samples 
  tot_svg_result_perm <- NULL
  if (perm == TRUE){
    
    message("Permutation test begin....")
    tot_svg_result_perm <- plyr::llply(tot_list, function(a){
      
      if (svg_method == "SPARK-X"){
        
        exp_mat_s <- t(a[[1]])
      } else {
        
        exp_mat_s <- a[[1]]
      }
      loc_info_s <- a[[2]]
      perm_tot <- data.frame()
      for (i in 1: 10){
        
        index_perm <- sample(1: ncol(exp_mat_s))
        exp_mat_perm <- exp_mat_s[, index_perm]
        if (svg_method == "SPARK-X"){
          
          sparkx_perm <- SparkX.svg.test(exp_mat_perm, loc_info_s)
          colnames(sparkx_perm) <- c("Gene", "Pval", "adjPval")
          sparkx_perm$Gene <- paste0(sparkx_perm$Gene, "_", i)
          perm_tot <- rbind(perm_tot, sparkx_perm)
        }
        if (svg_method == "SPARK"){
          
          spark_perm <- Spark.svg.test(exp_mat_perm, as.data.frame(loc_info_s))
          colnames(spark_perm) <- c("Gene", "Pval", "adjPval")
          spark_perm$Gene <- paste0(spark_perm$Gene, "_", i)
          perm_tot <- rbind(perm_tot, spark_perm)
        }
        if (svg_method == "SpatialDE"){
          
          spatialDE_perm <- SpatialDE.svg.test(exp_mat_perm, loc_info_s)
          colnames(spatialDE_perm) <- c("Gene", "Pval", "adjPval")
          spatialDE_perm$Gene <- paste0(spatialDE_perm$Gene, "_", i)
          perm_tot <- rbind(perm_tot, spatialDE_perm)
        }
      }
      return(perm_tot)
    })
    message("Permutation test end.")
  }
  
  ## estimate pattern
  tot_pattern_result <- NULL
  if (pattern_est == TRUE){
    
    tot_pattern_result <- plyr::alply(c(1: sample_size), 1, function(a){
      if (svg_method == "SPARK-X"){
        
        exp_mat_s <- t(tot_list[[a]][[1]])
      } else {
        
        exp_mat_s <- tot_list[[a]][[1]]
      }
      loc_info_s <- tot_list[[a]][[2]]
      svg_result_s <- tot_svg_result[[a]]
      pattern_dat_s <- pattern.define(exp_mat = exp_mat_s,
                                      loc_info = loc_info_s,
                                      svg_result = svg_result_s,
                                      pval_thr = pval_thr,
                                      clus_num = clus_num) %>%
        reshape2::melt(.,
                       id.vars = c("cell", "x", "y"),
                       variable.name = "patternType",
                       value.name = "patternValue")
      return(pattern_dat_s)
      
    })
  } 
  
  ## output parameters
  svg_result_path <- paste0(out_path, "/svg_call_", svg_method,
                            "_result.RData")
  write.table(c(svg_result_path, svg_method,
                sample_size, pattern_est),
              file = paste0(out_path, "/svg_call_file.txt"),
              col.names = F, row.names = F, quote = F)
  
  ## output data
  save(tot_svg_result, tot_svg_result_perm, tot_pattern_result,
       file = svg_result_path)
  
  return(0)
}

# Function 7: post function
svg.post <- function(out_path = NULL     ## out path 
){
  
  ## load
  call_file <- read.table(paste0(out_path, "/svg_call_file.txt"))[, 1]
  svg_result_file <- call_file[1]
  svg_method <- call_file[2]
  load(svg_result_file)
  sample_size <- length(tot_svg_result)
  
  ## process svg
  sig_svg <- plyr::llply(tot_svg_result, function(a){
    
    svg_s <- a[order(a$Pval), ]
    svg_s <- a[a$adjPval < 0.05, 1]
    if (length(svg_s) >= 20000){
      return(svg_s[1: 20000])
    } else {
      return(svg_s)
    }
  }) %>% unlist() %>% unique()
  # all_gene <- tot_svg_result[[1]][, 1]
  # sig_svg_dat <- data.frame(Gene = all_gene, 
  #                           Sel = ifelse(all_gene %in% sig_svg, 1, 0))
  
  # output
  if (!file.exists(paste0(out_path, "/svg_result"))) {
    system(paste0("mkdir '", out_path, "/svg_result'"))
  }
  result_out <- plyr::aaply(c(1: sample_size), 1, function(a){
    write.csv(tot_svg_result[[a]], 
              file = paste0(out_path, "/svg_result/svg_", 
                            svg_method, "_result_s", a, ".csv"), 
              row.names = F, quote = F)
    return(a)
  })
  
  write.csv(sig_svg, row.names = F, quote = F, 
            file = paste0(out_path, "/svg_result/svg_", 
                          svg_method, "_result_sig.csv"))
  write.table(paste0(out_path, "/svg_result/svg_", 
                     svg_method, "_result_sig.csv"), 
              file = paste0(out_path, "/svg_post_file.txt"), 
              col.names = F, row.names = F, quote = F)
  
  return(0)
}
