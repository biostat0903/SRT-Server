#! /usr/bin/env Rscript
# Perform over representation analysis
# up-stream procedure: de.R, svg.R
# down-stream procedure: ora_plt.R

# Load packages 
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(dplyr)
library(plyr)
library(tibble)
library(Matrix)
library(bigreadr)

# Function 1: gsea check
ora.check <- function(data_path1 = NULL,              ## String: data path of svg
                      data_path2 = NULL,              ## String: data path of de
                      method = "svg",                 ## String: up-stream methods, including svg and de
                      species = "Mouse",              ## String: species to analysis, including"Human" or "Mouse" (required)******
                      pathway_db = NULL,              ## String: pathway databases, {Go}
                      out_path = NULL                 ## String: output path of ora procedure
){
  
  ## check 
  post_file <- ifelse(method == "svg", 
                      paste0(data_path1, "/svg_post_file.txt"), 
                      paste0(data_path2, "/de_post_file.txt"))
  gene_sig <- ifelse(read.table(post_file)[1, 1] == "0", 
                     read.table(post_file)[2, 1],
                     read.table(post_file)[1, 1])
  ## output
  write.table(c(gene_sig, species, pathway_db, method), 
              file = paste0(out_path, "/ora_check_file.txt"), 
              row.names = F, quote = F, col.names = F)
  
  return(0)
}

# Function 2: ora function
ora.func <- function(gene_id, 
                     species, 
                     species_db, 
                     pathway_db
){

  ## GO analysis
  if (pathway_db == "GO" | pathway_db == "ALL"){
    
    ora_res_GO <- enrichGO(gene = gene_id[, 1],
                           keyType = 'SYMBOL',
                           OrgDb = species_db,
                           ont = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 1,
                           readable = TRUE)@result
  }
  ## KEGG analysis
  if (pathway_db == "KEGG" | pathway_db == "ALL"){

    species_str <- ifelse(species == "Mouse", "mmu", "hsa")
    ora_res_KEGG <- enrichKEGG(gene = gene_id[, 2],
                               organism = species_str,
                               keyType = "ncbi-geneid",
                               pvalueCutoff = 1)@result
  }
  ## WikiPathways analysis
  if (pathway_db == "WikiPathways" | pathway_db == "ALL"){

    species_str <- ifelse(species == "Mouse",
                          "Mus musculus", "Homo sapiens")
    ora_res_WikiPathways <- enrichWP(gene = gene_id[, 2],
                                     organism = species_str)@result
  }
  ## Reactome pathway analysis
  if (pathway_db == "Reactome" | pathway_db == "ALL"){
    
    species_str <- ifelse(species == "Mouse", 
                          "mouse", "human")
    ora_res_Reactome <- enrichPathway(gene = gene_id[, 2], 
                                      organism = species_str,
                                      readable = TRUE)@result
  }
  ## Disease Oncology pathway analysis
  ora_res_DO <- NULL
  if (species == "Human"){
    
    if (pathway_db == "DO" | pathway_db == "ALL"){
      
      ora_res_DO <- enrichDO(gene = gene_id[, 2], 
                             readable = TRUE)@result
    }
  }
  ## ALL pathway summary
  if(pathway_db == "ALL"){
    
    ### GO 
    db_info <- ifelse(ora_res_GO$ONTOLOGY == "BP", "GO Biological Process", 
                      ifelse(ora_res_GO$ONTOLOGY == "CC", "GO Cell Component", "GO Molecular Function"))
    ora_res_GO$ONTOLOGY <- NULL
    ora_res_GO$PathwayInfo <- db_info
    ### KEGG
    ora_res_KEGG$PathwayInfo <- "KEGG"
    ### Reactome
    ora_res_Reactome$PathwayInfo <- "Reactome"
    ### WikiPathways
    ora_res_WikiPathways$PathwayInfo <- "WikiPathways"
    ### DO
    if (!is.null(ora_res_DO)){
      
      ora_res_DO$PathwayInfo <- "DO"
    }
    # 
    ora_res <- rbind(ora_res_GO,  
                     ora_res_KEGG, ora_res_WikiPathways, 
                     ora_res_Reactome, ora_res_DO)
  } else {
    
    ora_res <- eval(parse(text = paste0("ora_res_", pathway_db)))
  }
  return(ora_res)
}

# Function 3: ora call function
ora.call <- function(out_path
){

  ## load
  check_file <- read.table(paste0(out_path, "/ora_check_file.txt"))[, 1]
  gene_dat <- fread2(check_file[1])
  if(ncol(gene_dat) > 1){
    
    gene_name <- gene_dat[gene_dat$p_adj < 0.05, 1] %>% unique()
  } else {
    
    gene_name <- gene_dat[, 1]
  }
  species <- check_file[2]
  pathway_db <- check_file[3]
  method <- check_file[4]
  
  ## transfer gene name
  if(check_file[2] == "Mouse"){
    
    species_db <- org.Mm.eg.db
  }else{
    
    species_db <- org.Hs.eg.db
  }
  gene_id <- bitr(gene_name, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = species_db)
  ora_res <- ora.func(gene_id, species, species_db, pathway_db)
  ora_res_clx <- NULL
  # if (method == "deg"){
  # 
  #   cl_num <- length(unique(gene_dat$cluster))
  #   ora_res_clx <- alply(c(1: cl_num), 1, function(cl){
  #     cat(cl, "\n")
  #     gene_dat_cl <- gene_dat[gene_dat$p_adj < 0.05 & gene_dat$cluster == cl, 1]
  #     gene_id_cl <- gene_id[gene_id[, 1] %in% gene_dat_cl, ]
  #     ora_res_cl <- try(ora.func(gene_id_cl, species, species_db, pathway_db), 
  #                       silent = T)
  #     if (inherits(ora_res_cl, "try-error") == F){
  #       
  #       ora_res_cl <- cbind.data.frame(cluster = cl, ora_res_cl)
  #     } 
  #     return(ora_res_cl)
  #   }) %>% do.call("rbind", .)
  # }
  ## output
  save(ora_res, ora_res_clx, 
       file = paste0(out_path, "/ora_call_result.RData"))
  write.table(paste0(out_path, "/ora_call_result.RData"), 
              file = paste0(out_path, "/ora_call_file.txt"), 
              row.names = F, col.names = F, quote = F)
  return(0)
}

# Function 3: ora post file
ora.post <- function(out_path
){
  
  ## load
  load(paste0(out_path, "/ora_call_result.RData"))
  check_file <- read.table(paste0(out_path, "/ora_check_file.txt"))[, 1]
  method <- check_file[4]
  ## output
  result_dir <- paste0(out_path, "/ora_result")
  if (!file.exists(result_dir)) {
    
    system(paste0("mkdir -p ", result_dir))
  }
  write.table(ora_res, 
              file = paste0(result_dir, "/ora_pathway_", method, ".txt"), 
              row.names = F, quote = F, sep = "\t")
  if (method == "deg"){
    
    write.table(ora_res_clx, 
                file = paste0(result_dir, "/ora_pathway_", method, "_cluster.txt"), 
                row.names = F, quote = F, sep = "\t")
  }
  write.table(paste0(result_dir, "/ora_pathway_", method, ".txt"), 
              file = paste0(out_path, "/ora_post_file.txt"), 
              row.names = F, quote = F, col.names = F)
  return(0)
}
