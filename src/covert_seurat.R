#
library(optparse)
library(Seurat)
library(SeuratDisk)
library(glue)
## Input parameters
args_list = list(
  make_option("--spec", type="character", default=NULL,
              help="INPUT: species",
              metavar="character"),
  make_option("--ref", type="character", default=NULL,
              help="INPUT: reference name",
              metavar="character"))
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)
refx = opt$ref
spec = opt$spec

#
setwd("/public/home/biostat03/project/stwebProject/02_data/reference_data/SRT-Server/scRNA-seq/")
# spec <- "Human"
# all_ref <- list.files(paste0("./", spec, "/"))
# for (refx in all_ref) {
# set path
sc_count_file <- paste0("./", spec, "/", refx, "/sc_count.RData")
sc_meta_file <- paste0("./", spec, "/", refx, "/sc_meta.RData")
CARD_path <- paste0("./", spec, "/", refx, "/CARD")
c2l_path <- paste0("./", spec, "/", refx, "/cell2location")
tg_path <- paste0("./", spec, "/", refx, "/tangram")
# load data
if (file.exists(sc_count_file) &
    file.exists(sc_meta_file)) {
  load(sc_count_file)
  load(sc_meta_file)
  # create Seurat Object
  seurat_obj <- CreateSeuratObject(counts = sc_count, 
                                   meta.data = sc_meta)
  SaveH5Seurat(seurat_obj, filename = paste0("./", spec, "/", refx, "/sc.h5Seurat"))
  # convert to h5ad
  if (!file.exists(tg_path)) {
    system(glue("mkdir {tg_path}"))
  }
  Convert(paste0("./", spec, "/", refx, "/sc.h5Seurat"), 
          dest = paste0(tg_path, "/sc.h5ad"))
  # move to CARD
  if (!file.exists(CARD_path)) {
    system(glue("mkdir {CARD_path}"))
  }
  system(glue("mv {sc_count_file} {CARD_path}/"))
  system(glue("mv {sc_meta_file} {CARD_path}/"))
  
} else {
  message("sc files not exist!")
}
# }

