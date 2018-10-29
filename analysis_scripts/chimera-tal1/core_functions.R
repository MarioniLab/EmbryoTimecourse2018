# get in the embryo functions

source("/nfs/research1/marioni/jonny/embryos/scripts/core_functions.R")


load_data = function(normalise = TRUE){
  require(scran)
  require(scater)
  require(SingleCellExperiment)
  require(Matrix)
  
  # counts = readMM("/nfs/research1/marioni/jonny/chimera-tal1/data/raw_counts.mtx")
  counts = readRDS("/nfs/research1/marioni/jonny/chimera-tal1/data/raw_counts.rds")
  genes <<- read.table("/nfs/research1/marioni/jonny/chimera-tal1/data/genes.tsv", stringsAsFactors = F)
  meta <<- read.table("/nfs/research1/marioni/jonny/chimera-tal1/data/meta.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  
  rownames(counts) = genes[,1] #ensembl
  colnames(counts) = meta$cell

  sce <<- SingleCellExperiment(assays = list("counts" = as(counts, "dgCMatrix")))
  
  if(normalise){
    sfs <<- read.table("/nfs/research1/marioni/jonny/chimera-tal1/data/sizefactors.tab", stringsAsFactors = F)[,1]
    sizeFactors(sce) = sfs
    sce <<- normalize(sce)
  }
  
  invisible(0)
}


