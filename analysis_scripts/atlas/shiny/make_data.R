set.seed(42)

source("/nfs/research1/marioni/jonny/embryos/scripts/core_functions.R")
load_data(remove_doublets = TRUE, remove_stripped = TRUE, load_corrected = TRUE)


library(Matrix)
library(irlba)
library(scran)
library(Rtsne)
library(igraph)
library(biomaRt)
library(rhdf5)
library(rgexf)

library(BiocParallel)
mcparam = SnowParam(workers = 16)

translate_ensembl = function(ensembl_gene_vec, gene_df = genes){
  return(gene_df[match(ensembl_gene_vec, gene_df[,1]) ,2])
}

translate_mgi = function(mgi_gene_vec, gene_df = genes){
  return(gene_df[match(mgi_gene_vec, gene_df[,2]) ,1])
}



# HDF5 COUNTS ####
# library(HDF5Array)
# 
# count_hdf5 = "/nfs/research1/marioni/jonny/embryos/scripts/shiny/counts.hdf5"
# 
# if(file.exists(count_hdf5)){
#   file.remove(count_hdf5)
# }
# writeHDF5Array(x = DelayedArray(t(logcounts(sce))), file = count_hdf5, name = "logcounts", verbose = TRUE)

# GRAPH FOR LOCAL VIS
# graph = buildKNNGraph(sce)


# tSNEs
set.seed(42)
pca = corrected$all
tsne = Rtsne(pca, pca = FALSE)$Y

tsnes = lapply(unique(meta$stage), function(x){
  pca = corrected$stage[[x]]
  return(Rtsne(pca, pca = FALSE)$Y)
})

tsnes_ts = lapply(unique(meta$theiler), function(x){
  pca = corrected$theiler[[x]]
  return(Rtsne(pca, pca = FALSE)$Y)
})


tsnes = c(tsnes, list(all = tsne), tsnes_ts)

names(tsnes) = c(unique(meta$stage), "all", unique(meta$theiler))


# UMAPs

do_umap = function(pca){
  write.table(pca, file = "/nfs/research1/marioni/jonny/embryos/scripts/shiny/scanpy_input.tab", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  system("python3 /nfs/research1/marioni/jonny/embryos/scripts/shiny/umap.py")
  umap = read.table("/nfs/research1/marioni/jonny/embryos/scripts/shiny/scanpy_output.tab", sep = "\t", header = FALSE)
  return(as.matrix(umap))
}

umap = read.table("/nfs/research1/marioni/jonny/embryos/scripts/vis/umap/umap.tab", header = FALSE)

umap_stage = lapply(unique(meta$stage), function(x){
  pca = corrected$stage[[x]]
  return(do_umap(pca))
})

umap_ts = lapply(unique(meta$theiler), function(x){
  pca = corrected$theiler[[x]]
  return(do_umap(pca))
})

umaps = c(umap_stage, list(all = umap), umap_ts)

names(umaps) = c(unique(meta$stage), "all", unique(meta$theiler))


save(tsnes,
     umaps,
     meta,
     genes,
     file = "/nfs/research1/marioni/jonny/embryos/scripts/shiny/data.RData")


# CELL TYPE MARKERS ####


markers_celltype = findMarkers(sce,
                               clusters = meta$celltype,
                               block = meta$sample,
                               pval.type = "all",
                               direction = "up")



# SAVE RESULTS ####


save(tsnes,
     umaps,
     meta,
     genes,
     markers_celltype,
     file = "/nfs/research1/marioni/jonny/embryos/scripts/shiny/data.RData")


saveRDS(genes, file = "/nfs/research1/marioni/jonny/embryos/scripts/shiny/genes.rds")
saveRDS(meta, file = "/nfs/research1/marioni/jonny/embryos/scripts/shiny/meta.rds")
