celltype_colours = c("Epiblast" = "#635547",
                     "Primitive Streak" = "#DABE99",
                     "Caudal epiblast" = "#9e6762",
                     
                     "PGC" = "#FACB12",
                     
                     "Anterior Primitive Streak" = "#c19f70",
                     "Notochord" = "#0F4A9C",
                     "Def. endoderm" = "#F397C0",
                     "Gut" = "#EF5A9D",
                     
                     "Nascent mesoderm" = "#C594BF",
                     "Mixed mesoderm" = "#DFCDE4",
                     "Intermediate mesoderm" = "#139992",
                     "Caudal Mesoderm" = "#3F84AA",
                     "Paraxial mesoderm" = "#8DB5CE",
                     "Somitic mesoderm" = "#005579",
                     "Pharyngeal mesoderm" = "#C9EBFB",
                     "Cardiomyocytes" = "#B51D8D",
                     "Allantois" = "#532C8A",
                     "ExE mesoderm" = "#8870ad",
                     "Mesenchyme" = "#cc7818",
                     
                     "Haematoendothelial progenitors" = "#FBBE92",
                     "Endothelium" = "#ff891c",
                     "Blood progenitors 1" = "#f9decf",
                     "Blood progenitors 2" = "#c9a997",
                     "Erythroid1" = "#C72228",
                     "Erythroid2" = "#f79083",
                     "Erythroid3" = "#EF4E22",
                     
                     "NMP" = "#8EC792",
                     
                     "Rostral neurectoderm" = "#65A83E",
                     "Caudal neurectoderm" = "#354E23",
                     "Neural crest" = "#C3C388",
                     "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                     "Spinal cord" = "#CDE088",
                     
                     "Surface ectoderm" = "#f7f79e",
                     
                     "Visceral endoderm" = "#F6BFCB",
                     "ExE endoderm" = "#7F6874",
                     "ExE ectoderm" = "#989898",
                     "Parietal endoderm" = "#1A1A1A"
                     
)

haem_colours = c(
  "Mes1"= "#c4a6b2",
  "Mes2"= "#ca728c",
  
  "Cardiomyocytes" =  "#B51D8D",  
  
  "BP1" = "#6460c5",
  "BP2" = "#96b8e4",
  "Haem3"= "#02f9ff",
  "BP3" = "#07499f",
  "BP4" = "#036ef9",
  
  "Haem1"= "#bb22a7",
  "Haem2" = "#f695e9",
  "Haem4" = "#4c4a81",
  
  "EC1"= "#006737",
  
  "EC2" = "#5eba46",
  "EC3" = "#818068",
  "EC4"="#d6de22",
  "EC5"="#5f7238",
  "EC6"="#2ab572",
  "EC7"="#000000",
  "EC8"="#a0cb3b",
  
  "Ery1"="#f67a58",
  "Ery2" ="#a26724",
  "Ery3"="#cdaf7f",
  "Ery4"= "#625218",
  
  "My" = "#c62127",
  "Mk"= "#f6931d"
)

stage_colours = c("E6.5" = "#D53E4F",
                  "E6.75" = "#F46D43",
                  "E7.0" = "#FDAE61",
                  "E7.25" = "#FEE08B",
                  "E7.5" = "#FFFFBF",
                  "E7.75" = "#E6F598",
                  "E8.0" = "#ABDDA4",
                  "E8.25" = "#66C2A5",
                  "E8.5" = "#3288BD",
                  "mixed_gastrulation" = "#A9A9A9")

stage_labels = names(stage_colours)
names(stage_labels) = names(stage_colours)
stage_labels[10] = "Mixed"

# ROUTINELY USED PACKAGES
library(irlba)
library(cowplot)


load_data = function(normalise = TRUE, remove_doublets = FALSE, remove_stripped = FALSE, load_corrected = FALSE){
  
  if(load_corrected & (!remove_doublets | !remove_stripped)){
    message("Using corrected PCs, also removing doublets + stripped now.")
    remove_doublets = TRUE
    remove_stripped = TRUE
  }
  
  require(scran)
  require(scater)
  require(SingleCellExperiment)
  require(Matrix)
  
  # counts = readMM("/nfs/research1/marioni/jonny/embryos/data/raw_counts.mtx")
  counts = readRDS("/nfs/research1/marioni/jonny/embryos/data/raw_counts.rds")
  genes = read.table("/nfs/research1/marioni/jonny/embryos/data/genes.tsv", stringsAsFactors = F)
  meta = read.table("/nfs/research1/marioni/jonny/embryos/data/meta.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
  
  rownames(counts) = genes[,1] #ensembl
  colnames(counts) = meta$cell

  sce = SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs = read.table("/nfs/research1/marioni/jonny/embryos/data/sizefactors.tab", stringsAsFactors = F)[,1]
    sizeFactors(sce) = sfs
    sce = scater::normalize(sce)
  }
  
  if(remove_doublets){
    sce = scater::normalize(sce[,!meta$doublet])
    meta = meta[!meta$doublet,]
  }
  
  if(remove_stripped){
    sce = scater::normalize(sce[,!meta$stripped])
    meta = meta[!meta$stripped, ]
  }
  
  if(load_corrected){
    corrected = readRDS("/nfs/research1/marioni/jonny/embryos/data/corrected_pcas.rds")
    assign("corrected", corrected, envir = .GlobalEnv)
    
  }
  

  
  
  assign("genes", genes, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  assign("sce", sce, envir = .GlobalEnv)
  
  
  invisible(0)
}


#removed yellow from position 2 ("#FFFF00")
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour", "Publication",
                 manual_pal(values = c(
                   "#000000", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                   "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                   "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                   "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                   "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                   "#372101", "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                   "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                   "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                   "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                   "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                   "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                   "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                   "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")), ...)
}

scale_color_Publication = scale_colour_Publication

#removed yellow from position 2 ("#FFFF00")
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill", "Publication",
                 manual_pal(values = c(
                   "#000000", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                   "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                   "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                   "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                   "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                   "#372101", "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                   "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                   "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                   "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                   "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                   "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                   "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                   "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")), ...)
}


#ensure counts has columns names for the cells
#match timepoints,samples to the count table
#timepoint_order, sample_order should contain each sample/timepoint ONCE, in correct order for correction
doBatchCorrect = function(counts, timepoints, samples, timepoint_order, sample_order, npc = 50, pc_override = NULL, BPPARAM = SerialParam()){
  require(scran)
  require(irlba)
  require(BiocParallel)
  
  if(!is.null(pc_override)){
    pca = pc_override
  } else {
    pca = prcomp_irlba(t(counts), n = npc)$x
    rownames(pca) = colnames(counts)
  }
  
  if(length(unique(samples)) == 1){
    return(pca)
  }
  
  #create nested list
  pc_list = lapply(unique(timepoints), function(tp){
    sub_pc = pca[timepoints == tp, , drop = FALSE]
    sub_samp = samples[timepoints == tp]
    list = lapply(unique(sub_samp), function(samp){
      sub_pc[sub_samp == samp, , drop = FALSE]
    })
    names(list) = unique(sub_samp)
    return(list)
  })
  
  names(pc_list) = unique(timepoints)
  
  #arrange to match timepoint order
  pc_list = pc_list[order(match(names(pc_list), timepoint_order))]
  pc_list = lapply(pc_list, function(x){
    x[order(match(names(x), sample_order))]
  })
  
  #perform corrections within list elements (i.e. within stages)
  correct_list = lapply(pc_list, function(x){
    if(length(x) > 1){
      return(do.call(fastMNN, c(x, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected)
    } else {
      return(x[[1]])
    }
  })
  
  #perform correction over list
  if(length(correct_list)>1){
    correct = do.call(fastMNN, c(correct_list, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected
  } else {
    correct = correct_list[[1]]
  }
  
  correct = correct[match(colnames(counts), rownames(correct)),]
  
  return(correct)
  
}

#gene_df is the cellranger gene table
getHVGs = function(sce, min.mean = 1e-3, gene_df = genes){
  require(biomaRt)
  trend = scran::trendVar(sce, use.spikes = FALSE, loess.args = list(span = 0.05))
  decomp = scran::decomposeVar(sce, fit = trend)
  decomp = decomp[decomp$mean > min.mean,]
  
  #exclude sex genes
  xist = "ENSMUSG00000086503"
  # mouse_ensembl = useMart("ensembl")
  # mouse_ensembl = useDataset("mmusculus_gene_ensembl", mart = mouse_ensembl)
  # gene_map = getBM(attributes=c("ensembl_gene_id", "chromosome_name"), filters = "ensembl_gene_id", values = rownames(decomp), mart = mouse_ensembl)
  # ychr = gene_map[gene_map[,2] == "Y", 1]
  ychr = read.table("/nfs/research1/marioni/jonny/embryos/data/ygenes.tab", stringsAsFactors = FALSE)[,1]
  other = c("tomato-td") #for the chimera
  decomp = decomp[!rownames(decomp) %in% c(xist, ychr, other),]
  
  decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.05])
}



#use ensembl gene IDs
enrich_genes_topGO = function(hits, universe, n.out = 100){
  require(topGO)
  require(org.Mm.eg.db)
  
  hits = as.character(hits)
  universe = as.character(universe)
  
  input = as.numeric(universe %in% hits)
  names(input) = universe
  
  topgo = new("topGOdata",
              ontology = "BP",
              geneSelectionFun = function(x){x==1},
              allGenes = input,
              # nodeSize = 10,
              annot = annFUN.org,
              mapping = "org.Mm.eg.db",
              ID = "ensembl")
  
  ks.test = suppressMessages(
    runTest(topgo, algorithm = "elim", statistic = "ks")
    )
  out = GenTable(topgo, ks.test, topNodes = n.out)
  return(out)
}

#use ensembl gene IDs
enrich_genes_goana = function(hits, universe, warn_missing = TRUE){
  require(limma)
  require(org.Mm.eg.db)
  
  db = org.Mm.egENSEMBL
  mappedgenes = mappedkeys(db)
  mapping = unlist(as.list(db[mappedgenes]))
  
  hits = as.character(hits)
  universe = as.character(universe)
  
  hits_missing = hits[which(!hits %in% mapping)]
  universe_missing = universe[which(!universe %in% mapping)]
  
  if(length(universe)>0 & warn_missing){
    warning("Some genes do not map to Entrez identifiers")
  }
  
  hits_sub = hits[which(hits %in% mapping)]
  universe_sub = universe[which(universe %in% mapping)]
  
  hits_entrez = names(mapping)[match(hits_sub, mapping)]
  universe_entrez = names(mapping)[match(universe_sub, mapping)]
 
  out = goana(de = hits_entrez,
              universe = universe_entrez,
              species = "Mm")
  
  out = out[order(out$P.DE),]
  
  return(list(result = out,
              excluded = list(hits = hits_missing,
                              universe = universe_missing)))
   
}

#MAPPING FUNCTIONS

getmode <- function(v, dist) {
  tab = table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied = names(tab)[tab == max(tab)]
    sub = dist[v %in% tied]
    names(sub) = v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

mnnMap = function(atlas_pca, atlas_meta, map_pca, map_meta, k_map = 10){
  require(kmknn)
  require(scran)
  
  correct = fastMNN(atlas_pca, map_pca, pc.input = TRUE)$corrected
  atlas = 1:nrow(atlas_pca)
  correct_atlas = correct[atlas,]
  correct_map = correct[-atlas,]
  
  knns = kmknn::queryKNN(correct_atlas, correct_map, k = k_map, get.index = TRUE, get.distance = FALSE)
  
  #get closest k matching cells
  k.mapped = t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  celltypes = t(apply(k.mapped, 1, function(x) atlas_meta$celltype[match(x, atlas_meta$cell)]))
  stages = t(apply(k.mapped, 1, function(x) atlas_meta$stage[match(x, atlas_meta$cell)]))
  celltype.mapped = apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
  stage.mapped = apply(stages, 1, function(x) getmode(x, 1:length(x)))
  
  out = lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped = k.mapped[x,],
         celltype.mapped = celltype.mapped[x],
         stage.mapped = stage.mapped[x])
  })
  
  names(out) = map_meta$cell
  
  return(out)
  
}

#meta MUST have a cell column, celltype column and a stage column, spelled exactly like that
mapWrap = function(atlas_sce, atlas_meta, map_sce, map_meta, k = 10, return.list = FALSE){
  require(scran)
  #prevent duplicate rownames
  colnames(map_sce) = paste0("map_", colnames(map_sce))
  map_meta$cell = paste0("map_", map_meta$cell)
  
  big_sce = scater::normalize(cbind(atlas_sce, map_sce))
  hvgs = getHVGs(big_sce)
  big_pca = prcomp_irlba(t(logcounts(big_sce[hvgs,])), n = 50)$x
  rownames(big_pca) = colnames(big_sce) 
  atlas_pca = big_pca[1:ncol(atlas_sce),]
  map_pca = big_pca[-(1:ncol(atlas_sce)),]
  
  #correct the atlas first
  order_df = atlas_meta[!duplicated(atlas_meta$sample), c("stage", "sample")]
  order_df$ncells = sapply(order_df$sample, function(x) sum(atlas_meta$sample == x))
  order_df$stage = factor(order_df$stage, 
                          levels = rev(c("E8.5", 
                                         "E8.25", 
                                         "E8.0", 
                                         "E7.75", 
                                         "E7.5", 
                                         "E7.25", 
                                         "mixed_gastrulation", 
                                         "E7.0", 
                                         "E6.75", 
                                         "E6.5")))
  order_df = order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
  order_df$stage = as.character(order_df$stage)
  
  set.seed(42)
  atlas_corrected = doBatchCorrect(counts = logcounts(atlas_sce[hvgs,]), 
                                   timepoints = atlas_meta$stage, 
                                   samples = atlas_meta$sample, 
                                   timepoint_order = order_df$stage, 
                                   sample_order = order_df$sample, 
                                   pc_override = atlas_pca)
  
  mapping = mnnMap(atlas_pca = atlas_corrected,
                   atlas_meta = atlas_meta,
                   map_pca = map_pca,
                   map_meta = map_meta)
  if(return.list){
    return(mapping)
  }
  
  cn = substr(names(mapping), 5, nchar(names(mapping))) # remove 
  ct = sapply(mapping, function(x) x$celltype.mapped)
  st = sapply(mapping, function(x) x$stage.mapped)
  closest = sapply(mapping, function(x) x$cells.mapped[1])
  
  out = data.frame(cell = cn, celltype.mapped = ct, stage.mapped = st, closest.cell = closest)
  return(out)
  
}
