#############################################
#title: "10X embryo data - Characterising the haemato-endothelial landscape"
#author: "Blanca Pijuan-Sala"
#date: "09 September 2018"
#############################################
wd <- "/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/My_papers_submitted/20180601_PijuanSalaEtAl_Gastrulation10X/02.Review_01/SCRIPTS/"


library(anSeq)##You can find the anSeq package here: https://github.com/BPijuanSala/anSeq
############################################
# FUNCTIONS
############################################



heatmapRedYelBlue <- c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090",
                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")

palette <- colorRampPalette(rev(heatmapRedYelBlue))
bluePal <- c("#BFBFBF","#6495ED","#000000")
redPal <- c("gray","#ff6666","#e60000","#990000")
all_colours = c(
  "Allantois" = "#532C8A",#[32] "Allantois" 
  "Anterior Primitive Streak" = "#c19f70",
  "Blood progenitors 1" = "#f9decf",
  "Blood progenitors 2" = "#c9a997",
  "Cardiomyocytes" =  "#B51D8D",#[34] "Cardiomyocytes"     
  "Caudal epiblast" = "#9e6762",
  "Caudal Mesoderm" = "#3F84AA",
  #"Caudal Primitive Streak epiblast"= "#702f1f",
  "Def. endoderm" = "#F397C0",#[24] "Def. endoderm"    
  "Nascent mesoderm" =  "#C594BF",#[7] "Early mixed mesoderm"                        
  "Mixed mesoderm" =  "#DFCDE4",#[26] "Early ExE mesoderm"    
  
  "Endothelium" =  "#eda450",#[20] "Endothelium"                                 
  "Epiblast" =  "#635547",#[1] "Epiblast"  
  "Erythroid1" =  "#C72228",#[15] "Erythroid 1"                                 
  "Erythroid2" =  "#EF4E22",#[37] "Erythroid2"     
  "Erythroid3" =  "#f77b59",
  "ExE ectoderm" =  "#989898",#[30] "ExE ectoderm 1"                              
  
  "ExE endoderm" = "#7F6874",#[5] "ExE endoderm"                                
  "ExE mesoderm" =  "#8870ad",#[12] "ExE mesoderm" 
  
  
  "Rostral neurectoderm" =  "#65A83E",#[8] "Forebrain"   
  "Forebrain/Midbrain/Hindbrain" = "#647a4f",
  
  
  "Gut" =  "#EF5A9D",#[19] "Foregut"                                     
  "Haematoendothelial progenitors" =  "#FBBE92",#[9] "Hemato-endothelial progenitors"              
  "Caudal neurectoderm"= "#354E23",
  
  
  "Intermediate mesoderm" =  "#139992",#[31] "Intermediate mesoderm"                       
  "Neural crest"= "#C3C388",
  
  "NMP" =  "#8EC792",#[14] "NMPs"                                         
  "Notochord" =  "#0F4A9C",#[21] "Notochord"                                   
  "Paraxial mesoderm" =  "#8DB5CE",#[33] "Late paraxial mesoderm (presomitic mesoderm)"
  "Parietal endoderm" =  "#1A1A1A",#[29] "Parietal endoderm"                           
  "PGC" =  "#FACB12",#[25] "PGC"                                         
  
  "Pharyngeal mesoderm" =  "#C9EBFB",#[13] "Late mixed mesoderm"                         
  "Primitive Streak" =  "#DABE99",#[2] "Primitive Streak"    
  "Mesenchyme" = "#ed8f84",
  "Somitic mesoderm" =  "#005579",#[16] "Somites"                                     
  "Spinal cord" =  "#CDE088",#[38] "Spinal cord"                                 
  "Surface ectoderm" = "#BBDCA8",#[22] "Surface ectoderm"                            
  
  
  "Visceral endoderm" = "#F6BFCB"#[3] "Visceral endoderm"   
  

)


spectralPal = c(
  'E6.5'="#D53E4F",
  'E6.75'="#F46D43",
  'E7.0'="#FDAE61",
  'E7.5'="#FFFFBF",
  'E7.25'="#FEE08B",
  'E7.75'="#E6F598",
  'E8.0'="#ABDDA4",
  'E8.5'="#3288BD",
  'E8.25'="#66C2A5",
  'mixed_gastrulation'= "#A9A9A9"  
  
)

all_colours_sub = c(
  "Mes1"= "#c4a6b2",#
  "Mes2"= "#ca728c",#
  
  
  "BP1" = "#6460c5",#
  "BP2" = "#96b8e4",#
  "BP3" = "#07499f",#changed
  "BP4" = "#036ef9",#changed
  #"BP6"  = "#03bbf9",
  
  "Haem1"= "#bb22a7",
  "Haem2" = "#f695e9",
  "Haem3"= "#02f9ff",#changed
  
  "Haem4" = "#4c4a81",#changed
  
  "EC1"= "#006737",#
  
  "EC2" = "#5eba46",#
  "EC3" = "#818068",#
  "EC4"="#d6de22",#
  "EC5"="#5f7238",#
  "EC6"="#2ab572",#
  "EC7"="#000000",#
  "EC8"="#a0cb3b",#
  
  "Ery1"="#f67a58",#
  "Ery2" ="#a26724",#
  "Ery3"="#cdaf7f",#
  "Ery4"= "#625218",#
  
  "My" = "#c62127" ,#
  "Mk"= "#f6931d"#
  
  
)


plotGeneLevelsLocal <- function(data, x, y, gene, cols=c("#BFBFBF","#6495ED","#000000"),
                                xlab="x",ylab="y",titlePlot=gene,cexType=1,ybsub=0.1){
  
  redRamp <- colorRampPalette(cols)
  df <- data.frame(x = x, y = y, exp = data[gene,])
  df <- df[order(df$exp,decreasing=F),]
  dfsub <- df[df$exp>0,]
  
  interval <- findInterval(dfsub$exp, seq(min(dfsub$exp), 
                                          max(dfsub$exp), 
                                          (max(dfsub$exp)-min(dfsub$exp))/10))
  
  interval[interval==0]<-1
  colorsPlot <- redRamp(11)[interval]
  
  par(mar=c(8,4,8,4),xpd=NA)
  plot(df$x, df$y, col=cols[1], pch=20, cex=cexType,
       xlab="", ylab="", main=titlePlot, axes=F, cex.main=1.5)
  box(bty="l")
  points(dfsub$x, dfsub$y, pch=20, cex=cexType, 
         col=colorsPlot)
  
  
}
####====================
# Load data
####====================

###Read metadata with the cell type annotation
meta <- read.table(paste0(wd,"DATA/metadata_mergedClustering_Celltypes_20180910.tab"), header = T, stringsAsFactors = F, sep = "\t")    
rownames(meta) <- meta$index
colnames(meta)[1] <- 'cell'

###Load the counts of haemato-endothelial related clusters 

counts <- as.matrix(Matrix::readMM(file=paste0(wd,"DATA/20180908_countsBlood.mtx")))

colnames(counts) <- as.character(read.table(paste0(wd,"DATA/20180908_cells_countsBlood.tab"),header=T)[,1])
rownames(counts) <-  read.table(paste0(wd,"DATA/20180908_genes_countsBlood.tab"),header=TRUE)[,1]

###Read metadata for the haemato-endothelial related clusters 

metaSub <- meta[colnames(counts),]
rownames(metaSub) <- metaSub$cell

metaSub <- read.table(paste0(wd,"DATA/metadata_Hematoendothelium_subclustering.txt"),
                      header = T, stringsAsFactors = F, sep = "\t")    
rownames(metaSub) <- metaSub$cell



## Force directed graph has been performed in Gephi v.0.9.2
gephi <- read.table(paste0(wd,"DATA/blood_graph_coords_v7_20180908_mod.gdf"),
                    sep=',')
gephi <- gephi[order(gephi$V1),]

metaSub$gephiX <- gephi$V5
metaSub$gephiY <- gephi$V6


###################=================================
# gephi representation
###################=================================
par(mfrow=c(1,1))

plot(metaSub$gephiX,metaSub$gephiY,
     col=all_colours[as.character(metaSub$celltype_new)],
     pch=20, cex=0.35,xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
box(bty="l")





plot(metaSub$gephiX,metaSub$gephiY,
     col=spectralPal[as.character(metaSub$stage)],
     pch=20, cex=0.35,xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
box(bty="l")




genestoPlot <-c("Kdr","Lmo2","Tal1",
                "Etv2","Itga2b","Hoxa11",
                "Hoxa13","Cxcr4","Hbb-bh1","Hbb-y",
                "Hbb-bt",  "Hbb-bs",  "Hbb-bh3", "Hbb-bh2", "Hbb-bh1", "Hbb-bh0", "Hbb-y",
                "Alox5ap","Spi1","Nfe2","Lyve1",
                "Mpl","Csf1r","Ltc4s","Pf4","Fcgr3","Meis1","Cdh5","Pecam1",
                "F2rl2",  "Gp9","Sox17","Cdx4","Hlf","Hoxa3","Hoxa5","Hoxa7",
                "Hoxa9",
                "Hoxa10",
                "Sox17","Efna1","Runx1") 

countsPlot <- counts[,rownames(metaSub)]


for (i in 1:length(genestoPlot)){
  geneName <-genestoPlot[i]
  tissue <- "blood"
  
  gene <- anSeq::getGeneID(geneName)$geneIDs[1]
  if (gene %in% rownames(countsPlot)){
    if (sum(countsPlot[gene,])>0){

      plotGeneLevelsLocal(countsPlot,x=metaSub$gephiX,metaSub$gephiY,
                          gene=gene,ybsub=0.25,cexType = 1,
                          titlePlot=geneName,
                          
                          cols=bluePal)

    }
  }
}


cols=bluePal
redRamp <- colorRampPalette(cols)

plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='',ylim=c(0,0.5),xlim=c(0,2))

xl <- 0.2; yb <- 0; xr <- 1.5; yt <- 0.2
rect(xleft=head(seq(xl,xr,(xr-xl)/10),-1), ybottom = yb, 
     xright = tail(seq(xl,xr,(xr-xl)/10),-1), ytop = yt,
     col=redRamp(10), border=redRamp(10), lwd=0.5, xpd=NA)


####====================
# Batch correction
####====================

#functions
getHVGs = function(counts, min.mean = 1e-3, gene_df = genes){
  require(biomaRt)
  trend = scran::trendVar(counts,
                          loess.args = list(span = 0.05))
  decomp = scran::decomposeVar(counts, fit = trend)
  decomp = decomp[decomp$mean > min.mean,]
  #exclude sex genes
  xist = "ENSMUSG00000086503"
  mouse_ensembl = useMart("ensembl")
  mouse_ensembl = useDataset("mmusculus_gene_ensembl", mart = mouse_ensembl)
  # gene_map = getBM(attributes=c("ensembl_gene_id", "chromosome_name"), filters = "ensembl_gene_id", values = rownames(decomp), mart = mouse_ensembl)
  # ychr = gene_map[gene_map[,2] == "Y", 1]
  ychr = read.table(paste0(wd,"DATA/mmusculus_Y_genes.txt"), stringsAsFactors = FALSE)[,1]
  #other = c("tomato-td") #for the chimera
  decomp = decomp[!rownames(decomp) %in% c(xist, ychr),]
  decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.05])
}

#ensure counts has columns names for the cells
#match timepoints,samples to the count table
#timepoint_order, sample_order should contain each sample/timepoint ONCE, in correct order for correction
doBatchCorrect = function(counts, timepoints, samples, timepoint_order, sample_order, npc = 50, BPPARAM = SerialParam()){
  require(scran)
  require(irlba)
  require(BiocParallel)
  #compute pca
  pca = prcomp_irlba(t(counts), n = npc)$x
  rownames(pca) = colnames(counts)
  if(length(unique(samples)) == 1){
    return(pca)
  }
  #create nested list
  pc_list = lapply(unique(timepoints), function(tp){
    sub_pc = pca[timepoints == tp, ]
    sub_samp = samples[timepoints == tp]
    list = lapply(unique(sub_samp), function(samp){
      sub_pc[sub_samp == samp, ]
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



hvgs = getHVGs(counts)


#get order: oldest to youngest; most cells to least cells
order_df = metaSub[!duplicated(metaSub$sample), c("stage", "sample")]
order_df$ncells = sapply(order_df$sample, function(x) sum(metaSub$sample == x))
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



cat("Batch correction...")

all_correct = doBatchCorrect(counts = counts[rownames(counts) %in% hvgs,],
                             timepoints = metaSub$stage,
                             samples = metaSub$sample,
                             timepoint_order = order_df$stage,
                             sample_order = order_df$sample,
                             npc = 50)
dim(all_correct)


####====================
# Louvain clustering
####====================
#First subclustering of the data 

PCAbatchCorrected <- as.matrix(Matrix::readMM(paste0(wd,"DATA/20180908_countsBlood_PCAcorrected.mtx")))
rownames(PCAbatchCorrected)<-rownames(metaSub)


graph <- scran::buildKNNGraph(as.matrix(t(PCAbatchCorrected)),k=10)


library(igraph)
library(moduleColor)

louv <- igraph::cluster_louvain(graph, weights = NULL)

clustersLouv <- igraph::membership(louv)
names(clustersLouv) <-rownames(PCAbatchCorrected)

library(moduleColor)
clustColLouv<-labels2colors(clustersLouv)
names(clustColLouv)<-names(clustersLouv)


metaSub$subclust <- metaSub$celltype_new
metaSub[names(clustColLouv),'subclust'] <-clustColLouv


metaSub$subclustCol = metaSub$subclust
colorssubClust <- c("turquoise"="BP1",
                    "pink"="Mes2", 
                    "tan" ="BP2",
                    "green"="Mes1",  
                    "brown" ="Haem2",   
                    "purple"="Haem1" , 
                    "black" = "EC3",
                    "yellow" ="EC2",
                    "greenyellow" = "Haem3",
                    "red" = "EC1", 
                    "blue"= "Ery1", 
                    "salmon"="Ery2",
                    "magenta"="Ery3", 
                    "cyan"="Ery4")
for (i in 1:length(colorssubClust)){
  col <- names(colorssubClust)[i]
  lab <- colorssubClust[i]
  
  metaSub[,"subclust"][metaSub$subclust==col] <- lab
  
}




plot(metaSub$gephiX,metaSub$gephiY,col=all_colours_sub[metaSub$subclust],
     
     pch=20, cex=0.35,xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
box(bty="l")





####====================
# Louvain clustering for heterogeneous cluster containing My and Mk
####====================

cellsBranch <- rownames(metaSub)[metaSub$subclust%in%c("Haem3")]

PCAbatchCorrectedHemendo <- PCAbatchCorrected[cellsBranch,]

graph <- scran::buildKNNGraph(as.matrix(t(PCAbatchCorrectedHemendo)),k=15)


library(igraph)
library(moduleColor)

louv <- igraph::cluster_louvain(graph, weights = NULL)

clustersLouv <- igraph::membership(louv)
names(clustersLouv) <-rownames(PCAbatchCorrectedHemendo)

library(moduleColor)
clustColLouv<-labels2colors(clustersLouv)
names(clustColLouv)<-names(clustersLouv)

par(mfrow=c(1,1))

plot(metaSub$gephiX,metaSub$gephiY,pch=20,col="gray",
     main="",xlab="",ylab="")
points(metaSub[names(clustColLouv),"gephiX"],metaSub[names(clustColLouv),"gephiY"],pch=20,
       col=clustColLouv
       ,
       main="",xlab="",ylab="")

par(mfrow=c(2,4))

for (i in unique(clustColLouv)){
  namesClust <- names(clustColLouv)[clustColLouv==i]
  plot(metaSub$gephiX,metaSub$gephiY,pch=20,col="gray",
       main="",xlab="",ylab="")
  points(metaSub[namesClust,"gephiX"],metaSub[namesClust,"gephiY"],pch=20,
         col=i
         ,
         main="",xlab="",ylab="")
  
}

table(clustColLouv)
#blue     brown     green       red turquoise    yellow 
#96       102        32        56       139        54 

#Red Macrophages
#green Mk


metaSub$subclust2 <- metaSub$subclust
metaSub[names(clustColLouv),"subclust2"]<-clustColLouv
metaSub$subclust2Col <- metaSub$subclust2
colorSub2 <- c(
  'green'="Mk",
  'yellow'="BP4",
  'blue'="BP3",
  'red'="My",
  'turquoise'="Haem3",
  'brown' ="Haem4"
  
)


for (i in 1:length(colorSub2)){
  col <- names(colorSub2)[i]
  lab <- colorSub2[i]
  
  metaSub[,"subclust2"][metaSub$subclust2==col] <- lab
  
}


par(mfrow=c(1,1))
plot(metaSub$gephiX,metaSub$gephiY,col=all_colours_sub[metaSub$subclust2],
     
     pch=20, cex=0.35,xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
box(bty="l")




####====================
# Louvain clustering for original EC3 (heterogeneous)
####====================

cellsBranch <- rownames(metaSub)[metaSub$subclust%in%c("EC3")]

PCAbatchCorrectedHemendo <- PCAbatchCorrected[cellsBranch,]

graph <- scran::buildKNNGraph(as.matrix(t(PCAbatchCorrectedHemendo)),k=30)


library(igraph)
library(moduleColor)

louv <- igraph::cluster_louvain(graph, weights = NULL)

clustersLouv <- igraph::membership(louv)
names(clustersLouv) <-rownames(PCAbatchCorrectedHemendo)

library(moduleColor)
clustColLouv<-labels2colors(clustersLouv)
names(clustColLouv)<-names(clustersLouv)

par(mfrow=c(1,1))

plot(metaSub$gephiX,metaSub$gephiY,pch=20,col="gray",
     main="",xlab="",ylab="")
points(metaSub[names(clustColLouv),"gephiX"],metaSub[names(clustColLouv),"gephiY"],pch=20,
       col=clustColLouv
       ,
       main="",xlab="",ylab="")

par(mfrow=c(2,4))

for (i in unique(clustColLouv)){
  namesClust <- names(clustColLouv)[clustColLouv==i]
  plot(metaSub$gephiX,metaSub$gephiY,pch=20,col="gray",
       main="",xlab="",ylab="")
  points(metaSub[namesClust,"gephiX"],metaSub[namesClust,"gephiY"],pch=20,
         col=i
         ,
         main="",xlab="",ylab="")
  
}

table(clustColLouv)
#blue     brown     green       red turquoise    yellow 
#338       170       164       174       197       169 




metaSub$subclust3 <- metaSub$subclust2
metaSub[names(clustColLouv),"subclust3"]<-clustColLouv
metaSub$subclust3Col <- metaSub$subclust3
colorSub3 <- c(
  'blue'="EC3",
  
  'green'="EC4",
  'yellow'="EC5",
  'red'="EC8",
  'brown' ="EC6",
  'turquoise'="EC7"
  
)


for (i in 1:length(colorSub3)){
  col <- names(colorSub3)[i]
  lab <- colorSub3[i]
  
  metaSub[,"subclust3"][metaSub$subclust3==col] <- lab
  
}



par(mfrow=c(1,1))
plot(metaSub$gephiX,metaSub$gephiY,col=all_colours_sub[as.character(metaSub$subclust3)],
     
     pch=20, cex=0.35,xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
box(bty="l")







####====================
# Estimation of My and Mk prevalence
####====================


meta <- read.table(paste0(wd,"DATA/metadata_mergedClustering_Celltypes_20180910.tab"), header = T, stringsAsFactors = F, sep = "\t")    
rownames(meta) <- meta$index
colnames(meta)[1] <- 'cell'


metaSub <- read.table(file=paste0(wd,"DATA/metadata_subclusters_bloodLineage_v3_20180910.tab"),
                      sep="\t",header=T)
rownames(metaSub) <- metaSub$cell



cellsMy <- rownames(metaSub)[metaSub$subclust3=="My"]
cellsMk <- rownames(metaSub)[metaSub$subclust3=="Mk"]

tab <- table(meta[cellsMy,"stage"])
numStage <- table(meta$stage)[c('E8.0','E8.25','E8.5')]
mean(tab/numStage)*100

tab <- table(meta[cellsMk,"stage"])
numStage <- table(meta$stage)[c('E8.25','E8.5')]
mean(tab/numStage)*100

####====================
#  Markers of My and Mk
####====================

HVG <- as.character(read.table(paste0(wd,"DATA/20180908_HVGs_countsBlood.tab"),header=T)[,1])

groups <- metaSub$subclust3
names(groups)<-rownames(metaSub)
markers <- anSeq::findMarkers(counts[HVG,],groups=groups,
                              topGeneThres=50)

for (i in 1:length(markers)){
  cat(dim(i),"\n")
  markers[[i]]$geneNames <- anSeq::getGeneName(rownames(markers[[i]]))[[2]]
}


logFCThres <- 2.5
signqVal <- 0.05

MySign <-  rownames(markers[["My"]])[(markers[["My"]]$logFC>logFCThres & markers[["My"]]$fdr < signqVal)]
MySignName <- anSeq::getGeneName(MySign)[[2]]


MkSign <-  rownames(markers[["Mk"]])[(markers[["Mk"]]$logFC>logFCThres & markers[["Mk"]]$fdr < signqVal)]
MkSignName <- anSeq::getGeneName(MkSign)[[2]]

EC7Sign <-  rownames(markers[["EC7"]])[(markers[["EC7"]]$logFC>logFCThres & markers[["EC7"]]$fdr < signqVal)]
EC7SignName <- anSeq::getGeneName(EC7Sign)[[2]]


BP4Sign <-  rownames(markers[["BP4"]])[(markers[["BP4"]]$logFC>logFCThres & markers[["BP4"]]$fdr < signqVal)]
BP4SignName <- anSeq::getGeneName(BP4Sign)[[2]]


Haem4Sign <-  rownames(markers[["Haem4"]])[(markers[["Haem4"]]$logFC>logFCThres & markers[["Haem4"]]$fdr < signqVal)]
Haem4SignName <- anSeq::getGeneName(Haem4Sign)[[2]]

Myunique <- setdiff(MySign,MkSign)
Mkunique <- setdiff(MkSign,MySign)
EC7unique <- setdiff(setdiff(EC7Sign,MySign),Mkunique)
BP4unique <- setdiff(setdiff(setdiff(BP4Sign,Myunique),Mkunique),EC7unique)
Haem4unique <- setdiff(setdiff(setdiff(setdiff(Haem4Sign,Myunique),Mkunique),EC7unique),BP4unique)

Myunique <- Myunique[Myunique%in%rownames(counts)]
Mkunique <- Mkunique[Mkunique%in%rownames(counts)]
EC7unique <- EC7unique[EC7unique%in%rownames(counts)]
BP4unique <- BP4unique[BP4unique%in%rownames(counts)]
Haem4unique <- Haem4unique[Haem4unique%in%rownames(counts)]


cols0 <- c(rep(all_colours_sub["Mk"],length(Mkunique)),
           rep(all_colours_sub["My"],length(Myunique)),
           rep(all_colours_sub["EC7"],length(EC7unique)),
           rep(all_colours_sub["BP4"],length(BP4unique)),
           rep(all_colours_sub["Haem4"],length(Haem4unique)))


labelsCells0 <- as.character(metaSub$subclust3)[as.character(metaSub$subclust3)%in%c("My","Mk",
                                                                                     'EC7','BP4',
                                                                                     "Haem4")]

names(labelsCells0) <- rownames(metaSub)[as.character(metaSub$subclust3)%in%c("My","Mk",'EC7',"BP4",'Haem4')]


labelsCells0 <- sort(labelsCells0)
labelsCells <- c()

orderHeat <- c("EC7","Haem4","My","BP4","Mk")
for (i in orderHeat){
  labelsCells <- c(labelsCells,labelsCells0[labelsCells0==i])
}


genesSel0 <- cols0
names(genesSel0)<-c(Mkunique,Myunique,EC7unique,BP4unique,Haem4unique)


genesSel <- genesSel0

cellsSel <- names(labelsCells) 

heatCounts <- counts[names(genesSel),cellsSel]
cols<-genesSel

heatCountsVal <- heatCounts[apply(heatCounts,1,function(x){sum(x)>0}),]

heatCountsStd <- t(apply(heatCounts, 1, function(x) x/max(x)))  # standarise


g <- anSeq::getGeneName(rownames(heatCounts))[[2]]
c <- all_colours_sub[labelsCells]


gplots::heatmap.2(heatCountsStd, trace="none", 
                  col=palette, 
                  Colv = F, Rowv = T, 
                  ColSideColors = c, 
                  #RowSideColors = cols, 
                  dendrogram = "none", density.info = 'none',
                  labRow = g,key=FALSE,labCol="")



####====================
#  Plot genes selected for heatmap
####====================


genesSel1 <- c("Gap43","Icam2","Anxa12","Gpr182","Cldn5","Igf1","Oit3","Gimap1","Fxyd5","Cd34","Selenop","Plvap","Lyve1",
               "Gp5","Thbs1","Gp1bb","Treml1","Gp9","F2rl2","Sla","Slc35d3","Gnaz","Itgb3","Dok2","Mpl","Rgs18",
               "Bin2","Lat","Lyn","Timp3","Clec1b","Rab27b","Plek","Pf4","Mef2c","Mfng","Nrros","Gimap5",
               "Hcls1","Lyz2","Fcgr3","Ncf2","Celf2","Coro1a","Tyrobp","Ptpn7","Spi1","Alox5ap","Fcer1g")

heatCounts <- counts[anSeq::getGeneID(genesSel1)[[2]],cellsSel]


heatCountsVal <- heatCounts[apply(heatCounts,1,function(x){sum(x)>0}),]

heatCountsStd <- t(apply(heatCounts, 1, function(x) x/max(x)))  # standarise
g <- anSeq::getGeneName(rownames(heatCounts))[[2]]
c <- all_colours_sub[labelsCells]




gplots::heatmap.2(heatCountsStd, trace="none", 
                  col=palette, 
                  Colv = F, Rowv = T, 
                  ColSideColors = c, 
                  #RowSideColors = cols, 
                  dendrogram = "none", density.info = 'none',
                  labRow = g,key=FALSE,labCol="")







####====================
#  Microglial genes plot
####====================

genesPlot <- c("Cx3cr1","Adgre1","Tmem119","Fcgr3","Kit","Ptprc","Csf1r")
cols <- c(rep("turquoise",3),rep("gold",4))


labelsCells <- as.character(metaSub$subclust3)[as.character(metaSub$subclust3)%in%c("My")]
names(labelsCells) <- rownames(metaSub)[as.character(metaSub$subclust2)%in%c("My")]
labelsCells <- sort(labelsCells)

cellsSel <- names(labelsCells) 

heatCounts <- counts[anSeq::getGeneID(genesPlot)[[2]],cellsSel]



heatCountsVal <- heatCounts[apply(heatCounts,1,function(x){sum(x)>0}),]

heatCountsStd <- t(apply(heatCounts, 1, function(x) x/(max(x)+0.00001)))  # standarise

g <- getGeneName(rownames(heatCounts))[[2]]
c <- all_colours_sub[labelsCells]
pdf(paste0(wd,"PhD_BPS32/release6/plots/blood_lineage/Embryo10Xv6_microglia_EMP.pdf"),
    width=10,height=3)
gplots::heatmap.2(heatCounts, trace="none", 
                  col=palette, 
                  Colv = F, Rowv = F, 
                  ColSideColors = c, 
                  RowSideColors = cols, 
                  dendrogram = "none", density.info = 'none', labRow = g,key=FALSE,labCol="",
                  margins=c(5,10))
dev.off()





####====================
# Fraction EC3 in stages
####====================
sum(apply(table(metaSub$stage,metaSub$subclust3),2,function(x){x/sum(x)})[c("E7.75","E8.0","E8.25"),"EC3"])
apply(table(metaSub$stage,metaSub$subclust3),2,function(x){x/sum(x)})[,"EC3"]
