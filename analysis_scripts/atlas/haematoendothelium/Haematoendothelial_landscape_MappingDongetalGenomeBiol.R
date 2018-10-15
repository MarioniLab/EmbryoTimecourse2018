#############################################
#title: "10X embryo data - Mapping to Dong et al Genome Biology"
#author: "Blanca Pijuan-Sala"
#date: "19 July 2018"
#############################################
wd <- "/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/Experiments/"



library(Matrix)
library(anSeq)##You can find the anSeq package here: https://github.com/BPijuanSala/anSeq
library(proSeq)##You can find the proSeq package here: https://github.com/BPijuanSala/anSeq

####################
# Palettes
##################

heatmapRedYelBlue <- c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090",
                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")
palette <- colorRampPalette(rev(heatmapRedYelBlue))


DongHCcolor <- c(
  "Definitive erythroid cells"="orange",
  "Primitive erythrocytes"="brown",
  "EMPs1"="blue",
  "EMPs2"="cyan",
  "Primitive macrophages"="purple",
  "Myeloid_progenitors2"="black",
  "Myeloid progenitors"="grey",
  "EMPs3"="turquoise"
  
)


DongHCcolor2 <- c(
  "Definitive erythroid cells"="orange",
  "Primitive erythrocytes"="brown",
  "EMPs1"="blue",
  "EMPs2"="blue",
  "Primitive macrophages"="purple",
  "Myeloid_progenitors2"="black",
  "Myeloid progenitors"="black",
  "EMPs3"="blue"
  
)
bluePal <- c("#BFBFBF","#6495ED","#000000")
redPal <- c("gray","#ff6666","#e60000","#990000")



###------------------------------------------------------
##Load data 10X data for haemato-endothelial lineage
###------------------------------------------------------

counts <- as.matrix(Matrix::readMM(file=paste0(wd,"PhD_BPS32/release6/data/20180908_countsBlood.mtx")))

colnames(counts) <- as.character(read.table(paste0(wd,"PhD_BPS32/release6/data/20180908_cells_countsBlood.tab"),header=T)[,1])
rownames(counts) <-  read.table(paste0(wd,"PhD_BPS32/release6/data/20180908_genes_countsBlood.tab"),header=TRUE)[,1]


metaSub <- read.table(file=paste0(wd,"PhD_BPS32/release6/data/metadata_subclusters_bloodLineage_v3_20180910.tab"),
                      sep="\t",header=T)
rownames(metaSub) <- metaSub$cell


###------------------------------------------------------
##Subset myeloid cells
###------------------------------------------------------

cellsMy <- rownames(metaSub)[as.character(metaSub$subclust3)%in%c("My")]

countsMy <- counts[,cellsMy]

metaMy <- metaSub[cellsMy,]


###------------------------------------------------------
##Load data from Dong et al Genome Biology 2018 
###------------------------------------------------------

countsDong <- read.table(file=paste0(wd,"PhD_BPS50/data/Dongetal_GenomeBiol_2018/counts_Dong_TPM.txt"),
              header=T,sep="\t")
metaDong <- read.table(file=paste0(wd,"PhD_BPS50/data/Dongetal_GenomeBiol_2018/meta_DongEtAl.txt"),header=TRUE,
            sep="\t")


genesIDDong <- anSeq::getGeneID(rownames(countsDong))[[1]]
countsHCF1Dong <- countsDong[rownames(countsDong) %in% genesIDDong$Associated.Gene.Name,]

genesIDF1Dong <- genesIDDong[genesIDDong$Associated.Gene.Name%in%rownames(countsHCF1Dong),]
genesIDF1Dong<- genesIDF1Dong[match(genesIDF1Dong$Associated.Gene.Name,rownames(countsHCF1Dong)),]
dupDong<-duplicated(genesIDF1Dong$Associated.Gene.Name)
genesIDF2Dong <- genesIDF1Dong[!dupDong,]

countsHCF2Dong <- countsHCF1Dong
countsHCF2Dong <- countsHCF2Dong[as.character(genesIDF2Dong$Associated.Gene.Name),]
rownames(countsHCF2Dong) <- as.character(genesIDF2Dong[,"Gene.ID"])

###------------------------------------------------------
##Take Haematopoietic cluster (HC)
###------------------------------------------------------

metaHC <- metaDong[metaDong$Cluster=="HC",]

countsHCF2 <- log2((countsHCF2Dong[,rownames(metaHC)]/10) +1)

###------------------------------------------------------
##Find HVG
###------------------------------------------------------

HVG <- proSeq::findHVGMatrix(as.matrix(countsHCF2),UseSpike = F)

###------------------------------------------------------
##Compute PCA to find structure
###------------------------------------------------------


countsPCA <- prcomp(t(countsHCF2[HVG,]))

countsPCAx <- countsPCA$x[,1:50]
rownames(countsPCAx)<-colnames(countsHCF2)

plot(countsPCAx[,1],countsPCAx[,2],pch=20,col=metaHC$Organ)



###------------------------------------------------------
## Cluster cells
###------------------------------------------------------

graph <- scran::buildKNNGraph(as.matrix(t(countsPCAx)),k=10)


library(igraph)
library(moduleColor)

louv <- igraph::cluster_louvain(graph, weights = NULL)

clustersLouv <- igraph::membership(louv)
names(clustersLouv) <-rownames(countsPCAx)

library(moduleColor)
clustColLouv<-labels2colors(clustersLouv)
names(clustColLouv)<-names(clustersLouv)
#pdf(paste0(wd,"PhD_BPS50/plots/PCA_clustLouv_Dong.pdf"),
#    width=7.5,height=8)
plot(countsPCAx[,1],countsPCAx[,2],pch=20,col=clustColLouv,
     main="Stage",xlab="PC1",ylab="PC2")
#dev.off()

metaHC$clustColLouv <- as.character(clustColLouv[rownames(metaHC)])


###------------------------------------------------------
## Take genes they use to classify their populations in Fig. 8
###------------------------------------------------------

genestoPlot <-c("Csf1r","Cx3cr1","Cst3","Tyrobp","Cstc","Lcp1","B2m","Unc93b1","Cyth4","Capza2",
                "Itga2b","Elf1","Dapp1","Rgs18","Tagln2","Rasgrp2","Plek","Spns1","Adgrg1","Srgn","Muc13",
                "Npm1","Nolc1","Hspd1","Casp3","Ncl","Snhg5","Gm15915","Ddx21","Pgm1","Zmynd8",
                "Rgcc","Cd24a","Ubac1","Cyb5r3","Gid8","Icam4","Nxpe2","Rhd",
                
                "Hba-x","Hba-y","Hbb-bh1","Nudt4","Bpgm","Alas2","Sepw1","Slc39a8","H19","Ube2l6",
                "Kit",
                "Sox6"
) 

genesSelectedFirst <- anSeq::getGeneID(genestoPlot)[[2]]
genesSelectedFirst<- genesSelectedFirst[genesSelectedFirst%in%rownames(countsMy)]
genesSelectedFirst<- genesSelectedFirst[genesSelectedFirst%in%rownames(countsHCF2)]

labelsCells <- clustColLouv
labelsCells <- sort(labelsCells)
cellsSel <- names(labelsCells)

heatCounts <- countsHCF2[genesSelectedFirst,cellsSel]
heatCountsVal <- heatCounts[apply(heatCounts,1,function(x){sum(x)>0}),]

heatCountsStd <- t(apply(heatCountsVal, 1, function(x) x/max(x)))  # standarise

#cols <- terrain.colors(n=2)
#clust.genesFirst0 <- clust.genesFirst[match(row.names(osc), names(clust.genesFirst))]
#cols <- cols[clust.genesFirst0]
g <- anSeq::getGeneName(rownames(heatCountsStd))[[2]]
c <- labelsCells

pdf(paste0(wd,"PhD_BPS32/release6/plots/blood_lineage/Embryo10Xv6_Dongetal_HC_markerGenespaper.pdf"),
    width=9,height=9)

heat <- gplots::heatmap.2(as.matrix(heatCountsVal), trace="none", 
                          col=palette, 
                          Colv = FALSE, Rowv = F, 
                          ColSideColors = c, 
                          #RowSideColors = cols, 
                          dendrogram = "none", density.info = 'none', labRow = g,key=FALSE,labCol="")
dev.off()



###------------------------------------------------------
## Annotate clusters based on their criteria
###------------------------------------------------------

metaHC$ann <- as.character(clustColLouv[rownames(metaHC)])
s1 <- gsub(".*_E","E",rownames(metaHC))
s2 <- gsub("_emb.*","",s1)
metaHC$stage <- s2


labelsHC <- c(
  "black"="Definitive erythroid cells",
  "blue"="Primitive erythrocytes",
  "brown"="EMPs1",
  "green"="EMPs2",
  "pink"="Primitive macrophages",
  "red"="Myeloid_progenitors2",
  "turquoise"="Myeloid progenitors",
  "yellow"="EMPs3"
  
)
metaHC$colClustLouv <- metaHC$ann
for (i in 1:length(metaHC$ann)){
  metaHC[i,"ann"]<-labelsHC[metaHC[i,"ann"]]
}


heat <- gplots::heatmap.2(as.matrix(heatCountsVal), trace="none", 
                          col=palette, 
                          Colv = FALSE, Rowv = F, 
                          ColSideColors = DongHCcolor[labelsHC[labelsCells]], 
                          #RowSideColors = cols, 
                          dendrogram = "none", density.info = 'none', labRow = g,key=FALSE,labCol="")


colorCol <-DongHCcolor2[labelsHC[labelsCells]]
names(colorCol) <- names(labelsCells)

colorCol <- sort(colorCol)


orderCol<- c("purple","black","blue","orange","brown")
colorCol2 <- c()
for (i in orderCol){
  colorCol2 <-c(colorCol2,colorCol[colorCol==i])
}

pdf(paste0(wd,"PhD_BPS32/release6/plots/blood_lineage/Embryo10Xv6_Dongetal_HC_markerGenespaper.pdf"),
    width=5,height=6)
heat <- gplots::heatmap.2(as.matrix(heatCounts[,names(colorCol2)]), 
                          trace="none", 
                          col=palette, 
                          Colv = FALSE, Rowv = F, 
                          ColSideColors = as.character(colorCol2), 
                          #RowSideColors = cols, 
                          dendrogram = "none", density.info = 'none', labRow = g,key=FALSE,labCol="")
dev.off()






###------------------------------------------------------
## Map 10X myeloid cells to their HC landscape
###------------------------------------------------------


commonGenes <- HVG[HVG%in%rownames(countsMy)]
#A. Log norm
correl <- cor(countsMy[commonGenes,],
              countsHCF2[commonGenes,],method="spearman")
#correl<-t(correl)
#dist_correl <- sqrt(0.5*((1-correl)))

TissueLabel <- metaHC$ann
names(TissueLabel)<-rownames(metaHC)
#names(TissueLabel) <- cellsDong
k=5
dist_correl <- sqrt(0.5*((1-correl)))


numbersCellTypes <- matrix(0L,nrow=nrow(dist_correl),ncol=10)
colnames(numbersCellTypes) <- c(as.character(unique(metaHC$ann)),
                             "significant","assignation")
numbersCellTypes[,"significant"] <- as.logical(numbersCellTypes[,"significant"])
for (i in 1:nrow(dist_correl)){
  x <- dist_correl[i,]
  y <- x
  names(y) <- colnames(dist_correl)
  z <- sort(y,decreasing=F)[1:k]
  
  c <- as.character(TissueLabel[names(z)])
  numbersCellTypes[i,names(table(c))] <- table(c)
  if (length(table(c))>1){
    numbersCellTypes[i,"significant"] <- sort(table(c),decreasing=TRUE)[1] > (sort(table(c),decreasing=TRUE)[2])
    
  } else { numbersCellTypes[i,"significant"] <- TRUE }
  
  numbersCellTypes[i,"assignation"] <- names(sort(table(c),decreasing=TRUE))[1]
  
}

rownames(numbersCellTypes) <- rownames(correl)

table(numbersCellTypes[,"significantFreq"])
numbersCellTypes<-as.data.frame(numbersCellTypes)

table(numbersCellTypes$assignation)


numbersCellTypes1 <- numbersCellTypes[numbersCellTypes$significant==T,]


tab1<-t(table(metaSub[rownames(numbersCellTypes1),'subclust3'],
              as.character(numbersCellTypes1$assignation)))[,"My"]
#Myeloid progenitors  Myeloid_progenitors2 Primitive macrophages 
#30                    21                     1
tab3<- c((30+21),1)/(30+21+1)


pdf(paste0(wd,"PhD_BPS32/release6/plots/blood_lineage/Embryo10Xv6_Dongetal_barplot_assignation_onlyMy.pdf"),
    width=3,height=5)
barplot(tab3,
        col=DongHCcolor2[names(tab1)[c(1,3)]],ylim=c(0,1))

dev.off()
