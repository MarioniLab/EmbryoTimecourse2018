#############################################
#title: E8.25 YS, EP, AL with Smart-seq2 
#author: "Blanca Pijuan-Sala"
#date: "28 September 2017"
#############################################


wd <- "/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/My_papers_submitted/20180601_PijuanSalaEtAl_Gastrulation10X/02.Review_01/SCRIPTS/"

library(Biobase)
library(ggplot2)
library(Matrix)
library(knitr)
library(scran)
library(igraph)
library(Rtsne)
library(proSeq)#You can find the proSeq package here: https://github.com/BPijuanSala/proSeq
###------------------------------------------------------
##Load data Smartseq2 data
###------------------------------------------------------

counts = read.table(paste0(wd, "DATA/Smartseq2_YSEPAL/PhD_BPS13_run2_countsRaw_ensembl92.txt"),sep="\t",header=TRUE)

colnames(counts)<-sub("SLX-","SLX.",gsub(".H.*","",colnames(counts)))
colnames(counts)<-sub("_",".",colnames(counts))

colnames(counts)<-sub(".i5","_i5",colnames(counts))

load(paste0(wd,"DATA/Smartseq2_YSEPAL/PhD_BPS13_run2_metadata.rda"))
colnames(metadata)<-gsub('\\.','',colnames(metadata))

rownames(metadata) <- as.character(paste0(metadata$CRIidentifier,".",metadata$CIannotationofindex))
rownames(metadata)<-sub("-",".",rownames(metadata))
rownames(metadata)<-sub(".i5","_i5",rownames(metadata))


###################


countsfeat = counts[grep("^__",rownames(counts)),]
rownames(countsfeat) = c("no.feature","ambiguous","lowQuality","unaligned","align.not.unique")
countsRaw = counts[grep("^__",rownames(counts),invert = TRUE),]



spikes = substr(rownames(countsRaw),1,4)=="ERCC"
names(spikes)=rownames(countsRaw)




metadata = metadata[rownames(metadata)%in%colnames(countsRaw),]
countsRaw <- countsRaw[,colnames(countsRaw)%in%rownames(metadata)]
metadata <- metadata[colnames(countsRaw),]
countsfeat1 <- countsfeat[,colnames(countsRaw)]

data <- new("RNAseq", countsRaw=as.matrix(countsRaw), metadata=as.matrix(metadata),SpikeIn=spikes,CountsFeat=as.matrix(countsfeat1))


###------------------------------------------------------
## Do CellQC
###------------------------------------------------------

annTable <- read.table(paste0(wd,"DATA/geneTable_mart_export_mm10_ensembl92.txt"),header=T,sep="\t")
colnames(annTable)<-colnames(geneTable)
cellQC_data <- runCellQC(data,checkCountsFeat = TRUE,plotting="tiff",annTable = annTable,
                         outputPlots = paste0(wd,"DATA/Smartseq2_YSEPAL/QCplots_"),maxMapmit=0.10,
                         minGenesExpr = 4000,minMapreads = 50000,maxMapSpike=0.5)



cellsFailQC(data)<-c(cellQC_data)


#These are the numbers for the cells that have failed

table(metadata[names(cellQC_data[cellQC_data==TRUE]),"Celltypegeneral"])
#Allantois EmbryoProper      YolkSac 
#14            7           17 



#cat("cells passed QC")
table(metadata[names(cellQC_data[cellQC_data==FALSE]),"Celltypegeneral"])
#Allantois EmbryoProper      YolkSac 
#82           89           79 


metadata$passCellQC <- !cellQC_data


xlsx::write.xlsx(metadata,paste0(wd,"DATA/Smartseq2_YSEPAL/PhD_BPS13_run2_metadata_processed.xls"))


###------------------------------------------------------
## Normalise data
###------------------------------------------------------


NormOutput <- normalise(data)
countsNorm(data) <- NormOutput$countsNorm

CellsSizeFac(data) <- list(sizeFactorsGenes = NormOutput$sizeFactorsGenes, 
                             sizeFactorsSpikeIn = NormOutput$sizeFactorsSpikeIn)







sce_loc = SingleCellExperiment(assays = list("counts" = as.matrix(countsRaw)[,names(NormOutput$sizeFactorsGenes)]))                 

sizeFactors(sce_loc) = NormOutput$sizeFactorsGenes



#save(sce_loc,file=paste0(wd,"DATA/Smartseq2_YSEPAL/sce_countsYSEPAL.rda"))








