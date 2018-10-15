#############################################
#title: "10X embryo data - Tal1-/- chimeras"
#author: "Blanca Pijuan-Sala"
#date: "25 September 2018"
#############################################

wd <- "/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/My_papers_submitted/20180601_PijuanSalaEtAl_Gastrulation10X/02.Review_01/SCRIPTS/"
library(bglab)
library(anSeq)#You can find the anSeq package here: https://github.com/BPijuanSala/anSeq



############################################
# FUNCTIONS
############################################


smart_decomp = function(counts, stage, sample, loess.span = 0.05, assay = "logcounts"){
  require(scran)
  
  method = switch( (length(unique(sample)) > 1) + 1,
                   "base",
                   switch((length(unique(stage)) > 1) + 1,
                          "sample",
                          "stage_sample"))
  
  #fit.base = trendVar(sce, use.spikes=FALSE, loess.args = list(span = loess.span), 
  #                    assay.type = assay)
  fit.base <- scran::trendVar(as.matrix(counts), loess.args = list(span = loess.span))
  decomp.all = decomposeVar(as.matrix(counts), fit.base)
  
  if(method %in% c("sample", "stage_sample")){
    
    fit.sample <- scran::trendVar(as.matrix(counts), loess.args = list(span = loess.span),
                                  design = model.matrix(~ factor(sample)))
    
    decomp.sample = decomposeVar(as.matrix(counts), fit.sample)
  }
  
  if(method == "stage_sample"){
    fit.stage <- scran::trendVar(as.matrix(counts), loess.args = list(span = loess.span),
                                 design = model.matrix(~ factor(stage)))
    
    
    decomp.stage = decomposeVar(as.matrix(counts), fit.stage)
  }
  
  
  if(method == "base"){
    decomp.new = decomp.all
  } else if(method == "sample"){
    decomp.new = decomp.sample
  } else {
    decomp.new = data.frame(mean = decomp.all$mean, 
                            total = decomp.all$total + decomp.sample$total - decomp.stage$total,
                            row.names = rownames(decomp.all),
                            stringsAsFactors = FALSE)
    
    #repeat the trend fitting like scran does
    #i.e. using log variance, so the loess can never fit a variance < 0
    fitvar = log(decomp.new$total[decomp.new$total>0])
    fitmean = decomp.new$mean[decomp.new$total>0]
    new_fit = loess(fitvar ~ fitmean, span = 0.05)
    decomp.new$tech = exp(predict(new_fit, data.frame(fitmean = decomp.new$mean)))
    decomp.new$tech[is.na(decomp.new$tech)] = 0
    
    #test via scran
    decomp.new$bio = decomp.new$total - decomp.new$tech
    mm = model.matrix(~stage + sample) #used both of these covariates
    decomp.new$p.value = testVar(decomp.new$total, decomp.new$tech, 
                                 df = nrow(mm) - ncol(mm), test = "chisq")
    decomp.new$FDR = p.adjust(decomp.new$p.value, method = "fdr")
  }
  
  if(method %in% c("sample", "stage_sample")){
    #add the variance ratios to data frame
    #that is, the ratio of the variance of the mean between samples
    #to the mean of the variance within samples
    samplesDiscard <- names(table(sample))[table(sample)==1]
    samplesToTake <- setdiff(names(table(sample)),samplesDiscard)
    #countsTest <- counts[,sample != samplesDiscard]
    #cellsForRatio <- tablerownames(metaEpi)[metaEpi$sample!=9]
    
    var_within = sapply(unique(as.character(samplesToTake)), 
                        function(x){
                          rowVars(as.matrix(counts[,sample == x]))
                          
                        })
    
    means = sapply(unique(samplesToTake), function(x){
      Matrix::rowMeans(as.matrix(counts[,sample == x]))
    })
    mean_var_within = rowMeans(var_within)
    var_mean = rowVars(means)
    decomp.new$ratio = rowVars(means)/c(rowMeans(var_within)+0.001)
  }
  
  
  return(decomp.new)
}


choose_hvgs = function(decomp, fdr.pval = 0.01, nhvg = NULL, min.mean = 0.5e-3, 
                       max.mean = NULL, use.inflection = TRUE, return.plots=FALSE){
  require(biomaRt)
  require(ggplot2)
  require(gam)
  
  #retain for plot
  decomp.plot = as.data.frame(decomp)
  
  #remove unexpressed/lowly expressed
  decomp = decomp[decomp$mean > min.mean,]
  
  ratio_present = "ratio" %in% names(decomp)
  #override max.mean if inflection specified
  if(use.inflection & ratio_present){
    keep = !is.na(decomp$ratio)
    inflection_df = data.frame(gene = rownames(decomp)[keep], mean = decomp$mean[keep], ratio = decomp$ratio[keep])
    model_df = data.frame(logratio = log10(inflection_df$ratio),
                          logmean = log10(inflection_df$mean))
    #Fit a GAM model, see result in plots$gam_fit
    model = gam(data = model_df, formula = logratio ~ logmean + s(logmean, df = 10))
    #predict an evenly spread number of points across the data
    points = seq(from = min(log10(inflection_df$mean)),
                 to = max(log10(inflection_df$mean)),
                 length.out = 1e4)
    vals = predict(model, data.frame(logmean = points))
    #take the second derivative i.e. how is the gradient changing?
    twodiff = diff(diff(vals, 1), 1)
    names(twodiff) = points[-c(1:2)]
    #find the peak i.e. the steepest changing part on the peak at the end
    peak = which.max(twodiff)
    #low identifies when we start the ascent up to this peak
    #didn't use third derivative as it's unstable 
    # low = max(which(twodiff[seq_len(peak)]< 0)) #derivative hits 0
    low = max(which(twodiff[seq_len(peak)] < (max(twodiff) * 0.1) )) #derivative hits 10%
    #redefine max mean. [low] will be more conservative than [peak].
    max.mean = 10^(as.numeric(names(twodiff)[low]))
  } else if (use.inflection & !ratio_present){
    max.mean = as.numeric(quantile(decomp.plot$mean, 1-0.073)) #0.0932 from HVG script, mean value between stage hvg calls
  }
  
  #remove highly expressed genes
  if(!is.null(max.mean)){
    decomp = decomp[decomp$mean < max.mean,]
  }
  
  
  
  sex_genes <- as.vector(read.table(paste0(wd,"DATA/mmusculus_Y_genes.txt"),header=FALSE)[,1])
  
  decomp = decomp[!rownames(decomp) %in% sex_genes,]
  
  #re-correct FDR (after thresholding)
  decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
  
  #select HVGs
  if(all(is.null(fdr.pval), is.null(nhvg)) | (!is.null(fdr.pval) & !is.null(nhvg))){
    stop("Please use exactly one of fdr.pval or nhvg")
  }
  if(!is.null(fdr.pval)){
    decomp$select = decomp$FDR <= fdr.pval
  } else {
    decomp$select = FALSE
    decomp$select[order(decomp$p.value, -decomp$bio, decreasing = FALSE)[seq_len(nhvg)]] = TRUE
  }
  
  hvgs = rownames(decomp)[decomp$select]
  
  if(return.plots){
    decomp.plot$select = rownames(decomp.plot) %in% hvgs
    nhvg = sum(decomp.plot$select)
    
    chosen_genes = ggplot(decomp.plot, aes(x = mean, y = total)) +
      geom_point(aes(col = select), size= 0.8) +
      geom_line(aes(y = tech), col = "cornflowerblue") +
      scale_x_log10(breaks = c(1e-3, 1e-1, 1e1), labels = c("0.001", "0.1", "10")) + 
      scale_y_log10(breaks = c(1e-4, 1e-2, 1e0), labels = c("0.0001", "0.01", "1")) +
      theme_bw() +
      scale_color_manual(values = c("TRUE" = "coral", "FALSE" = "grey30")) +
      theme(legend.position = "none") +
      labs(x = "log10-count mean", y = "log10-count variance") +
      ggtitle(paste0("nHVG = ", nhvg))
    
    plots = list(chosen_genes = chosen_genes)
    
    if(use.inflection & ratio_present){
      
      gam_fit = ggplot(model_df, aes(x = logmean, y = logratio)) +
        geom_point(colour = "grey80") +
        geom_line(data = data.frame(X = points, Y = vals), mapping = aes(x = X, y = Y), col = "red") +
        theme_bw() +
        labs(x = "log10 mean(log-counts)", y = "log10 ratio") +
        geom_vline(mapping = aes(xintercept = log10(max.mean)), col = "coral", lty = "twodash")
      
      twodiff = ggplot(data.frame(x = points[-c(1:2)], y = twodiff), mapping = aes(x =x, y=y)) +
        geom_line() +
        labs(x = "log10 mean(log-counts)", y= "Second derivative") +
        theme_bw() +
        geom_vline(mapping = aes(xintercept = log10(max.mean)), col = "coral", lty = "twodash")
      
      plots = list(chosen_genes = chosen_genes, gam_fit = gam_fit, twodiff = twodiff)
    }
    return(list(hvgs = hvgs, plots = plots))
  }
  
  return(hvgs)
}



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


####################
# Palettes
##################

heatmapRedYelBlue <- c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090",
                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")

heatmapOrangeGrad <- rev(c("#fef0d9","#fdd49e","#fdbb84","#fc8d59","#e34a33",
                           "#b30000"))


heatmapColsBW <- rev(c("#ffffff",
                       "#f0f0f0",
                       "#d9d9d9",
                       "#bdbdbd",
                       "#969696",
                       "#737373",
                       "#525252",
                       "#252525"))

palette <- colorRampPalette(rev(heatmapRedYelBlue))


bluePal <- c("#BFBFBF","#6495ED","#000000")
redPal <- c("gray","#ff6666","#e60000","#990000")
all_colours = c(
  "Allantois" = "#532C8A",#[32] "Allantois" 
  "Anterior Primitive Streak" = "#c19f70",
  "Blood progenitors 1" = "#f9decf",
  "Blood progenitors 2" = "#c9a997",
  "Cardiomyocytes" =  "#B51D8D",
  "Caudal epiblast" = "#9e6762",
  "Caudal Mesoderm" = "#3F84AA",
  "Def. endoderm" = "#F397C0",
  "Nascent mesoderm" =  "#C594BF",                    
  "Mixed mesoderm" =  "#DFCDE4",
  
  "Endothelium" =  "#eda450",                           
  "Epiblast" =  "#635547",
  "Erythroid1" =  "#C72228",                          
  "Erythroid2" =  "#EF4E22",  
  "Erythroid3" =  "#f77b59",
  "ExE ectoderm" =  "#989898",                      
  
  "ExE endoderm" = "#7F6874",                          
  "ExE mesoderm" =  "#8870ad",
  
  
  "Rostral neurectoderm" =  "#65A83E",
  "Forebrain/Midbrain/Hindbrain" = "#647a4f",
  
  
  "Gut" =  "#EF5A9D",                               
  "Haematoendothelial progenitors" =  "#FBBE92",
  "Caudal neurectoderm"= "#354E23",
  
  
  "Intermediate mesoderm" =  "#139992",               
  "Neural crest"= "#C3C388",
  
  "NMP" =  "#8EC792",                                        
  "Notochord" =  "#0F4A9C",                                 
  "Paraxial mesoderm" =  "#8DB5CE",
  "Parietal endoderm" =  "#1A1A1A",    
  "PGC" =  "#FACB12",                                     
  
  "Pharyngeal mesoderm" =  "#C9EBFB",             
  "Primitive Streak" =  "#DABE99",
  "Mesenchyme" = "#ed8f84",
  "Somitic mesoderm" =  "#005579",                                     
  "Spinal cord" =  "#CDE088",                               
  "Surface ectoderm" = "#BBDCA8",                       
  
  
  "Visceral endoderm" = "#F6BFCB"
  
  
)

all_colours_sub = c(
  "Mes1"= "#c4a6b2",#
  "Mes2"= "#ca728c",#
  
  "Cardiomyocytes" =  "#B51D8D",  
  
  "BP1" = "#6460c5",#
  "BP2" = "#96b8e4",#
  "Haem3"= "#02f9ff",#changed
  "BP3" = "#07499f",#changed
  "BP4" = "#036ef9",#changed
  #"BP6"  = "#03bbf9",
  
  "Haem1"= "#bb22a7",
  "Haem2" = "#f695e9",
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

stagePal <- c(
  "E6.5" ="#E3FCFA",
  "E6.75" = "#C1ECEF",
  "E7.0" ="#A3D4E3",
  "E7.25" ="#86B8D6",
  "E7.5"="#6C98CA",
  "E7.75"="#5476BE",
  "E8.0"="#3E52B1",
  "E8.25"="#2B2DA5",
  "E8.5"="#2B1999",
  "mixed_gastrulation"="gray"
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


#=====================================
# Setting up data
#=====================================

meta <- read.table(paste0(wd,"DATA/metadata_mergedClustering_Celltypes_20180910.tab"), header = T, stringsAsFactors = F, sep = "\t")    
rownames(meta) <- meta$index
meta$cell <- meta$index
metaSub <- read.table(file=paste0(wd,"DATA/metadata_subclusters_bloodLineage_v3_20180910.tab"),
                      sep="\t",header=T)
rownames(metaSub) <- metaSub$cell
metaSub$celltype_new <- as.character(metaSub$celltype_new)



metaTal1 <- read.table(paste0(wd,"DATA/chimera-tal1/meta.tab"), header = T, stringsAsFactors = F, sep = "\t")    
rownames(metaTal1) <- metaTal1$cell

metaTal1Sub <- metaTal1[as.character(metaTal1$celltype.mapped)%in%as.character(unique(metaSub$celltype_new)),]


table(metaTal1Sub$tomato)

closestCells <- read.table(paste0(wd,"DATA/chimera-tal1/blood_lineage_bloodlineage_nearestneighbours.tab"), header = F, stringsAsFactors = F, sep = "\t")    
colnames(closestCells)<- c("chimera",paste0("atlas",seq(1:10)))
rownames(closestCells)<- closestCells$chimera


metaTal1Sub <- cbind(metaTal1Sub,closestCells[rownames(metaTal1Sub),])
metaTal1Sub$subclust3 <- metaSub[as.character(metaTal1Sub$atlas1),'subclust3']



countsTal1 <- readRDS(paste0(wd,"DATA/chimera-tal1/raw_counts.rds"))

sceTal1 <- SingleCellExperiment::SingleCellExperiment(assays=c("counts"=countsTal1))
SingleCellExperiment::sizeFactors(sceTal1)<-as.numeric(as.character(read.table(paste0(wd,"DATA/chimera-tal1/sizefactors.tab"))[,1])) 

sceTal1 <- scater::normalize(sceTal1)
rownames(sceTal1)<- (as.character(read.table(paste0(wd,"DATA/chimera-tal1/genes.tsv"))[,1])) 
colnames(sceTal1) <- rownames(metaTal1)


counts <- as.matrix(Matrix::readMM(file=paste0(wd,"PhD_BPS32/release6/data/20180908_countsBlood.mtx")))

colnames(counts) <- as.character(read.table(paste0(wd,"DATA/20180908_cells_countsBlood.tab"),header=T)[,1])
rownames(counts) <-  read.table(paste0(wd,"DATA/20180908_genes_countsBlood.tab"),header=TRUE)[,1]


metaSub$celltype_new <- as.character(metaSub$celltype_new)


countsCM <- as.matrix(Matrix::readMM(file=paste0(wd,"DATA/CM_10X_sample/countsData_norm_CM.mtx")))
colnames(countsCM) <- as.character(read.table(paste0(wd,"DATA/CM_10X_sample/cells_CM.txt"),header=F)[,1])
rownames(countsCM) <-  read.table(paste0(wd,"DATA/CM_10X_sample//genes_CM.txt"),header=F)[,1]

#=====================================
# Plot landscape
#=====================================


plot(metaSub$gephiX,metaSub$gephiY,
     col="gray",
     pch=20, cex=0.2,xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
points(metaSub[unique(as.character(closestCells$atlas1)),'gephiX'],
       metaSub[unique(as.character(closestCells$atlas1)),'gephiY'],
       col=all_colours_sub[as.character(metaSub[unique(as.character(closestCells$atlas1)),'subclust3'])],
       pch=20)
#box(bty="l")
#rect(xleft=-4000, ybottom=-2500, xright= 1500, ytop=1500,lty=4,lwd=2)




cellsChim<-rownames(metaTal1Sub)[metaTal1Sub$tomato]


plot(metaSub$gephiX,metaSub$gephiY,
     col="gray",
     pch=20, cex=0.2,xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
points(metaSub[as.character(closestCells[cellsChim,"atlas1"]),'gephiX'],
       metaSub[as.character(closestCells[cellsChim,"atlas1"]),'gephiY'],
       col="red",pch=20,cex=0.4)
#box(bty="l")



cellsWT<-rownames(metaTal1Sub)[!metaTal1Sub$tomato]

plot(metaSub$gephiX,metaSub$gephiY,
     col="gray",
     pch=20, cex=0.2,xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
points(metaSub[as.character(closestCells[cellsWT,"atlas1"]),'gephiX'],
       metaSub[as.character(closestCells[cellsWT,"atlas1"]),'gephiY'],
       col="black",pch=20,cex=0.4)


metaZoomed0 <- metaSub[(metaSub$gephiX> -4000)&(metaSub$gephiY> -2500)&(metaSub$gephiX < 1500)&(metaSub$gephiY < 1500),]





plot(metaZoomed0$gephiX,metaZoomed0$gephiY,
     col="gray",
     pch=20, cex=1,xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
points(metaZoomed0[as.character(closestCells[cellsWT,"atlas1"]),'gephiX'],
       metaZoomed0[as.character(closestCells[cellsWT,"atlas1"]),'gephiY'],
       col="black",pch=20,cex=1)
#box()


plot(metaZoomed0$gephiX,metaZoomed0$gephiY,
     col="gray",
     pch=20, cex=1,xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
points(metaZoomed0[as.character(closestCells[cellsChim,"atlas1"]),'gephiX'],
       metaZoomed0[as.character(closestCells[cellsChim,"atlas1"]),'gephiY'],
       col="red",pch=20,cex=1)
#box()




#=====================================
# Plot heatmap Itga2b and Cbfa2t3 in EC6 and 7
#=====================================


cells2nd <- rownames(metaTal1Sub)[metaTal1Sub$atlas1 %in% rownames(metaSub)[metaSub$subclust3%in%c("EC6","EC7")]]


cols<-c("Cdh5","Pecam1","Itga2b","Cbfa2t3")


heatCounts <- SingleCellExperiment::logcounts(sceTal1)[anSeq::getGeneID(cols)[[2]],cells2nd]

heatCountsVal <- heatCounts[apply(heatCounts,1,function(x){sum(x)>0}),]

heatCountsStd <- t(apply(heatCounts, 1, function(x) x/max(x)))  # standarise


g <- anSeq::getGeneName(rownames(heatCounts))[[2]]

colChim <- c(
  "TRUE"="red",
  "FALSE"="black"
)
c <- colChim[as.character(metaTal1[cells2nd,"tomato"])]
names(c)<- cells2nd
d <- order(metaTal1Sub[cells2nd,"subclust3"],metaTal1Sub[cells2nd,"tomato"])
ord <- c[d]


bglab::heatmap.minus(as.matrix(heatCountsStd[,d]), trace="none", 
                     col=palette, 
                     Colv = F, Rowv = F, 
                     ColSideColors = cbind(all_colours_sub[as.character(metaTal1Sub[names(ord),"subclust3"])],ord), 
                     #RowSideColors = cols, 
                     dendrogram = "none", density.info = 'none',axis.cex=1,
                     labRow = g,key=FALSE,labCol="",margins=c(5,15)
                     
)






###==================================
# Contribution
###===================================

require(edgeR)
require(ggplot2)

meta_tal1 = metaTal1

#Get counts of the Tal1 celltypes by sample
tab_tal1 = table(meta_tal1$celltype.mapped, meta_tal1$tomato)

#exclude cells to which the Tal1 KO cannot contribute
tab_tal1 = tab_tal1[!rownames(tab_tal1) %in% c("ExE endoderm", "Visceral endoderm", "ExE ectoderm", "Parietal endoderm", "Doublet", "Stripped"),]
tab_tal1 = tab_tal1[!grepl("Eryth", rownames(tab_tal1)),]
tab_tal1 = tab_tal1[!grepl("Blood", rownames(tab_tal1)),]
#exclude technical cells
tab_tal1 = tab_tal1[!rownames(tab_tal1) %in% c("Stripped", "Doublet"),]



tab_tal1Sub = table(metaTal1Sub$subclust3, metaTal1Sub$tomato)

#exclude cells to which the Tal1 KO cannot contribute
tab_tal1Sub = tab_tal1Sub[!rownames(tab_tal1Sub) %in% c("ExE endoderm", "Visceral endoderm", "ExE ectoderm", "Parietal endoderm", "Doublet", "Stripped"),]
tab_tal1Sub = tab_tal1Sub[!grepl("Ery", rownames(tab_tal1Sub)),]
tab_tal1Sub = tab_tal1Sub[!grepl("BP", rownames(tab_tal1Sub)),]
#exclude technical cells
tab_tal1Sub = tab_tal1Sub[!rownames(tab_tal1Sub) %in% c("Stripped", "Doublet"),]
tab_tal1Sub <- tab_tal1Sub[apply(tab_tal1Sub,1,min) > 0,]


#Now to make the plot

freq_small = sweep(tab_tal1Sub, 2, colSums(tab_tal1), "/")
freq_small = freq_small[apply(tab_tal1Sub,1,max) > 20,]

fc_tal1 = log2(freq_small[,2]) - log2(freq_small[,1])


par(mfrow=c(1,1))
barplot(fc_tal1,col=all_colours_sub[names(fc_tal1)],las=2,ylim=c(-5,5),axes=F)
axis(side = 2, at = c(-4,-2,0,2,4),las=2)
dev.off()




#=====================================
# EC3 exploration only
#=====================================


sceTal1EC3 <- sceTal1[,rownames(metaTal1Sub)[as.character(metaTal1Sub$subclust3)=="EC3" & metaTal1Sub$tomato==T]]
Tal1EC3 <- as.matrix(SingleCellExperiment::logcounts(sceTal1EC3))
colnames(Tal1EC3)<-paste0("chim_",colnames(Tal1EC3))
#rm(sceTal1)
Tal1EC3 <- Tal1EC3[rownames(Tal1EC3)!="tomato-td",]


sceTal1EC3F <- sceTal1[,rownames(metaTal1Sub)[as.character(metaTal1Sub$subclust3)=="EC3" & metaTal1Sub$tomato==F]]
Tal1EC3F <- as.matrix(SingleCellExperiment::logcounts(sceTal1EC3F))
colnames(Tal1EC3F)<-paste0("chim_",colnames(Tal1EC3F))
#rm(sceTal1)
Tal1EC3F <- Tal1EC3F[rownames(Tal1EC3F)!="tomato-td",]


countsEC3 <- counts[,unique(metaTal1Sub$atlas1[as.character(metaTal1Sub$subclust3)=="EC3" & metaTal1Sub$tomato==T])]


countsEC3all <- cbind(Tal1EC3,Tal1EC3F, countsEC3,countsCM)
metaEC3all <- data.frame(
  cell = colnames(countsEC3all),
  chimera = c(rep("chimKO",ncol(Tal1EC3)),rep("chimWT",ncol(Tal1EC3F)),rep("atlas",ncol(countsEC3)),rep("CM",ncol(countsCM))),
  stage = c(metaTal1Sub[gsub("chim_","",colnames(Tal1EC3)),"stage.mapped"],metaTal1Sub[gsub("chim_","",colnames(Tal1EC3F)),"stage.mapped"],metaSub[colnames(countsEC3),"stage"],meta[colnames(countsCM),"stage"]),
  sample = c(paste0("chim",metaTal1Sub[gsub("chim_","",colnames(Tal1EC3)),"sample"]),paste0("chim",metaTal1Sub[gsub("chim_","",colnames(Tal1EC3F)),"sample"]),metaSub[colnames(countsEC3),"sample"],meta[colnames(countsCM),"sample"])
)


colorsClust <- c(
  "atlas" = "#eda450",
  "chimKO"="red",
  "CM" = "#B51D8D",
  "chimWT"="black"
)

#HVG for 10X subset
decomp <- smart_decomp(countsEC3all,stage=metaEC3all$stage,
                       sample=metaEC3all$sample)
hvg <- choose_hvgs(decomp,return.plots=TRUE,use.inflection = FALSE,max.mean = 1)

HVGFilt <- hvg$hvgs


###########Comparison KO vs atlas endothelium

groups <- metaEC3all$chimera[metaEC3all$chimera %in% c("atlas","chimKO")]
names(groups)<-as.character(metaEC3all$cell)[metaEC3all$chimera %in% c("atlas","chimKO")]
markerschimEC <- anSeq::findMarkers(countsEC3all[,names(groups)],groups=groups,groupTarget = 'chimKO')
markerschimEC$geneNames <- anSeq::getGeneName(rownames(markerschimEC))[[2]]
#save(markerschimEC,file=paste0(wd,"PhD_BPS46/release2/rda/markerschimEC_2.rda"))
#load(paste0(wd,"PhD_BPS46/release2/rda/markerschimEC_2.rda"))

###########Comparison atlas CM vs atlas endothelium

groups <- metaEC3all$chimera[metaEC3all$chimera %in% c("CM","atlas")]
names(groups)<-as.character(metaEC3all$cell)[metaEC3all$chimera %in% c("CM","atlas")]
markerschimCMat<- anSeq::findMarkers(countsEC3all[,names(groups)],groups=groups,groupTarget = 'atlas')
markerschimCMat$geneNames <- anSeq::getGeneName(rownames(markerschimCMat))[[2]]

#save(markerschimCMat,file=paste0(wd,"PhD_BPS46/release2/rda/markerschimCMat.rda"))
#load(paste0(wd,"PhD_BPS46/release2/rda/markerschimCMat.rda"))


###########Comparison KO vs WT

groups <- metaEC3all$chimera[metaEC3all$chimera %in% c("chimKO","chimWT")]
names(groups)<-as.character(metaEC3all$cell)[metaEC3all$chimera %in% c("chimKO","chimWT")]
markerschimOnly<- anSeq::findMarkers(countsEC3all[,names(groups)],groups=groups,groupTarget = 'chimKO')
markerschimOnly$geneNames <- anSeq::getGeneName(rownames(markerschimOnly))[[2]]

#save(markerschimOnly,file=paste0(wd,"PhD_BPS46/release2/rda/markerschimOnly.rda"))
#load(paste0(wd,"PhD_BPS46/release2/rda/markerschimOnly.rda"))


genesSelEC <- rownames(markerschimCMat)[markerschimCMat$logFC>2.9]
genesSelCM <- rownames(markerschimCMat)[markerschimCMat$logFC<(-2.9)]

genesSelECchim <- rownames(markerschimEC)[markerschimEC$logFC>=2]
genesSelCHIMchim <- rownames(markerschimOnly)[markerschimOnly$logFC>(2)]




genesSel2 <- unique(c(genesSelECchim,genesSelEC,genesSelCM,genesSelCHIMchim))
genesSel2 <- genesSel2[genesSel2!=anSeq::getGeneID("Tal1")[[2]]]


heatCounts <- countsEC3all[unique(genesSel2),]
cols<-rep("gray",length(unique(genesSel2)))
cols[unique(genesSel2)%in%genesSelCM]<-"#B51D8D"
cols[unique(genesSel2)%in%genesSelEC]<-"#eda450"
names(cols)<-genesSel2
cols <- sort(cols)

g <- anSeq::getGeneName(names(cols))[[2]]
c <- colorsClust[as.character(metaEC3all$chimera)]
names(c)<-metaEC3all$cell

den <- anSeq::leaf_reordering(heatCounts)

heatCountsVal <- heatCounts[apply(heatCounts,1,function(x){sum(x)>0}),]

heatCountsStd <- t(apply(heatCounts, 1, function(x) x/max(x)))  # standarise


c2 <- c[den$labels[den$order]]
c3 <- c()

colorsClust <- colorsClust[c("CM","chimKO","chimWT","atlas")]
for (i in colorsClust){
  c3 <- c(c3,c2[c2 == i])
}


gplots::heatmap.2(heatCountsStd[names(cols),names(c3)], trace="none", 
                  col=palette, 
                  Colv = F, Rowv = T, 
                  ColSideColors = c3, 
                  RowSideColors = cols, 
                  #dendrogram = "row", 
                  density.info = 'none',
                  labRow = g,key=FALSE,labCol="",margins = c(5,8))




#Selection of representative genes


genesSel3 <- c("Gata5","Myl4","Tnnc1","Nexn","Nkx2-5","Myl7","Tnnt2","Tnni","Acta2","Mef2c","Prrx2",
               "Sox17","Ctla2a","Flt1","Pecam1","Cd34","Cdh5","Fli1","Etv2","Lmo2","Esam","Kdr",
               "Slc24a5","Rpgrip1","Dpysl2","Impac1","Pcolce","Rab3il1","Tdo2","Plagl1")



heatCounts <- countsEC3all[anSeq::getGeneID(genesSel3)[[2]],]
cols<-rep("gray",nrow(heatCounts))
cols[(rownames(heatCounts))%in%c(genesSelCM,anSeq::getGeneID("Myl7")[[2]])]<-"#B51D8D"
cols[(rownames(heatCounts))%in%genesSelEC]<-"#eda450"
names(cols)<-rownames(heatCounts)
cols <- sort(cols)

g <- anSeq::getGeneName(names(cols))[[2]]
c <- colorsClust[as.character(metaEC3all$chimera)]
names(c)<-metaEC3all$cell

den <- anSeq::leaf_reordering(heatCounts)

heatCountsVal <- heatCounts[apply(heatCounts,1,function(x){sum(x)>0}),]

heatCountsStd <- t(apply(heatCounts, 1, function(x) x/max(x)))  # standarise



gplots::heatmap.2(heatCountsStd[names(cols),names(c3)], trace="none", 
                  col=palette, 
                  Colv = F, Rowv = F, 
                  ColSideColors = c3, 
                  RowSideColors = cols, 
                  #dendrogram = "row", 
                  density.info = 'none',
                  labRow = g,key=FALSE,labCol="",margins = c(5,8))







