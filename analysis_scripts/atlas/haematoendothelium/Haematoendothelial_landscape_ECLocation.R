#############################################
#title: "10X embryo data - Mapping EC location"
#author: "Blanca Pijuan-Sala"
#date: "09 September 2018"
#############################################
wd <- "/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/My_papers_submitted/20180601_PijuanSalaEtAl_Gastrulation10X/02.Review_01/SCRIPTS/"

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


############################################
# palettes
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




###------------------------------------------------------
##Load data 10X data for haemato-endothelial lineage
###------------------------------------------------------


metaSub <- read.table(file=paste0(wd,"DATA/metadata_subclusters_bloodLineage_v3_20180910.tab"),
                      sep="\t",header=T)
rownames(metaSub)<-metaSub$cell
counts <- as.matrix(Matrix::readMM(file=paste0(wd,"DATA/20180908_countsBlood.mtx")))

colnames(counts) <- as.character(read.table(paste0(wd,"DATA/20180908_cells_countsBlood.tab"),header=T)[,1])
rownames(counts) <-  read.table(paste0(wd,"DATA/20180908_genes_countsBlood.tab"),header=TRUE)[,1]

#######================
# CLASSIFY ENDOTHELIUM CLUSTER WITH NEAREST NEIGHBOURS INTO YS, EP, AL
#######================
cellsSel <- rownames(metaSub)[metaSub$subclust3%in%paste0("EC",1:8)]

countsendoClass <- counts[,cellsSel]
metaendoClass  <- metaSub[cellsSel,]
plot(metaSub$gephiX,metaSub$gephiY,
     col="gray",
     pch=20, cex=0.35,xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
points(metaendoClass$gephiX,metaendoClass$gephiY,pch=20, cex=0.35,col="red")
box(bty="l")


metaSmart <- xlsx::read.xlsx(paste0(wd,"DATA/Smartseq2_YSEPAL/PhD_BPS13_run2_metadata_processed.xls"),
                             sheetIndex = 1,header=TRUE)
rownames(metaSmart)=metaSmart[,1]

cellsWT <- rownames(metaSmart)[as.logical(metaSmart$passCellQC)]
metaSmartWT <- metaSmart[cellsWT,]


load(file=paste0(wd,"DATA/Smartseq2_YSEPAL/sce_countsYSEPAL.rda"))
sce_loc <- scater::normalize(sce_loc)
countsSmartWT <- as.matrix(SingleCellExperiment::logcounts(sce_loc))
countsSmartWTfilt <- countsSmartWT[,apply(countsSmartWT,2, function(x) sum(x)>0)]


countsSmartWTNormLog <- countsSmartWT



#HVG for 10X subset
decomp <- smart_decomp(countsendoClass,stage=metaendoClass$stage,
                       sample=metaendoClass$sample)
hvg <- choose_hvgs(decomp,return.plots=TRUE,use.inflection = FALSE,max.mean = 1)

HVGFilt <- hvg$hvgs

HVGinSmart <-HVGFilt[HVGFilt%in%rownames(countsSmartWT)]




cellsSelTrain <- rownames(metaSmartWT)


correl <- cor(countsendoClass[HVGinSmart,],
              countsSmartWTNormLog[HVGinSmart,cellsSelTrain],method="spearman")

dist_correl <- sqrt(0.5*((1-correl)))

TissueLabel <- metaSmartWT$Celltypegeneral
names(TissueLabel) <- rownames(metaSmartWT)
k=5
numbersEPYSAL <- matrix(0L,nrow=nrow(dist_correl),ncol=5)
colnames(numbersEPYSAL) <- c(as.character(unique(metaSmartWT$Celltypegeneral)),
                             "significant","assignation")
numbersEPYSAL[,"significant"] <- as.logical(numbersEPYSAL[,"significant"])
for (i in 1:nrow(dist_correl)){
  x <- dist_correl[i,]
  y <- x
  names(y) <- colnames(dist_correl)
  z <- sort(y,decreasing=F)[1:k]
  
  c <- as.character(TissueLabel[names(z)])
  numbersEPYSAL[i,names(table(c))] <- table(c)
  if (length(table(c))>1){
    numbersEPYSAL[i,"significant"] <- sort(table(c),decreasing=TRUE)[1] > (sort(table(c),decreasing=TRUE)[2])
    
  } else { numbersEPYSAL[i,"significant"] <- TRUE }
  
  numbersEPYSAL[i,"assignation"] <- names(sort(table(c),decreasing=TRUE))[1]
  
}

rownames(numbersEPYSAL) <- rownames(dist_correl)
numbersEPYSAL[1,"significant"]<- "TRUE"
numbersEPYSAL[,"significant"] <- as.logical(numbersEPYSAL[,"significant"])

table(numbersEPYSAL[,"significant"])

numbersEPYSAL[numbersEPYSAL[,"significant"]==FALSE,]

colorCell <- c("#f58231","#0082c8","#ffe119")
names(colorCell) <- as.character(unique(TissueLabel))

#tiff(paste0(wd,"PhD_BPS32/release4/bloodJourneys/plots/Embryo10Xv4_Preselection7_bloodJourneys_endotheliumAssignation_all.tiff"),
#     width=7.5,height=8,units = "in",res=300)
plot(metaSub$gephiX,metaSub$gephiY,pch=20,
     col="gray",main="",xlab="",ylab="",cex=0.5,cex.main=2,xaxt='n',yaxt='n',axes=F)
box(bty="l")

for (i in 1:length(unique(TissueLabel))){
  Code <- as.character(unique(TissueLabel)[i])
  cellsCode0 <- rownames(numbersEPYSAL)[numbersEPYSAL[,"assignation"]==Code]
  cellsCode <- cellsCode0[as.logical(numbersEPYSAL[cellsCode0,"significant"])]
  plot(metaSub$gephiX,metaSub$gephiY,pch=20,
       col="gray",main="",xlab="",ylab="",cex=0.5,cex.main=2,xaxt='n',yaxt='n',axes=F)
  box(bty="l")
  points(metaSub[cellsCode,"gephiX"],metaSub[cellsCode,"gephiY"],pch=20,
         col=colorCell[as.character(Code)],xlab="GEPHI1",ylab="GEPHI2",cex=0.5)
}
#dev.off()

numbersEPYSALsign <- numbersEPYSAL[as.logical(numbersEPYSAL[,'significant']),]
metaSub$Location <- as.character(rep("None",nrow(metaSub)))
metaSub[rownames(numbersEPYSALsign),"Location"] <- unname(numbersEPYSALsign[,"assignation"])




tab <- table(metaSub$subclust3,metaSub$Location)
tab1 <- tab[,c("Allantois","EmbryoProper","YolkSac")]
tab2 <- tab1[apply(tab1,1,function(x){sum(x)>0}),]
tab3<- t(apply(tab2,1,function(x){x/sum(x)}))



colorLocation <- c("#f58231","#0082c8","#ffe119")
names(colorLocation) <- as.character(unique(TissueLabel))

tab3 <- tab3[c(paste0("EC",1:8)),]




par(mfrow=c(1,1))
barplot(t(tab3),col=colorLocation[colnames(tab3)],ylim=c(0,1),axes = FALSE,las=2)
axis(side = 2, at = c(0,0.25,0.5,0.75,1),las=2,labels=c(0,NA,0.5,NA,1))





