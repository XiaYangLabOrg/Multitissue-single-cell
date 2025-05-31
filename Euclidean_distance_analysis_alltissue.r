#Eucdistance-allcombined
library(Seurat)
library(dplyr)
library(plyr)
library(parallel)
library(doSNOW)
library(reshape)
library(ggplot2)
library(ggrepel)
rm(list = ls())

rm(list=ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/SI_combined_Set2_final.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SI.rda")

seuratObject=dropEST.combined.filtered
DefaultAssay(seuratObject) <- "RNA"
#metaData1=NULL
#metaDataSubset=NULL
cellTypeMetaData="celltype"
nPerm=5000
do.par=TRUE
num.cores= 4
comparisonMetaData="data.diseasestatus"
group1="Fruc"
group2="Control"
#distanceScorePlot=NULL
#fcPlot=NULL
#distanceSave=NULL
seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)

#Number of permutations
numPerm = nPerm
if(do.par==TRUE){
  #auto detect number of available cores
  if(num.cores == "auto"){
    cores=detectCores()
    #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
    cores = max(1,cores[1]-1)
    #setup and start parallel processing cluster
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else{
    #if a number of cores for parallel processing is specified, use this number of cores
    if(is.numeric(num.cores)){
      cores = num.cores
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
  }
}
cellTypePermutationDistances = list()
cellTypeDistances = c()
allcelltypes = levels(seuratObject@active.ident)
availcelltypes <- NULL
for(i in 1:length(allcelltypes)){
  cellType = allcelltypes[i]
  print(cellType)
  #subset to cell type of interest
  cellTypeSubset = subset(seuratObject,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
  
  
  group1Length = length(which(cellTypeSubset@active.ident==group1))
  group2length = length(which(cellTypeSubset@active.ident==group2))
  
  #Only continue if we have at least 3 cells of the cell type in each group
  if(group1Length > 9 & group2length > 9){
    availcelltypes <- c(availcelltypes, cellType)
    #Label group1 and group2 cells
    group1Cells = WhichCells(cellTypeSubset,idents = group1)
    group2Cells = WhichCells(cellTypeSubset,idents = group2)
    
    #Get the average expression of group1 and group2 cells
    cellTypeSubset <- ScaleData(cellTypeSubset, vars.to.regress = NULL)
    
    AverageExpall <- rowMeans(as.matrix(cellTypeSubset@assays$RNA@counts))
    AverageExpall <- sort(AverageExpall, decreasing = T)
    AverageExpall <- AverageExpall[AverageExpall > quantile(AverageExpall,probs = 0.9)]
    #AverageExpall <- AverageExpall[-c(1:20)]
    selectedgenes <- names(AverageExpall)
    selectedgenes <- selectedgenes[1:1000]
    #compute in non-log space
    #avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes))
    
    avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes, use.scale = T))
    avgExp = avgExp$RNA
    #top and bottom filtration
    #calculate the distance between the group1 and group2 cells
    cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
    
    cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
    
    #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
    allCells = c(group1Cells,group2Cells)
    
    #set up a progress bar
    pb <- txtProgressBar(max = numPerm, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
      set.seed(j)
      permSubset = cellTypeSubset
      sampleCells = sample(allCells,length(allCells),replace = FALSE)
      #place the cells into either group1 or group2 (same size as original group1 and group2)
      group1Synth = sampleCells[1:length(group1Cells)]
      group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
      
      permStatus = permSubset@active.ident
      permStatus[group1Synth] = group1
      permStatus[group2Synth] = group2
      permSubset$permStatus = permStatus
      permSubset = SetIdent(permSubset, value = "permStatus")
      
      #compute in non-log space
      #permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes))
      
      permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes, use.scale = T))
      permAvgExp = permAvgExp$RNA
      
      permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
      permCellTypeDistance
    }
    cellTypePermutationDistances[[cellType]] = permutationDistances
    close(pb)
  }
}

stopCluster(cl)

names(cellTypeDistances) = availcelltypes

cellTypeDistancePvals = c()
for(cellType in levels(seuratObject@active.ident)){
  numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
  cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
}
names(cellTypeDistancePvals) = availcelltypes


melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
#availcelltypes = unique(melted.cellTypePermutationDistances$L1)

melted.cellTypePvals = rep(cellTypeDistancePvals[availcelltypes],each=numPerm)
melted.cellTypeDistances = rep(cellTypeDistances[availcelltypes],each=numPerm)

melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)

#Bonferroni correction
allDistancePvals = c(cellTypeDistancePvals[availcelltypes])
# all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm

#Combine neuronal subtype and cell types to plot
colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
melted.combinedPermutationDistances = melted.cellTypePermutationDistances
melted.combinedPermutationDistances$pvals = allDistancePvals
melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))


#########
#Plot
ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
  geom_density(adjust = 1) +
  facet_wrap( ~ L1, ncol=5,scales="free") +
  xlab("Euclidian Distance") +
  ylab("Density") +
  # Add line
  geom_vline(aes(xintercept=distances),
             color="blue", linetype="dashed", size=1) +
  annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)

pdf(file = paste0("./Plots/EuclideanDistanceFructose/Fructoseplot1.pdf"),height = 16,width = 16)
print(p1)
dev.off()

#########
#get the log FC of the euclidean distance versus the median of the background distribution

medianBackground = c()
for(cellType in levels(seuratObject@active.ident)){
  medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
}
names(medianBackground) = availcelltypes

logFC = log10(cellTypeDistances[availcelltypes]) - log10(medianBackground)
pval = -log10(allDistancePvals)
toPlot = as.data.frame(cbind(logFC,pval))

p2 <- ggplot(toPlot, aes(pval, logFC)) +
  geom_point(color = 'black') +
  theme_classic(base_size = 10) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 

pdf(file = paste0("./Plots/EuclideanDistanceFructose/Fructoseplot2.pdf"))
print(p2)
dev.off()

save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Fructosedistanceresultzscore.rda"))

rm(list=ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined")
#load("SI_Set1_DropEST.Object1.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SI.rda")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents =  c("HFHS","Control"))


seuratObject=dropEST.combined.filtered
DefaultAssay(seuratObject) <- "RNA"
#metaData1=NULL
#metaDataSubset=NULL
cellTypeMetaData="celltype"
nPerm=5000
do.par=TRUE
num.cores= 4
comparisonMetaData="data.diseasestatus"
group1="HFHS"
group2="Control"
#distanceScorePlot=NULL
#fcPlot=NULL
#distanceSave=NULL
seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)

#Number of permutations
numPerm = nPerm
if(do.par==TRUE){
  #auto detect number of available cores
  if(num.cores == "auto"){
    cores=detectCores()
    #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
    cores = max(1,cores[1]-1)
    #setup and start parallel processing cluster
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else{
    #if a number of cores for parallel processing is specified, use this number of cores
    if(is.numeric(num.cores)){
      cores = num.cores
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
  }
}
cellTypePermutationDistances = list()
cellTypeDistances = c()
allcelltypes = levels(seuratObject@active.ident)
availcelltypes <- NULL
for(i in 1:length(allcelltypes)){
  cellType = allcelltypes[i]
  print(cellType)
  #subset to cell type of interest
  cellTypeSubset = subset(seuratObject,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
  
  group1Length = length(which(cellTypeSubset@active.ident==group1))
  group2length = length(which(cellTypeSubset@active.ident==group2))
  
  #Only continue if we have at least 3 cells of the cell type in each group
  if(group1Length > 9 & group2length > 9){
    availcelltypes <- c(availcelltypes, cellType)
    #Label group1 and group2 cells
    cellTypeSubset <- ScaleData(cellTypeSubset, vars.to.regress = NULL)
    
    group1Cells = WhichCells(cellTypeSubset,idents = group1)
    group2Cells = WhichCells(cellTypeSubset,idents = group2)
    
    #Seletcing gene subsets for distance analysis
    AverageExpall <- rowMeans(as.matrix(cellTypeSubset@assays$RNA@counts))
    AverageExpall <- sort(AverageExpall, decreasing = T)
    AverageExpall <- AverageExpall[AverageExpall > quantile(AverageExpall,probs = 0.9)]
    #AverageExpall <- AverageExpall[-c(1:20)]
    selectedgenes <- names(AverageExpall)
    #Get the average expression of group1 and group2 cells
    #compute in non-log space
    #avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes))
    avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes,use.scale = T))
    
    avgExp = avgExp$RNA
    
    #calculate the distance between the group1 and group2 cells
    cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
    
    cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
    
    #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
    allCells = c(group1Cells,group2Cells)
    
    #set up a progress bar
    pb <- txtProgressBar(max = numPerm, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
      set.seed(j)
      permSubset = cellTypeSubset
      sampleCells = sample(allCells,length(allCells),replace = FALSE)
      #place the cells into either group1 or group2 (same size as original group1 and group2)
      group1Synth = sampleCells[1:length(group1Cells)]
      group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
      
      permStatus = permSubset@active.ident
      permStatus[group1Synth] = group1
      permStatus[group2Synth] = group2
      permSubset$permStatus = permStatus
      permSubset = SetIdent(permSubset, value = "permStatus")
      
      #compute in non-log space
      #permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes))
      permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes, use.scale = T))
      
      permAvgExp = permAvgExp$RNA
      
      permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
      permCellTypeDistance
    }
    cellTypePermutationDistances[[cellType]] = permutationDistances
    close(pb)
  }
}

stopCluster(cl)

names(cellTypeDistances) = availcelltypes

cellTypeDistancePvals = c()
for(cellType in levels(seuratObject@active.ident)){
  numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
  cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
}
names(cellTypeDistancePvals) = availcelltypes


melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
#availcelltypes = unique(melted.cellTypePermutationDistances$L1)

melted.cellTypePvals = rep(cellTypeDistancePvals[availcelltypes],each=numPerm)
melted.cellTypeDistances = rep(cellTypeDistances[availcelltypes],each=numPerm)

melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)

#Bonferroni correction
allDistancePvals = c(cellTypeDistancePvals[availcelltypes])
# all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm

#Combine neuronal subtype and cell types to plot
colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
melted.combinedPermutationDistances = melted.cellTypePermutationDistances
melted.combinedPermutationDistances$pvals = allDistancePvals
melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))


#########
#Plot
ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
  geom_density(adjust = 1) +
  facet_wrap( ~ L1, ncol=5,scales="free") +
  xlab("Euclidian Distance") +
  ylab("Density") +
  # Add line
  geom_vline(aes(xintercept=distances),
             color="blue", linetype="dashed", size=1) +
  annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)

pdf(file = paste0("./Plots/EuclideanDistanceHFHS/HFHSplot1.pdf"),height = 16,width = 16)
print(p1)
dev.off()

#########
#get the log FC of the euclidean distance versus the median of the background distribution

medianBackground = c()
for(cellType in levels(seuratObject@active.ident)){
  medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
}
names(medianBackground) = availcelltypes

logFC = log10(cellTypeDistances[availcelltypes]) - log10(medianBackground)
pval = -log10(allDistancePvals)
toPlot = as.data.frame(cbind(logFC,pval))

p2 <- ggplot(toPlot, aes(pval, logFC)) +
  geom_point(color = 'black') +
  theme_classic(base_size = 10) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 

pdf(file = paste0("./Plots/EuclideanDistanceHFHS/HFHSplot2.pdf"))
print(p2)
dev.off()

#Save incase we need to alter the figures, don't want to have to rerun everything
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresult.rda"))
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresult5000.rda"))
save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresultscore.rda"))


#SVF

rm(list=ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined/SVF_Final_Set1only.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SVF.rda")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set1")
#SI only Set2, Set1 is not making a lot sense with lowerr quality
#dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents =  c("Fruc","Control"))


seuratObject=dropEST.combined.filtered
DefaultAssay(seuratObject) <- "RNA"
#metaData1=NULL
#metaDataSubset=NULL
cellTypeMetaData="celltype"
nPerm=5000
do.par=TRUE
num.cores= 4
comparisonMetaData="data.diseasestatus"
group1="Fruc"
group2="Control"
#distanceScorePlot=NULL
#fcPlot=NULL
#distanceSave=NULL
seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)

#Number of permutations
numPerm = nPerm
if(do.par==TRUE){
  #auto detect number of available cores
  if(num.cores == "auto"){
    cores=detectCores()
    #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
    cores = max(1,cores[1]-1)
    #setup and start parallel processing cluster
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else{
    #if a number of cores for parallel processing is specified, use this number of cores
    if(is.numeric(num.cores)){
      cores = num.cores
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
  }
}
cellTypePermutationDistances = list()
cellTypeDistances = c()
allcelltypes = levels(seuratObject@active.ident)
availcelltypes <- NULL
for(i in 1:length(allcelltypes)){
  cellType = allcelltypes[i]
  print(cellType)
  #subset to cell type of interest
  cellTypeSubset = subset(seuratObject,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
  
  
  group1Length = length(which(cellTypeSubset@active.ident==group1))
  group2length = length(which(cellTypeSubset@active.ident==group2))
  
  #Only continue if we have at least 3 cells of the cell type in each group
  if(group1Length > 9 & group2length > 9){
    availcelltypes <- c(availcelltypes, cellType)
    #Label group1 and group2 cells
    group1Cells = WhichCells(cellTypeSubset,idents = group1)
    group2Cells = WhichCells(cellTypeSubset,idents = group2)
    
    #Get the average expression of group1 and group2 cells
    cellTypeSubset <- ScaleData(cellTypeSubset, vars.to.regress = NULL)
    
    AverageExpall <- rowMeans(as.matrix(cellTypeSubset@assays$RNA@counts))
    AverageExpall <- sort(AverageExpall, decreasing = T)
    AverageExpall <- AverageExpall[AverageExpall > quantile(AverageExpall,probs = 0.9)]
    #AverageExpall <- AverageExpall[-c(1:20)]
    selectedgenes <- names(AverageExpall)
    selectedgenes <- selectedgenes[1:1000]
    #compute in non-log space
    #avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes))
    
    avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes, use.scale = T))
    avgExp = avgExp$RNA
    #top and bottom filtration
    #calculate the distance between the group1 and group2 cells
    cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
    
    cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
    
    #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
    allCells = c(group1Cells,group2Cells)
    
    #set up a progress bar
    pb <- txtProgressBar(max = numPerm, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
      set.seed(j)
      permSubset = cellTypeSubset
      sampleCells = sample(allCells,length(allCells),replace = FALSE)
      #place the cells into either group1 or group2 (same size as original group1 and group2)
      group1Synth = sampleCells[1:length(group1Cells)]
      group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
      
      permStatus = permSubset@active.ident
      permStatus[group1Synth] = group1
      permStatus[group2Synth] = group2
      permSubset$permStatus = permStatus
      permSubset = SetIdent(permSubset, value = "permStatus")
      
      #compute in non-log space
      #permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes))
      
      permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes, use.scale = T))
      permAvgExp = permAvgExp$RNA
      
      permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
      permCellTypeDistance
    }
    cellTypePermutationDistances[[cellType]] = permutationDistances
    close(pb)
  }
}

stopCluster(cl)

names(cellTypeDistances) = availcelltypes

cellTypeDistancePvals = c()
for(cellType in levels(seuratObject@active.ident)){
  numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
  cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
}
names(cellTypeDistancePvals) = availcelltypes


melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
#availcelltypes = unique(melted.cellTypePermutationDistances$L1)

melted.cellTypePvals = rep(cellTypeDistancePvals[availcelltypes],each=numPerm)
melted.cellTypeDistances = rep(cellTypeDistances[availcelltypes],each=numPerm)

melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)

#Bonferroni correction
allDistancePvals = c(cellTypeDistancePvals[availcelltypes])
# all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm

#Combine neuronal subtype and cell types to plot
colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
melted.combinedPermutationDistances = melted.cellTypePermutationDistances
melted.combinedPermutationDistances$pvals = allDistancePvals
melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))


#########
#Plot
ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
  geom_density(adjust = 1) +
  facet_wrap( ~ L1, ncol=5,scales="free") +
  xlab("Euclidian Distance") +
  ylab("Density") +
  # Add line
  geom_vline(aes(xintercept=distances),
             color="blue", linetype="dashed", size=1) +
  annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)

pdf(file = paste0("./Plots/EuclideanDistanceFructose/Fructoseplot1.pdf"),height = 16,width = 16)
print(p1)
dev.off()

#########
#get the log FC of the euclidean distance versus the median of the background distribution

medianBackground = c()
for(cellType in levels(seuratObject@active.ident)){
  medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
}
names(medianBackground) = availcelltypes

logFC = log10(cellTypeDistances[availcelltypes]) - log10(medianBackground)
pval = -log10(allDistancePvals)
toPlot = as.data.frame(cbind(logFC,pval))

p2 <- ggplot(toPlot, aes(pval, logFC)) +
  geom_point(color = 'black') +
  theme_classic(base_size = 10) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 

pdf(file = paste0("./Plots/EuclideanDistanceFructose/Fructoseplot2.pdf"))
print(p2)
dev.off()

#Save incase we need to alter the figures, don't want to have to rerun everything
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Fructosedistanceresult.rda"))
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Fructosedistanceresult5000.rda"))
save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Fructosedistanceresultzscore.rda"))


rm(list=ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined/SVF_Final_Set1only.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SVF.rda")
#Focusing on Set1 only, except SVF since almost no batch effect in SVF
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set1")
#dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents =  c("HFHS","Control"))

seuratObject=dropEST.combined.filtered
DefaultAssay(seuratObject) <- "RNA"
#metaData1=NULL
#metaDataSubset=NULL
cellTypeMetaData="celltype"
nPerm=5000
do.par=TRUE
num.cores= 4
comparisonMetaData="data.diseasestatus"
group1="HFHS"
group2="Control"
#distanceScorePlot=NULL
#fcPlot=NULL
#distanceSave=NULL
seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)

#Number of permutations
numPerm = nPerm
if(do.par==TRUE){
  #auto detect number of available cores
  if(num.cores == "auto"){
    cores=detectCores()
    #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
    cores = max(1,cores[1]-1)
    #setup and start parallel processing cluster
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else{
    #if a number of cores for parallel processing is specified, use this number of cores
    if(is.numeric(num.cores)){
      cores = num.cores
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
  }
}
cellTypePermutationDistances = list()
cellTypeDistances = c()
allcelltypes = levels(seuratObject@active.ident)
availcelltypes <- NULL
for(i in 1:length(allcelltypes)){
  cellType = allcelltypes[i]
  print(cellType)
  #subset to cell type of interest
  cellTypeSubset = subset(seuratObject,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
  
  group1Length = length(which(cellTypeSubset@active.ident==group1))
  group2length = length(which(cellTypeSubset@active.ident==group2))
  
  #Only continue if we have at least 3 cells of the cell type in each group
  if(group1Length > 9 & group2length > 9){
    availcelltypes <- c(availcelltypes, cellType)
    #Label group1 and group2 cells
    cellTypeSubset <- ScaleData(cellTypeSubset, vars.to.regress = NULL)
    
    group1Cells = WhichCells(cellTypeSubset,idents = group1)
    group2Cells = WhichCells(cellTypeSubset,idents = group2)
    
    #Seletcing gene subsets for distance analysis
    AverageExpall <- rowMeans(as.matrix(cellTypeSubset@assays$RNA@counts))
    AverageExpall <- sort(AverageExpall, decreasing = T)
    AverageExpall <- AverageExpall[AverageExpall > quantile(AverageExpall,probs = 0.9)]
    #AverageExpall <- AverageExpall[-c(1:20)]
    selectedgenes <- names(AverageExpall)
    #Get the average expression of group1 and group2 cells
    #compute in non-log space
    #avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes))
    avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes,use.scale = T))
    
    avgExp = avgExp$RNA
    
    #calculate the distance between the group1 and group2 cells
    cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
    
    cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
    
    #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
    allCells = c(group1Cells,group2Cells)
    
    #set up a progress bar
    pb <- txtProgressBar(max = numPerm, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
      set.seed(j)
      permSubset = cellTypeSubset
      sampleCells = sample(allCells,length(allCells),replace = FALSE)
      #place the cells into either group1 or group2 (same size as original group1 and group2)
      group1Synth = sampleCells[1:length(group1Cells)]
      group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
      
      permStatus = permSubset@active.ident
      permStatus[group1Synth] = group1
      permStatus[group2Synth] = group2
      permSubset$permStatus = permStatus
      permSubset = SetIdent(permSubset, value = "permStatus")
      
      #compute in non-log space
      #permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes))
      permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes, use.scale = T))
      
      permAvgExp = permAvgExp$RNA
      
      permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
      permCellTypeDistance
    }
    cellTypePermutationDistances[[cellType]] = permutationDistances
    close(pb)
  }
}

stopCluster(cl)

names(cellTypeDistances) = availcelltypes

cellTypeDistancePvals = c()
for(cellType in levels(seuratObject@active.ident)){
  numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
  cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
}
names(cellTypeDistancePvals) = availcelltypes


melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
#availcelltypes = unique(melted.cellTypePermutationDistances$L1)

melted.cellTypePvals = rep(cellTypeDistancePvals[availcelltypes],each=numPerm)
melted.cellTypeDistances = rep(cellTypeDistances[availcelltypes],each=numPerm)

melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)

#Bonferroni correction
allDistancePvals = c(cellTypeDistancePvals[availcelltypes])
# all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm

#Combine neuronal subtype and cell types to plot
colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
melted.combinedPermutationDistances = melted.cellTypePermutationDistances
melted.combinedPermutationDistances$pvals = allDistancePvals
melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))


#########
#Plot
ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
  geom_density(adjust = 1) +
  facet_wrap( ~ L1, ncol=5,scales="free") +
  xlab("Euclidian Distance") +
  ylab("Density") +
  # Add line
  geom_vline(aes(xintercept=distances),
             color="blue", linetype="dashed", size=1) +
  annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)

pdf(file = paste0("./Plots/EuclideanDistanceHFHS/HFHSplot1.pdf"),height = 16,width = 16)
print(p1)
dev.off()

#########
#get the log FC of the euclidean distance versus the median of the background distribution

medianBackground = c()
for(cellType in levels(seuratObject@active.ident)){
  medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
}
names(medianBackground) = availcelltypes

logFC = log10(cellTypeDistances[availcelltypes]) - log10(medianBackground)
pval = -log10(allDistancePvals)
toPlot = as.data.frame(cbind(logFC,pval))

p2 <- ggplot(toPlot, aes(pval, logFC)) +
  geom_point(color = 'black') +
  theme_classic(base_size = 10) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 

pdf(file = paste0("./Plots/EuclideanDistanceHFHS/HFHSplot2.pdf"))
print(p2)
dev.off()

#Save incase we need to alter the figures, don't want to have to rerun everything
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresult.rda"))
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresult5000.rda"))
save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresultscore.rda"))

#liver combined, only use first batch
rm(list=ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Combined_allliver.rda")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "Sinusoidal endothelial cells","1" = "Periportal hepatocytes","2" = "Periportal hepatocytes","3" = "Kupffer", "4" = "Pericentral hepatocytes","5" = "Periportal hepatocytes", "6" = "NKT cells",
#                                                    "7" = "Periportal hepatocytes","8" = "Classical dendritic cells","9" = "B cells","10" = "Hepatic stellate cells","11" = "Cholangiocytes","12" = "Dividing cells","13" = "Plasmacytoid dendritic cells"))
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/liver.rda")
#Focusing on Set1 only, 
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set1")
#SI only Set2, Set1 is not making a lot sense with lowerr quality
#dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents =  c("Fruc","Control"))


seuratObject=dropEST.combined.filtered
DefaultAssay(seuratObject) <- "RNA"
#metaData1=NULL
#metaDataSubset=NULL
cellTypeMetaData="celltype"
nPerm=5000
do.par=TRUE
num.cores= 4
comparisonMetaData="data.diseasestatus"
group1="Fruc"
group2="Control"
#distanceScorePlot=NULL
#fcPlot=NULL
#distanceSave=NULL
seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)

#Number of permutations
numPerm = nPerm
if(do.par==TRUE){
  #auto detect number of available cores
  if(num.cores == "auto"){
    cores=detectCores()
    #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
    cores = max(1,cores[1]-1)
    #setup and start parallel processing cluster
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else{
    #if a number of cores for parallel processing is specified, use this number of cores
    if(is.numeric(num.cores)){
      cores = num.cores
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
  }
}
cellTypePermutationDistances = list()
cellTypeDistances = c()
allcelltypes = levels(seuratObject@active.ident)
availcelltypes <- NULL
for(i in 1:length(allcelltypes)){
  cellType = allcelltypes[i]
  print(cellType)
  #subset to cell type of interest
  cellTypeSubset = subset(seuratObject,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
  
  
  group1Length = length(which(cellTypeSubset@active.ident==group1))
  group2length = length(which(cellTypeSubset@active.ident==group2))
  
  #Only continue if we have at least 3 cells of the cell type in each group
  if(group1Length > 9 & group2length > 9){
    availcelltypes <- c(availcelltypes, cellType)
    #Label group1 and group2 cells
    group1Cells = WhichCells(cellTypeSubset,idents = group1)
    group2Cells = WhichCells(cellTypeSubset,idents = group2)
    
    #Get the average expression of group1 and group2 cells
    cellTypeSubset <- ScaleData(cellTypeSubset, vars.to.regress = NULL)
    
    AverageExpall <- rowMeans(as.matrix(cellTypeSubset@assays$RNA@counts))
    AverageExpall <- sort(AverageExpall, decreasing = T)
    AverageExpall <- AverageExpall[AverageExpall > quantile(AverageExpall,probs = 0.9)]
    #AverageExpall <- AverageExpall[-c(1:20)]
    selectedgenes <- names(AverageExpall)
    selectedgenes <- selectedgenes[1:1000]
    #compute in non-log space
    #avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes))
    
    avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes, use.scale = T))
    avgExp = avgExp$RNA
    #top and bottom filtration
    #calculate the distance between the group1 and group2 cells
    cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
    
    cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
    
    #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
    allCells = c(group1Cells,group2Cells)
    
    #set up a progress bar
    pb <- txtProgressBar(max = numPerm, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
      set.seed(j)
      permSubset = cellTypeSubset
      sampleCells = sample(allCells,length(allCells),replace = FALSE)
      #place the cells into either group1 or group2 (same size as original group1 and group2)
      group1Synth = sampleCells[1:length(group1Cells)]
      group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
      
      permStatus = permSubset@active.ident
      permStatus[group1Synth] = group1
      permStatus[group2Synth] = group2
      permSubset$permStatus = permStatus
      permSubset = SetIdent(permSubset, value = "permStatus")
      
      #compute in non-log space
      #permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes))
      
      permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes, use.scale = T))
      permAvgExp = permAvgExp$RNA
      
      permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
      permCellTypeDistance
    }
    cellTypePermutationDistances[[cellType]] = permutationDistances
    close(pb)
  }
}

stopCluster(cl)

names(cellTypeDistances) = availcelltypes

cellTypeDistancePvals = c()
for(cellType in levels(seuratObject@active.ident)){
  numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
  cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
}
names(cellTypeDistancePvals) = availcelltypes


melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
#availcelltypes = unique(melted.cellTypePermutationDistances$L1)

melted.cellTypePvals = rep(cellTypeDistancePvals[availcelltypes],each=numPerm)
melted.cellTypeDistances = rep(cellTypeDistances[availcelltypes],each=numPerm)

melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)

#Bonferroni correction
allDistancePvals = c(cellTypeDistancePvals[availcelltypes])
# all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm

#Combine neuronal subtype and cell types to plot
colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
melted.combinedPermutationDistances = melted.cellTypePermutationDistances
melted.combinedPermutationDistances$pvals = allDistancePvals
melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))


#########
#Plot
ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
  geom_density(adjust = 1) +
  facet_wrap( ~ L1, ncol=5,scales="free") +
  xlab("Euclidian Distance") +
  ylab("Density") +
  # Add line
  geom_vline(aes(xintercept=distances),
             color="blue", linetype="dashed", size=1) +
  annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)

pdf(file = paste0("./Plots/EuclideanDistanceFructose/Fructoseplot1.pdf"),height = 16,width = 16)
print(p1)
dev.off()

#########
#get the log FC of the euclidean distance versus the median of the background distribution

medianBackground = c()
for(cellType in levels(seuratObject@active.ident)){
  medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
}
names(medianBackground) = availcelltypes

logFC = log10(cellTypeDistances[availcelltypes]) - log10(medianBackground)
pval = -log10(allDistancePvals)
toPlot = as.data.frame(cbind(logFC,pval))

p2 <- ggplot(toPlot, aes(pval, logFC)) +
  geom_point(color = 'black') +
  theme_classic(base_size = 10) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 

pdf(file = paste0("./Plots/EuclideanDistanceFructose/Fructoseplot2.pdf"))
print(p2)
dev.off()

#Save incase we need to alter the figures, don't want to have to rerun everything
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Fructosedistanceresult.rda"))
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Fructosedistanceresult5000.rda"))
save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Set1_Fructosedistanceresultzscore.rda"))


#NPC_liver Set2
rm(list=ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Combined_allliver.rda")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "Sinusoidal endothelial cells","1" = "Periportal hepatocytes","2" = "Periportal hepatocytes","3" = "Kupffer", "4" = "Pericentral hepatocytes","5" = "Periportal hepatocytes", "6" = "NKT cells",
#                                                    "7" = "Periportal hepatocytes","8" = "Classical dendritic cells","9" = "B cells","10" = "Hepatic stellate cells","11" = "Cholangiocytes","12" = "Dividing cells","13" = "Plasmacytoid dendritic cells"))
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/liver.rda")
#Focusing on Set1 only, 
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
#SI only Set2, Set1 is not making a lot sense with lowerr quality
#dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents =  c("Fruc","Control"))


seuratObject=dropEST.combined.filtered
DefaultAssay(seuratObject) <- "RNA"
#metaData1=NULL
#metaDataSubset=NULL
cellTypeMetaData="celltype"
nPerm=5000
do.par=TRUE
num.cores= 4
comparisonMetaData="data.diseasestatus"
group1="Fruc"
group2="Control"
#distanceScorePlot=NULL
#fcPlot=NULL
#distanceSave=NULL
seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)

#Number of permutations
numPerm = nPerm
if(do.par==TRUE){
  #auto detect number of available cores
  if(num.cores == "auto"){
    cores=detectCores()
    #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
    cores = max(1,cores[1]-1)
    #setup and start parallel processing cluster
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else{
    #if a number of cores for parallel processing is specified, use this number of cores
    if(is.numeric(num.cores)){
      cores = num.cores
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
  }
}
cellTypePermutationDistances = list()
cellTypeDistances = c()
allcelltypes = levels(seuratObject@active.ident)
availcelltypes <- NULL
for(i in 1:length(allcelltypes)){
  cellType = allcelltypes[i]
  print(cellType)
  #subset to cell type of interest
  cellTypeSubset = subset(seuratObject,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
  
  
  group1Length = length(which(cellTypeSubset@active.ident==group1))
  group2length = length(which(cellTypeSubset@active.ident==group2))
  
  #Only continue if we have at least 3 cells of the cell type in each group
  if(group1Length > 9 & group2length > 9){
    availcelltypes <- c(availcelltypes, cellType)
    #Label group1 and group2 cells
    group1Cells = WhichCells(cellTypeSubset,idents = group1)
    group2Cells = WhichCells(cellTypeSubset,idents = group2)
    
    #Get the average expression of group1 and group2 cells
    cellTypeSubset <- ScaleData(cellTypeSubset, vars.to.regress = NULL)
    
    AverageExpall <- rowMeans(as.matrix(cellTypeSubset@assays$RNA@counts))
    AverageExpall <- sort(AverageExpall, decreasing = T)
    AverageExpall <- AverageExpall[AverageExpall > quantile(AverageExpall,probs = 0.9)]
    #AverageExpall <- AverageExpall[-c(1:20)]
    selectedgenes <- names(AverageExpall)
    selectedgenes <- selectedgenes[1:1000]
    #compute in non-log space
    #avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes))
    
    avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes, use.scale = T))
    avgExp = avgExp$RNA
    #top and bottom filtration
    #calculate the distance between the group1 and group2 cells
    cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
    
    cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
    
    #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
    allCells = c(group1Cells,group2Cells)
    
    #set up a progress bar
    pb <- txtProgressBar(max = numPerm, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
      set.seed(j)
      permSubset = cellTypeSubset
      sampleCells = sample(allCells,length(allCells),replace = FALSE)
      #place the cells into either group1 or group2 (same size as original group1 and group2)
      group1Synth = sampleCells[1:length(group1Cells)]
      group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
      
      permStatus = permSubset@active.ident
      permStatus[group1Synth] = group1
      permStatus[group2Synth] = group2
      permSubset$permStatus = permStatus
      permSubset = SetIdent(permSubset, value = "permStatus")
      
      #compute in non-log space
      #permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes))
      
      permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes, use.scale = T))
      permAvgExp = permAvgExp$RNA
      
      permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
      permCellTypeDistance
    }
    cellTypePermutationDistances[[cellType]] = permutationDistances
    close(pb)
  }
}

stopCluster(cl)

names(cellTypeDistances) = availcelltypes

cellTypeDistancePvals = c()
for(cellType in levels(seuratObject@active.ident)){
  numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
  cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
}
names(cellTypeDistancePvals) = availcelltypes


melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
#availcelltypes = unique(melted.cellTypePermutationDistances$L1)

melted.cellTypePvals = rep(cellTypeDistancePvals[availcelltypes],each=numPerm)
melted.cellTypeDistances = rep(cellTypeDistances[availcelltypes],each=numPerm)

melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)

#Bonferroni correction
allDistancePvals = c(cellTypeDistancePvals[availcelltypes])
# all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm

#Combine neuronal subtype and cell types to plot
colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
melted.combinedPermutationDistances = melted.cellTypePermutationDistances
melted.combinedPermutationDistances$pvals = allDistancePvals
melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))


#########
#Plot
ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
  geom_density(adjust = 1) +
  facet_wrap( ~ L1, ncol=5,scales="free") +
  xlab("Euclidian Distance") +
  ylab("Density") +
  # Add line
  geom_vline(aes(xintercept=distances),
             color="blue", linetype="dashed", size=1) +
  annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)

pdf(file = paste0("./Plots/EuclideanDistanceFructose/Fructoseplot1.pdf"),height = 16,width = 16)
print(p1)
dev.off()

#########
#get the log FC of the euclidean distance versus the median of the background distribution

medianBackground = c()
for(cellType in levels(seuratObject@active.ident)){
  medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
}
names(medianBackground) = availcelltypes

logFC = log10(cellTypeDistances[availcelltypes]) - log10(medianBackground)
pval = -log10(allDistancePvals)
toPlot = as.data.frame(cbind(logFC,pval))

p2 <- ggplot(toPlot, aes(pval, logFC)) +
  geom_point(color = 'black') +
  theme_classic(base_size = 10) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 

pdf(file = paste0("./Plots/EuclideanDistanceFructose/Fructoseplot2.pdf"))
print(p2)
dev.off()

#Save incase we need to alter the figures, don't want to have to rerun everything
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Fructosedistanceresult.rda"))
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Fructosedistanceresult5000.rda"))
save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Set2_Fructosedistanceresultzscore.rda"))




#NPC liver HFHS Set1
rm(list = ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Combined_allliver.rda")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "Sinusoidal endothelial cells","1" = "Periportal hepatocytes","2" = "Periportal hepatocytes","3" = "Kupffer", "4" = "Pericentral hepatocytes","5" = "Periportal hepatocytes", "6" = "NKT cells",
#                                                    "7" = "Periportal hepatocytes","8" = "Classical dendritic cells","9" = "B cells","10" = "Hepatic stellate cells","11" = "Cholangiocytes","12" = "Dividing cells","13" = "Plasmacytoid dendritic cells"))
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/liver.rda")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set1")
#dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents =  c("HFHS","Control"))

seuratObject=dropEST.combined.filtered
DefaultAssay(seuratObject) <- "RNA"
#metaData1=NULL
#metaDataSubset=NULL
cellTypeMetaData="celltype"
nPerm=5000
do.par=TRUE
num.cores= 4
comparisonMetaData="data.diseasestatus"
group1="HFHS"
group2="Control"
#distanceScorePlot=NULL
#fcPlot=NULL
#distanceSave=NULL
seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)

#Number of permutations
numPerm = nPerm
if(do.par==TRUE){
  #auto detect number of available cores
  if(num.cores == "auto"){
    cores=detectCores()
    #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
    cores = max(1,cores[1]-1)
    #setup and start parallel processing cluster
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else{
    #if a number of cores for parallel processing is specified, use this number of cores
    if(is.numeric(num.cores)){
      cores = num.cores
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
  }
}
cellTypePermutationDistances = list()
cellTypeDistances = c()
allcelltypes = levels(seuratObject@active.ident)
availcelltypes <- NULL
for(i in 1:length(allcelltypes)){
  cellType = allcelltypes[i]
  print(cellType)
  #subset to cell type of interest
  cellTypeSubset = subset(seuratObject,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
  
  group1Length = length(which(cellTypeSubset@active.ident==group1))
  group2length = length(which(cellTypeSubset@active.ident==group2))
  
  #Only continue if we have at least 3 cells of the cell type in each group
  if(group1Length > 9 & group2length > 9){
    availcelltypes <- c(availcelltypes, cellType)
    #Label group1 and group2 cells
    cellTypeSubset <- ScaleData(cellTypeSubset, vars.to.regress = NULL)
    
    group1Cells = WhichCells(cellTypeSubset,idents = group1)
    group2Cells = WhichCells(cellTypeSubset,idents = group2)
    
    #Seletcing gene subsets for distance analysis
    AverageExpall <- rowMeans(as.matrix(cellTypeSubset@assays$RNA@counts))
    AverageExpall <- sort(AverageExpall, decreasing = T)
    AverageExpall <- AverageExpall[AverageExpall > quantile(AverageExpall,probs = 0.9)]
    #AverageExpall <- AverageExpall[-c(1:20)]
    selectedgenes <- names(AverageExpall)
    #Get the average expression of group1 and group2 cells
    #compute in non-log space
    #avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes))
    avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes,use.scale = T))
    
    avgExp = avgExp$RNA
    
    #calculate the distance between the group1 and group2 cells
    cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
    
    cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
    
    #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
    allCells = c(group1Cells,group2Cells)
    
    #set up a progress bar
    pb <- txtProgressBar(max = numPerm, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
      set.seed(j)
      permSubset = cellTypeSubset
      sampleCells = sample(allCells,length(allCells),replace = FALSE)
      #place the cells into either group1 or group2 (same size as original group1 and group2)
      group1Synth = sampleCells[1:length(group1Cells)]
      group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
      
      permStatus = permSubset@active.ident
      permStatus[group1Synth] = group1
      permStatus[group2Synth] = group2
      permSubset$permStatus = permStatus
      permSubset = SetIdent(permSubset, value = "permStatus")
      
      #compute in non-log space
      #permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes))
      permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes, use.scale = T))
      
      permAvgExp = permAvgExp$RNA
      
      permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
      permCellTypeDistance
    }
    cellTypePermutationDistances[[cellType]] = permutationDistances
    close(pb)
  }
}

stopCluster(cl)

names(cellTypeDistances) = availcelltypes

cellTypeDistancePvals = c()
for(cellType in levels(seuratObject@active.ident)){
  numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
  cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
}
names(cellTypeDistancePvals) = availcelltypes


melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
#availcelltypes = unique(melted.cellTypePermutationDistances$L1)

melted.cellTypePvals = rep(cellTypeDistancePvals[availcelltypes],each=numPerm)
melted.cellTypeDistances = rep(cellTypeDistances[availcelltypes],each=numPerm)

melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)

#Bonferroni correction
allDistancePvals = c(cellTypeDistancePvals[availcelltypes])
# all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm

#Combine neuronal subtype and cell types to plot
colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
melted.combinedPermutationDistances = melted.cellTypePermutationDistances
melted.combinedPermutationDistances$pvals = allDistancePvals
melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))


#########
#Plot
ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
  geom_density(adjust = 1) +
  facet_wrap( ~ L1, ncol=5,scales="free") +
  xlab("Euclidian Distance") +
  ylab("Density") +
  # Add line
  geom_vline(aes(xintercept=distances),
             color="blue", linetype="dashed", size=1) +
  annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)

pdf(file = paste0("./Plots/EuclideanDistanceHFHS/HFHSplot1.pdf"),height = 16,width = 16)
print(p1)
dev.off()

#########
#get the log FC of the euclidean distance versus the median of the background distribution

medianBackground = c()
for(cellType in levels(seuratObject@active.ident)){
  medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
}
names(medianBackground) = availcelltypes

logFC = log10(cellTypeDistances[availcelltypes]) - log10(medianBackground)
pval = -log10(allDistancePvals)
toPlot = as.data.frame(cbind(logFC,pval))

p2 <- ggplot(toPlot, aes(pval, logFC)) +
  geom_point(color = 'black') +
  theme_classic(base_size = 10) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 

pdf(file = paste0("./Plots/EuclideanDistanceHFHS/HFHSplot2.pdf"))
print(p2)
dev.off()

#Save incase we need to alter the figures, don't want to have to rerun everything
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresult.rda"))
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresult5000.rda"))
save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Set1_HFHSdistanceresultscore.rda"))


#NPC_liver HFHS Set2
rm(list = ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Combined_allliver.rda")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "Sinusoidal endothelial cells","1" = "Periportal hepatocytes","2" = "Periportal hepatocytes","3" = "Kupffer", "4" = "Pericentral hepatocytes","5" = "Periportal hepatocytes", "6" = "NKT cells",
#                                                    "7" = "Periportal hepatocytes","8" = "Classical dendritic cells","9" = "B cells","10" = "Hepatic stellate cells","11" = "Cholangiocytes","12" = "Dividing cells","13" = "Plasmacytoid dendritic cells"))
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/liver.rda")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
#dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents =  c("HFHS","Control"))

seuratObject=dropEST.combined.filtered
DefaultAssay(seuratObject) <- "RNA"
#metaData1=NULL
#metaDataSubset=NULL
cellTypeMetaData="celltype"
nPerm=5000
do.par=TRUE
num.cores= 4
comparisonMetaData="data.diseasestatus"
group1="HFHS"
group2="Control"
#distanceScorePlot=NULL
#fcPlot=NULL
#distanceSave=NULL
seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)

#Number of permutations
numPerm = nPerm
if(do.par==TRUE){
  #auto detect number of available cores
  if(num.cores == "auto"){
    cores=detectCores()
    #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
    cores = max(1,cores[1]-1)
    #setup and start parallel processing cluster
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else{
    #if a number of cores for parallel processing is specified, use this number of cores
    if(is.numeric(num.cores)){
      cores = num.cores
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
  }
}
cellTypePermutationDistances = list()
cellTypeDistances = c()
allcelltypes = levels(seuratObject@active.ident)
availcelltypes <- NULL
for(i in 1:length(allcelltypes)){
  cellType = allcelltypes[i]
  print(cellType)
  #subset to cell type of interest
  cellTypeSubset = subset(seuratObject,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
  
  group1Length = length(which(cellTypeSubset@active.ident==group1))
  group2length = length(which(cellTypeSubset@active.ident==group2))
  
  #Only continue if we have at least 3 cells of the cell type in each group
  if(group1Length > 9 & group2length > 9){
    availcelltypes <- c(availcelltypes, cellType)
    #Label group1 and group2 cells
    cellTypeSubset <- ScaleData(cellTypeSubset, vars.to.regress = NULL)
    
    group1Cells = WhichCells(cellTypeSubset,idents = group1)
    group2Cells = WhichCells(cellTypeSubset,idents = group2)
    
    #Seletcing gene subsets for distance analysis
    AverageExpall <- rowMeans(as.matrix(cellTypeSubset@assays$RNA@counts))
    AverageExpall <- sort(AverageExpall, decreasing = T)
    AverageExpall <- AverageExpall[AverageExpall > quantile(AverageExpall,probs = 0.9)]
    #AverageExpall <- AverageExpall[-c(1:20)]
    selectedgenes <- names(AverageExpall)
    #Get the average expression of group1 and group2 cells
    #compute in non-log space
    #avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes))
    avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes,use.scale = T))
    
    avgExp = avgExp$RNA
    
    #calculate the distance between the group1 and group2 cells
    cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
    
    cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
    
    #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
    allCells = c(group1Cells,group2Cells)
    
    #set up a progress bar
    pb <- txtProgressBar(max = numPerm, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
      set.seed(j)
      permSubset = cellTypeSubset
      sampleCells = sample(allCells,length(allCells),replace = FALSE)
      #place the cells into either group1 or group2 (same size as original group1 and group2)
      group1Synth = sampleCells[1:length(group1Cells)]
      group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
      
      permStatus = permSubset@active.ident
      permStatus[group1Synth] = group1
      permStatus[group2Synth] = group2
      permSubset$permStatus = permStatus
      permSubset = SetIdent(permSubset, value = "permStatus")
      
      #compute in non-log space
      #permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes))
      permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes, use.scale = T))
      
      permAvgExp = permAvgExp$RNA
      
      permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
      permCellTypeDistance
    }
    cellTypePermutationDistances[[cellType]] = permutationDistances
    close(pb)
  }
}

stopCluster(cl)

names(cellTypeDistances) = availcelltypes

cellTypeDistancePvals = c()
for(cellType in levels(seuratObject@active.ident)){
  numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
  cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
}
names(cellTypeDistancePvals) = availcelltypes


melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
#availcelltypes = unique(melted.cellTypePermutationDistances$L1)

melted.cellTypePvals = rep(cellTypeDistancePvals[availcelltypes],each=numPerm)
melted.cellTypeDistances = rep(cellTypeDistances[availcelltypes],each=numPerm)

melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)

#Bonferroni correction
allDistancePvals = c(cellTypeDistancePvals[availcelltypes])
# all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm

#Combine neuronal subtype and cell types to plot
colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
melted.combinedPermutationDistances = melted.cellTypePermutationDistances
melted.combinedPermutationDistances$pvals = allDistancePvals
melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))


#########
#Plot
ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
  geom_density(adjust = 1) +
  facet_wrap( ~ L1, ncol=5,scales="free") +
  xlab("Euclidian Distance") +
  ylab("Density") +
  # Add line
  geom_vline(aes(xintercept=distances),
             color="blue", linetype="dashed", size=1) +
  annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)

pdf(file = paste0("./Plots/EuclideanDistanceHFHS/HFHSplot1.pdf"),height = 16,width = 16)
print(p1)
dev.off()

#########
#get the log FC of the euclidean distance versus the median of the background distribution

medianBackground = c()
for(cellType in levels(seuratObject@active.ident)){
  medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
}
names(medianBackground) = availcelltypes

logFC = log10(cellTypeDistances[availcelltypes]) - log10(medianBackground)
pval = -log10(allDistancePvals)
toPlot = as.data.frame(cbind(logFC,pval))

p2 <- ggplot(toPlot, aes(pval, logFC)) +
  geom_point(color = 'black') +
  theme_classic(base_size = 10) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 

pdf(file = paste0("./Plots/EuclideanDistanceHFHS/HFHSplot2.pdf"))
print(p2)
dev.off()

#Save incase we need to alter the figures, don't want to have to rerun everything
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresult.rda"))
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresult5000.rda"))
save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Set2_HFHSdistanceresultscore.rda"))




rm(list = ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Hyp_combined_Set1_final.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/Hyp.rda")
#Focusing on Set1 only, 
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set1")
#SI only Set2, Set1 is not making a lot sense with lowerr quality
#dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents =  c("Fruc","Control"))


seuratObject=dropEST.combined.filtered
DefaultAssay(seuratObject) <- "RNA"
#metaData1=NULL
#metaDataSubset=NULL
cellTypeMetaData="celltype"
nPerm=5000
do.par=TRUE
num.cores= 4
comparisonMetaData="data.diseasestatus"
group1="Fruc"
group2="Control"
#distanceScorePlot=NULL
#fcPlot=NULL
#distanceSave=NULL
seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)

#Number of permutations
numPerm = nPerm
if(do.par==TRUE){
  #auto detect number of available cores
  if(num.cores == "auto"){
    cores=detectCores()
    #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
    cores = max(1,cores[1]-1)
    #setup and start parallel processing cluster
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else{
    #if a number of cores for parallel processing is specified, use this number of cores
    if(is.numeric(num.cores)){
      cores = num.cores
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
  }
}
cellTypePermutationDistances = list()
cellTypeDistances = c()
allcelltypes = levels(seuratObject@active.ident)
availcelltypes <- NULL
for(i in 1:length(allcelltypes)){
  cellType = allcelltypes[i]
  print(cellType)
  #subset to cell type of interest
  cellTypeSubset = subset(seuratObject,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
  
  
  group1Length = length(which(cellTypeSubset@active.ident==group1))
  group2length = length(which(cellTypeSubset@active.ident==group2))
  
  #Only continue if we have at least 3 cells of the cell type in each group
  if(group1Length > 9 & group2length > 9){
    availcelltypes <- c(availcelltypes, cellType)
    #Label group1 and group2 cells
    group1Cells = WhichCells(cellTypeSubset,idents = group1)
    group2Cells = WhichCells(cellTypeSubset,idents = group2)
    
    #Get the average expression of group1 and group2 cells
    cellTypeSubset <- ScaleData(cellTypeSubset, vars.to.regress = NULL)
    
    AverageExpall <- rowMeans(as.matrix(cellTypeSubset@assays$RNA@counts))
    AverageExpall <- sort(AverageExpall, decreasing = T)
    AverageExpall <- AverageExpall[AverageExpall > quantile(AverageExpall,probs = 0.9)]
    #AverageExpall <- AverageExpall[-c(1:20)]
    selectedgenes <- names(AverageExpall)
    selectedgenes <- selectedgenes[1:1000]
    #compute in non-log space
    #avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes))
    
    avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes, use.scale = T))
    avgExp = avgExp$RNA
    #top and bottom filtration
    #calculate the distance between the group1 and group2 cells
    cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
    
    cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
    
    #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
    allCells = c(group1Cells,group2Cells)
    
    #set up a progress bar
    pb <- txtProgressBar(max = numPerm, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
      set.seed(j)
      permSubset = cellTypeSubset
      sampleCells = sample(allCells,length(allCells),replace = FALSE)
      #place the cells into either group1 or group2 (same size as original group1 and group2)
      group1Synth = sampleCells[1:length(group1Cells)]
      group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
      
      permStatus = permSubset@active.ident
      permStatus[group1Synth] = group1
      permStatus[group2Synth] = group2
      permSubset$permStatus = permStatus
      permSubset = SetIdent(permSubset, value = "permStatus")
      
      #compute in non-log space
      #permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes))
      
      permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes, use.scale = T))
      permAvgExp = permAvgExp$RNA
      
      permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
      permCellTypeDistance
    }
    cellTypePermutationDistances[[cellType]] = permutationDistances
    close(pb)
  }
}

stopCluster(cl)

names(cellTypeDistances) = availcelltypes

cellTypeDistancePvals = c()
for(cellType in levels(seuratObject@active.ident)){
  numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
  cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
}
names(cellTypeDistancePvals) = availcelltypes


melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
#availcelltypes = unique(melted.cellTypePermutationDistances$L1)

melted.cellTypePvals = rep(cellTypeDistancePvals[availcelltypes],each=numPerm)
melted.cellTypeDistances = rep(cellTypeDistances[availcelltypes],each=numPerm)

melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)

#Bonferroni correction
allDistancePvals = c(cellTypeDistancePvals[availcelltypes])
# all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm

#Combine neuronal subtype and cell types to plot
colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
melted.combinedPermutationDistances = melted.cellTypePermutationDistances
melted.combinedPermutationDistances$pvals = allDistancePvals
melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))


#########
#Plot
ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
  geom_density(adjust = 1) +
  facet_wrap( ~ L1, ncol=5,scales="free") +
  xlab("Euclidian Distance") +
  ylab("Density") +
  # Add line
  geom_vline(aes(xintercept=distances),
             color="blue", linetype="dashed", size=1) +
  annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)

pdf(file = paste0("./Plots/EuclideanDistanceFructose/Fructoseplot1.pdf"),height = 16,width = 16)
print(p1)
dev.off()

#########
#get the log FC of the euclidean distance versus the median of the background distribution

medianBackground = c()
for(cellType in levels(seuratObject@active.ident)){
  medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
}
names(medianBackground) = availcelltypes

logFC = log10(cellTypeDistances[availcelltypes]) - log10(medianBackground)
pval = -log10(allDistancePvals)
toPlot = as.data.frame(cbind(logFC,pval))

p2 <- ggplot(toPlot, aes(pval, logFC)) +
  geom_point(color = 'black') +
  theme_classic(base_size = 10) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 

pdf(file = paste0("./Plots/EuclideanDistanceFructose/Fructoseplot2.pdf"))
print(p2)
dev.off()

#Save incase we need to alter the figures, don't want to have to rerun everything
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Fructosedistanceresult.rda"))
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Fructosedistanceresult5000.rda"))
save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Fructosedistanceresultzscore.rda"))
rm(list = ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Hyp_combined_Set1_final.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/Hyp.rda")
#Focusing on Set1 only, except SVF since almost no batch effect in SVF
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set1")
#dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents =  c("HFHS","Control"))

seuratObject=dropEST.combined.filtered
DefaultAssay(seuratObject) <- "RNA"
#metaData1=NULL
#metaDataSubset=NULL
cellTypeMetaData="celltype"
nPerm=5000
do.par=TRUE
num.cores= 4
comparisonMetaData="data.diseasestatus"
group1="HFHS"
group2="Control"
#distanceScorePlot=NULL
#fcPlot=NULL
#distanceSave=NULL
seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)

#Number of permutations
numPerm = nPerm
if(do.par==TRUE){
  #auto detect number of available cores
  if(num.cores == "auto"){
    cores=detectCores()
    #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
    cores = max(1,cores[1]-1)
    #setup and start parallel processing cluster
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else{
    #if a number of cores for parallel processing is specified, use this number of cores
    if(is.numeric(num.cores)){
      cores = num.cores
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
  }
}
cellTypePermutationDistances = list()
cellTypeDistances = c()
allcelltypes = levels(seuratObject@active.ident)
availcelltypes <- NULL
for(i in 1:length(allcelltypes)){
  cellType = allcelltypes[i]
  print(cellType)
  #subset to cell type of interest
  cellTypeSubset = subset(seuratObject,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
  
  group1Length = length(which(cellTypeSubset@active.ident==group1))
  group2length = length(which(cellTypeSubset@active.ident==group2))
  
  #Only continue if we have at least 3 cells of the cell type in each group
  if(group1Length > 9 & group2length > 9){
    availcelltypes <- c(availcelltypes, cellType)
    #Label group1 and group2 cells
    cellTypeSubset <- ScaleData(cellTypeSubset, vars.to.regress = NULL)
    
    group1Cells = WhichCells(cellTypeSubset,idents = group1)
    group2Cells = WhichCells(cellTypeSubset,idents = group2)
    
    #Seletcing gene subsets for distance analysis
    AverageExpall <- rowMeans(as.matrix(cellTypeSubset@assays$RNA@counts))
    AverageExpall <- sort(AverageExpall, decreasing = T)
    AverageExpall <- AverageExpall[AverageExpall > quantile(AverageExpall,probs = 0.9)]
    #AverageExpall <- AverageExpall[-c(1:20)]
    selectedgenes <- names(AverageExpall)
    #Get the average expression of group1 and group2 cells
    #compute in non-log space
    #avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes))
    avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes,use.scale = T))
    
    avgExp = avgExp$RNA
    
    #calculate the distance between the group1 and group2 cells
    cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
    
    cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
    
    #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
    allCells = c(group1Cells,group2Cells)
    
    #set up a progress bar
    pb <- txtProgressBar(max = numPerm, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
      set.seed(j)
      permSubset = cellTypeSubset
      sampleCells = sample(allCells,length(allCells),replace = FALSE)
      #place the cells into either group1 or group2 (same size as original group1 and group2)
      group1Synth = sampleCells[1:length(group1Cells)]
      group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
      
      permStatus = permSubset@active.ident
      permStatus[group1Synth] = group1
      permStatus[group2Synth] = group2
      permSubset$permStatus = permStatus
      permSubset = SetIdent(permSubset, value = "permStatus")
      
      #compute in non-log space
      #permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes))
      permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes, use.scale = T))
      
      permAvgExp = permAvgExp$RNA
      
      permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
      permCellTypeDistance
    }
    cellTypePermutationDistances[[cellType]] = permutationDistances
    close(pb)
  }
}

stopCluster(cl)

names(cellTypeDistances) = availcelltypes

cellTypeDistancePvals = c()
for(cellType in levels(seuratObject@active.ident)){
  numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
  cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
}
names(cellTypeDistancePvals) = availcelltypes


melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
#availcelltypes = unique(melted.cellTypePermutationDistances$L1)

melted.cellTypePvals = rep(cellTypeDistancePvals[availcelltypes],each=numPerm)
melted.cellTypeDistances = rep(cellTypeDistances[availcelltypes],each=numPerm)

melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)

#Bonferroni correction
allDistancePvals = c(cellTypeDistancePvals[availcelltypes])
# all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm

#Combine neuronal subtype and cell types to plot
colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
melted.combinedPermutationDistances = melted.cellTypePermutationDistances
melted.combinedPermutationDistances$pvals = allDistancePvals
melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))


#########
#Plot
ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
  geom_density(adjust = 1) +
  facet_wrap( ~ L1, ncol=5,scales="free") +
  xlab("Euclidian Distance") +
  ylab("Density") +
  # Add line
  geom_vline(aes(xintercept=distances),
             color="blue", linetype="dashed", size=1) +
  annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)

pdf(file = paste0("./Plots/EuclideanDistanceHFHS/HFHSplot1.pdf"),height = 16,width = 16)
print(p1)
dev.off()

#########
#get the log FC of the euclidean distance versus the median of the background distribution

medianBackground = c()
for(cellType in levels(seuratObject@active.ident)){
  medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
}
names(medianBackground) = availcelltypes

logFC = log10(cellTypeDistances[availcelltypes]) - log10(medianBackground)
pval = -log10(allDistancePvals)
toPlot = as.data.frame(cbind(logFC,pval))

p2 <- ggplot(toPlot, aes(pval, logFC)) +
  geom_point(color = 'black') +
  theme_classic(base_size = 10) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 

pdf(file = paste0("./Plots/EuclideanDistanceHFHS/HFHSplot2.pdf"))
print(p2)
dev.off()

#Save incase we need to alter the figures, don't want to have to rerun everything
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresult.rda"))
#save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresult5000.rda"))
save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresultscore.rda"))



rm(list = ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Neuroncelltype.rda")
dropEST.combined.filtered <- CellTypeSubset
dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.1.5,
                                                  c("0" = "General neurons","1" = "General neurons","2" = "General neurons","3" = "Trh_Gabaergic", 
                                                    "4" = "General neurons","5" = "Grp_Glutamatergic","6" = "Prok2_Gabaergic","7" = "Agrp_Gabaergic",
                                                    "8" = "Pomc_Glutamatergic","9" = "Sst_Bcl11b_Gabaergic","10" = "Prlr_Glutamatergic","11" = "Sim1_Glutamatergic","12" = "Ghrh_Gabaergic",
                                                    "13" = "Crabp1_Gabaergic","14" = "Kl_Gabaergic","15" = "Avp_Glutamatergic","16" = "Lhx1_Gabaergic",
                                                    "17" = "Lhx6_Gabaergic","18" = "Prph_Histaminergic","19" = "Cx3cr1_Gabaergic","20" = "Ndnf_Glutamatergic"))



#Focusing on Set1 only, 
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set1")
#SI only Set2, Set1 is not making a lot sense with lowerr quality
#dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents =  c("Fruc","Control"))


seuratObject=dropEST.combined.filtered
DefaultAssay(seuratObject) <- "RNA"
#metaData1=NULL
#metaDataSubset=NULL
cellTypeMetaData="celltype"
nPerm=5000
do.par=TRUE
num.cores= 4
comparisonMetaData="data.diseasestatus"
group1="Fruc"
group2="Control"
#distanceScorePlot=NULL
#fcPlot=NULL
#distanceSave=NULL
seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)

#Number of permutations
numPerm = nPerm
if(do.par==TRUE){
  #auto detect number of available cores
  if(num.cores == "auto"){
    cores=detectCores()
    #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
    cores = max(1,cores[1]-1)
    #setup and start parallel processing cluster
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else{
    #if a number of cores for parallel processing is specified, use this number of cores
    if(is.numeric(num.cores)){
      cores = num.cores
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
  }
}
cellTypePermutationDistances = list()
cellTypeDistances = c()
allcelltypes = levels(seuratObject@active.ident)
availcelltypes <- NULL
for(i in 1:length(allcelltypes)){
  cellType = allcelltypes[i]
  print(cellType)
  #subset to cell type of interest
  cellTypeSubset = subset(seuratObject,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
  
  
  group1Length = length(which(cellTypeSubset@active.ident==group1))
  group2length = length(which(cellTypeSubset@active.ident==group2))
  
  #Only continue if we have at least 3 cells of the cell type in each group
  if(group1Length > 9 & group2length > 9){
    availcelltypes <- c(availcelltypes, cellType)
    #Label group1 and group2 cells
    group1Cells = WhichCells(cellTypeSubset,idents = group1)
    group2Cells = WhichCells(cellTypeSubset,idents = group2)
    
    #Get the average expression of group1 and group2 cells
    cellTypeSubset <- ScaleData(cellTypeSubset, vars.to.regress = NULL)
    
    AverageExpall <- rowMeans(as.matrix(cellTypeSubset@assays$RNA@counts))
    AverageExpall <- sort(AverageExpall, decreasing = T)
    AverageExpall <- AverageExpall[AverageExpall > quantile(AverageExpall,probs = 0.9)]
    #AverageExpall <- AverageExpall[-c(1:20)]
    selectedgenes <- names(AverageExpall)
    selectedgenes <- selectedgenes[1:1000]
    #compute in non-log space
    #avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes))
    
    avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes, use.scale = T))
    avgExp = avgExp$RNA
    #top and bottom filtration
    #calculate the distance between the group1 and group2 cells
    cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
    
    cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
    
    #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
    allCells = c(group1Cells,group2Cells)
    
    #set up a progress bar
    pb <- txtProgressBar(max = numPerm, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
      set.seed(j)
      permSubset = cellTypeSubset
      sampleCells = sample(allCells,length(allCells),replace = FALSE)
      #place the cells into either group1 or group2 (same size as original group1 and group2)
      group1Synth = sampleCells[1:length(group1Cells)]
      group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
      
      permStatus = permSubset@active.ident
      permStatus[group1Synth] = group1
      permStatus[group2Synth] = group2
      permSubset$permStatus = permStatus
      permSubset = SetIdent(permSubset, value = "permStatus")
      
      #compute in non-log space
      #permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes))
      
      permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes, use.scale = T))
      permAvgExp = permAvgExp$RNA
      
      permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
      permCellTypeDistance
    }
    cellTypePermutationDistances[[cellType]] = permutationDistances
    close(pb)
  }
}

stopCluster(cl)

names(cellTypeDistances) = availcelltypes

cellTypeDistancePvals = c()
for(cellType in levels(seuratObject@active.ident)){
  numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
  cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
}
names(cellTypeDistancePvals) = availcelltypes


melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
#availcelltypes = unique(melted.cellTypePermutationDistances$L1)

melted.cellTypePvals = rep(cellTypeDistancePvals[availcelltypes],each=numPerm)
melted.cellTypeDistances = rep(cellTypeDistances[availcelltypes],each=numPerm)

melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)

#Bonferroni correction
allDistancePvals = c(cellTypeDistancePvals[availcelltypes])
# all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm

#Combine neuronal subtype and cell types to plot
colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
melted.combinedPermutationDistances = melted.cellTypePermutationDistances
melted.combinedPermutationDistances$pvals = allDistancePvals
melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))


#########
#Plot
ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
  geom_density(adjust = 1) +
  facet_wrap( ~ L1, ncol=5,scales="free") +
  xlab("Euclidian Distance") +
  ylab("Density") +
  # Add line
  geom_vline(aes(xintercept=distances),
             color="blue", linetype="dashed", size=1) +
  annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)

pdf(file = paste0("./Plots/EuclideanDistanceFructose/Fructoseplot1.pdf"),height = 16,width = 16)
print(p1)
dev.off()

#########
#get the log FC of the euclidean distance versus the median of the background distribution

medianBackground = c()
for(cellType in levels(seuratObject@active.ident)){
  medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
}
names(medianBackground) = availcelltypes

logFC = log10(cellTypeDistances[availcelltypes]) - log10(medianBackground)
pval = -log10(allDistancePvals)
toPlot = as.data.frame(cbind(logFC,pval))

p2 <- ggplot(toPlot, aes(pval, logFC)) +
  geom_point(color = 'black') +
  theme_classic(base_size = 10) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 

pdf(file = paste0("./Plots/EuclideanDistanceFructose/Fructoseplot2.pdf"))
print(p2)
dev.off()

save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./Fructosedistanceresult_subsetzscore.rda"))
rm(list = ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Neuroncelltype.rda")
dropEST.combined.filtered <- CellTypeSubset
dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.1.5,
                                                  c("0" = "General neurons","1" = "General neurons","2" = "General neurons","3" = "Trh_Gabaergic", 
                                                    "4" = "General neurons","5" = "Grp_Glutamatergic","6" = "Prok2_Gabaergic","7" = "Agrp_Gabaergic",
                                                    "8" = "Pomc_Glutamatergic","9" = "Sst_Bcl11b_Gabaergic","10" = "Prlr_Glutamatergic","11" = "Sim1_Glutamatergic","12" = "Ghrh_Gabaergic",
                                                    "13" = "Crabp1_Gabaergic","14" = "Kl_Gabaergic","15" = "Avp_Glutamatergic","16" = "Lhx1_Gabaergic",
                                                    "17" = "Lhx6_Gabaergic","18" = "Prph_Histaminergic","19" = "Cx3cr1_Gabaergic","20" = "Ndnf_Glutamatergic"))

#Focusing on Set1 only, except SVF since almost no batch effect in SVF
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set1")
#dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents  = "Set2")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents =  c("HFHS","Control"))

seuratObject=dropEST.combined.filtered
DefaultAssay(seuratObject) <- "RNA"
#metaData1=NULL
#metaDataSubset=NULL
cellTypeMetaData="celltype"
nPerm=5000
do.par=TRUE
num.cores= 4
comparisonMetaData="data.diseasestatus"
group1="HFHS"
group2="Control"
#distanceScorePlot=NULL
#fcPlot=NULL
#distanceSave=NULL
seuratObject = SetIdent(seuratObject, value = cellTypeMetaData)

#Number of permutations
numPerm = nPerm
if(do.par==TRUE){
  #auto detect number of available cores
  if(num.cores == "auto"){
    cores=detectCores()
    #set the number of cores for processing to the number of available cores minus 1 or 1 (if only 1 available core)
    cores = max(1,cores[1]-1)
    #setup and start parallel processing cluster
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  } else{
    #if a number of cores for parallel processing is specified, use this number of cores
    if(is.numeric(num.cores)){
      cores = num.cores
      #setup and start parallel processing cluster
      cl <- makeCluster(cores)
      registerDoSNOW(cl)
    }
  }
}
cellTypePermutationDistances = list()
cellTypeDistances = c()
allcelltypes = levels(seuratObject@active.ident)
availcelltypes <- NULL
for(i in 1:length(allcelltypes)){
  cellType = allcelltypes[i]
  print(cellType)
  #subset to cell type of interest
  cellTypeSubset = subset(seuratObject,idents = cellType)
  cellTypeSubset = SetIdent(cellTypeSubset, value = comparisonMetaData)
  
  group1Length = length(which(cellTypeSubset@active.ident==group1))
  group2length = length(which(cellTypeSubset@active.ident==group2))
  
  #Only continue if we have at least 3 cells of the cell type in each group
  if(group1Length > 9 & group2length > 9){
    availcelltypes <- c(availcelltypes, cellType)
    #Label group1 and group2 cells
    cellTypeSubset <- ScaleData(cellTypeSubset, vars.to.regress = NULL)
    
    group1Cells = WhichCells(cellTypeSubset,idents = group1)
    group2Cells = WhichCells(cellTypeSubset,idents = group2)
    
    #Seletcing gene subsets for distance analysis
    AverageExpall <- rowMeans(as.matrix(cellTypeSubset@assays$RNA@counts))
    AverageExpall <- sort(AverageExpall, decreasing = T)
    AverageExpall <- AverageExpall[AverageExpall > quantile(AverageExpall,probs = 0.9)]
    #AverageExpall <- AverageExpall[-c(1:20)]
    selectedgenes <- names(AverageExpall)
    #Get the average expression of group1 and group2 cells
    #compute in non-log space
    #avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes))
    avgExp = suppressWarnings(AverageExpression(cellTypeSubset, assays = "RNA", features = selectedgenes,use.scale = T))
    
    avgExp = avgExp$RNA
    
    #calculate the distance between the group1 and group2 cells
    cellTypeDistance = sqrt(sum((avgExp[,1] - avgExp[,2]) ^ 2))
    
    cellTypeDistances = c(cellTypeDistances,cellTypeDistance)
    
    #Generate null distribution -- randomly pull cells from both groups, to make a group the size of group1 and a group the size of group2
    allCells = c(group1Cells,group2Cells)
    
    #set up a progress bar
    pb <- txtProgressBar(max = numPerm, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    permutationDistances <- foreach(j=1:numPerm, .combine=c, .options.snow=opts, .packages='Seurat') %dopar% {
      set.seed(j)
      permSubset = cellTypeSubset
      sampleCells = sample(allCells,length(allCells),replace = FALSE)
      #place the cells into either group1 or group2 (same size as original group1 and group2)
      group1Synth = sampleCells[1:length(group1Cells)]
      group2Synth = sampleCells[(length(group1Cells)+1):(length(group1Cells)+length(group2Cells))]
      
      permStatus = permSubset@active.ident
      permStatus[group1Synth] = group1
      permStatus[group2Synth] = group2
      permSubset$permStatus = permStatus
      permSubset = SetIdent(permSubset, value = "permStatus")
      
      #compute in non-log space
      #permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes))
      permAvgExp = suppressWarnings(AverageExpression(permSubset, assays = "RNA", features = selectedgenes, use.scale = T))
      
      permAvgExp = permAvgExp$RNA
      
      permCellTypeDistance = sqrt(sum((permAvgExp[,1] - permAvgExp[,2]) ^ 2))
      permCellTypeDistance
    }
    cellTypePermutationDistances[[cellType]] = permutationDistances
    close(pb)
  }
}

stopCluster(cl)

names(cellTypeDistances) = availcelltypes

cellTypeDistancePvals = c()
for(cellType in levels(seuratObject@active.ident)){
  numGreater = length(which(cellTypePermutationDistances[[cellType]]>cellTypeDistances[cellType]))
  cellTypeDistancePvals = c(cellTypeDistancePvals,numGreater/numPerm)
}
names(cellTypeDistancePvals) = availcelltypes


melted.cellTypePermutationDistances = melt(cellTypePermutationDistances)
#availcelltypes = unique(melted.cellTypePermutationDistances$L1)

melted.cellTypePvals = rep(cellTypeDistancePvals[availcelltypes],each=numPerm)
melted.cellTypeDistances = rep(cellTypeDistances[availcelltypes],each=numPerm)

melted.cellTypePermutationDistances = cbind(melted.cellTypePermutationDistances,melted.cellTypePvals,melted.cellTypeDistances)

#Bonferroni correction
allDistancePvals = c(cellTypeDistancePvals[availcelltypes])
# all_distance_pvals = p.adjust(all_distance_pvals,method = "bonferroni")
allDistancePvals[which(allDistancePvals == 0)] <- 1/numPerm

#Combine neuronal subtype and cell types to plot
colnames(melted.cellTypePermutationDistances) <- c("value","L1","pvals","distances")
melted.combinedPermutationDistances = melted.cellTypePermutationDistances
melted.combinedPermutationDistances$pvals = allDistancePvals
melted.combinedPermutationDistances$L1 <- factor(melted.combinedPermutationDistances$L1, levels=c(levels(seuratObject@active.ident)))


#########
#Plot
ifelse(!dir.exists(file.path("./Plots/EuclideanDistance/")), dir.create(file.path("./Plots/EuclideanDistance/"),recursive = T), FALSE)
p1 <- ggplot(melted.combinedPermutationDistances, aes(value)) +
  geom_density(adjust = 1) +
  facet_wrap( ~ L1, ncol=5,scales="free") +
  xlab("Euclidian Distance") +
  ylab("Density") +
  # Add line
  geom_vline(aes(xintercept=distances),
             color="blue", linetype="dashed", size=1) +
  annotate(geom = "text", -Inf , Inf, hjust = -0.75, vjust = 2 , label=paste0("P adj: ",sprintf("%0.3f", round(allDistancePvals, digits = 3))), color = "red", size = 4)

pdf(file = paste0("./Plots/EuclideanDistanceHFHS/HFHSplot1.pdf"),height = 16,width = 16)
print(p1)
dev.off()

#########
#get the log FC of the euclidean distance versus the median of the background distribution

medianBackground = c()
for(cellType in levels(seuratObject@active.ident)){
  medianBackground = c(medianBackground,median(cellTypePermutationDistances[[cellType]]))
}
names(medianBackground) = availcelltypes

logFC = log10(cellTypeDistances[availcelltypes]) - log10(medianBackground)
pval = -log10(allDistancePvals)
toPlot = as.data.frame(cbind(logFC,pval))

p2 <- ggplot(toPlot, aes(pval, logFC)) +
  geom_point(color = 'black') +
  theme_classic(base_size = 10) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = rownames(toPlot)), size = 3.5) 

pdf(file = paste0("./Plots/EuclideanDistanceHFHS/HFHSplot2.pdf"))
print(p2)
dev.off()
save(cellTypeDistances,cellTypePermutationDistances,numPerm,melted.combinedPermutationDistances,allDistancePvals,toPlot,file = paste0("./HFHSdistanceresult_subsetzscore.rda"))
