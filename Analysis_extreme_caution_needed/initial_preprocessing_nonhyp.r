#dropseq analysis
#CCA integration based on sample condition (if some sample is of low numbers -> integrating on sample level will do worse)
#Seurat label transfer function
library(Seurat)
mincells <- 3           #min number of cells that a gene has to be detected in to be kept
mingenes <- 300         #min number of genes to keep the cell, minimum genes changed from 400
#not used in current version
#minUMIs <- 700          #min number of UMIs to keep the cell
minPercentMT <- -Inf    #min mitochondrial percentage to keep a cell
maxgenes <- 3000        #max number of genes to keep a cell
#maxPercentMT <- 0.1     #max mitochondrial percentage to keep a cell
maxPercentMT <- 0.7 #for liver cells only
maxUMIs <- 5000         #max UMIs to keep a cell (< 2* max genes )
maxPercentRibo <- 0.1 #max ribosomal percentage to keep a cell
maxPercentPred <- 0.05  #max predicted genes percentage to keep a cell
#Pseudogenes filtration? LncRNA
#Ideal cell number, drop seq: 2-4% of total input cells after filtration
#10X, 50-60%

#NPC specific setting
maxPercentMT <- 0.15     #max mitochondrial percentage to keep a cell
minUMIs <- 700          #min number of UMIs to keep the cell
mingenes <- 200         #
maxPercentRibo <- 0.2
#Read in the data
tissue = "Hip"
#Get dropEST matrices directories
#"/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/Hep/F4_Hep_HFHS_6.1"
#combined batch liver
directories = c("/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/Hep/F1_Hep_Control_2.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/Hep/F2_Hep_Fruc_3.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/Hep/F5_Hep_Fruc_4.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Control/Hep/F1_Hep_Control_2.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/Hep/F2_Hep_Fruc_3.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/Hep/F3_Hep_Fruc_3.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/Hep/F4_Hep_HFHS_6.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/Hep/F5_Hep_HFHS_6.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/Hep/F6_Hep_HFHS_5.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/Hep/F4_Hep_HFHS_6.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/Hep/F2_Hep_HFHS_6.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Mixed/Mixed/F3_Hep_Control_2.2_SVF_HFHS_6.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/SVF/F1_SVF_HFHS_6.3")

#All NPC sets and SI combined
directories = c("/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/NPC/F4_NPC_Control_2.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/NPC/F2_NPC_Fruc_3.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/NPC/F3_NPC_Fruc_4.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Mixed/Mixed/F5_SI_HFHS_5.1_NPC_HFHS_6.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/SI/F1_SI_Fruc_4.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/SI/F5_SI_Control_1.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Control/NPC/F1_NPC_Control_1.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Control/NPC/F2_NPC_Control_2.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/NPC/F1_NPC_Fruc_3.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/NPC/F3_NPC_Fruc_3.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/NPC/F6_NPC_Fruc_4.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/NPC/F4_NPC_HFHS_6.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/NPC/F5_NPC_HFHS_5.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/NPC/F6_NPC_HFHS_6.3"
)

#SVF combined 
directories = c("/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/SVF/F2_SVF_Control_2.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/SVF/F5_SVF_Control_1.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/SVF/F4_SVF_Fruc_3.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/SVF/F4_SVF_Fruc_4.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Control/SVF/F1_SVF_Control_1.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Control/SVF/F2_SVF_Control_2.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/SVF/F3_SVF_Fruc_3.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/SVF/F4_SVF_Fruc_3.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/SVF/F5_SVF_HFHS_6.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/SVF/F1_SVF_HFHS_6.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Mixed/Mixed/F3_SVF_HFHS_6.1_Hep_Control_2.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/Hep/F1_Hep_Control_2.1")

directories = c("/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/SVF/F2_SVF_Control_2.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/SVF/F5_SVF_Control_1.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/SVF/F4_SVF_Fruc_3.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/SVF/F4_SVF_Fruc_4.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Control/SVF/F1_SVF_Control_1.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Control/SVF/F2_SVF_Control_2.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/SVF/F3_SVF_Fruc_3.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/SVF/F4_SVF_Fruc_3.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/SVF/F5_SVF_HFHS_6.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/SVF/F1_SVF_HFHS_6.3")

for(sample in directories){
  tissue <- unlist(strsplit(sample,"/"))[10]
  splitSample <- unlist(strsplit(sample,"/"))[11]
  batch <- unlist(strsplit(sample,"/"))[8]
  #get meta data
  
  splitSample = unlist(strsplit(splitSample,"_"))
  #tissue = splitSample[2]
  condition = splitSample[3]
  sampleName = paste0(batch,splitSample[4])
  print(sampleName)
  
  
  #Old sample loading
  #dropEST.data <- read.delim(sample,sep = "\t", header = T, row.names = 1)
  #dropEST.seurat <- CreateSeuratObject(counts = dropEST.data, project = sampleName)
  #
  #load the 10X data
  dropEST.data <- Read10X(data.dir = sample)
  dropEST.data@Dimnames[[2]] = paste0(sampleName,"_",dropEST.data@Dimnames[[2]])
  dropEST.seurat <- CreateSeuratObject(counts = dropEST.data, project = sampleName)
  rm(dropEST.data)
  
  #define meta data
  dropEST.seurat@meta.data$data.diseasestatus = condition
  dropEST.seurat@meta.data$data.tissue = tissue
  dropEST.seurat@meta.data$batch = batch
  dropEST.seurat@meta.data$set_treatment = paste(tissue,condition,batch,sep = "_")
  #merge seurat objects
  if(!exists("dropEST.combined")){
    dropEST.combined <- dropEST.seurat
    firstSampleName = sampleName
    firstSample = TRUE
  } else{
    if(firstSample==TRUE){
      dropEST.combined <- merge(x = dropEST.combined, y = dropEST.seurat, project = "7Day24HrHip")
      firstSample = FALSE
    } else{
      dropEST.combined <- merge(x = dropEST.combined, y = dropEST.seurat, project = "7Day24HrHip")
    }
  }
  #clean up environment
  rm(dropEST.seurat)
}
mito.features <- grep(pattern = "^mt-", x = rownames(x = dropEST.combined), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts")[mito.features,])/Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts"))

# get ribosome %
ribo.features <- grep(pattern = "^Rps", x = rownames(x = dropEST.combined), value = TRUE)
ribo.features <- c(ribo.features, grep(pattern = "^Rpl", x = rownames(x = dropEST.combined), value = TRUE))
percent.ribo <- Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts")[ribo.features,])/Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts"))

# get predicted genes %
pred.features <- grep(pattern = "^Gm1", x = rownames(x = dropEST.combined), value = TRUE)
pred.features <- c(pred.features,grep(pattern = "^Gm2", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm3", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm4", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm5", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm6", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm7", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm8", x = rownames(x = dropEST.combined), value = TRUE))
pred.features <- c(pred.features,grep(pattern = "^Gm9", x = rownames(x = dropEST.combined), value = TRUE))
percent.pred <- Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts")[pred.features,])/Matrix::colSums(x = GetAssayData(object = dropEST.combined, slot = "counts"))

# Add % mito, ribo and pred to the meta data
dropEST.combined[["percent.mito"]] <- percent.mito
dropEST.combined[["percent.ribo"]] <- percent.ribo
dropEST.combined[["percent.pred"]] <- percent.pred

#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1") #29905 -> 9900, 3 sample (one seems failed), 1485 in min gene > 300
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set2") #12000 -> 7000, 3 sample, 3915 in min gene > 300
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/NPC_set1") #13832 -> 3498, 3 sample, 958  in min gene > 300
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/NPC_set2") #11973 -> 337, 5 sample,  discontinued due to extremely high mitochondrial content..
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_set1") #12960 -> 9164, 4 sample (one seems failed), 4183 in min gene > 300
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_set2") #16634 -> 12335, 4 sample (one seems failed), 6919 in min gene > 300
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Liver_combined")
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined")
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_NPC")
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/old_NPC_set1")

pdf(file="NumGenesNumUMIPercentMito_PreFilter.pdf",height=10,width=20)
VlnPlot(object = dropEST.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo", "percent.pred"),pt.size = 0.15)
dev.off()

dropEST.combined = subset(x = dropEST.combined, subset = nFeature_RNA > mingenes & nFeature_RNA < maxgenes & percent.mito < maxPercentMT 
                          & percent.ribo < maxPercentRibo & percent.pred < maxPercentPred & nCount_RNA < maxUMIs)
#dropEST.combined = subset(x = dropEST.combined, subset = nFeature_RNA > mingenes & nFeature_RNA < maxgenes & percent.mito < maxPercentMT 
#                           & percent.pred < maxPercentPred & nCount_RNA < maxUMIs)

dropEST.combined.filtered = dropEST.combined

table(dropEST.combined.filtered@active.ident)

rm(dropEST.combined,animal,condition,directories,firstSample,firstSampleName,maxgenes,maxPercentMT,maxUMIs,mincells,mingenes,minPercentMT,
   minUMIs,maxPercentPred,maxPercentRibo,mito.features,percent.mito,pred.features,percent.pred,ribo.features,percent.ribo,
   sample,sampleName,splitSample,timePoint,tissue)
save(dropEST.combined.filtered, file = "filtered_raw.rda")

dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "orig.ident")
pdf(file="PostFilter.pdf",height=10,width=20)
VlnPlot(object = dropEST.combined.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo", "percent.pred"),pt.size = 0.15)
dev.off()

dropEST.combined.filtered <- NormalizeData(object = dropEST.combined.filtered, normalization.method = "LogNormalize",scale.factor = 10000)

dropEST.combined.filtered <- FindVariableFeatures(object = dropEST.combined.filtered)
length(dropEST.combined.filtered@assays$RNA@var.features) #2000


dropEST.combined.filtered <- ScaleData(object = dropEST.combined.filtered, features = rownames(x = dropEST.combined.filtered),
                                       vars.to.regress = c("nCount_RNA","percent.mito","percent.ribo"))

#Perform PCA on the scaled data (uses the highly var genes)
dropEST.combined.filtered <- RunPCA(object = dropEST.combined.filtered, verbose = T, npcs = 75, ndims.print = 1:5, nfeatures.print = 10)

#Perform ICA on the scaled data (uses the highly var genes)
dropEST.combined.filtered <- RunICA(object = dropEST.combined.filtered, verbose = T, nics = 75, ndims.print = 1:5, nfeatures.print = 10)

ifelse(!dir.exists(file.path("./PCA")), dir.create(file.path("./PCA")), FALSE)
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "condition")
pdf(file="./PCA/Hip_PCA_HighlyVarGenes.pdf")
DimPlot(object = dropEST.combined.filtered, dims = c(1,2), reduction = "pca")
dev.off()


#Plot Heatmaps of PCs
#pdf(file="./PCA/Hip_PCA_Heatmaps.pdf",height=20,width=20)
#DimHeatmap(object = dropEST.combined.filtered, dims = 1:12, balanced = TRUE, cells = 100, reduction = "pca")
#dev.off()


#Plot ICA Plots
ifelse(!dir.exists(file.path("./ICA")), dir.create(file.path("./ICA")), FALSE)
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "condition")
pdf(file="./ICA/Hip_ICA_HighlyVarGenes.pdf")
DimPlot(object = dropEST.combined.filtered, dims = c(1,2), reduction = "ica")
dev.off()

#Plot Heatmaps of ICs
#pdf(file="./ICA/Hip_ICA_Heatmaps.pdf",height=20,width=20)
#DimHeatmap(object = dropEST.combined.filtered, dims = 1:12, balanced = TRUE, cells = 100, reduction = "ica")
#dev.off()

#dropEST.combined.filtered <- JackStraw(object = dropEST.combined.filtered, reduction = "pca", num.replicate = 50, 
#                                       verbose = TRUE, dims = 50)
#dropEST.combined.filtered = ScoreJackStraw(object = dropEST.combined.filtered, dims = 1:50, reduction = "pca")

#pdf(file="ackStraw.pdf",height=45,width=10) #50 Sig PCs
#JackStrawPlot(object = dropEST.combined.filtered, dims = 1:50)
#dev.off()
#ElbowPlot(dropEST.combined.filtered) #10PCs

#CCA integration of each condition
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered,value="set_treatment")
seurat.list = list()
for(i in 1:length(levels(dropEST.combined.filtered@active.ident))){
  seurat.list[[i]] = subset(dropEST.combined.filtered,idents=levels(dropEST.combined.filtered@active.ident)[i])
}
for (i in 1:length(x = seurat.list)) {
  seurat.list[[i]] <- NormalizeData(object = seurat.list[[i]], verbose = FALSE)
  seurat.list[[i]] <- FindVariableFeatures(object = seurat.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:30)
seuratIntegrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:30)
library(ggplot2)
library(cowplot)
DefaultAssay(object = seuratIntegrated) <- "integrated"
seuratIntegrated <- ScaleData(object = seuratIntegrated, verbose = FALSE)
seuratIntegrated <- RunPCA(object = seuratIntegrated, npcs = 30, verbose = FALSE)
seuratIntegrated <- RunUMAP(object = seuratIntegrated, reduction = "pca", dims = 1:30)
seuratIntegrated <- RunTSNE(object = seuratIntegrated, reduction = "pca", dims = 1:30)
seuratIntegrated <- FindNeighbors(object = seuratIntegrated, reduction = "pca", dims = 1:30, k.param = 25)
seuratIntegrated <- FindClusters(object = seuratIntegrated, resolution = seq(0.5,4,by=0.5), verbose = T, reduction = "pca")
ifelse(!dir.exists(file.path(paste0("./Plots/ccaTSNEcondition"))), dir.create(file.path(paste0("./Plots/ccaTSNEcondition")),recursive = T), FALSE)
ifelse(!dir.exists(file.path(paste0("./Plots/ccaUMAPcondition"))), dir.create(file.path(paste0("./Plots/ccaUMAPcondition")),recursive = T), FALSE)

tissue <- "SVF"
pdf(file=paste0("./Plots/ccaUMAPcondition/",tissue,".pdf"))
print(DimPlot(object = seuratIntegrated, reduction = "umap", group.by = "set_treatment"))
dev.off()
pdf(file=paste0("./Plots/ccaTSNEcondition/",tissue,".pdf"))
print(DimPlot(object = seuratIntegrated, reduction = "tsne", group.by = "set_treatment"))
dev.off()

#For hepatocytes
dropEST.combined.filtered <- seuratIntegrated
#
#Demultiplexing of mixed cells in NPC and SVF
#####
seuratIntegrated <- SetIdent(seuratIntegrated,value="set_treatment")
selectedcells <- FeatureLocator(DimPlot(object = seuratIntegrated, do.identify = TRUE,  reduction = "umap",pt.size = 0.5))
selectedcells <- FeatureLocator(DimPlot(object = seuratIntegrated, do.identify = TRUE,  reduction = "tsne",pt.size = 0.5))

CellTypeSubset = subset(seuratIntegrated, cells = selectedcells,idents = c("Mixed_HFHS_Set1","NPC_Control_Set1","NPC_Control_Set2","NPC_Fruc_Set1",
                                                                           "NPC_Fruc_Set2","NPC_HFHS_Set2")) #For NPC UMAP 
CellTypeSubset = subset(seuratIntegrated, cells = selectedcells,idents = c("Mixed_HFHS_Set1","SVF_Control_Set1","SVF_Control_Set2","SVF_Fruc_Set1",
                                                                           "SVF_Fruc_Set2","SVF_HFHS_Set1","SVF_HFHS_Set2")) #For SVF UMAP 

CellTypeSubset = subset(seuratIntegrated, cells = selectedcells,idents = c("Hep_Control_Set1","Mixed_Control_Set1","Hep_Fruc_Set1","Hep_HFHS_Set1")) #For SVF UMAP 

dropEST.combined.filtered <- CellTypeSubset
rm(CellTypeSubset)

dropEST.combined.filtered  <- ScaleData(object = dropEST.combined.filtered , verbose = FALSE)
dropEST.combined.filtered  <- RunPCA(object = dropEST.combined.filtered , npcs = 30, verbose = FALSE)
dropEST.combined.filtered  <- FindNeighbors(object = dropEST.combined.filtered , reduction = "pca", dims = 1:30, k.param = 25)
dropEST.combined.filtered  <- FindClusters(object = dropEST.combined.filtered , resolution = seq(0.5,4,by=0.5), verbose = T, reduction = "pca")
dropEST.combined.filtered  <- RunUMAP(object = dropEST.combined.filtered , reduction = "pca", dims = 1:30)
dropEST.combined.filtered  <- RunTSNE(object = dropEST.combined.filtered , reduction = "pca", dims = 1:30)
######


resolutions <- c("0.5","1","1.5","2","2.5","3","3.5","4")

dir.create(file.path("./SampleTSNE"))
dir.create(file.path("./SampleUMAP"))
#dir.create(file.path("./SampleTSNE_allcombine"))
#dir.create(file.path("./SampleUMAP_allcombine"))
for(reso in resolutions){
  pdf(file=paste0("./SampleTSNE/",reso,".pdf"))
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = paste0("integrated_snn_res.",reso))
  print(DimPlot(dropEST.combined.filtered, label = T, reduction = "tsne"))
  dev.off()
  pdf(file=paste0("./SampleUMAP/",reso,".pdf"))
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = paste0("integrated_snn_res.",reso))
  print(DimPlot(dropEST.combined.filtered, label = T, reduction = "umap"))
  dev.off()
}

c("CD4Tcells","Cd14Moncytes","Bcells","Cd8Tcells","FCGR3Amonocytes","NkCells","DendriticCells1","DendriticCells2","Megakaryocytes")

pdf("Marker_set1.pdf",width = 10, height = 10)
FeaturePlot(object =dropEST.combined.filtered, features = c("Il7r","Cd14","Ms4a1","Cd8a","Ms4a7","Nkg7","Fcer1a","Cst3","Ppbp"), reduction = "umap", min.cutoff = 0,pt.size =0.0000000001)
dev.off()
pdf("Marker_set2.pdf",width = 10, height = 10)
FeaturePlot(object =dropEST.combined.filtered, features = c("Cd3e","Fcgr4","Lyz2","Lyz1"), reduction = "umap", min.cutoff = 0,pt.size =0.0000000001)
dev.off()
pdf("Marker_set3.pdf",width = 10, height = 10)
FeaturePlot(object =dropEST.combined.filtered, features = c("Postn","Siglech","Fabp1","Fbp1","Alb","Clec4f","Mmrn2","Lyve1","Tek"), reduction = "umap", min.cutoff = 0,pt.size =0.0000000001)
dev.off()

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
pdf(file="batch_SampleUMAP.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.5, reduction = "umap")
dev.off()
pdf(file="batch_SampleTSNE.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.5, reduction = "tsne")
dev.off()

#dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
#pdf(file="batch_noCCA_SampleUMAP.pdf")
#DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.5, reduction = "umap")
#dev.off()

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
pdf(file="TreatmentUMAP.pdf")
DimPlot(dropEST.combined.filtered,do.label = F, reduction = "umap")
dev.off()
pdf(file="TreatmentTSNE.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.1, reduction = "tsne")
dev.off()

save(dropEST.combined.filtered, file = "HEP_Set1_DropEST.Object1.rda")
pdf("liver_zonation.pdf", width = 10,height = 6)
FeaturePlot(object =dropEST.combined.filtered, features = c("Cyp2f2","Cyp2e1"), reduction = "umap", min.cutoff = 0)
dev.off()
FeaturePlot(object =dropEST.combined.filtered, features = c("Glul","Ass1"), reduction = "tsne", min.cutoff = 0)

#removing HFHS samples
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
CellTypeSubset = subset(dropEST.combined.filtered, idents = c("Fruc","Control"))
CellTypeSubset <- SetIdent(CellTypeSubset, value = "set_treatment")
DimPlot(CellTypeSubset,do.label = F, reduction = "umap")

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "orig.ident")
DimPlot(dropEST.combined.filtered,do.label = F, reduction = "umap")

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "set_treatment")


library(scMCA)
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")
test <- exp(AverageExpression(dropEST.combined.filtered)$integrated)
mca_result <- scMCA(scdata = test, numbers_plot = 3)
scMCA_vis(mca_result)
########################
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1") #29905 -> 9900, 3 sample (one seems failed), 1485 in min gene > 300
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set2") #12000 -> 7000, 3 sample, 3915 in min gene > 300
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/NPC_set1") #13832 -> 3498, 3 sample, 958  in min gene > 300

#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/NPC_set2") #11973 -> 337, 5 sample,  discontinued due to extremely high mitochondrial content..
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_set1") #12960 -> 9164, 4 sample (one seems failed), 4183 in min gene > 300
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_set2") #16634 -> 12335, 4 sample (one seems failed), 6919 in min gene > 300
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Liver_combined")
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined")
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_NPC")
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/old_NPC_set1")

library(Seurat)
#re-install reticulate fu**er after mac update
#ibrary(reticulate)
#py_install("umap-learn")

load("HEP_Set1_DropEST.Object1.rda")
#dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "RNA_snn_res.1.5")
#CellTypeSubset = subset(dropEST.combined.filtered, idents = c("0","9")) #For SI_NPC combined set
#CellTypeSubset = subset(dropEST.combined.filtered, idents = as.character(c(0:7,9:16))) #For SVF UMAP 
#dropEST.combined.filtered <- CellTypeSubset
#CellTypeSubset = subset(dropEST.combined.filtered, cells = selectedcells) #For SVF UMAP 
#dropEST.combined.filtered <- CellTypeSubset


dropEST.combined.filtered <- FindNeighbors(object = dropEST.combined.filtered, reduction = "ica", dims = 1:50, k.param = 25)
dropEST.combined.filtered <- FindClusters(object = dropEST.combined.filtered, resolution = c(0.5,1,1.5,2,2.5,3,3.5,4), verbose = T, reduction = "ica")

dropEST.combined.filtered <- RunUMAP(object = dropEST.combined.filtered, reduction = "ica", dims = 1:10)
dropEST.combined.filtered <- RunTSNE(object = dropEST.combined.filtered, reduction = "ica", dims = 1:10)

#Pericentral cells
FeaturePlot(object =dropEST.combined.filtered, features = c("Glul"), reduction = "umap", min.cutoff = 0, pt.size = 1)
FeaturePlot(object =dropEST.combined.filtered, features = c("Cyp2e1"), reduction = "umap", min.cutoff = 0, pt.size = 1)
FeaturePlot(object =dropEST.combined.filtered, features = c("Ass1"), reduction = "umap", min.cutoff = 0, pt.size = 1)
FeaturePlot(object =dropEST.combined.filtered, features = c("Asl"), reduction = "umap", min.cutoff = 0, pt.size = 1)
FeaturePlot(object =dropEST.combined.filtered, features = c("Cyp2f2"), reduction = "umap", min.cutoff = 0, pt.size = 1)
FeaturePlot(object =dropEST.combined.filtered, features = c("nCount_RNA"), reduction = "umap", min.cutoff = 0, pt.size = 1)

#Pericentral cells
FeaturePlot(object =dropEST.combined.filtered, features = c("Axin2","Cyp1a2","Gstm3","Psmd4"), reduction = "tsne", min.cutoff = 0, pt.size = 1)
#Non monotonic
FeaturePlot(object =dropEST.combined.filtered, features = c("Hamp","Igfbp2","Cyp8b1","Mup3"), reduction = "tsne", min.cutoff = 0, pt.size =1)
#Periportal
FeaturePlot(object =dropEST.combined.filtered, features = c("Arg1","Pck1","C2","Sdhd"), reduction = "tsne", min.cutoff = 0, pt.size = 1)
#NPC cells
pdf("NPCmarkers.pdf",width = 10, height = 12)
FeaturePlot(object =dropEST.combined.filtered, features = c("Ccl5","Spp1","Dcn","Csf1r","Kdr","Cd79b"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
dev.off()
FeaturePlot(object =dropEST.combined.filtered, features = c("Apoc3","Mup3","Tek"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
FeaturePlot(object =dropEST.combined.filtered, features = c("Ccr2","Cd5l","C1qb"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
FeaturePlot(object =dropEST.combined.filtered, features = c("Col1a1","Cxcl12","Pdgfra"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
FeaturePlot(object =dropEST.combined.filtered, features = c("Bmp2","Kdr","Tek"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)

FeaturePlot(object =dropEST.combined.filtered, features = c("Hbb"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
#Kupffer cell
pdf("Kupffer cell.pdf")
FeaturePlot(object =dropEST.combined.filtered, features = c("Irf7","Spic","Clec4f"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
dev.off()
#sinusoidal cell
pdf("sinosoidal endothelial cells.pdf", width = 16, height = 16)
FeaturePlot(object =dropEST.combined.filtered, features = c("Ushbp1", "Myf6", "Oit3", "Il1a", "F8", "Bmp2", "C1qtnf1", "Mmrn2", "Pcdh12", "Dpp4"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
dev.off()
#hepatic satellite cells (HSC)
FeaturePlot(object =dropEST.combined.filtered, features = c("Reln"), reduction = "umap", min.cutoff = 0, pt.size = 1)
#Cholangiocytes
pdf("Cholangiocytes.pdf")
FeaturePlot(object =dropEST.combined.filtered, features = c("Krt19",	"Krt7",	"Spp1"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
dev.off()
#committed Adipose progenitor cells
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")
VlnPlot(dropEST.combined.filtered,features = c("Cd34"))

#CD34+, CD29+, CD13+, CD73+, CD31-, CD45-
pdf("ADSC_marker.pdf",width = 10, height = 12)
FeaturePlot(object =dropEST.combined.filtered, features = c("Cd34","Itgb1","Anpep","Nt5e","Pecam1","Ptprc"), reduction = "tsne", min.cutoff = 0, pt.size = 0.001)
dev.off()
#EPC endothelial precursor cells
#CD34+, CD31+,CD133+, CD146+,CD45-
pdf("EPC_marker.pdf",width = 10, height = 10)
FeaturePlot(object =dropEST.combined.filtered, features = c("Cd34","Pecam1","Prom1","Mcam","Ptprc"), reduction = "tsne", min.cutoff = 0, pt.size = 0.001)
dev.off()
#Endothelial Cells CD31+,F8+, CD45-, none
FeaturePlot(object =dropEST.combined.filtered, features = c("Pecam1","F8","Ptprc"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
#T cell reg
#CD4+, CD25+, Foxp3+, CD8+
FeaturePlot(object =dropEST.combined.filtered, features = c("Cd4","Il2ra","Foxp3","Cd8a"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
#Macrophage
#CD45, CD14, CD34, CD206
FeaturePlot(object =dropEST.combined.filtered, features = c("Cd34","Ptprc","Cd14","Mrc1"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
#smooth muscle
FeaturePlot(object =dropEST.combined.filtered, features = c("Acta2"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
#Pre-adipocyte
#CD34+,  CD31-, CD45- 
FeaturePlot(object =dropEST.combined.filtered, features = c("Cd34","Pecam1","Ptprc"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
#Pericytes
#CD29+, CD34-, CD45-
FeaturePlot(object =dropEST.combined.filtered, features = c("Itgb1","Cd34","Ptprc"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)

#brown tissue progenitor
FeaturePlot(object =dropEST.combined.filtered, features = c("Myf5","Pax5","Myod1"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)

resolutions <- c("0.5","1","1.5","2","2.5","3","3.5","4")

dir.create(file.path("./SampleTSNE"))
dir.create(file.path("./SampleUMAP"))
#dir.create(file.path("./SampleTSNE_allcombine"))
#dir.create(file.path("./SampleUMAP_allcombine"))
for(reso in resolutions){
  pdf(file=paste0("./SampleTSNE/",reso,".pdf"))
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = paste0("RNA_snn_res.",reso))
  print(DimPlot(dropEST.combined.filtered, label = T, reduction = "tsne"))
  dev.off()
  pdf(file=paste0("./SampleUMAP/",reso,".pdf"))
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = paste0("RNA_snn_res.",reso))
  print(DimPlot(dropEST.combined.filtered, label = T, reduction = "umap"))
  dev.off()
}


#For SI and intestine combined only, select cells by using interactive function
#dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.tissue")
#pdf(file="TissueUMAP.pdf")
#DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.5, reduction = "umap")
#dev.off()
#pdf(file="TissueTSNE.pdf")
#DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.5, reduction = "tsne")
#dev.off()
#selectedcells <- FeatureLocator(DimPlot(object = dropEST.combined.filtered, do.identify = TRUE,  reduction = "tsne",pt.size = 0.5))
#selectedcells <- FeatureLocator(DimPlot(object = dropEST.combined.filtered, do.identify = TRUE,  reduction = "umap",pt.size = 0.5))

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
pdf(file="TreatmentUMAP.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.5, reduction = "umap")
dev.off()
pdf(file="TreatmentTSNE.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.5, reduction = "tsne")
dev.off()


dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "orig.ident")
pdf(file="SampleUMAP.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.5, reduction = "umap")
dev.off()
pdf(file="SampleTSNE.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.5, reduction = "tsne")
dev.off()


dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
CellTypeSubset = subset(dropEST.combined.filtered, idents = c("Fruc","Control"))
CellTypeSubset = SetIdent(CellTypeSubset, value = "batch")
DimPlot(CellTypeSubset,do.label = F,pt.size = 0.5, reduction = "umap")

pdf(file="batch_SampleUMAP.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.05, reduction = "umap")
dev.off()
pdf(file="batch_SampleTSNE.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.5, reduction = "tsne")
dev.off()

#Marker matching
library(scMCA)
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "RNA_snn_res.0.5")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")

test <- AverageExpression(dropEST.combined.filtered)$RNA
test <- exp(AverageExpression(dropEST.combined.filtered)$integrated)
mca_result <- scMCA(scdata = test, numbers_plot = 3)
scMCA_vis(mca_result)

#Marker Feature plot for mouse immune cells

FeaturePlot(object =dropEST.combined.filtered, features = c("Ms4a1","Cd3e","Cd14","Fcer1a","Fcgr4","Lyz2","Lyz1","Ppbp","Cd8a"), reduction = "umap", min.cutoff = 0)

pdf("Marker_set1.pdf",width = 10, height = 10)
FeaturePlot(object =dropEST.combined.filtered, features = c("Il7r","Cd14","Ms4a1","Cd8a","Ms4a7","Nkg7","Fcer1a","Cst3","Ppbp","Fcgr3"), reduction = "umap", min.cutoff = 0)
c("CD4Tcells","Cd14Moncytes","Bcells","Cd8Tcells","FCGR3Amonocytes","NkCells","DendriticCells1","DendriticCells2","Megakaryocytes","Dendritic cell 3")
dev.off()
pdf("Marker_set2.pdf",width = 10, height = 10)
FeaturePlot(object =dropEST.combined.filtered, features = c("Cd3e","Fcgr4","Lyz2","Lyz1"), reduction = "umap", min.cutoff = 0)
dev.off()
#Dendritic cell
#FeaturePlot(object =dropEST.combined.filtered, features = c("Itgax","Cd14","Fcgr3"), reduction = "tsne")

#FeaturePlot(object =dropEST.combined.filtered, features = c("Acvrl1","Cd8a"), reduction = "tsne")
#FeaturePlot(object =dropEST.combined.filtered, features = c("Basp1", "Ccl5", "Fscn1", "S100a4", "Tcf7"), reduction = "tsne")
FeaturePlot(object =dropEST.combined.filtered, features = c("Itgax","Cd14","Fcgr3"), reduction = "tsne")
FeaturePlot(object =dropEST.combined.filtered, features = c("Bst2", "Ly6d"), reduction = "tsne")
FeaturePlot(object =dropEST.combined.filtered, features = c("Hbb-bs"), reduction = "tsne", min.cutoff = 0)

#Monocyte Ly6c+
FeaturePlot(object =dropEST.combined.filtered, features = c("Ccr2", "Lyz2", "S100a4"), reduction = "tsne")

c("Tcell","Monocyte","CD14monocyte","CD14Monocyte")

FeaturePlot(object =dropEST.combined.filtered, features = c("Nkg7","Gzma","Cd3d","Il7r","Cst3"))
FeaturePlot(object =dropEST.combined.filtered, features = c("S100a8","Cd79a"))
FeaturePlot(object =dropEST.combined.filtered, features = c("Pi16"))

FeaturePlot(object =dropEST.combined.filtered, features = c("Apoc3","Mup3"), min.cutoff = 0)
FeaturePlot(object =dropEST.combined.filtered, features = c("Col1a2"), min.cutoff = 0)
FeaturePlot(object =dropEST.combined.filtered, features = c("Acta2"), min.cutoff = 0)

FeaturePlot(object =dropEST.combined.filtered, features = c("Zfp423","Mmp2","Eng","Hsp47","Fsp1","Cdh5"), min.cutoff = 0)
FeaturePlot(object =dropEST.combined.filtered, features = c("Itgam","mac1", "Ly-71"), min.cutoff = 0)
#dendritic cell marker
FeaturePlot(object =dropEST.combined.filtered, features = c("Itgax","Cd14","Fcgr3"), min.cutoff = 0)


#SVF combined cell indetification
#SVF, use cluster c(1,4,5,6, c(0,2,3,10:15))
#1. Macrophage
#c(0,2,3,7,10:15): Stromal/Mesechymal cells/Muscle cells
#4 Macrophage
#5 Endotheelial or Mural cells
#6 Epithelial cells, mesenchymal alveolar niche cells
#8,9 Ductal cells
#Hep set 1, hepatpcyte confirmed in both group
#Hep set 1, hepatpcyte confirmed in both group (with slightly higher score)
#NPC set 1, endotheial cell and macrophages (Erythroblast?)
###NPC cells are not well segregated, 0 is more to endothelial cell and 1 is more to macrophages
#Findmarkers
#Original identity comparison (for homogeneous liver cells)
#For other cells, subcluster analysis
#SVF, use cluster c(1,4,5,6, c(0,2,3,10:15))
#for NPC set1, set 0,1

###DEG analysis
library(Seurat)
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1") #29905 -> 9900, 3 sample (one seems failed), 1485 in min gene > 300
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set2") #12000 -> 7000, 3 sample, 3915 in min gene > 300
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/NPC_set1") #13832 -> 3498, 3 sample, 958  in min gene > 300

#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/NPC_set2") #11973 -> 337, 5 sample,  discontinued due to extremely high mitochondrial content..
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_set1") #12960 -> 9164, 4 sample (one seems failed), 4183 in min gene > 300
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_set2") #16634 -> 12335, 4 sample (one seems failed), 6919 in min gene > 300
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Liver_combined")
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined")
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_NPC")
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/old_NPC_set1")

load("HEP_Set1_DropEST.Object1.rda")

#re-install reticulate fu**er after mac update
#ibrary(reticulate)
#py_install("umap-learn")


#result for integrated space ±2000
#######
Degs = list()
topDEGs <- list()
DEGnum <- NULL
pathwaynum <- NULL
#dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "RNA_snn_res.1")
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")

table(dropEST.combined.filtered@active.ident) 
alllevels <- levels(dropEST.combined.filtered@active.ident)
#Full name space

for(i in 1:length(alllevels)){
  CellType <- alllevels[i]
  CellTypeSubset = subset(dropEST.combined.filtered, idents = CellType)
  
  CellTypeSubset = SetIdent(CellTypeSubset, value = "data.diseasestatus")
  DefaultAssay(CellTypeSubset) <- "RNA"
  siggenes <- NULL
  try({Degs[[CellType]] = FindMarkers(CellTypeSubset, ident.1 = "Fruc", ident.2 = "Control")})
  #Degs[[CellType]] = FindMarkers(CellTypeSubset, ident.1 = "water", ident.2 = "BPA10")
  siggenes <- rownames(Degs[[CellType]][Degs[[CellType]]$p_val_adj < 0.05,])
  if(length(grep("mt|Rps|Rpl",siggenes)) > 0){
    siggenes <- siggenes[-grep("mt|Rps|Rpl",siggenes)]
  }
  
  if(length(siggenes) == 0){
    pathwaynum[i] <- NA
    DEGnum[i] <- NA
    next
  }
  DEGnum[i] <- length(siggenes)
  topDEGs[[CellType]] <- siggenes
  write.csv(Degs[[CellType]][Degs[[CellType]]$p_val_adj < 0.05,],file = paste0(getwd(),"/",CellType,".csv"))
  
  #write.csv(siggenes,file = paste0(CellType,".csv"))
  #write.csv(siggenes,file = paste0(CellType,"_BPA100.csv"))
  
  rm("drugframe")
  res <- enrichr(siggenes,databases = c("KEGG_2016","GO_Biological_Process_2017"))
  for(j in 1:length(res)){
    currentframe <- res[[j]]
    currentframe <- currentframe[currentframe$Adjusted.P.value < 0.05,]
    if(nrow(currentframe)  == 0){next
    }else{
      if(!exists("drugframe")){
        drugframe <- currentframe
      }else{
        drugframe <- rbind.data.frame(drugframe,currentframe)
      }
    }
  }
  if(exists("drugframe")){
    write.csv(drugframe,file = paste0(getwd(),"/",CellType,"DEGpathways.csv"))
    pathwaynum[i] <- nrow(drugframe)
  }else{
    pathwaynum[i] <- NA
  }
}

finalframe <- data.frame(subpopulation = names(table(dropEST.combined.filtered@active.ident)),
                         cellnumber = as.numeric(table(dropEST.combined.filtered@active.ident) ),
                         DEGnum = DEGnum,
                         pathwaynum = pathwaynum)
write.csv(finalframe,"summary_DEG.csv")
#########

#combined results for full feature space
Degs = list()
topDEGs <- list()
DEGnum <- NULL
pathwaynum <- NULL
#dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "RNA_snn_res.1")
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")

table(dropEST.combined.filtered@active.ident) 
alllevels <- levels(dropEST.combined.filtered@active.ident)
#Full name space
DefaultAssay(dropEST.combined.filtered) <- "RNA"

for(i in 1:length(alllevels)){
  CellType <- alllevels[i]
  CellTypeSubset = subset(dropEST.combined.filtered, idents = CellType)
  
  CellTypeSubset = SetIdent(CellTypeSubset, value = "data.diseasestatus")
  DefaultAssay(CellTypeSubset) <- "RNA"
  siggenes <- NULL
  try({Degs[[CellType]] = FindMarkers(CellTypeSubset, ident.1 = "Fruc", ident.2 = "Control")})
  #Degs[[CellType]] = FindMarkers(CellTypeSubset, ident.1 = "water", ident.2 = "BPA10")
  siggenes <- rownames(Degs[[CellType]][Degs[[CellType]]$p_val_adj < 0.05,])
  if(length(grep("mt|Rps|Rpl",siggenes)) > 0){
    siggenes <- siggenes[-grep("mt|Rps|Rpl",siggenes)]
  }
  
  if(length(siggenes) == 0){
    pathwaynum[i] <- NA
    DEGnum[i] <- NA
    next
  }
  DEGnum[i] <- length(siggenes)
  topDEGs[[CellType]] <- siggenes
  write.csv(Degs[[CellType]][Degs[[CellType]]$p_val_adj < 0.05,],file = paste0(getwd(),"/",CellType,"full.csv"))
  
  #write.csv(siggenes,file = paste0(CellType,".csv"))
  #write.csv(siggenes,file = paste0(CellType,"_BPA100.csv"))
  
  rm("drugframe")
  res <- enrichr(siggenes,databases = c("KEGG_2016","GO_Biological_Process_2017"))
  for(j in 1:length(res)){
    currentframe <- res[[j]]
    currentframe <- currentframe[currentframe$Adjusted.P.value < 0.05,]
    if(nrow(currentframe)  == 0){next
    }else{
      if(!exists("drugframe")){
        drugframe <- currentframe
      }else{
        drugframe <- rbind.data.frame(drugframe,currentframe)
      }
    }
  }
  if(exists("drugframe")){
    write.csv(drugframe,file = paste0(getwd(),"/",CellType,"DEGpathways_full.csv"))
    pathwaynum[i] <- nrow(drugframe)
  }else{
    pathwaynum[i] <- NA
  }
}

finalframe <- data.frame(subpopulation = names(table(dropEST.combined.filtered@active.ident)),
                         cellnumber = as.numeric(table(dropEST.combined.filtered@active.ident) ),
                         DEGnum = DEGnum,
                         pathwaynum = pathwaynum)
write.csv(finalframe,"summary_DEG_full.csv")

pdf("UMIfeature.pdf", width = 6, height = 6)
FeaturePlot(dropEST.combined.filtered,"nCount_RNA",reduction = "umap")
dev.off()
