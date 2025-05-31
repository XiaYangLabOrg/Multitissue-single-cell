#NPC+liver
library(Seurat)
mincells <- 3           #min number of cells that a gene has to be detected in to be kept
mingenes <- 300         #min number of genes to keep the cell, minimum genes changed from 400
#not used in current version
#minUMIs <- 700          #min number of UMIs to keep the cell
minPercentMT <- -Inf    #min mitochondrial percentage to keep a cell
maxgenes <- 3000        #max number of genes to keep a cell
maxPercentMT <- 0.15     #max mitochondrial percentage to keep a cell
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

directories = c("/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/NPC/F4_NPC_Control_2.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/NPC/F2_NPC_Fruc_3.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/NPC/F3_NPC_Fruc_4.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Mixed/Mixed/F5_SI_HFHS_5.1_NPC_HFHS_6.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/SI/F1_SI_Fruc_4.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/SI/F5_SI_Control_1.1")


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

setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/NPC_set1")

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


ifelse(!dir.exists(file.path("./ICA")), dir.create(file.path("./ICA")), FALSE)
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "condition")
pdf(file="./ICA/Hip_ICA_HighlyVarGenes.pdf")
DimPlot(object = dropEST.combined.filtered, dims = c(1,2), reduction = "ica")
dev.off()
dropEST.combined.filtered  <- FindNeighbors(object = dropEST.combined.filtered , reduction = "pca", dims = 1:30, k.param = 25)
dropEST.combined.filtered  <- FindClusters(object = dropEST.combined.filtered , resolution = seq(0.5,4,by=0.5), verbose = T, reduction = "pca")
dropEST.combined.filtered  <- RunUMAP(object = dropEST.combined.filtered , reduction = "pca", dims = 1:30)
dropEST.combined.filtered  <- RunTSNE(object = dropEST.combined.filtered , reduction = "pca", dims = 1:30)

seuratIntegrated <- dropEST.combined.filtered
print(DimPlot(object = seuratIntegrated, reduction = "umap", group.by = "set_treatment"))
seuratIntegrated <- SetIdent(seuratIntegrated,value="set_treatment")
selectedcells <- CellSelector(DimPlot(object = seuratIntegrated, do.identify = TRUE,  reduction = "umap",pt.size = 0.5))

CellTypeSubset = subset(seuratIntegrated, cells = selectedcells,idents = c("Mixed_HFHS_Set1","NPC_Control_Set1","NPC_Fruc_Set1"))
CellTypeSubset <- SetIdent(CellTypeSubset, value = "orig.ident")

#Add liver cellsmincells <- 3           #min number of cells that a gene has to be detected in to be kept
mingenes <- 300         #min number of genes to keep the cell, minimum genes changed from 400
#not used in current version
#minUMIs <- 700          #min number of UMIs to keep the cell
minPercentMT <- -Inf    #min mitochondrial percentage to keep a cell
maxgenes <- 3000        #max number of genes to keep a cell
maxPercentMT <- 0.1     #max mitochondrial percentage to keep a cell
maxUMIs <- 5000         #max UMIs to keep a cell (< 2* max genes )
maxPercentRibo <- 0.1 #max ribosomal percentage to keep a cell
maxPercentPred <- 0.05  #max predicted genes percentage to keep a cell



directories = c("/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/Hep/F1_Hep_Control_2.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/Hep/F2_Hep_Fruc_3.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/Hep/F5_Hep_Fruc_4.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/Hep/F4_Hep_HFHS_6.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/Hep/F2_Hep_HFHS_6.1")
#directories = c("/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/NPC/F4_NPC_Control_2.1",
#                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/NPC/F2_NPC_Fruc_3.1",
#                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/NPC/F3_NPC_Fruc_4.3",
#                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/Hep/F1_Hep_Control_2.1",
#              "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/Hep/F2_Hep_Fruc_3.1",
#               "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/Hep/F5_Hep_Fruc_4.3",
#                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/Hep/F4_Hep_HFHS_6.3",
#               "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/Hep/F2_Hep_HFHS_6.1")

tissue = "Hip"
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

dropEST.combined = subset(x = dropEST.combined, subset = nFeature_RNA > mingenes & nFeature_RNA < maxgenes & percent.mito < maxPercentMT 
                          & percent.ribo < maxPercentRibo & percent.pred < maxPercentPred & nCount_RNA < maxUMIs)
#dropEST.combined = subset(x = dropEST.combined, subset = nFeature_RNA > mingenes & nFeature_RNA < maxgenes & percent.mito < maxPercentMT 
#                           & percent.pred < maxPercentPred & nCount_RNA < maxUMIs)

dropEST.combined.filtered = dropEST.combined

table(dropEST.combined.filtered@active.ident)


rm(dropEST.combined,animal,condition,directories,firstSample,firstSampleName,maxgenes,maxPercentMT,maxUMIs,mincells,mingenes,minPercentMT,
   minUMIs,maxPercentPred,maxPercentRibo,mito.features,percent.mito,pred.features,percent.pred,ribo.features,percent.ribo,
   sample,sampleName,splitSample,timePoint,tissue)
dropEST.combined.filtered <- merge(x = dropEST.combined.filtered,y = CellTypeSubset, project = "NPCproject")

table(dropEST.combined.filtered@meta.data$set_treatment)

dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "orig.ident")


dropEST.combined.filtered <- NormalizeData(object = dropEST.combined.filtered, normalization.method = "LogNormalize",scale.factor = 10000)
dropEST.combined.filtered <- FindVariableFeatures(object = dropEST.combined.filtered)
length(dropEST.combined.filtered@assays$RNA@var.features) #2000


dropEST.combined.filtered <- ScaleData(object = dropEST.combined.filtered, features = rownames(x = dropEST.combined.filtered),
                                       vars.to.regress = c("nCount_RNA","percent.mito","percent.ribo"))
dropEST.combined.filtered <- ScaleData(object = dropEST.combined.filtered, features = rownames(x = dropEST.combined.filtered),
                                       vars.to.regress = NULL)

#Perform PCA on the scaled data (uses the highly var genes)
dropEST.combined.filtered <- RunPCA(object = dropEST.combined.filtered, verbose = T, npcs = 75, ndims.print = 1:5, nfeatures.print = 10)

#Perform ICA on the scaled data (uses the highly var genes)
dropEST.combined.filtered <- RunICA(object = dropEST.combined.filtered, verbose = T, nics = 75, ndims.print = 1:5, nfeatures.print = 10)


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

setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER")
ifelse(!dir.exists(file.path(paste0("./Plots/ccaTSNEcondition"))), dir.create(file.path(paste0("./Plots/ccaTSNEcondition")),recursive = T), FALSE)
ifelse(!dir.exists(file.path(paste0("./Plots/ccaUMAPcondition"))), dir.create(file.path(paste0("./Plots/ccaUMAPcondition")),recursive = T), FALSE)

tissue <- "NPC_HEP"
pdf(file=paste0("./Plots/ccaUMAPcondition/",tissue,".pdf"))
print(DimPlot(object = seuratIntegrated, reduction = "umap", group.by = "set_treatment"))
dev.off()
pdf(file=paste0("./Plots/ccaTSNEcondition/",tissue,".pdf"))
print(DimPlot(object = seuratIntegrated, reduction = "tsne", group.by = "set_treatment"))
dev.off()

dropEST.combined.filtered <- seuratIntegrated
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
pdf(file="TreatmentUMAP.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.1, reduction = "umap")
dev.off()
pdf(file="TreatmentTSNE.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.1, reduction = "tsne")
dev.off()


#Kupffer cell
pdf("Kupffer cell.pdf")
FeaturePlot(object =dropEST.combined.filtered, features = c("Irf7","Spic","Clec4f"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
dev.off()
#sinusoidal cell
pdf("sinosoidal endothelial cells.pdf", width = 16, height = 16)
FeaturePlot(object =dropEST.combined.filtered, features = c("Ushbp1", "Myf6", "Oit3", "Il1a", "F8", "Bmp2", "C1qtnf1", "Mmrn2", "Pcdh12", "Dpp4"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
dev.off()
#hepatic satellite cells (HSC)
pdf("hepatic satellite cells.pdf", width = 16, height = 16)
FeaturePlot(object =dropEST.combined.filtered, features = c("Reln"), reduction = "umap", min.cutoff = 0, pt.size = 1)
dev.off()
#Cholangiocytes
pdf("Cholangiocytes.pdf")
FeaturePlot(object =dropEST.combined.filtered, features = c("Krt19",	"Krt7",	"Spp1","Epcam"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
dev.off()
pdf("Hepatocyte.pdf")
FeaturePlot(object =dropEST.combined.filtered, features = c("Fabp1","Alb","Fbp1"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
dev.off()
pdf("liver_zonation.pdf", width = 10,height = 6)
FeaturePlot(object =dropEST.combined.filtered, features = c("Cyp2f2","Cyp2e1"), reduction = "umap", min.cutoff = 0)
dev.off()

pdf("UMI.pdf", width = 10,height = 6)
FeaturePlot(object =dropEST.combined.filtered, features = c("nCount_RNA"), reduction = "umap", min.cutoff = 0)
dev.off()

pdf("Marker_set1.pdf",width = 10, height = 10)
FeaturePlot(object =dropEST.combined.filtered, features = c("Il7r","Cd14","Ms4a1","Cd8a","Ms4a7","Nkg7","Fcer1a","Cst3","Ppbp"), reduction = "umap", min.cutoff = 0,pt.size =0.0000000001)
dev.off()
pdf("Marker_set2.pdf",width = 10, height = 10)
FeaturePlot(object =dropEST.combined.filtered, features = c("Cd3e","Fcgr4","Lyz2","Lyz1"), reduction = "umap", min.cutoff = 0,pt.size =0.0000000001)
dev.off()
pdf("Marker_set3.pdf",width = 10, height = 10)
FeaturePlot(object =dropEST.combined.filtered, features = c("Postn","Siglech","Fabp1","Fbp1","Alb","Clec4f","Mmrn2","Lyve1","Tek"), reduction = "umap", min.cutoff = 0,pt.size =0.0000000001)
dev.off()


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
FeaturePlot(object =dropEST.combined.filtered, features = c("Siglech"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)
#Myeloid prohenitor cells
FeaturePlot(object =dropEST.combined.filtered, features = c("Flt3"), reduction = "umap", min.cutoff = 0, pt.size = 0.001)

dropEST.combined.filtered2 <- SetIdent(dropEST.combined.filtered2, value = "data.diseasestatus")
pdf(file="TreatmentUMAP_all.pdf")
DimPlot(dropEST.combined.filtered,do.label = F, reduction = "umap")
dev.off()

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "set_treatment")
save(dropEST.combined.filtered, file = "NPC_liver_combined.rda")

dropEST.combined.filtered2 <- subset(dropEST.combined.filtered ,idents = c("NPC_Control_Set1","NPC_Fruc_Set1","Hep_Control_Set1","Hep_Fruc_Set1"))
dropEST.combined.filtered2 <- SetIdent(dropEST.combined.filtered2, value = "data.diseasestatus")

pdf(file="TreatmentUMAP2.pdf")
DimPlot(dropEST.combined.filtered2,do.label = F, reduction = "umap",cols = c("black","orange"))
dev.off()
pdf(file="TreatmentTSNE2.pdf")
DimPlot(dropEST.combined.filtered2,do.label = F,pt.size = 0.1, reduction = "tsne")
dev.off()
DefaultAssay(dropEST.combined.filtered) <- "RNA"
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")
res <- FindMarkers(dropEST.combined.filtered, ident.1 = 10)

dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
                                                  c("0" = "LowUMI_HEP","1" = "SEC","2" = "Kupffer","3" = "NKT", "4" = "CV_Hep","5" = "immature_macrophage", "6" = "PN_Hep",
                                                    "7" = "Hepatic_satellite","8" = "Bcell","9" = "Proliferating_NK","10" = "Unk10","11" = "PDC","12" = "Cholangiocytes"))
save(dropEST.combined.filtered, file = "NPC_liver.rda")
