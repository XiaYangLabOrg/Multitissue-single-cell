#dropseq analysis
#CCA integration based on sample condition (if some sample is of low numbers -> integrating on sample level will do worse)
#Seurat label transfer function
library(metap)
library(Seurat)
mincells <- 3           #min number of cells that a gene has to be detected in to be kept
mingenes <- 300         #min number of genes to keep the cell, minimum genes changed from 400
#not used in current version
#minUMIs <- 700          #min number of UMIs to keep the cell
minPercentMT <- -Inf    #min mitochondrial percentage to keep a cell
maxgenes <- 3000        #max number of genes to keep a cell
maxPercentMT <- 0.1     #max mitochondrial percentage to keep a cell
maxUMIs <- 5000         #max UMIs to keep a cell (< 2* max genes )
maxPercentRibo <- 0.1 #max ribosomal percentage to keep a cell
maxPercentPred <- 0.05  #max predicted genes percentage to keep a cell

#hypothalamus analysis full
tissue = "Hip"
#hypothalamus combined
directories = c("/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/Hyp/F3_Hyp_Fruc_3.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/Hyp/F1_Hyp_Fruc_3.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/Hyp/F4_Hyp_Control_1.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/Hyp/F2_Hyp_Control_1.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/Hyp/F3_Hyp_HFHS_5.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/Hyp/F1_Hyp_HFHS_6.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/Hyp/F1_Hyp_HFHS_5.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/Hyp/F6_Hyp_HFHS_6.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/Hyp/F3_Hyp_Fruc_4.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/Hyp/F4_Hyp_Fruc_4.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Control/Hyp/F5_Hyp_Control_1.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Control/Hyp/F2_Hyp_Control_2.3")

directories = c("/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/Hyp/F3_Hyp_Fruc_3.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Fructose/Hyp/F1_Hyp_Fruc_3.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/Hyp/F4_Hyp_Control_1.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/Control/Hyp/F2_Hyp_Control_1.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/Hyp/F3_Hyp_HFHS_5.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set1/HFHS/Hyp/F1_Hyp_HFHS_6.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/Hyp/F1_Hyp_HFHS_5.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/HFHS/Hyp/F6_Hyp_HFHS_6.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/Hyp/F3_Hyp_Fruc_4.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Fructose/Hyp/F4_Hyp_Fruc_4.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Control/Hyp/F5_Hyp_Control_1.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set2/Control/Hyp/F2_Hyp_Control_2.3",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set3/Control/Hyp/F2_Hyp_Control_2.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set3/Control/Hyp/F2_Hyp_Control_2.4",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set3/Control/Hyp/F2_Hyp_Control_2.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set3/Fructose/Hyp/F1_Hyp_Fruc_1.1",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set3/Fructose/Hyp/F1_Hyp_Fruc_1.2",
                "/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Set3/Fructose/Hyp/F2_Hyp_Fruc_2.x")

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

setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/batch3")

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
rm(dropEST.combined.filtered)
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

tissue <- "HYP"
pdf(file=paste0("./Plots/ccaUMAPcondition/",tissue,".pdf"))
print(DimPlot(object = seuratIntegrated, reduction = "umap", group.by = "set_treatment"))
dev.off()
pdf(file=paste0("./Plots/ccaTSNEcondition/",tissue,".pdf"))
print(DimPlot(object = seuratIntegrated, reduction = "tsne", group.by = "set_treatment"))
dev.off()

dropEST.combined.filtered <- seuratIntegrated  



#####separating neuronal cells
########
seuratIntegrated <- SetIdent(seuratIntegrated,value="set_treatment")
selectedcells <- FeatureLocator(DimPlot(object = seuratIntegrated, do.identify = TRUE,  reduction = "umap",pt.size = 0.5))

CellTypeSubset = subset(seuratIntegrated, cells = selectedcells) 

DefaultAssay(CellTypeSubset) <- "RNA" 
CellTypeSubset = SetIdent(CellTypeSubset,value="set_treatment")
seurat.list = list()
for(i in 1:length(levels(CellTypeSubset@active.ident))){
  seurat.list[[i]] = subset(CellTypeSubset,idents=levels(CellTypeSubset@active.ident)[i])
  cat(ncol(seurat.list[[i]])," ")
}

for (i in 1:length(x = seurat.list)) {
  seurat.list[[i]] <- NormalizeData(object = seurat.list[[i]], verbose = FALSE)
  seurat.list[[i]] <- FindVariableFeatures(object = seurat.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

#This works for small dataset integration
k.filter <- min(200, min(sapply(seurat.list, ncol)))
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:30, k.filter = k.filter)

CellTypeSubset <- IntegrateData(anchorset = seurat.anchors, dims = 1:30)

DefaultAssay(object = CellTypeSubset) <- "integrated"

CellTypeSubset  <- ScaleData(object = CellTypeSubset , verbose = FALSE)
CellTypeSubset  <- RunPCA(object = CellTypeSubset , npcs = 30, verbose = FALSE)
CellTypeSubset  <- FindNeighbors(object = CellTypeSubset , reduction = "pca", dims = 1:30, k.param = 25)
CellTypeSubset  <- FindClusters(object = CellTypeSubset , resolution = seq(0.5,4,by=0.5), verbose = T, reduction = "pca")
CellTypeSubset  <- RunUMAP(object = CellTypeSubset , reduction = "pca", dims = 1:30)
CellTypeSubset  <- RunTSNE(object = CellTypeSubset , reduction = "pca", dims = 1:30)

pdf("UMIcountplot_subset.pdf")
FeaturePlot(CellTypeSubset,"nCount_RNA",reduction = "tsne")
dev.off()

CellTypeSubset <- SetIdent(CellTypeSubset, value = "batch")
pdf(file="batch_SampleUMAP_subset.pdf")
DimPlot(CellTypeSubset,do.label = F,pt.size = 1, reduction = "umap")
dev.off()
pdf(file="batch_SampleTSNE_subset.pdf")
DimPlot(CellTypeSubset,do.label = F,pt.size = 1, reduction = "tsne")
dev.off()

CellTypeSubset <- SetIdent(CellTypeSubset, value = "data.diseasestatus")
pdf(file="TreatmentUMAP_subset.pdf")
DimPlot(CellTypeSubset,do.label = F,pt.size = 0.1, reduction = "umap")
dev.off()
pdf(file="TreatmentTSNE_subset.pdf")
DimPlot(CellTypeSubset,do.label = F,pt.size = 0.1, reduction = "tsne")
dev.off()

CellTypeSubset <- SetIdent(CellTypeSubset, value = "orig.ident")
pdf(file="sampleUMAP_subset.pdf")
DimPlot(CellTypeSubset,do.label = F,pt.size = 0.1, reduction = "umap")
dev.off()
pdf(file="sampleTSNE_subset.pdf")
DimPlot(CellTypeSubset,do.label = F,pt.size = 0.1, reduction = "tsne")
dev.off()

CellTypeSubset <- SetIdent(CellTypeSubset, value = "set_treatment")
pdf(file="sample_setUMAP_subset.pdf")
DimPlot(CellTypeSubset,do.label = F,pt.size = 0.1, reduction = "umap")
dev.off()
pdf(file="sample_setTSNE_subset.pdf")
DimPlot(CellTypeSubset,do.label = F,pt.size = 0.1, reduction = "tsne")
dev.off()
table(CellTypeSubset@meta.data$set_treatment)
save(CellTypeSubset,file = "Neuroncelltype.rda")
#######


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


dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "batch")
pdf(file="batch_SampleUMAP.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.01, reduction = "umap")
dev.off()
pdf(file="batch_SampleTSNE.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.01, reduction = "tsne")
dev.off()

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
pdf(file="TreatmentUMAP.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.1, reduction = "umap")
dev.off()
pdf(file="TreatmentTSNE.pdf")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.1, reduction = "tsne")
dev.off()

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "orig.ident")
DimPlot(dropEST.combined.filtered,do.label = F,pt.size = 0.1, reduction = "tsne")

save(dropEST.combined.filtered, file = "SI_Set1_DropEST.Object1_3.rda")

#identification of cell types
library(Seurat)
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
load("SI_Set1_DropEST.Object1.rda")

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")
DefaultAssay(dropEST.combined.filtered) <- "RNA"

#cluster identification
ifelse(!dir.exists(file.path("./IdentifyClusters")), dir.create(file.path("./IdentifyClusters")), FALSE)
ifelse(!dir.exists(file.path("./IdentifyClusters","./MarkerGenes")), dir.create(file.path("./IdentifyClusters","./MarkerGenes")), FALSE)
ifelse(!dir.exists(file.path("./IdentifyClusters","./MarkerGenes","./Hyp")), dir.create(file.path("./IdentifyClusters","./MarkerGenes","./Hyp")), FALSE)

cell_type_genes = c("Snap25","Syt1","Olig1","Mobp","Fyn","Pdgfra","Agt","Aldoc","Aqp4","Gfap","Rax","Ccdc153","Flt1","Aif1","Cx3cr1","Mrc1","Acta2","Col1a1","Lum","Folr1")
cell_type_genes_description = c("Pan.Neuronal1","Pan.Neuronal2","Pan.Oligodendrocytes","Mylenating.Oligodendrocytes","Newly.Formed.Oligo",
                                "Oligo.PCs","Astrocytes","Astrocytes2","Astrocytes3","Astrocytes4","Tanycyte","Ependymal","Endothelial","Perivascular.Macrophages.Microglia",
                                "Microglia","Macrophages","Mural.Pericytes.VSMCs","Vascular.Leptomeningeal.Cells1","Vascular.Leptomeningeal.Cells2","ChoroidPlexusEpethelial")

for(i in 1:length(cell_type_genes)){
  pdf(file=paste0("./IdentifyClusters/MarkerGenes/Hyp/FeaturePlot_",cell_type_genes_description[i],"MarkerGene.pdf"))
  print(FeaturePlot(dropEST.combined.filtered,cell_type_genes[i],cols = c("lightgrey","blue"), reduction = "umap"))
  dev.off()
}
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.2.5")
for(i in 1:length(cell_type_genes)){
  p = VlnPlot(dropEST.combined.filtered,features = cell_type_genes[i],sort=T)
  pdf(file=paste0("./IdentifyClusters/MarkerGenes/Hyp/VlnPlot_",cell_type_genes_description[i],"MarkerGene.pdf"))
  print(p)
  dev.off()
}

#Neuronal subtype identification
library(Seurat)
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Neuroncelltype.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/batch3/Neuroncelltype.rda")

#Activation, inhibition
DefaultAssay(CellTypeSubset) -> "RNA"
FeaturePlot(CellTypeSubset,c("Gad1","Slc17a6"), reduction = "umap",min.cutoff = 0)
FeaturePlot(CellTypeSubset,c("Gad2","Nnat"), reduction = "umap",min.cutoff = 0)


resolutions <- c("0.5","1","1.5","2","2.5","3","3.5","4")

dir.create(file.path("./SampleTSNE_sub"))
dir.create(file.path("./SampleUMAP_sub"))

for(reso in resolutions){
  pdf(file=paste0("./SampleTSNE_sub/",reso,".pdf"))
  CellTypeSubset <- SetIdent(CellTypeSubset, value = paste0("integrated_snn_res.",reso))
  print(DimPlot(CellTypeSubset, label = T, reduction = "tsne"))
  dev.off()
  pdf(file=paste0("./SampleUMAP_sub/",reso,".pdf"))
  CellTypeSubset <- SetIdent(CellTypeSubset, value = paste0("integrated_snn_res.",reso))
  print(DimPlot(CellTypeSubset, label = T, reduction = "umap"))
  dev.off()
}



#vasopressin and Oxytocin for Supraoptic nucleus
pdf("Supraoptic nucleus.pdf",width = 6,height = 6)
FeaturePlot(CellTypeSubset,c("Avp","Oxt"), reduction = "umap",min.cutoff = 0)
dev.off()

#TRH for Paraventricular nucleus
pdf("Paraventricular nucleus.pdf",width = 6,height = 6)
FeaturePlot(CellTypeSubset,c("Trh","Crh","Sst"), reduction = "umap",min.cutoff = 0)
dev.off()

#Lateral hypothalamus, non-specific expression
pdf("Lateral hypothalamus.pdf",width = 6,height = 6)

FeaturePlot(CellTypeSubset,c("Hcrt"), reduction = "umap",min.cutoff = 0)
dev.off()

#Gabanergic neuron
pdf("Gabanergic.pdf",width = 6,height = 6)
FeaturePlot(CellTypeSubset,c("Gad1","Abat","Slc6a1"), reduction = "umap",min.cutoff = 0)
dev.off()

#Growth hormone–releasing hormone
pdf("Ghrh.pdf",width = 6,height = 6)

FeaturePlot(CellTypeSubset,c("Ghrh"), reduction = "umap",min.cutoff = 0)
dev.off()

#Medial preoptic nucleus, Sexually related
pdf("Medial preoptic nucleus Sexually related.pdf",width = 6,height = 6)
FeaturePlot(CellTypeSubset,c("Gnrh1","Gnrh2"), reduction = "umap",min.cutoff = 0)
dev.off()

#dopaminergic neuron
pdf("Dopaminergic.pdf",width = 6,height = 6)

FeaturePlot(CellTypeSubset,c("Ddc","Th"), reduction = "umap",min.cutoff = 0)
dev.off()

#arousal related neuron
pdf("arousal related neuron.pdf",width = 6,height = 6)

FeaturePlot(CellTypeSubset,c("Esr1"), reduction = "umap",min.cutoff = 0)
dev.off()
#Thermorregulation? Galantin?
pdf("Thermorregulation.pdf",width = 6,height = 6)

FeaturePlot(CellTypeSubset,c("Gal"), reduction = "umap",min.cutoff = 0)
dev.off()
#Circadian genes, not working
pdf("Circadian.pdf",width = 6,height = 6)

FeaturePlot(CellTypeSubset,c("Per1","Per2","Per3","Cry1","Cry2"), reduction = "umap",min.cutoff = 0)
dev.off()
#NPY genes
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4066284/
pdf("NPY.pdf",width = 6,height = 6)

FeaturePlot(CellTypeSubset,c("Npy","Npy2r","Agrp"), reduction = "umap",min.cutoff = 0)
dev.off()

#Prolactin genes-preoptic hypothalamus
#https://www.ncbi.nlm.nih.gov/pubmed/10668647
pdf("Prolactin.pdf",width = 6,height = 6)

FeaturePlot(CellTypeSubset,c("Prlr"), reduction = "umap",min.cutoff = 0)
dev.off()

#pan-glutamatergic marker
pdf("pan-glutamatergic.pdf",width = 6,height = 6)

FeaturePlot(CellTypeSubset,c("Slc32a1"), reduction = "umap",min.cutoff = 0)
dev.off()
#Prph-posterior hypothalamus
pdf("Prph-posterior hypothalamu.pdf",width = 6,height = 6)

FeaturePlot(CellTypeSubset,c("Prph"), reduction = "umap",min.cutoff = 0)
dev.off()
#Grp https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2646528/
pdf("Grp.pdf",width = 6,height = 6)

FeaturePlot(CellTypeSubset,c("Grp"), reduction = "umap",min.cutoff = 0)
dev.off()

CellTypeSubset <- SetIdent(CellTypeSubset,value = "integrated_snn_res.0.5")
Degs = FindMarkers(CellTypeSubset , ident.1 = "12")

Degs_9 = FindMarkers(CellTypeSubset , ident.1 = "9")

Degs_3 = FindMarkers(CellTypeSubset , ident.1 = "3")
Degs_7 = FindMarkers(CellTypeSubset , ident.1 = "7")
Degs_4 = FindMarkers(CellTypeSubset , ident.1 = "4")
Degs_11 = FindMarkers(CellTypeSubset , ident.1 = "11")
Degs_8 = FindMarkers(CellTypeSubset , ident.1 = "8")
Degs_2 = FindMarkers(CellTypeSubset , ident.1 = "2")


#DEG analysis-full

library(Seurat)
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/SI_Set1_DropEST.Object1.rda")

load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Neuroncelltype.rda")
dropEST.combined.filtered <- CellTypeSubset
#single
Degs = list()
topDEGs <- list()
pathwaynum <- NULL
DEGnum <- NULL
#dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "RNA_snn_res.1")
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")

table(dropEST.combined.filtered@active.ident) 
alllevels <- levels(dropEST.combined.filtered@active.ident)

for(i  in 1:length(alllevels)){
  CellType <- alllevels[i]
  CellTypeSubset = subset(dropEST.combined.filtered, idents = CellType)
  
  CellTypeSubset = SetIdent(CellTypeSubset, value = "data.diseasestatus")
  siggenes <- NULL
  try({Degs[[CellType]] = FindMarkers(CellTypeSubset, ident.1 = "Fruc", ident.2 = "Control")})
  #Degs[[CellType]] = FindMarkers(CellTypeSubset, ident.1 = "water", ident.2 = "BPA10")
  siggenes <- rownames(Degs[[CellType]][Degs[[CellType]]$p_val_adj < 0.05,])
  if(length(grep("mt|Rps|Rpl",siggenes)) > 0){
    siggenes <- siggenes[-grep("mt|Rps|Rpl",siggenes)]
  }
  
  if(length(siggenes) == 0){
    DEGnum[i] <- NA
    pathwaynum[i] <- NA
    next
  }
  DEGnum[i] <- length(siggenes)
  topDEGs[[CellType]] <- siggenes
  #write.csv(Degs[[CellType]][Degs[[CellType]]$p_val_adj < 0.05,],file = paste0(getwd(),"/",CellType,".csv"))
  write.csv(Degs[[CellType]][Degs[[CellType]]$p_val_adj < 0.05,],file = paste0(getwd(),"/sub",CellType,".csv"))
  
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
    
    #write.csv(drugframe,file = paste0(getwd(),"/",CellType,"DEGpathways.csv"))
    write.csv(drugframe,file = paste0(getwd(),"/sub",CellType,"DEGpathways.csv"))
    pathwaynum[i] <- nrow(drugframe)
  }else{
    pathwaynum[i] <- NA
  }
}

finalframe <- data.frame(subpopulation = names(table(dropEST.combined.filtered@active.ident)),
                         cellnumber = as.numeric(table(dropEST.combined.filtered@active.ident) ),
                         DEGnum = DEGnum,
                         pathwaynum = pathwaynum)
#write.csv(finalframe,"allcells_summary_DEG.csv")
write.csv(finalframe,"subcells_summary_DEG.csv")

##Testing on full DEG space
library(Seurat)
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/SI_Set1_DropEST.Object1.rda")

load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Neuroncelltype.rda")
dropEST.combined.filtered <- CellTypeSubset
#combined results for full feature space
Degs = list()
topDEGs <- list()
DEGnum <- NULL
pathwaynum <- NULL
#dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "RNA_snn_res.1")
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")

table(dropEST.combined.filtered@active.ident) 
alllevels <- levels(dropEST.combined.filtered@active.ident)
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
  #write.csv(Degs[[CellType]][Degs[[CellType]]$p_val_adj < 0.05,],file = paste0(getwd(),"/",CellType,"full.csv"))
  write.csv(Degs[[CellType]][Degs[[CellType]]$p_val_adj < 0.05,],file = paste0(getwd(),"/neuron",CellType,"full.csv"))
  
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
    #write.csv(drugframe,file = paste0(getwd(),"/",CellType,"DEGpathways_full.csv"))
    write.csv(drugframe,file = paste0(getwd(),"/",CellType,"neuron_DEGpathways_full.csv"))
    pathwaynum[i] <- nrow(drugframe)
  }else{
    pathwaynum[i] <- NA
  }
}

finalframe <- data.frame(subpopulation = names(table(dropEST.combined.filtered@active.ident)),
                         cellnumber = as.numeric(table(dropEST.combined.filtered@active.ident) ),
                         DEGnum = DEGnum,
                         pathwaynum = pathwaynum)
#write.csv(finalframe,"summary_DEG_full.csv")
write.csv(finalframe,"neuron_summary_DEG_full.csv")


