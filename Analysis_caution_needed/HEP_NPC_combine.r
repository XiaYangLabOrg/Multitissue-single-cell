

#rm(list=ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_NPC")
load("HEP_Set1_DropEST.Object1.rda")
dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
                                                  c("0" = "SEC","1" = "Kupffer","2" = "NKT","3" = "SEC","4" = "Bcell","5" = "Macrphage_Lyz2","6" = "Kupffer", "8" = "Hepatocyte","7" = "Macrphage_Flt3", "10" = "SEC_Gm26870",
                                                    "11" = "Dendritic cells","9" = "Hepatic Stellate","13" = "Cholangiocytes","12" = "Neutrophil_Top2a","14" = "Neutrophil_S100a9"))
dropEST.combined.filtered[["tissue_set_treatment"]]  <- paste0("NPC_",dropEST.combined.filtered$set_treatment)
dropEST.combined.filtered_NPC <- dropEST.combined.filtered

setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1") #29905 -> 9900, 3 sample (one seems failed), 1485 in min gene > 300
load("HEP_Set1_DropEST.Object1.rda")
dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
                                                  c("0" = "lowUMI_PN1","1" = "Perinodal","2" = "Centralvein","3" = "LowUMI_PN2", "4" = "sinusoidal cells"))


dropEST.combined.filtered <- merge(dropEST.combined.filtered, y =dropEST.combined.filtered_NPC, add.cell.ids = c("HEP", "LIVER"), project = "PBMC12K")

DefaultAssay(dropEST.combined.filtered ) <- "RNA"
dropEST.combined.filtered <- NormalizeData(object = dropEST.combined.filtered, normalization.method = "LogNormalize",scale.factor = 10000)

dropEST.combined.filtered <- FindVariableFeatures(object = dropEST.combined.filtered)
length(dropEST.combined.filtered@assays$RNA@var.features) #2000


#dropEST.combined.filtered <- ScaleData(object = dropEST.combined.filtered, features = rownames(x = dropEST.combined.filtered),
#                                       vars.to.regress = c("nCount_RNA","percent.mito","percent.ribo"))
#dropEST.combined.filtered <- RunPCA(object = dropEST.combined.filtered, verbose = T, npcs = 75, ndims.print = 1:5, nfeatures.print = 10)
#dropEST.combined.filtered <- RunICA(object = dropEST.combined.filtered, verbose = T, nics = 75, ndims.print = 1:5, nfeatures.print = 10)

#dropEST.combined.filtered = SetIdent(dropEST.combined.filtered,value="set_treatment")
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered,value="orig.ident")
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

DefaultAssay(object = seuratIntegrated) <- "integrated"
seuratIntegrated <- ScaleData(object = seuratIntegrated, verbose = FALSE)
seuratIntegrated <- RunPCA(object = seuratIntegrated, npcs = 30, verbose = FALSE)
seuratIntegrated <- RunUMAP(object = seuratIntegrated, reduction = "pca", dims = 1:30)
seuratIntegrated <- RunTSNE(object = seuratIntegrated, reduction = "pca", dims = 1:30)
seuratIntegrated <- FindNeighbors(object = seuratIntegrated, reduction = "pca", dims = 1:30, k.param = 25)
seuratIntegrated <- FindClusters(object = seuratIntegrated, resolution = seq(0.5,4,by=0.5), verbose = T, reduction = "pca")

print(DimPlot(object = seuratIntegrated, reduction = "tsne", group.by = "set_treatment"),label = T)

print(DimPlot(object = seuratIntegrated, reduction = "tsne", group.by = "celltype",label = T))
print(DimPlot(object = seuratIntegrated, reduction = "tsne", group.by = "integrated_snn_res.0.5",label = T))
print(DimPlot(object = seuratIntegrated, reduction = "tsne", group.by = "data.diseasestatus",label = T))
print(DimPlot(object = seuratIntegrated, reduction = "tsne", group.by = "batch",label = T))

dropEST.combined.filtered <- seuratIntegrated
save(dropEST.combined.filtered, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Combined_from_dataset_liver.rda")
