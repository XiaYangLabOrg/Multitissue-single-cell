#working on other doublet finder results
mincells <- 3           #min number of cells that a gene has to be detected in to be kept
mingenes <- 300         #min number of genes to keep the cell, minimum genes changed from 400
#not used in current version
#minUMIs <- 700          #min number of UMIs to keep the cell
minPercentMT <- -Inf    #min mitochondrial percentage to keep a cell
maxgenes <- 3000        #max number of genes to keep a cell
maxPercentMT <- 0.1     #max mitochondrial percentage to keep a cell
#maxPercentMT <- 0.7 #for liver cells only
maxUMIs <- 5000         #max UMIs to keep a cell (< 2* max genes )
maxPercentRibo <- 0.1 #max ribosomal percentage to keep a cell
maxPercentPred <- 0.05  #max predicted genes percentage to keep a cell

library(Seurat)
library(SingleCellExperiment)
library(DoubletFinder)
library(Seurat)


alldirectory <- list.dirs("/Users/Tsai_Lab/Desktop/Box Sync/FructoseDropSeq/Data/Final",recursive = F)
alldirectory <- alldirectory[-grep("HEP",alldirectory)]
rm("allmerged")
alldirectory <- alldirectory[-c(1:which(alldirectory %in% directory))]

for(directory in alldirectory){
  samplename <-  unlist(strsplit(directory,"/"))[9]
  organ <- unlist(strsplit(samplename, "_"))[1]
  samplename <- gsub("_SI_mixed|_liver_mixed","",samplename)
  if(organ %in% "NPC"){
    maxPercentMT <- 0.15     #max mitochondrial percentage to keep a cell
    mingenes <- 200         #
    maxPercentRibo <- 0.2
  }else{
    mingenes <- 300 
    maxgenes <- 3000        #max number of genes to keep a cell
    maxPercentMT <- 0.1     #max mitochondrial percentage to keep a cell
    #maxPercentMT <- 0.7 #for liver cells only
    maxUMIs <- 5000         #max UMIs to keep a cell (< 2* max genes )
    maxPercentRibo <- 0.1 #max ribosomal percentage to keep a cell
    maxPercentPred <- 0.05  #max predicted genes percentage to keep a cell
    
  }
  dropEST.data <- Read10X(data.dir = directory)
  dropEST.seurat <- CreateSeuratObject(counts = dropEST.data, project = samplename)
  
  mito.features <- grep(pattern = "^mt-", x = rownames(x = dropEST.seurat), value = TRUE)
  percent.mito <- Matrix::colSums(x = GetAssayData(object = dropEST.seurat, slot = "counts")[mito.features,])/Matrix::colSums(x = GetAssayData(object = dropEST.seurat, slot = "counts"))
  
  # get ribosome %
  ribo.features <- grep(pattern = "^Rps", x = rownames(x = dropEST.seurat), value = TRUE)
  ribo.features <- c(ribo.features, grep(pattern = "^Rpl", x = rownames(x = dropEST.seurat), value = TRUE))
  percent.ribo <- Matrix::colSums(x = GetAssayData(object = dropEST.seurat, slot = "counts")[ribo.features,])/Matrix::colSums(x = GetAssayData(object = dropEST.seurat, slot = "counts"))
  
  # get predicted genes %
  pred.features <- grep(pattern = "^Gm1", x = rownames(x = dropEST.seurat), value = TRUE)
  pred.features <- c(pred.features,grep(pattern = "^Gm2", x = rownames(x = dropEST.seurat), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = "^Gm3", x = rownames(x = dropEST.seurat), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = "^Gm4", x = rownames(x = dropEST.seurat), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = "^Gm5", x = rownames(x = dropEST.seurat), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = "^Gm6", x = rownames(x = dropEST.seurat), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = "^Gm7", x = rownames(x = dropEST.seurat), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = "^Gm8", x = rownames(x = dropEST.seurat), value = TRUE))
  pred.features <- c(pred.features,grep(pattern = "^Gm9", x = rownames(x = dropEST.seurat), value = TRUE))
  percent.pred <- Matrix::colSums(x = GetAssayData(object = dropEST.seurat, slot = "counts")[pred.features,])/Matrix::colSums(x = GetAssayData(object = dropEST.seurat, slot = "counts"))
  
  # Add % mito, ribo and pred to the meta data
  dropEST.seurat[["percent.mito"]] <- percent.mito
  dropEST.seurat[["percent.ribo"]] <- percent.ribo
  dropEST.seurat[["percent.pred"]] <- percent.pred
  dropEST.seurat = subset(x = dropEST.seurat, subset = nFeature_RNA > mingenes & nFeature_RNA < maxgenes & percent.mito < maxPercentMT 
                          & percent.ribo < maxPercentRibo & percent.pred < maxPercentPred & nCount_RNA < maxUMIs)
  
  
  
  dropEST.seurat <- NormalizeData(dropEST.seurat)
  dropEST.seurat <- FindVariableFeatures(dropEST.seurat, selection.method = "vst", nfeatures = 2000)
  dropEST.seurat <- ScaleData(dropEST.seurat)
  dropEST.seurat <- RunPCA(dropEST.seurat, npcs = 30)
  dropEST.seurat <- RunUMAP(dropEST.seurat, dims = 1:10)
  
  dropEST.combined.filtered <- dropEST.seurat
  rm(dropEST.seurat, dropEST.data)
  
  stdv <- dropEST.combined.filtered[["pca"]]@stdev
  sum.stdv <- sum(dropEST.combined.filtered[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  sweep.res.list_kidney <-  paramSweep_v3(dropEST.combined.filtered, PCs = 1:10, sct = FALSE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn <- find.pK(sweep.stats_kidney)
  
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  
  annotations <- dropEST.combined.filtered@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(dropEST.combined.filtered@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  mouse.sample <- doubletFinder_v3(seu = dropEST.combined.filtered, 
                                   PCs = 1:min.pc, 
                                   pK = optimal.pk,
                                   nExp = nExp.poi.adj)
  
  metadata <- mouse.sample@meta.data
  colnames(metadata)[5] <- "doublet_finder"
  classificationcol <- grep("classification",colnames(metadata))
  
  result <- data.frame(name = rownames(metadata),doublet = metadata$doublet_finder,prediction = metadata[,classificationcol], file = samplename)
  if(!exists("allmerged")){
    allmerged <- result
  }else{
    allmerged <- rbind(allmerged, result)
  }
}
allmerged_doublet_finder_nonhep <- allmerged
save(allmerged_doublet_finder_nonhep, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Doubletfinder_nonhep.rda")

rm(list = setdiff(ls(), lsf.str()))
firstup <- function(x) {
  result <- NULL
  for(i in 1:length(x)){
    current <- x[i]
    current <- tolower(current)
    if(startsWith(current,"mt-")){
      substr(current, 4, 4) <- toupper(substr(current, 4, 4))
    }else{
      substr(current, 1, 1) <- toupper(substr(current, 1, 1))
    }
    result[i] <- current
  }
  result
}
sampleconverter <- function(string){
  number <- as.numeric(substr(string,5,5))
  for(i in 1:length(number)){
    if(number[i] %in% c(3,4)){number[i] <- number[i] - 2
    }else if(number[i] %in% c(5,6,7)){number[i] <- number[i] - 4
    }else if(number[i] %in% 8){number[i] <-  number[i] - 5}
  }
  return(number)
}

library(dplyr)
library(Seurat)
library(plyr)
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Doubletfinder_nonhep.rda")

load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Doubletfinder.rda")
allmerged <- rbind.data.frame(allmerged_doublet_finder,allmerged_doublet_finder_nonhep)

#liver

#rm(list = setdiff(ls(), lsf.str()))
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/portal_data")
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Combined_allliver.rda")
dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
                                                  c("0" = "Sinusoidal endothelial cells","1" = "Periportal hepatocytes","2" = "Periportal hepatocytes","3" = "Kupffer", "4" = "Pericentral hepatocytes","5" = "Periportal hepatocytes", "6" = "NKT cells",
                                                    "7" = "Periportal hepatocytes","8" = "Classical dendritic cells","9" = "B cells","10" = "Hepatic stellate cells","11" = "Cholangiocytes","12" = "Dividing cells","13" = "Plasmacytoid dendritic cells"))

DefaultAssay(dropEST.combined.filtered) <- "RNA"
table(dropEST.combined.filtered$orig.ident, dropEST.combined.filtered$tissue)
#Making changes to have sample names sorted
dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set12.1"] <- "Set11.1"
dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set16.1"] <- "Set15.1"
dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set22.1" & dropEST.combined.filtered$tissue %in% "HEP"] <- "Set21.1"
dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set26.3" & dropEST.combined.filtered$tissue %in% "NPC"] <- "Set27.1"
dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set24.2"] <- "Set28.2"
dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set23.3"] <- "Set24.3"

dropEST.combined.filtered[["samplename"]] <- paste0(dropEST.combined.filtered$tissue,"_",dropEST.combined.filtered$data.diseasestatus,"_",dropEST.combined.filtered$batch,"_", sampleconverter(dropEST.combined.filtered$orig.ident))

test <- dropEST.combined.filtered@meta.data
rownames(test) <- gsub("^_","",rownames(test))

test$decontname <-  paste0(sapply(strsplit(rownames(test),"_"), `[[`, 3),"_",test$samplename )
allmerged$decontname <- paste0(allmerged$name,"_",allmerged$file)
test7 <- left_join(test,allmerged)
test7$prediction[is.na(test7$prediction)] <- "Singlet"
test8 <- aggregate(test7[,c("nCount_RNA")], list(test7$integrated_snn_res.0.5),median, na.rm = T)
colnames(test8) <- c("integrated_snn_res.0.5","Avg_UMI")
test7 <- left_join(test7,test8)

dropEST.combined.filtered[["doublet_finder"]] <- test7$prediction
dropEST.combined.filtered[["Avg_UMI"]] <- test7$Avg_UMI
#FeaturePlot(dropEST.combined.filtered,reduction = "tsne", "doublet_finder")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value =  "celltype")
DimPlot(dropEST.combined.filtered,reduction = "tsne")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value =  "doublet_finder")
DimPlot(dropEST.combined.filtered,reduction = "tsne")
FeaturePlot(dropEST.combined.filtered,"nCount_RNA", reduction = "tsne")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value =  "integrated_snn_res.0.5")
DimPlot(dropEST.combined.filtered,reduction = "tsne")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value =  "Avg_UMI")
DimPlot(dropEST.combined.filtered,reduction = "tsne", label = T)

#to remove
test7 <- dropEST.combined.filtered@meta.data
to_remove <- which(test7$nCount_RNA < 1000 & test7$data.tissue %in% "Hep")
to_remove <- union(to_remove, which(test7$doublet_finder %in% "Doublet"))
to_remove <- rownames(dropEST.combined.filtered@meta.data)[to_remove]
dropEST.combined.filtered2 <-  dropEST.combined.filtered[,!colnames(dropEST.combined.filtered) %in% to_remove]
dropEST.combined.filtered2 <- SetIdent(dropEST.combined.filtered2, value =  "celltype")
DimPlot(dropEST.combined.filtered2,reduction = "tsne")
length(to_remove)/nrow(test7)

save(dropEST.combined.filtered2, file = "/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/liver.rda")
#SI
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/SI_combined_Set2_final.rda")
dropEST.combined.filtered[["samplename"]] <- paste0(dropEST.combined.filtered$data.diseasestatus,"_", sampleconverter(dropEST.combined.filtered$orig.ident))
DefaultAssay(dropEST.combined.filtered) <- "RNA"
test <- dropEST.combined.filtered@meta.data
rownames(test) <- gsub("^_","",rownames(test))

test$decontname <-  paste0(sapply(strsplit(rownames(test),"_"), `[[`, 2),"_SI_",test$samplename )
allmerged$decontname <- paste0(allmerged$name,"_",allmerged$file)
test7 <- left_join(test,allmerged)
test8 <- aggregate(test7[,c("nCount_RNA")], list(test7$integrated_snn_res.0.5),median, na.rm = T)
colnames(test8) <- c("integrated_snn_res.0.5","Avg_UMI")
test7 <- left_join(test7,test8)

dropEST.combined.filtered[["doublet_finder"]] <- test7$prediction
dropEST.combined.filtered[["Avg_UMI"]] <- test7$Avg_UMI

#FeaturePlot(dropEST.combined.filtered,reduction = "tsne", "doublet_finder")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value =  "celltype")
DimPlot(dropEST.combined.filtered,reduction = "tsne")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value =  "doublet_finder")
DimPlot(dropEST.combined.filtered,reduction = "tsne")
FeaturePlot(dropEST.combined.filtered,"nCount_RNA", reduction = "tsne")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value =  "Avg_UMI")
DimPlot(dropEST.combined.filtered,reduction = "tsne", label = T)

#Remove
test7 <- dropEST.combined.filtered@meta.data
to_remove <-  which(test7$doublet_finder %in% "Doublet")
to_remove <- rownames(dropEST.combined.filtered@meta.data)[to_remove]
dropEST.combined.filtered2 <-  dropEST.combined.filtered[,!colnames(dropEST.combined.filtered) %in% to_remove]
dropEST.combined.filtered2 <- SetIdent(dropEST.combined.filtered2, value =  "celltype")
DimPlot(dropEST.combined.filtered2,reduction = "tsne")
length(to_remove)/nrow(test7)
save(dropEST.combined.filtered2, file = "/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SI.rda")

#SVF
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined/SVF_Final_Set1only.rda")
dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set16.1"] <- "Set15.1"
dropEST.combined.filtered[["samplename"]] <- paste0(dropEST.combined.filtered$data.diseasestatus,"_", sampleconverter(dropEST.combined.filtered$orig.ident))

DefaultAssay(dropEST.combined.filtered) <- "RNA"
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "samplename")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents = c("Control_1","Control_2","Fruc_2","HFHS_1","HFHS_2"))

test <- dropEST.combined.filtered@meta.data
rownames(test) <- gsub("^_","",rownames(test))

test$decontname <-  paste0(sapply(strsplit(rownames(test),"_"), `[[`, 2),"_SVF_",test$samplename )
allmerged$decontname <- paste0(allmerged$name,"_",allmerged$file)
test7 <- left_join(test,allmerged)
dropEST.combined.filtered[["doublet_finder"]] <- test7$prediction
#FeaturePlot(dropEST.combined.filtered,reduction = "tsne", "doublet_finder")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value =  "celltype")
DimPlot(dropEST.combined.filtered,reduction = "tsne")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value =  "doublet_finder")
DimPlot(dropEST.combined.filtered,reduction = "tsne")


test7 <- dropEST.combined.filtered@meta.data
to_remove <-  which(test7$doublet_finder %in% "Doublet")
to_remove <- rownames(dropEST.combined.filtered@meta.data)[to_remove]
dropEST.combined.filtered2 <-  dropEST.combined.filtered[,!colnames(dropEST.combined.filtered) %in% to_remove]
dropEST.combined.filtered2 <- SetIdent(dropEST.combined.filtered2, value =  "celltype")
DimPlot(dropEST.combined.filtered2,reduction = "tsne")
length(to_remove)/nrow(test7)
length(to_remove)
save(dropEST.combined.filtered2, file = "/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SVF.rda")


#Hyp

load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Hyp_combined_Set1_final.rda")
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/portal_data")
dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set11.3"] <- "Set12.3"
dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set13.3"] <- "Set14.3"

dropEST.combined.filtered[["samplename"]] <- paste0(dropEST.combined.filtered$data.diseasestatus,"_", sampleconverter(dropEST.combined.filtered$orig.ident))

DefaultAssay(dropEST.combined.filtered) <- "RNA"

test <- dropEST.combined.filtered@meta.data
rownames(test) <- gsub("^_","",rownames(test))

test$decontname <-  paste0(sapply(strsplit(rownames(test),"_"), `[[`, 2),"_Hyp_",test$samplename )
allmerged$decontname <- paste0(allmerged$name,"_",allmerged$file)
test7 <- left_join(test,allmerged)
dropEST.combined.filtered[["doublet_finder"]] <- test7$prediction
#FeaturePlot(dropEST.combined.filtered,reduction = "tsne", "doublet_finder")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value =  "celltype")
DimPlot(dropEST.combined.filtered,reduction = "tsne")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value =  "doublet_finder")
DimPlot(dropEST.combined.filtered,reduction = "tsne")


test7 <- dropEST.combined.filtered@meta.data
to_remove <-  which(test7$doublet_finder %in% "Doublet")
to_remove <- rownames(dropEST.combined.filtered@meta.data)[to_remove]
dropEST.combined.filtered2 <-  dropEST.combined.filtered[,!colnames(dropEST.combined.filtered) %in% to_remove]
dropEST.combined.filtered2 <- SetIdent(dropEST.combined.filtered2, value =  "celltype")
DimPlot(dropEST.combined.filtered2,reduction = "tsne")
length(to_remove)/nrow(test7)
length(to_remove)

save(dropEST.combined.filtered2, file = "/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/Hyp.rda")
