library(Seurat)
library(plyr)
#library(dplyr)

#library(cerebroApp)
library(ggplot2)
library(gridExtra)
rm(list=ls())
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



#Cross DEG analysis
#Final All tSNE
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/SI_combined_Set2_final.rda")
#dropEST.combined.filtered[["samplename"]] <- paste0(dropEST.combined.filtered$data.diseasestatus,"_", sampleconverter(dropEST.combined.filtered$orig.ident))
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SI.rda")
intestine_frame <- unique(dropEST.combined.filtered@meta.data[,c("orig.ident","samplename")])

#DimPlot(dropEST.combined.filtered,group.by = "data.diseasestatus",reduction = "tsne",label = F, cols = c("Control"="darkslategray","Fruc" = "cyan3","HFHS" = "red"))+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=24))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("")
SI_Fig1 <- DimPlot(dropEST.combined.filtered,group.by = "data.diseasestatus",reduction = "tsne",label = F, 
                   cols = c("Control"="darkslategray","Fruc" = "cyan3","HFHS" = "hotpink"))+ 
  theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=3)))+ggtitle("")
#DimPlot(dropEST.combined.filtered,group.by = "samplename",reduction = "tsne",label = F)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=24))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("")
SI_Sup1 <- DimPlot(dropEST.combined.filtered,group.by = "samplename",reduction = "tsne",label = F)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=4)))+ggtitle("")
SI_Sup6 <- DimPlot(dropEST.combined.filtered,group.by = "celltype",reduction = "tsne",label = T, repel = T, label.size = 8)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=8)))+ggtitle("")

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "samplename")
SI_Sup2 <- VlnPlot(dropEST.combined.filtered,"nFeature_RNA",pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Gene count")+ theme(legend.position = "none")
SI_Sup3 <- VlnPlot(dropEST.combined.filtered,"nCount_RNA",pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("UMI count")+ theme(legend.position = "none")
SI_Sup4 <- VlnPlot(dropEST.combined.filtered,"percent.mito",pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Mitochondrial count percent")+ theme(legend.position = "none")
SI_Sup5 <- VlnPlot(dropEST.combined.filtered,"percent.ribo",pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Robisomal count percent")+ theme(legend.position = "none")
SI_Sup7 <- VlnPlot(dropEST.combined.filtered,"percent.pred",pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Robisomal count percent")+ theme(legend.position = "none")



DefaultAssay(dropEST.combined.filtered) <- "RNA"
cellTypeOrder = c("Distal enterocytes","Proximal enterocytes","T cells","Goblet","Macrophages","I cells","S cells","Tuft","Proliferative cells","B cells")
# cellTypesToRemove = c("Neurons","Cluster16","Cluster23")
dropEST.combined.filtered$celltype <- factor(dropEST.combined.filtered$celltype, levels = cellTypeOrder)
genesOrder = c("Fabp6","Pmp22","Lct","Ephx2","Gsta1","Gzma","Tff3","Lyz2","Cck","Chga","Dclk1","Top2a","Igj")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered,value = "celltype")

SI_Dot <- DotPlot(dropEST.combined.filtered, features = genesOrder,assay = "RNA",scale.by = "size", col.min = 0)+ylab("")+xlab("")+ theme(axis.text.x = element_text(size = 24, angle = 45, hjust = 1, face = "italic"),
                                                                                                                                          axis.text.y = element_text(size = 24),axis.title = element_text(size = 24),
                                                                                                                                          legend.text=element_text(size=20))+ guides(colour = guide_legend(override.aes = list(size=5)))


pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/SupFig1_SI.pdf",width = 24, height = 16)
grid.arrange(SI_Sup1,SI_Sup2,SI_Sup3,SI_Sup4,SI_Sup5,SI_Sup6,SI_Dot,layout_matrix = rbind(c(2,3,1,1),c(4,5,1,1),c(6,6,7,7),c(6,6,7,7)))
dev.off()

#SVF
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined/SVF_Final_Set1only.rda")
#dropEST.combined.filtered$celltype[dropEST.combined.filtered$celltype %in% "Macrophage_Il1b"] <- "M1 Macrophage"
#dropEST.combined.filtered$celltype[dropEST.combined.filtered$celltype %in% "Macrophage_Pf4"] <- "M2 Macrophage"
#save(dropEST.combined.filtered, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined/SVF_Final_Set1only.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SVF.rda")
dropEST.combined.filtered$celltype[dropEST.combined.filtered$celltype %in% "Macrophage_Pf4"]
#dropEST.combined.filtered$celltype <- gsub("APC_","",dropEST.combined.filtered$celltype )
DefaultAssay(dropEST.combined.filtered) <- "RNA"
#dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set16.1"] <- "Set15.1"
#dropEST.combined.filtered[["samplename"]] <- paste0(dropEST.combined.filtered$data.diseasestatus,"_", sampleconverter(dropEST.combined.filtered$orig.ident))
SVF_frame <- unique(dropEST.combined.filtered@meta.data[,c("orig.ident","samplename")])


dropEST.combined.filtered <- NormalizeData(object = dropEST.combined.filtered, normalization.method = "LogNormalize",scale.factor = 10000)
dropEST.combined.filtered <- FindVariableFeatures(object = dropEST.combined.filtered)
dropEST.combined.filtered <- ScaleData(object = dropEST.combined.filtered, features = rownames(x = dropEST.combined.filtered))
dropEST.combined.filtered <- RunPCA(object = dropEST.combined.filtered, verbose = T, npcs = 75, ndims.print = 1:5, nfeatures.print = 10)
#dropEST.combined.filtered  <- CellCycleScoring(dropEST.combined.filtered , s.features = firstup(cc.genes$s.genes), g2m.features = firstup(cc.genes$g2m.genes), set.ident = TRUE)
#head(dropEST.combined.filtered[[]])
#result <- as.data.frame(table(dropEST.combined.filtered$Phase,dropEST.combined.filtered$celltype,dropEST.combined.filtered$data.diseasestatus))
#detach(package:plyr)    
#test2 <- result %>%
#     group_by(Var2, Var3) %>%
#     mutate(freq = Freq / sum(Freq))

#RidgePlot(dropEST.combined.filtered , features = c( "Mki67"), ncol = 2, assay = "RNA", y.max)

#DimPlot(dropEST.combined.filtered,group.by = "data.diseasestatus",reduction = "tsne",label = F, cols = c("Control"="darkslategray","Fruc" = "cyan3","HFHS" = "red"))+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=24))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("")
SVF_Fig1 <- DimPlot(dropEST.combined.filtered,group.by = "data.diseasestatus",reduction = "tsne",label = F, 
                    cols = c("Control"="darkslategray","Fruc" = "cyan3","HFHS" = "hotpink"))+ 
  theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=3)))+ggtitle("")
#DimPlot(dropEST.combined.filtered,group.by = "samplename",reduction = "tsne",label = F)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=24))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("")
SVF_Sup1 <- DimPlot(dropEST.combined.filtered,group.by = "samplename",reduction = "tsne",label = F)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=4)))+ggtitle("")
SVF_Sup6 <- DimPlot(dropEST.combined.filtered,group.by = "celltype",reduction = "tsne",label = T, repel = T, label.size = 8)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=8)))+ggtitle("")

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "samplename")
SVF_Sup2 <- VlnPlot(dropEST.combined.filtered,"nFeature_RNA", pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Gene count")+ theme(legend.position = "none")
SVF_Sup3 <- VlnPlot(dropEST.combined.filtered,"nCount_RNA", pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("UMI count")+ theme(legend.position = "none")
SVF_Sup4 <- VlnPlot(dropEST.combined.filtered,"percent.mito", pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Mitochondrial count percent")+ theme(legend.position = "none")
SVF_Sup5 <- VlnPlot(dropEST.combined.filtered,"percent.ribo", pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Robisomal count percent")+ theme(legend.position = "none")
SVF_Sup7 <- VlnPlot(dropEST.combined.filtered,"percent.pred", pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Robisomal count percent")+ theme(legend.position = "none")


genesOrder = c("Pdgfra","Cd34","Pi16","Hsd11b1","Col4a2","Lyz2","Ccl6","Ccl8","Il1b","Mrc1","Tnfrsf1b","Agt","F3","Esam","Msln","Skap1","Top2a","Ms4a1","Mcpt4","Ccl21a")
#genesOrder = c("Cd68","Cd163","Mrc1","Cd80","Cd86","Fcgr1","Fcgr2b")
#genesOrder <- c("Lpxn", "Dhrs3", "Mical1", "Dnmt3a", "Jun", "Gab1", "P2ry1")
#genesOrder <- c("Il12a","Il12b","Irf9", "Irf7", "Ifi35", "Ifnar2", "Isg20", "Ifit2","Stat1")
#genesOrder <- c("Tgfb1","Il10","Ccl8")

which(genesOrder %in% rownames(dropEST.combined.filtered))

#cellTypeOrder = c("APC_Pi16","APC_Hsd11b1","APC_Agt","Macrophage_Il1b","Endothelial cells","Mesothelial","T cells","Macrophage_Pf4","Macrophage_Top2a","B cells","Mast cells","Unknown")
cellTypeOrder = c("APC_Pi16","APC_Hsd11b1","APC_Agt","M1 Macrophage","Endothelial cells","Mesothelial","T cells","M2 Macrophage","Macrophage_Top2a","B cells","Mast cells","Unknown")

dropEST.combined.filtered$celltype <- factor(dropEST.combined.filtered$celltype, levels = cellTypeOrder)


dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered,value = "celltype")
DefaultAssay(dropEST.combined.filtered) <- "RNA"
dropEST.combined.filtered[['integrated']] <- NULL
SVF_Dot <- DotPlot(dropEST.combined.filtered,features = genesOrder,scale.by = "size", col.min = 0)+ylab("")+xlab("")+ theme(axis.text.x = element_text(size = 24, angle = 45, hjust = 1, face = "italic"),
                                                                                                                                           axis.text.y = element_text(size = 24),axis.title = element_text(size = 24),
                                                                                                                                           legend.text=element_text(size=20))+ guides(colour = guide_legend(override.aes = list(size=5)))
pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/SupFig2_SVF.pdf",width = 26, height = 16)
grid.arrange(SVF_Sup1,SVF_Sup2,SVF_Sup3,SVF_Sup4,SVF_Sup5,SVF_Sup6,SVF_Dot,layout_matrix = rbind(c(2,3,1,1),c(4,5,1,1),c(6,6,7,7),c(6,6,7,7)))
dev.off()

#Liver
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Combined_allliver.rda")
#dropEST.combined.filtered[["celltype"]] = plyr::revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "Sinusoidal endothelial cells","1" = "Periportal hepatocytes","2" = "Periportal hepatocytes","3" = "Kupffer", "4" = "Pericentral hepatocytes","5" = "Periportal hepatocytes", "6" = "NKT cells",
#                                                    "7" = "Periportal hepatocytes","8" = "Classical dendritic cells","9" = "B cells","10" = "Hepatic stellate cells","11" = "Cholangiocytes","12" = "Dividing cells","13" = "Plasmacytoid dendritic cells"))


load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/liver.rda")
DefaultAssay(dropEST.combined.filtered) <- "RNA"
#dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set12.1"] <- "Set11.1"
#dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set16.1"] <- "Set15.1"
#dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set22.1" & dropEST.combined.filtered$tissue %in% "HEP"] <- "Set21.1"
#dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set26.3" & dropEST.combined.filtered$tissue %in% "NPC"] <- "Set27.1"
#dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set24.2"] <- "Set28.2"
#dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set23.3"] <- "Set24.3"

#dropEST.combined.filtered[["samplename"]] <- paste0(dropEST.combined.filtered$tissue,"_",dropEST.combined.filtered$data.diseasestatus,"_",dropEST.combined.filtered$batch,"_", sampleconverter(dropEST.combined.filtered$orig.ident))
Liver_frame <- unique(dropEST.combined.filtered@meta.data[,c("orig.ident","samplename")])

#DimPlot(dropEST.combined.filtered,group.by = "data.diseasestatus",reduction = "tsne",label = F, cols = c("Control"="darkslategray","Fruc" = "cyan3","HFHS" = "red"))+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=24))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("")
Liver_Fig1 <- DimPlot(dropEST.combined.filtered,group.by = "data.diseasestatus",reduction = "tsne",label = F, 
                      cols = c("Control"="darkslategray","Fruc" = "cyan3","HFHS" = "hotpink"))+ 
  theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=3)))+ggtitle("")
#DimPlot(dropEST.combined.filtered,group.by = "samplename",reduction = "tsne",label = F)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=24))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("")
Liver_Sup1 <- DimPlot(dropEST.combined.filtered,group.by = "samplename",reduction = "tsne",label = F)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=10))+ guides(colour = guide_legend(override.aes = list(size=4)))+ggtitle("")
Liver_Sup6 <- DimPlot(dropEST.combined.filtered,group.by = "celltype",reduction = "tsne",label = T, repel = T, label.size = 8)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=8)))+ggtitle("")

FeaturePlot(dropEST.combined.filtered,"nCount_RNA",reduction = "tsne")

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "samplename")
Liver_Sup2 <- VlnPlot(dropEST.combined.filtered,"nFeature_RNA", pt.size = 0)+ theme(axis.text = element_text(size = 10),axis.title = element_text(size = 8),legend.text=element_text(size=10))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Gene count")+ theme(legend.position = "none",plot.margin = unit(c(0,0.2,0,0.8), "cm"))
Liver_Sup3 <- VlnPlot(dropEST.combined.filtered,"nCount_RNA", pt.size = 0)+ theme(axis.text = element_text(size = 10),axis.title = element_text(size = 12),legend.text=element_text(size=10))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("UMI count")+ theme(legend.position = "none")
Liver_Sup4 <- VlnPlot(dropEST.combined.filtered,"percent.mito", pt.size = 0)+ theme(axis.text = element_text(size = 10),axis.title = element_text(size = 8),legend.text=element_text(size=10))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Mitochondrial count percent")+ theme(legend.position = "none",plot.margin = unit(c(0,0.2,0,0.8), "cm"))
Liver_Sup5 <- VlnPlot(dropEST.combined.filtered,"percent.ribo", pt.size = 0)+ theme(axis.text = element_text(size = 10),axis.title = element_text(size = 12),legend.text=element_text(size=10))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Robisomal count percent")+ theme(legend.position = "none")
Liver_Sup7 <- VlnPlot(dropEST.combined.filtered,"percent.pred", pt.size = 0)+ theme(axis.text = element_text(size = 10),axis.title = element_text(size = 12),legend.text=element_text(size=10))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Robisomal count percent")+ theme(legend.position = "none")


# cellTypesToRemove = c("Neurons","Cluster16","Cluster23")
genesOrder = c("Bmp2","Clec4f","Alb","Fabp1","Cyp2f2","Cyp2e1","Nkg7","Itgax","H2-Ab1","Ms4a1","Reln","Spp1","Top2a","Siglech","Ccr9")

cellTypeOrder = c("Sinusoidal endothelial cells","Kupffer","Periportal hepatocytes","Pericentral hepatocytes","NKT cells","Classical dendritic cells","B cells","Hepatic stellate cells","Cholangiocytes","Dividing cells","Plasmacytoid dendritic cells")
dropEST.combined.filtered$celltype <- factor(dropEST.combined.filtered$celltype, levels = cellTypeOrder)
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered,value = "celltype")
Liver_Dot <- DotPlot(dropEST.combined.filtered, features = genesOrder,assay = "RNA",scale.by = "size", col.min = 0)+ylab("")+xlab("")+ theme(axis.text.x = element_text(size = 24, angle = 45, hjust = 1, face = "italic"),
                                                                                                                                             axis.text.y = element_text(size = 24),axis.title = element_text(size = 24),
                                                                                                                                             legend.text=element_text(size=20))+ guides(colour = guide_legend(override.aes = list(size=5)))

pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/SupFig3_Liver.pdf",width = 30, height = 20)
grid.arrange(Liver_Sup1,Liver_Sup2,Liver_Sup3,Liver_Sup4,Liver_Sup5,Liver_Sup6,Liver_Dot,layout_matrix = rbind(c(2,3,1,1),c(4,5,1,1),c(6,6,7,7),c(6,6,7,7)))
dev.off()

#Hyp all
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Hyp_combined_Set1_final.rda")
#dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set11.3"] <- "Set12.3"
#dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set13.3"] <- "Set14.3"

load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/Hyp.rda")
DefaultAssay(dropEST.combined.filtered) <- "RNA"
#dropEST.combined.filtered[["samplename"]] <- paste0(dropEST.combined.filtered$data.diseasestatus,"_", sampleconverter(dropEST.combined.filtered$orig.ident))

Hyp_frame <- unique(dropEST.combined.filtered@meta.data[,c("orig.ident","samplename")])

#DimPlot(dropEST.combined.filtered,group.by = "data.diseasestatus",reduction = "tsne",label = F, cols = c("Control"="darkslategray","Fruc" = "cyan3","HFHS" = "red"))+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=24))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("")
Hyp_Fig1 <- DimPlot(dropEST.combined.filtered,group.by = "data.diseasestatus",reduction = "tsne",label = F,
                    cols = c("Control"="darkslategray","Fruc" = "cyan3","HFHS" = "hotpink"))+ 
  theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=3)))+ggtitle("")
#DimPlot(dropEST.combined.filtered,group.by = "samplename",reduction = "tsne",label = F)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=24))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("")
Hyp_Sup1 <- DimPlot(dropEST.combined.filtered,group.by = "samplename",reduction = "tsne",label = F)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=4)))+ggtitle("")
Hyp_Sup6 <-DimPlot(dropEST.combined.filtered,group.by = "celltype",reduction = "tsne",label = T, repel = T, label.size = 6)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=8)))+ggtitle("")

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "samplename")
Hyp_Sup2 <- VlnPlot(dropEST.combined.filtered,"nFeature_RNA", pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Gene count")+ theme(legend.position = "none")
Hyp_Sup3 <- VlnPlot(dropEST.combined.filtered,"nCount_RNA", pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("UMI count")+ theme(legend.position = "none")
Hyp_Sup4 <- VlnPlot(dropEST.combined.filtered,"percent.mito", pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Mitochondrial count percent")+ theme(legend.position = "none")
Hyp_Sup5 <- VlnPlot(dropEST.combined.filtered,"percent.ribo", pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Robisomal count percent")+ theme(legend.position = "none")
Hyp_Sup7 <- VlnPlot(dropEST.combined.filtered,"percent.pred", pt.size = 0)+ theme(axis.text = element_text(size = 14),axis.title = element_text(size = 16),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("") + xlab("") + ylab("Robisomal count percent")+ theme(legend.position = "none")



# cellTypesToRemove = c("Neurons","Cluster16","Cluster23")
genesOrder = c("Agt","Mobp","Cx3cr1","Nfkbiz","Flt1","Rax","Acta2","Snap25",
               "Syt1","Slc17a6","Pdgfrb","Kcnj8","Col1a1","Cck","Cyp2f2","Mrc1","Fyn",
               "Slc32a1","Kdr","Pdgfra","Ccdc153")

cellTypeOrder = c("Astrocytes","Myelinating oligodendrocytes","Microglial cells_Nfkbiz",
                  "Pericytes","Tanycytes","Muralpericytes_Acta2","Glutamatergic neurons","Pericyes_Kcnj8",
                  "Vascular Leptomenigeal cells","GABAergic neurons","Pars tuberalis",
                  "Macrophages","Newly formed oligodendrocytes","Endothelial cells",
                  "Oligodendrocytes precursor","General neurons","Ependymal cells",
                  "Microglial cells")
dropEST.combined.filtered$celltype <- factor(dropEST.combined.filtered$celltype, levels = cellTypeOrder)
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered,value = "celltype")

Hyp_Dot <- DotPlot(dropEST.combined.filtered, features = genesOrder,assay = "RNA",scale.by = "size", col.min = 0)+ylab("")+xlab("")+ theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1, face = "italic"),
                                                                                                                                           axis.text.y = element_text(size = 16),axis.title = element_text(size = 16),
                                                                                                                                           legend.text=element_text(size=12))+ guides(colour = guide_legend(override.aes = list(size=5)))
Hyp_Dot




#Hyp neurons
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Neuroncelltype.rda")
#dropEST.combined.filtered <- CellTypeSubset
#dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set11.3"] <- "Set12.3"
#dropEST.combined.filtered$orig.ident[dropEST.combined.filtered$orig.ident %in% "Set13.3"] <- "Set14.3"
#
#dropEST.combined.filtered[["celltype"]] =  plyr::revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.1.5,
#                                                  c("0" = "General neurons","1" = "General neurons","2" = "General neurons","3" = "Trh_GABAergic", 
#                                                    "4" = "General neurons","5" = "Grp_Glutamatergic","6" = "Prok2_GABAergic","7" = "Agrp_GABAergic",
  #                                                  "8" = "Pomc_Glutamatergic","9" = "Sst_Bcl11b_GABAergic","10" = "Prlr_Glutamatergic","11" = "Sim1_Glutamatergic","12" = "Ghrh_GABAergic",
  #                                                  "13" = "Crabp1_GABAergic","14" = "Kl_GABAergic","15" = "Avp_Glutamatergic","16" = "Lhx1_GABAergic",
  #                                                  "17" = "Lhx6_GABAergic","18" = "Prph_Histaminergic","19" = "Cx3cr1_GABAergic","20" = "Ndnf_Glutamatergic"))
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/Hyp_neurons.rda")
DefaultAssay(dropEST.combined.filtered) <- "RNA"
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered,value="batch")
dropEST.combined.filtered = subset(dropEST.combined.filtered,idents="Set1")
#dropEST.combined.filtered[["samplename"]] <- paste0(dropEST.combined.filtered$data.diseasestatus,"_", sampleconverter(dropEST.combined.filtered$orig.ident))

DimPlot(dropEST.combined.filtered,group.by = "data.diseasestatus",reduction = "tsne",label = F, cols = c("Control"="darkslategray","Fruc" = "cyan3","HFHS" = "red"))+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=24))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("")
Hypneuron_Fig1 <- DimPlot(dropEST.combined.filtered,group.by = "data.diseasestatus",reduction = "tsne",label = F, cols = c("Control"="darkslategray","Fruc" = "cyan3","HFHS" = "red"))+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=24))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("")
DimPlot(dropEST.combined.filtered,group.by = "samplename",reduction = "tsne",label = F)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=24))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("")
Hypneuron_Sup1 <- DimPlot(dropEST.combined.filtered,group.by = "samplename",reduction = "tsne",label = F)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=24))+ guides(colour = guide_legend(override.aes = list(size=5)))+ggtitle("")
Hypneuron_Sup6 <- DimPlot(dropEST.combined.filtered,group.by = "celltype",reduction = "tsne",label = T, repel = T, label.size = 8)+ theme(axis.text = element_text(size = 24),axis.title = element_text(size = 28),legend.text=element_text(size=14))+ guides(colour = guide_legend(override.aes = list(size=8)))+ggtitle("")

genesOrder <- c("Slc17a6","Slc32a1","Hdc","Trh","Grp","Prok2","Sst","Agrp","Pomc","Bcl11b","Prlr","Sim1","Ghrh","Crabp1","Kl","Avp","Lhx1","Lhx6","Prph","Cx3cr1","Ndnf")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered,value = "celltype")
Hypneuron_Dot <-  DotPlot(dropEST.combined.filtered, features = genesOrder,assay = "RNA",scale.by = "size", col.min = 0)+ylab("")+xlab("")+ theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1, face = "italic"),
                                                                                                                                                  axis.text.y = element_text(size = 16),axis.title = element_text(size = 16),
                                                                                                                                                  legend.text=element_text(size=12))+ guides(colour = guide_legend(override.aes = list(size=5)))
#pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Combined_figure1_tSNE.pdf",width = 10, height = 40)
#grid.arrange(SI_Fig1, SVF_Fig1, Liver_Fig1, Hyp_Fig1, Hypneuron_Fig1,ncol = 1)
#dev.off()
#pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Combined_figure1_Dot.pdf",width =16, height = 45)
#grid.arrange(SI_Dot,SVF_Dot, Liver_Dot, Hyp_Dot,Hypneuron_Dot,ncol = 1)
#dev.off()
pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/SupFig4_Hyp.pdf",width = 25, height = 25)
grid.arrange(Hyp_Sup1,Hyp_Sup2,Hyp_Sup3,Hyp_Sup4,Hyp_Sup5,Hyp_Sup6,Hyp_Dot,Hypneuron_Sup6,Hypneuron_Dot,layout_matrix = rbind(c(2,3,1,1),c(4,5,1,1),c(6,6,7,7),c(6,6,7,7),c(8,8,9,9),c(8,8,9,9)))
dev.off()

#png("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Combined_figure1_combined.png",width = 3000, height = 3000, pointsize = 8, res = 150)
png("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Combined_figure1_combined.png",width = 2000, height = 2000, pointsize = 8, res = 150)
#pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Combined_figure1_combined.pdf",width = 15, height = 15)
grid.arrange(SI_Fig1,SVF_Fig1,Liver_Fig1,Hyp_Fig1,nrow = 2)
dev.off()
save(SI_Fig1,SVF_Fig1,Liver_Fig1,Hyp_Fig1, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/tSNEFig2.rda")

pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Combined_figure1_combined.pdf",width = 23, height = 15)
grid.arrange(SI_Fig1,SVF_Fig1,Liver_Fig1,Hyp_Fig1,layout_matrix = rbind(c(1,2,3),c(4,NA,NA)))
dev.off()



#pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/New_Figure3.pdf",width =60, height = 15)
#grid.arrange(SI_Sup6,SVF_Sup6,Liver_Sup6,Hyp_Sup6,Hypneuron_Sup6,SI_Dot,SVF_Dot, Liver_Dot, Hyp_Dot,Hypneuron_Dot,nrow = 2)
pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/New_Figure3.pdf",width =40, height = 40)
grid.arrange(SI_Sup6,SVF_Sup6,Liver_Sup6,Hyp_Sup6,Hypneuron_Sup6,SI_Dot,SVF_Dot, Liver_Dot, Hyp_Dot,Hypneuron_Dot,layout_matrix = cbind(c(1,1,7,7,NA,NA,NA,NA),c(1,1,7,7,5,5,10,10),
                                                                                                                                        c(2,2,8,8,5,5,10,10),c(2,2,8,8,6,6,11,11),
                                                                                                                                        c(3,3,9,9,6,6,11,11), c(3,3,9,9,NA,NA,NA,NA)))
dev.off()

pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/New_Figure3_design.pdf",width =28, height = 50)
grid.arrange(SI_Sup6,SVF_Sup6,Liver_Sup6,Hyp_Sup6,Hypneuron_Sup6,SI_Dot,SVF_Dot, Liver_Dot, Hyp_Dot,Hypneuron_Dot,layout_matrix = cbind(c(1,2,3,4,5),c(6,7,8,9,10)))
dev.off()


pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/New_Figure.pdf",width =40, height = 40)
grid.arrange(SI_Sup6,SVF_Sup6,Liver_Sup6,Hyp_Sup6,Hypneuron_Sup6,SI_Dot,SVF_Dot, Liver_Dot, Hyp_Dot,Hypneuron_Dot,layout_matrix = cbind(c(1,1,7,7,NA,NA,NA,NA),c(1,1,7,7,5,5,10,10),
                                                                                                                                        c(2,2,8,8,5,5,10,10),c(2,2,8,8,6,6,11,11),
                                                                                                                                        c(3,3,9,9,6,6,11,11), c(3,3,9,9,NA,NA,NA,NA)))
dev.off()

pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/New_SupFig1.pdf",width =20, height = 20)
#grid.arrange(SI_Sup2,SI_Sup1,SVF_Sup2,SVF_Sup1,Liver_Sup2,Liver_Sup1,Hyp_Sup2,Hyp_Sup1,nrow = 4)
grid.arrange(SI_Sup2,SI_Sup3,SI_Sup4,SI_Sup5,SI_Sup1,SVF_Sup2,SVF_Sup3,SVF_Sup4,SVF_Sup5,SVF_Sup1,Liver_Sup2,Liver_Sup3,Liver_Sup4,Liver_Sup5,Liver_Sup1,Hyp_Sup2,Hyp_Sup3,Hyp_Sup4,Hyp_Sup5,Hyp_Sup1,
             layout_matrix = rbind(c(1,2,5,5),c(3,4,5,5),c(6,7,10,10),c(8,9,10,10),
                                   c(11,12,15,15),c(13,14,15,15),c(16,17,20,20),c(18,19,20,20)))

dev.off()

save(intestine_frame, Hyp_frame,Liver_frame,SVF_frame, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/sample_name_frame_for_upload.rda")

