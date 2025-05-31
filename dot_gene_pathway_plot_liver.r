#######Monocle  version!!!!!

library(metap)
library(Seurat)
library(dplyr)
library(plyr)
library(heatmap3)
library(RColorBrewer)
library(gplots)
library(ggplot2)
#Cross DEG analysis

#rm(list=ls())
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined")


#rm(list=ls())
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined")

rm(list=ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER")
#rm(list=ls())
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_NPC")

#rm(list=ls())
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")


#rm(list=ls())
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "lowUMI1","1" = "lowUMI2","3" = "Grp-PVN","2" = "Lateral-Dopamine", "4" = "Pmoc-arcuate","7" = "sst-neuron", "5" = "Npy-argp-arcuate",
#                                                    "6" = "PVN","8" = "Ghrh-arcuate","9" = "Prlr-preoptic","10" = "Supraoptic","11" = "Unk11","12" = "prph-Posterior"))

#rm(list=ls())
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1") #29905 -> 9900, 3 sample (one seems failed), 1485 in min gene > 300


#Set <- "Set1"
#Set <- "Set2" #SI
firstup <- function(x) {
  result <- NULL
  for(i in 1:length(x)){
    current <- x[i]
    #current <- tolower(current)
    if(startsWith(current,"mt-")){
      substr(current, 4, 4) <- toupper(substr(current, 4, 4))
    }else{
      substr(current, 1, 1) <- toupper(substr(current, 1, 1))
    }
    result[i] <- current
  }
  result
}

top5genelist <- function(DEGlist){
  finalgenelist <- NULL
  for(i in 1:length(DEGlist)){
    finalgenelist <- c(finalgenelist, DEGlist[[i]][1:3])
  }
  finalgenelist <- unique(finalgenelist[!is.na(finalgenelist)])
  
  return(finalgenelist)
}
#load( "Humaninefold.rda")
load( "HumaninefoldFruc_monocle.rda")
#load( "Humaninefold_subset.rda")

rownames(foldchangemerge_final) <- foldchangemerge_final$genename
foldchangemerge_final <- foldchangemerge_final[,!colnames(foldchangemerge_final) %in% "genename"]
foldchangemerge_final <- as.matrix(foldchangemerge_final)
foldchangemerge_final[is.na(foldchangemerge_final)] <- 0
#foldchangemerge_final <- foldchangemerge_final[,grep(Set,colnames(foldchangemerge_final))]
colnames(foldchangemerge_final) <- paste0(colnames(foldchangemerge_final),"_Fruc")


rownames(pvaluemerge_final) <- pvaluemerge_final$genename
pvaluemerge_final <- pvaluemerge_final[,!colnames(pvaluemerge_final) %in% "genename"]
pvaluemerge_final <- as.matrix(pvaluemerge_final)
pvaluemerge_final[is.na(pvaluemerge_final)] <- 1
pvaluemerge_final <- -log10(pvaluemerge_final)
pvaluemerge_final <- as.data.frame(pvaluemerge_final)
#pvaluemerge_final <- pvaluemerge_final[,grep(Set,colnames(pvaluemerge_final))]
colnames(pvaluemerge_final) <- paste0(colnames(pvaluemerge_final),"_Fruc")

topgenes_Fruc <- top5genelist(DEGlist)
foldchangemerge_final_Fruc <- foldchangemerge_final
pvaluemerge_final_Fruc <- pvaluemerge_final



load( "HumaninefoldHFHS_monocle.rda")
#load( "HumaninefoldHFHS_subset.rda")
rownames(foldchangemerge_final) <- foldchangemerge_final$genename
foldchangemerge_final <- foldchangemerge_final[,!colnames(foldchangemerge_final) %in% "genename"]
foldchangemerge_final <- as.matrix(foldchangemerge_final)
foldchangemerge_final[is.na(foldchangemerge_final)] <- 0
#foldchangemerge_final <- foldchangemerge_final[,grep(Set,colnames(foldchangemerge_final))]
colnames(foldchangemerge_final) <- paste0(colnames(foldchangemerge_final),"_HFHS")

rownames(pvaluemerge_final) <- pvaluemerge_final$genename
pvaluemerge_final <- pvaluemerge_final[,!colnames(pvaluemerge_final) %in% "genename"]
pvaluemerge_final <- as.matrix(pvaluemerge_final)
pvaluemerge_final[is.na(pvaluemerge_final)] <- 1
pvaluemerge_final <- -log10(pvaluemerge_final)
pvaluemerge_final <- as.data.frame(pvaluemerge_final)
#pvaluemerge_final <- pvaluemerge_final[,grep(Set,colnames(pvaluemerge_final))]
colnames(pvaluemerge_final) <- paste0(colnames(pvaluemerge_final),"_HFHS")

topgenes_HFHS <- top5genelist(DEGlist)
foldchangemerge_final_HFHS <- foldchangemerge_final
pvaluemerge_final_HFHS <- pvaluemerge_final





topgenes_shared <- union(topgenes_Fruc,topgenes_HFHS)
topgenes_shared <- intersect(topgenes_shared, rownames(foldchangemerge_final_Fruc))
topgenes_shared <- intersect(topgenes_shared, rownames(foldchangemerge_final_HFHS))

save(foldchangemerge_final_Fruc, foldchangemerge_final_HFHS, pvaluemerge_final_Fruc, pvaluemerge_final_HFHS, file = "Dotgenepackage.rda")
foldchangemerge_final <- cbind.data.frame(foldchangemerge_final_Fruc[topgenes_shared,],foldchangemerge_final_HFHS[topgenes_shared,])
foldchangemerge_final$gene <- rownames(foldchangemerge_final )
pvaluemerge_final <- cbind.data.frame(pvaluemerge_final_Fruc[topgenes_shared,],pvaluemerge_final_HFHS[topgenes_shared,])
pvaluemerge_final$gene <- rownames(pvaluemerge_final)

foldchangemerge_final_melt <- tidyr::gather(foldchangemerge_final, name, fold,-gene)
pvaluemerge_final_melt <-  tidyr::gather(pvaluemerge_final, name, FDR,-gene)
merged_final <- inner_join(foldchangemerge_final_melt,pvaluemerge_final_melt)
merged_final$FDR[merged_final$FDR > 3] <- 3
colnames(merged_final)[4] <- "FDR"

merged_final$fold[merged_final$fold > 5] <- 5
merged_final$fold[merged_final$fold < -5] <- -5
factorlevels <- colnames(foldchangemerge_final)
merged_final$gene <- factor(merged_final$gene,levels = rev(topgenes_shared))
merged_final$name <-  factor(merged_final$name ,levels = factorlevels)
merged_final <- merged_final[!merged_final$gene %in% "Malat1",]
p1 <- ggplot(merged_final, aes(y = gene,
                               x = name)) +        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = fold, 
                 size =FDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-5,5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+   ## color of the corresponding aes          ## to tune the size of circles
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ guides(colour=guide_legend(title="-log(fold change)"),
                                                                              size =guide_legend(title="-log(FDR)"))

p1



pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Liver/genedotplot.pdf",width = 10, height = 8)
print(p1)
dev.off()

newlevels <- unique(gsub("_Fruc|_HFHS","",gsub("oligodendrocytes precursor","Oligodendrocytes precursor",colnames(foldchangemerge_final))))
merged_final$treatment <- sapply(strsplit(as.character(merged_final$name),"_"), `[[`, 2)
merged_final$name <- factor(sapply(strsplit(as.character(merged_final$name),"_"), `[[`, 1),levels = newlevels)
p1 <- ggplot(merged_final, aes(y = gene,x = name)) +facet_grid(.~treatment, scales = "free", space = "free")+        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = fold, 
                 size = FDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-5,5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+   ## color of the corresponding aes          ## to tune the size of circles
  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                   strip.text.x = element_text(size = 20),
                   plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))+ guides(colour=guide_legend(title="-log(fold change)"),
                                                                              size =guide_legend(title="-log(FDR)"))

p1
#pdf("genedotplot_subset.pdf",width = 10, height = 8)
#print(p1)
#dev.off()

#cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(40)
#heatmap.2(as.matrix(foldchangemerge_final), scale = "none",col = cols, breaks = seq(-4,4,0.2),  trace = "none",density.info = "none", Rowv = F, Colv = F, margins = c(6, 10))
#heatmap3(as.matrix(foldchangemerge_final))


#building heatmap on top t genes across different signle cell types (top 5 in each cell)


#pvaluemerge_final2 <- pvaluemerge_final
#pvaluemerge_final2$score <- rowMeans(pvaluemerge_final > -log10(0.05))
#pvaluemerge_final2$score_set1 <- rowMeans(pvaluemerge_final[,grep("Set1",colnames(pvaluemerge_final))] > -log10(0.05))
#pvaluemerge_final2$score_set2 <- rowMeans(pvaluemerge_final[,grep("Set2",colnames(pvaluemerge_final))] > -log10(0.05))

#genelist <- top5genelist(foldchangemerge_final,pvaluemerge_final)



setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER")
#rm(list=ls())
#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_NPC")


#load("./pathwayinfo_Fruc.rda")
load( "HumaninefoldFruc_monocle.rda")
invalid <- NULL
for(i in 1:nrow(result_final)){
  currentgenes <- DEGlist[[result_final$celltype[i]]]
  if(length(currentgenes) < 20){
    invalid <- c(invalid, i)
  }
}
#load("./pathwayinfosubset_Fruc.rda")

#result_final$treatment <- "Fruc"
#result_final <- result_final[result_final$Set %in% Set,]
selectedpathways <- c("Translational initiation","Translation","Negative regulation of apoptotic process","Apoptosis","Triglyceride homeostasis",
                      "Fatty acid metabolism","Chylomicron assembly","Cholesterol homeostasis","High−density lipoprotein particle remodeling","Glycogen biosynthetic process","Gluconeogenesis",
                      "Carbon metabolism","PPAR signaling pathway","Peroxisome","Cellular protein metabolic process","Response to cytokine","Antigen processing and presentation")
result_final$celltype <- paste0(result_final$celltype,"_Fruc")
result_Fruc_original <- result_final
result_final$FDR[result_final$FDR > 3 ] <- 3
result_final$Term <- gsub(" (","",result_final$Term,fixed = T)
result_final$Term <- gsub(")","",result_final$Term,fixed = T)
result_final$Term <- gsub("GO:[0-9]+| hsa[0-9]+","",result_final$Term)
result_final$Term <- gsub(" Homo sapiens","",result_final$Term)
result_final$Term <- firstup(result_final$Term)
#result_final
result_Fruc <- result_final[result_final$Term %in% selectedpathways,]


#load("./pathwayinfo_HFHS.rda")
load( "HumaninefoldHFHS_monocle.rda")
#Filter out pathways for those with less than 10 genes
invalid <- NULL
for(i in 1:nrow(result_final)){
  currentgenes <- DEGlist[[result_final$celltype[i]]]
  if(length(currentgenes) < 20){
    invalid <- c(invalid, i)
  }
}
#load("pathwayinfo_subset_HFHS.rda")

#result_final$treatment <- "HFHS"
#result_final <- result_final[result_final$Set %in% Set,]
result_final$celltype <- paste0(result_final$celltype,"_HFHS")

result_HFHS_original <- result_final
result_final$FDR[result_final$FDR > 3 ] <- 3
result_final$Term <- gsub(" (","",result_final$Term,fixed = T)
result_final$Term <- gsub(")","",result_final$Term,fixed = T)
result_final$Term <- gsub("GO:[0-9]+| hsa[0-9]+","",result_final$Term)
result_final$Term <- gsub(" Homo sapiens","",result_final$Term)
result_final$Term <- firstup(result_final$Term)

#result_final
result_HFHS <- result_final[result_final$Term %in% selectedpathways,]

#result_HFHS <- result_final %>% dplyr::group_by(celltype) %>% dplyr::slice_min(P.value, n = 5)

result_final <- rbind.data.frame(result_Fruc,result_HFHS)
result_final$celltype <- factor(result_final$celltype,levels = factorlevels)


newlevels <- unique(gsub("_Fruc|_HFHS","",gsub("oligodendrocytes precursor","Oligodendrocytes precursor",colnames(foldchangemerge_final))))
result_final$treatment <- sapply(strsplit(as.character(result_final$celltype),"_"), `[[`, 2)
result_final$name <- factor(sapply(strsplit(as.character(result_final$celltype),"_"), `[[`, 1),levels = newlevels)
result_final$Term <- factor(result_final$Term, levels = rev(selectedpathways))
p2 <- ggplot(result_final, aes(y = Term,x = name)) +facet_grid(.~treatment, scales = "free", space = "free")+        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = fold, 
                 size = FDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-3,3),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+   ## color of the corresponding aes          ## to tune the size of circles
  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                   strip.text.x = element_text(size = 20),
                   plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))+ guides(colour=guide_legend(title="-log(fold change)"),
                                                                       size =guide_legend(title="-log(FDR)"))
p2

result_original <- rbind.data.frame(result_Fruc_original,result_HFHS_original)
write.csv(result_original, file = "/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Liver/pathway_final_info.csv")
p2
pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Liver/pathwaydotplot.pdf",width = 10, height = 10)
print(p2)
dev.off()

result_original2 <- result_original[,c("Term","fold","FDR","Genes","celltype")]
colnames(result_original2) <- c("name","fold","FDR","genevector","celltype")
write.csv(result_original2, file = "/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Liver/pathway_info_filtered.csv")

p1_liver <- p1
p2_liver <- p2
save(p1_liver, p2_liver, file = "dotgeneplot.rda")


#pdf("Pathwaydotplot_subset.pdf",width = 10, height = 10)
#print(p2)
#dev.off()
#closer look of top pathway gene components across different cell types

#trajectory analysis of specific cell type of interest?

#GWAS enrichment for further subclusters

