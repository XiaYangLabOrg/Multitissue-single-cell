options(stringsAsFactors = F)
#Dot heatmap selected pathway detailed gene view
library(metap)
library(Seurat)
library(heatmap3)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(dplyr)
library(rlang)
library(magrittr)
library(grid)
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


#Picking up dot gene pathway map
#SVF pathway
pathwayname <- "negative regulation of apoptotic process"
celltype <- c("APC_Hsd11b1_Fruc","APC_Pi16_Fruc","APC_Agt_Fruc","APCHsd11b1_HFHS","APC_Pi16_HFHS","APC_Agt_HFHS")
Tissue = "SVF"
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined/Dotgenepackage.rda")
setwd(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",Tissue))
result_final <-  read.csv("pathwayinfo_combined.csv")
result_final$genesize <- sapply(strsplit(result_final$genevector,";"),length)

result_final <- result_final[result_final$genesize  > 3,]
selectedpathways <- as.character(result_final$pathwayname[result_final$selected %in% "*"])
selectedpathways2 <- read.csv("./pathwayinfo_combined_filtered_selected.csv")
result_final <- inner_join(result_final,selectedpathways2)

result_final$FDR[result_final$FDR > 3 ] <- 3
result_final$fold[result_final$fold > 5 ] <- 5
result_final$fold[result_final$fold < -5 ] <- -5
result_final$celltype <- factor(result_final$celltype, levels = unique(result_final$celltype))
result_final <- result_final[order(result_final$group,result_final$pathwayname),]

result_final$pathwayname <- gsub(" (","",result_final$pathwayname,fixed = T)
result_final$pathwayname <- gsub(")","",result_final$pathwayname,fixed = T)
result_final$pathwayname <- gsub("GO:[0-9]+| hsa[0-9]+","",result_final$pathwayname)
result_final$pathwayname <- gsub(" Homo sapiens","",result_final$pathwayname)
result_final$pathwayname <- factor(result_final$pathwayname , levels = unique(result_final$pathwayname))
filtered <- result_final[result_final$pathwayname %in% pathwayname,]
filtered <- filtered[filtered$celltype %in% celltype,]
currentgenesets <- unique(firstup(unlist(strsplit(filtered$genevector,";"))))

#plotting with selected genes and gene sets
topgenes_shared <- currentgenesets
foldchangemerge_final <- cbind.data.frame(foldchangemerge_final_Fruc[topgenes_shared,],foldchangemerge_final_HFHS[topgenes_shared,])
foldchangemerge_final$gene <- rownames(foldchangemerge_final )
pvaluemerge_final <- cbind.data.frame(pvaluemerge_final_Fruc[topgenes_shared,],pvaluemerge_final_HFHS[topgenes_shared,])
pvaluemerge_final$gene <- rownames(pvaluemerge_final)

foldchangemerge_final_melt <- tidyr::gather(foldchangemerge_final, name, fold,-gene)
pvaluemerge_final_melt <-  tidyr::gather(pvaluemerge_final, name, FDR,-gene)
merged_final <- inner_join(foldchangemerge_final_melt,pvaluemerge_final_melt)
merged_final$FDR[merged_final$FDR > 3] <- 3
colnames(merged_final)[4] <- "FDR"
merged_final <- merged_final[merged_final$name %in% celltype,]
merged_final$gene <- factor(merged_final$gene,levels = rev(topgenes_shared))
merged_final$name <- factor(merged_final$name,levels = colnames(foldchangemerge_final))

newlevels <- unique(gsub("_Fruc|_HFHS","",gsub("oligodendrocytes precursor","Oligodendrocytes precursor",colnames(foldchangemerge_final))))
merged_final$name <- as.character(merged_final$name)
merged_final$name <- gsub("_Fruc",":Fruc",merged_final$name)
merged_final$name <- gsub("_HFHS",":HFHS",merged_final$name)
merged_final$treatment <- sapply(strsplit(as.character(merged_final$name),":"), `[[`, 2)
merged_final$name <- factor(sapply(strsplit(as.character(merged_final$name),":"), `[[`, 1),levels = newlevels)
#p1_SVF <- ggplot(merged_final, aes(y = gene,x = name)) +facet_wrap(.~treatment)+        ## global aes
#  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
#  geom_point(aes(colour = fold, 
#                 size = FDR))  +    ## geom_point for circle illusion
#  scale_color_gradientn(limits = c(-5,5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+   ## color of the corresponding aes          ## to tune the size of circles
#  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
#                   axis.text.x = element_text( angle = 45,hjust = 1,size = 18),
#                   strip.text.x = element_text(size = 20),
#                   plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))+ guides(colour=guide_legend(title="-log(fold change)"),
#                                                                       size =guide_legend(title="-log(FDR)")) 
merged_final <- merged_final[!merged_final$gene %in% "Malat1",]
p1_SVF <- ggplot(merged_final, aes(y = name,x = gene)) +facet_grid(treatment~.)+        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = fold, 
                 size = FDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-5,5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+   ## color of the corresponding aes          ## to tune the size of circles
  theme_bw()+theme(axis.text.y = element_text( size = 16),
                   axis.text.x = element_text( angle = 45, hjust = 1, size = 18,face = "italic"),
                   strip.text.y  = element_text(size = 20),
                   plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))+ guides(colour=guide_legend(title="-log(fold change)"),size =guide_legend(title="-log(FDR)")) 



pdf(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",pathwayname,"_",Tissue,".pdf"),width = 18, height = 5.5)
print(p1_SVF)
dev.off()


#Liver pathway
pathwayname <- "triglyceride homeostasis"
celltype <- c("Pericentral hepatocytes_Fruc","Periportal hepatocytes_Fruc","Pericentral hepatocytes_HFHS","Periportal hepatocytes_HFHS")
Tissue = "Liver"
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Dotgenepackage.rda")
setwd(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",Tissue))
result_final <-  read.csv("pathwayinfo_combined.csv")
result_final$genesize <- sapply(strsplit(result_final$genevector,";"),length)

result_final <- result_final[result_final$genesize  > 3,]
selectedpathways <- as.character(result_final$pathwayname[result_final$selected %in% "*"])
selectedpathways2 <- read.csv("./pathwayinfo_combined_filtered_selected.csv")
result_final <- inner_join(result_final,selectedpathways2)

result_final$FDR[result_final$FDR > 3 ] <- 3
result_final$fold[result_final$fold > 3 ] <- 3
result_final$fold[result_final$fold < -3 ] <- -3
result_final$celltype <- factor(result_final$celltype, levels = unique(result_final$celltype))
result_final <- result_final[order(result_final$group,result_final$pathwayname),]

result_final$pathwayname <- gsub(" (","",result_final$pathwayname,fixed = T)
result_final$pathwayname <- gsub(")","",result_final$pathwayname,fixed = T)
result_final$pathwayname <- gsub("GO:[0-9]+| hsa[0-9]+","",result_final$pathwayname)
result_final$pathwayname <- gsub(" Homo sapiens","",result_final$pathwayname)
result_final$pathwayname <- factor(result_final$pathwayname , levels = unique(result_final$pathwayname))
filtered <- result_final[result_final$pathwayname %in% pathwayname,]
filtered <- filtered[filtered$celltype %in% celltype,]
currentgenesets <- unique(firstup(unlist(strsplit(filtered$genevector,";"))))

#plotting with selected genes and gene sets
topgenes_shared <- currentgenesets
topgenes_shared <- intersect(topgenes_shared, rownames(foldchangemerge_final_Fruc))
topgenes_shared <- intersect(topgenes_shared, rownames(foldchangemerge_final_HFHS))

foldchangemerge_final <- cbind.data.frame(foldchangemerge_final_Fruc[topgenes_shared,],foldchangemerge_final_HFHS[topgenes_shared,])
foldchangemerge_final$gene <- rownames(foldchangemerge_final )
pvaluemerge_final <- cbind.data.frame(pvaluemerge_final_Fruc[topgenes_shared,],pvaluemerge_final_HFHS[topgenes_shared,])
pvaluemerge_final$gene <- rownames(pvaluemerge_final)

foldchangemerge_final_melt <- tidyr::gather(foldchangemerge_final, name, fold,-gene)
pvaluemerge_final_melt <-  tidyr::gather(pvaluemerge_final, name, FDR,-gene)
merged_final <- inner_join(foldchangemerge_final_melt,pvaluemerge_final_melt)
merged_final$FDR[merged_final$FDR > 3] <- 3
colnames(merged_final)[4] <- "FDR"
merged_final <- merged_final[merged_final$name %in% celltype,]
merged_final$gene <- factor(merged_final$gene,levels = rev(topgenes_shared))
merged_final$name <- factor(merged_final$name,levels = colnames(foldchangemerge_final))

newlevels <- unique(gsub("_Fruc|_HFHS","",gsub("oligodendrocytes precursor","Oligodendrocytes precursor",colnames(foldchangemerge_final))))
merged_final$name <- as.character(merged_final$name)
merged_final$name <- gsub("_Fruc",":Fruc",merged_final$name)
merged_final$name <- gsub("_HFHS",":HFHS",merged_final$name)
merged_final$treatment <- sapply(strsplit(as.character(merged_final$name),":"), `[[`, 2)
merged_final$name <- factor(sapply(strsplit(as.character(merged_final$name),":"), `[[`, 1),levels = newlevels)
#p1_liver <- ggplot(merged_final, aes(y = gene,x = name)) +facet_wrap(.~treatment)+        ## global aes
#  geom_tile(fill = "white") +  xlab("")+ylab("gene")+       ## to get the rect filled
#  geom_point(aes(colour = fold, 
#                 size = FDR))  +    ## geom_point for circle illusion
#  scale_color_gradientn(limits = c(-1,1),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+   ## color of the corresponding aes          ## to tune the size of circles
#  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
#                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
#                   strip.text.x = element_text(size = 20),
#                   plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))+ guides(colour=guide_legend(title="-log(fold change)"),
#                                                                       size =guide_legend(title="-log(FDR)"))

p1_liver <- ggplot(merged_final, aes(y = name,x = gene)) +facet_grid(treatment~.)+        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = fold, 
                 size = FDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-1,1),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+   ## color of the corresponding aes          ## to tune the size of circles
  theme_bw()+theme(axis.text.y = element_text( size = 16),
                   axis.text.x = element_text( angle = 45, hjust = 1, size = 18,face = "italic"),
                   strip.text.y  = element_text(size = 20),
                   plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))+ guides(colour=guide_legend(title="-log(fold change)"),size =guide_legend(title="-log(FDR)")) 



pdf(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",pathwayname,"_",Tissue,".pdf"),width = 10, height = 5.5)
print(p1_liver)
dev.off()



#Hypothalamus pathway
pathwayname <- "Oxytocin signaling pathway"
celltype <- c("GABAergic neurons_Fruc","Glutamatergic neurons_Fruc","GABAergic neurons_HFHS","Glutamatergic neurons_HFHS")
Tissue = "Hyp"
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Dotgenepackage.rda")
setwd(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",Tissue))
result_final <-  read.csv("pathwayinfo_combined.csv")
result_final$genesize <- sapply(strsplit(result_final$genevector,";"),length)

result_final <- result_final[result_final$genesize  > 3,]
selectedpathways <- as.character(result_final$pathwayname[result_final$selected %in% "*"])
selectedpathways2 <- read.csv("./pathwayinfo_combined_filtered_selected.csv")
result_final <- inner_join(result_final,selectedpathways2)

result_final$FDR[result_final$FDR > 3 ] <- 3
result_final$fold[result_final$fold > 5 ] <- 5
result_final$fold[result_final$fold < -5 ] <- -5
result_final$celltype <- factor(result_final$celltype, levels = unique(result_final$celltype))
result_final <- result_final[order(result_final$group,result_final$pathwayname),]

result_final$pathwayname <- gsub(" (","",result_final$pathwayname,fixed = T)
result_final$pathwayname <- gsub(")","",result_final$pathwayname,fixed = T)
result_final$pathwayname <- gsub("GO:[0-9]+| hsa[0-9]+","",result_final$pathwayname)
result_final$pathwayname <- gsub(" Homo sapiens","",result_final$pathwayname)
result_final$pathwayname <- factor(result_final$pathwayname , levels = unique(result_final$pathwayname))
filtered <- result_final[result_final$pathwayname %in% pathwayname,]
filtered <- filtered[filtered$celltype %in% celltype,]
currentgenesets <- unique(firstup(unlist(strsplit(filtered$genevector,";"))))

#plotting with selected genes and gene sets
topgenes_shared <- currentgenesets
topgenes_shared <- intersect(topgenes_shared, rownames(foldchangemerge_final_Fruc))
topgenes_shared <- intersect(topgenes_shared, rownames(foldchangemerge_final_HFHS))

foldchangemerge_final <- cbind.data.frame(foldchangemerge_final_Fruc[topgenes_shared,],foldchangemerge_final_HFHS[topgenes_shared,])
foldchangemerge_final$gene <- rownames(foldchangemerge_final )
pvaluemerge_final <- cbind.data.frame(pvaluemerge_final_Fruc[topgenes_shared,],pvaluemerge_final_HFHS[topgenes_shared,])
pvaluemerge_final$gene <- rownames(pvaluemerge_final)

foldchangemerge_final_melt <- tidyr::gather(foldchangemerge_final, name, fold,-gene)
pvaluemerge_final_melt <-  tidyr::gather(pvaluemerge_final, name, FDR,-gene)
merged_final <- inner_join(foldchangemerge_final_melt,pvaluemerge_final_melt)
merged_final$FDR[merged_final$FDR > 3] <- 3
colnames(merged_final)[4] <- "FDR"
merged_final <- merged_final[merged_final$name %in% celltype,]
merged_final$gene <- factor(merged_final$gene,levels = rev(topgenes_shared))
merged_final$name <- factor(merged_final$name,levels = colnames(foldchangemerge_final))



newlevels <- unique(gsub("_Fruc|_HFHS","",gsub("oligodendrocytes precursor","Oligodendrocytes precursor",colnames(foldchangemerge_final))))
merged_final$name <- as.character(merged_final$name)
merged_final$name <- gsub("_Fruc",":Fruc",merged_final$name)
merged_final$name <- gsub("_HFHS",":HFHS",merged_final$name)
merged_final$treatment <- sapply(strsplit(as.character(merged_final$name),":"), `[[`, 2)
merged_final$name <- factor(sapply(strsplit(as.character(merged_final$name),":"), `[[`, 1),levels = newlevels)
#p1_hyp <- ggplot(merged_final, aes(y = gene,x = name)) +facet_wrap(.~treatment)+        ## global aes
#  geom_tile(fill = "white") +  xlab("")+ylab("gene")+       ## to get the rect filled
# geom_point(aes(colour = fold, 
#                 size = FDR))  +    ## geom_point for circle illusion
#  scale_color_gradientn(limits = c(-5,5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+   ## color of the corresponding aes          ## to tune the size of circles
#  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
#                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
#                   strip.text.x = element_text(size = 20),
#                   plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))+ guides(colour=guide_legend(title="-log(fold change)"),
#                                                                       size =guide_legend(title="-log(FDR)"))
p1_hyp <- ggplot(merged_final, aes(y = name,x = gene)) +facet_grid(treatment~.)+        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = fold, 
                 size = FDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-5,5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+   ## color of the corresponding aes          ## to tune the size of circles
  theme_bw()+theme(axis.text.y = element_text( size = 16),
                   axis.text.x = element_text( angle = 45, hjust = 1, size = 18,face = "italic"),
                   strip.text.y  = element_text(size = 20),
                   plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))+ guides(colour=guide_legend(title="-log(fold change)"),size =guide_legend(title="-log(FDR)")) 



pdf(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",pathwayname,"_",Tissue,".pdf"),width = 20, height = 5.5)
print(p1_hyp)
dev.off()

library(gridExtra)
#pdf(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/combined_selected_pathwaygene.pdf"),width =14, height = 20)
#grid.arrange(p1_SVF,p1_liver,p1_hyp, layout_matrix = cbind(c(2,2,1,1,1,1),c(3,3,3,3,3,3)))
#dev.off()
save(p1_SVF,p1_liver,p1_hyp, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Dot_heatmap_selected_individual_pathwa.rda")
pdf(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/combined_selected_pathwaygene.pdf"),width =20, height = 14)
grid.arrange(p1_SVF,p1_liver,p1_hyp, layout_matrix = rbind(c(2,2,2,1,1,1,1),c(3,3,3,3,3,NA,NA)))
dev.off()





#picking up specific pathways - old method
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined/SVF_Final_Set1only.rda")
DefaultAssay(dropEST.combined.filtered) <- "RNA"
dropEST.combined.filtered <- ScaleData(object = dropEST.combined.filtered, features = rownames(x = dropEST.combined.filtered))
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "celltype")


pathwayname <- "negative regulation of apoptotic process"
status <- "Fruc"

filtered <- result_final[result_final$pathwayname %in% pathwayname,]
filtered <- filtered[grep(status,filtered$celltype),]
currentgenesets <- unique(firstup(unlist(strsplit(filtered$genevector,";"))))
object <- subset(dropEST.combined.filtered,idents= unique(gsub("_HFHS|_Fruc","",filtered$celltype)))
object <- SetIdent(object, value = "data.diseasestatus")
object <- subset(object, idents = c("Control",status))
object[["celltype_treatment"]] <- paste0(object$celltype, object$data.diseasestatus)
object <- SetIdent(object, value = "celltype_treatment")
object =  subset(object, downsample = 100)
object <- ScaleData(object = object, features = rownames(x = object))
cols.use <- list(data.diseasestatus=c('gray', 'blue'))
DoMultiBarHeatmap(object = object, assay = "RNA",slot = "scale.data",group.by="celltype",
                  additional.group.by="data.diseasestatus", additional.group.sort.by = "data.diseasestatus",features = currentgenesets, cols.use=cols.use)

status <- "HFHS"
filtered <- result_final[result_final$pathwayname %in% pathwayname,]
filtered <- filtered[grep(status,filtered$celltype),]
currentgenesets <- unique(firstup(unlist(strsplit(filtered$genevector,";"))))
object <- subset(dropEST.combined.filtered,idents= unique(gsub("_HFHS|_Fruc","",filtered$celltype)))
object <- SetIdent(object, value = "data.diseasestatus")
object <- subset(object, idents = c("Control",status))
object[["celltype_treatment"]] <- paste0(object$celltype, object$data.diseasestatus)
object <- SetIdent(object, value = "celltype_treatment")
object =  subset(object, downsample = 100)
object <- ScaleData(object = object, features = rownames(x = object))
cols.use <- list(data.diseasestatus=c('gray', 'red'))
DoMultiBarHeatmap(object = object, assay = "RNA",slot = "scale.data",group.by="celltype",
                  additional.group.by="data.diseasestatus",features = currentgenesets, cols.use=cols.use)



for(r in unique(result_final$group)){
  #currentgenesets <- firstup(unlist(strsplit(result_final$genevector[result_final$group %in% r],";")))
  currentgenesets <- unique(firstup(unlist(strsplit(result_final$genevector[result_final$pathwayname %in% pathwayname],";"))))
  matchingcelltypes <- s
  #genesOrder = c("Bmp2","Clec4f","Alb","Fabp1","Cyp2f2","Cyp2e1","Nkg7","Itgax","H2-Ab1","Ms4a1","Reln","Spp1","Top2a","Siglech","Ccr9")
  #dropEST.combined.filtered <- ScaleData(object = dropEST.combined.filtered, features = rownames(x = dropEST.combined.filtered))
  #DoHeatmap(subset(dropEST.combined.filtered, downsample = 100), features = genesOrder, size = 3,slot = "data")
  
}




DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = Matrix::t(x = GetAssayData(object = object, 
                                                                     slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])+ scale_fill_gradientn(colors = c("blue", "white", "red"))
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}