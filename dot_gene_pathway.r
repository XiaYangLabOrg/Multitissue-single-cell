#plot all top 5 genes
library(metap)
library(Seurat)
library(dplyr)
library(plyr)
library(heatmap3)
library(RColorBrewer)
library(gplots)
library(ggplot2)


top5genelist <- function(DEGlist){
  finalgenelist <- NULL
  for(i in 1:length(DEGlist)){
    finalgenelist <- c(finalgenelist, DEGlist[[i]][1:3])
  }
  finalgenelist <- unique(finalgenelist[!is.na(finalgenelist)])
  
  return(finalgenelist)
}

allpath <- c("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined",
             "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined",
             "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
allSet <- c("Set2","Set1","Set1")
alltissue <- c("SI","SVF","Hyp")


#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined")
#Set <- "Set2"

for(Tissue in 1:3){
  setwd(allpath[Tissue])
  Set <- allSet[Tissue]
  
  load( "Humaninefold.rda")
  #load( "Humaninefold_subset.rda")
  rownames(foldchangemerge_final) <- foldchangemerge_final$genename
  foldchangemerge_final <- foldchangemerge_final[,-1]
  foldchangemerge_final <- as.matrix(foldchangemerge_final)
  foldchangemerge_final[is.na(foldchangemerge_final)] <- 0
  foldchangemerge_final <- foldchangemerge_final[,grep(Set,colnames(foldchangemerge_final))]
  colnames(foldchangemerge_final) <- paste0(gsub(paste0(Set,"_"),"",colnames(foldchangemerge_final)),"_Fruc")
  
  
  rownames(pvaluemerge_final) <- pvaluemerge_final$genename
  pvaluemerge_final <- pvaluemerge_final[,-1]
  pvaluemerge_final <- as.matrix(pvaluemerge_final)
  pvaluemerge_final[is.na(pvaluemerge_final)] <- 1
  pvaluemerge_final <- -log10(pvaluemerge_final)
  pvaluemerge_final <- as.data.frame(pvaluemerge_final)
  pvaluemerge_final <- pvaluemerge_final[,grep(Set,colnames(pvaluemerge_final))]
  colnames(pvaluemerge_final) <- paste0(gsub(paste0(Set,"_"),"",colnames(pvaluemerge_final)),"_Fruc")
  
  topgenes_Fruc <- top5genelist(DEGlist)
  foldchangemerge_final_Fruc <- foldchangemerge_final
  pvaluemerge_final_Fruc <- pvaluemerge_final
  
  
  load( "HumaninefoldHFHS.rda")
  #load( "HumaninefoldHFHS_subset.rda")
  rownames(foldchangemerge_final) <- foldchangemerge_final$genename
  foldchangemerge_final <- foldchangemerge_final[,-1]
  foldchangemerge_final <- as.matrix(foldchangemerge_final)
  foldchangemerge_final[is.na(foldchangemerge_final)] <- 0
  foldchangemerge_final <- foldchangemerge_final[,grep(Set,colnames(foldchangemerge_final))]
  colnames(foldchangemerge_final) <- paste0(gsub(paste0(Set,"_"),"",colnames(foldchangemerge_final)),"_HFHS")
  
  rownames(pvaluemerge_final) <- pvaluemerge_final$genename
  pvaluemerge_final <- pvaluemerge_final[,-1]
  pvaluemerge_final <- as.matrix(pvaluemerge_final)
  pvaluemerge_final[is.na(pvaluemerge_final)] <- 1
  pvaluemerge_final <- -log10(pvaluemerge_final)
  pvaluemerge_final <- as.data.frame(pvaluemerge_final)
  pvaluemerge_final <- pvaluemerge_final[,grep(Set,colnames(pvaluemerge_final))]
  colnames(pvaluemerge_final) <- paste0(gsub(paste0(Set,"_"),"",colnames(pvaluemerge_final)),"_HFHS")
  
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
  merged_final$gene <- factor(merged_final$gene,levels = rev(topgenes_shared))
  merged_final$name <- gsub("oligodendrocytes precursor","Oligodendrocytes precursor",merged_final$name)
  
  merged_final$name <- factor(merged_final$name,levels =gsub("oligodendrocytes precursor","Oligodendrocytes precursor",colnames(foldchangemerge_final))  )
  merged_final$fold[merged_final$fold > 5] <- 5
  merged_final$fold[merged_final$fold < -5] <- -5
  merged_final <- merged_final[!merged_final$gene %in% "Malat1",]
  p1 <- ggplot(merged_final, aes(y = gene,
                                 x = name)) +        ## global aes
    geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
    geom_point(aes(colour = fold, 
                   size = FDR))  +    ## geom_point for circle illusion
    scale_color_gradientn(limits = c(-5,5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+   ## color of the corresponding aes          ## to tune the size of circles
    theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                     axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                     strip.text.x = element_text(size = 20))+ guides(colour=guide_legend(title="-log(fold change)"),
                                                                                size =guide_legend(title="-log(FDR)"))
  
  
  #p1
  #Another version- reduced name but with seperated label

  pdf(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",alltissue[Tissue],"/genedotplot.pdf"),width = 10, height = 8)
  print(p1)
  dev.off()
  assign(paste0("p1_",alltissue[Tissue]),p1)
  
  newlevels <- unique(gsub("_Fruc|_HFHS","",gsub("oligodendrocytes precursor","Oligodendrocytes precursor",colnames(foldchangemerge_final))))
  merged_final$name <- as.character(merged_final$name)
  merged_final$name <- gsub("_Fruc",":Fruc",merged_final$name)
  merged_final$name <- gsub("_HFHS",":HFHS",merged_final$name)
  merged_final$treatment <- sapply(strsplit(as.character(merged_final$name),":"), `[[`, 2)
  merged_final$name <- factor(sapply(strsplit(as.character(merged_final$name),":"), `[[`, 1),levels = newlevels)
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
  
  
  #p1
  assign(paste0("p1_",alltissue[Tissue]),p1)
  levelnames <- colnames(foldchangemerge_final)
  
  load("./pathwayinfo_Fruc.rda")
  load( "Humaninefold.rda")
  invalid <- NULL
  for(i in 1:nrow(result_final)){
    currentgenes <- DEGlist[[paste0(Set,"_",result_final$celltype[i])]]
    if(length(currentgenes) <= 20){
      invalid <- c(invalid, i)
    }
  }
  if(length(invalid) > 0){result_final <- result_final[-invalid,]}
  
  invalid <- NULL
  for(i in 1:nrow(result)){
    currentgenes <- DEGlist[[paste0(Set,"_",result$celltype[i])]]
    if(length(currentgenes) <= 20){
      invalid <- c(invalid, i)
    }
  }
  
  if(length(invalid) > 0){result_Fruc_original <- result[-invalid,]
  }else{result_Fruc_original <- result}
  result_Fruc_original$celltype <- paste0(result_Fruc_original$celltype,"_Fruc")
  
  #result_final <- result_final %>% 
  #  arrange(desc(FDR)) %>% 
  #  group_by(celltype) %>% slice(1:3)
  #load("./pathwayinfosubset_Fruc.rda")
  result_final$FDR[result_final$FDR > 3 ] <- 3
  result_final$fold[result_final$fold > 3 ] <- 3
  result_final$fold[result_final$fold < -3 ] <- -3
  #result_final$treatment <- "Fruc"
  result_final <- result_final[result_final$Set %in% Set,]
  result_final$celltype <- paste0(result_final$celltype,"_Fruc")
  result_Fruc <- result_final
  
  load("./pathwayinfo_HFHS.rda")
  load( "HumaninefoldHFHS.rda")
  invalid <- NULL
  for(i in 1:nrow(result_final)){
    currentgenes <- DEGlist[[paste0(Set,"_",result_final$celltype[i])]]
    if(length(currentgenes) <= 20){
      invalid <- c(invalid, i)
    }
  }
  if(length(invalid) > 0){result_final <- result_final[-invalid,]}
  
  
  invalid <- NULL
  for(i in 1:nrow(result)){
    currentgenes <- DEGlist[[paste0(Set,"_",result$celltype[i])]]
    if(length(currentgenes) <= 20){
      invalid <- c(invalid, i)
    }
  }
  if(length(invalid ) > 0){result_HFHS_original <- result[-invalid,]
  }else{result_HFHS_original <- result}
  
  result_HFHS_original$celltype <- paste0(result_HFHS_original$celltype,"_HFHS")
  #result_final <- result_final %>% 
  #  arrange(desc(FDR)) %>% 
  #  group_by(celltype) %>% slice(1:3)
  
  #load("./pathwayinfosubset_Fruc.rda")
  result_final$FDR[result_final$FDR > 3 ] <- 3
  result_final$fold[result_final$fold > 3 ] <- 3
  result_final$fold[result_final$fold < -3 ] <- -3
  #result_final$treatment <- "HFHS"
  result_final <- result_final[result_final$Set %in% Set,]
  result_final$celltype <- paste0(result_final$celltype,"_HFHS")
  result_HFHS <- result_final
  
  result_final <- rbind.data.frame(result_Fruc,result_HFHS)
  result_final$celltype <- factor(result_final$celltype,levels = levelnames)
  #intersect(unique(result_final$celltype),gsub("Set[0-9]_","",colnames(foldchangemerge_final)))
  result_final$pathwayname <- factor(result_final$pathwayname)
  
  result_original <- rbind.data.frame(result_Fruc_original,result_HFHS_original)
  write.csv(result_original, file = paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",alltissue[Tissue],"/pathwayinfo_combined.csv"))
  
  p2 <- ggplot(result_final, aes(y = pathwayname,x = celltype)) +        ## global aes
    geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
    geom_point(aes(colour = fold, 
                   size =FDR))  +    ## geom_point for circle illusion
    scale_color_gradientn(limits = c(-3.1,3.1),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+   ## color of the corresponding aes          ## to tune the size of circles
    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ guides(colour=guide_legend(title="-log(fold change)"),
                                                                                size =guide_legend(title="-log(FDR)"))
  
  p2
  assign( paste0("p2_",alltissue[Tissue]),p2)
  pdf(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",alltissue[Tissue],"/pathway_dotplot.pdf"),width = 14,height = 8)
  print(p2)
  dev.off()
  
}



library(gridExtra)
library(ggplot2)
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/dotgeneplot.rda")

pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Combined_figure3_genes.pdf",width = 16, height = 34)
grid.arrange(p1_SI,p1_SVF,p1_liver,p1_Hyp,heights = c(unit(0.5, "npc"),unit(0.4, "npc"),unit(0.65, "npc"),unit(0.75, "npc")),ncol = 1)
#grid.arrange(p1_SI,p1_liver,p1_SVF,p1_Hyp,heights = c(unit(0.5, "npc"),unit(0.65, "npc"),unit(0.4, "npc"),unit(0.65, "npc")),ncol = 2)
dev.off()

#pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Combined_figure3_pathway.pdf",width =16, height = 45)
#grid.arrange(p2_SI,p2_SVF,p2_liver,p2_Hyp,ncol = 1)
#dev.off()


