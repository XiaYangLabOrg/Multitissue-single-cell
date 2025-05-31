#plot selected  pathways
options(stringsAsFactors = F)
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

alltissue <- c("SI","SVF","Hyp","Liver")
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/")


#setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined")
#Set <- "Set2"

for(Tissue in 1:3){
  setwd(alltissue[Tissue])
  result_final <-  read.csv("pathwayinfo_combined_filtered_pathwaysize.csv")
  result_final$genesize <- sapply(strsplit(result_final$genevector,";"),length)
  
  ##writing down selected and filtered pathways
  result_final <- result_final[result_final$genesize  > 3,]
  #write.csv(result_final, file = "pathwayinfo_combined_filtered_pathwaysize.csv")
  #result_final$pathwayname <- firstup(result_final$pathwayname)
  selectedpathways <- as.character(result_final$pathwayname[result_final$selected %in% "*"])
  
  #write.csv(data.frame(pathwayname =  unique(selectedpathways), group = NA), file = "pathwayinfo_combined_filtered_selected.csv")
  #After this, further adding manual annotations on these pathways in order to further organize these pathways
  selectedpathways2 <- read.csv("./pathwayinfo_combined_filtered_selected.csv")
  result_final <- inner_join(result_final,selectedpathways2)
  #result_final <- result_final[result_final$pathwayname %in% selectedpathways,]
  
  result_final$FDR[result_final$FDR > 3 ] <- 3
  result_final$fold[result_final$fold > 3 ] <- 3
  result_final$fold[result_final$fold < -3 ] <- -3
  result_final$celltype <- factor(result_final$celltype, levels = unique(result_final$celltype))
  result_final <- result_final[order(result_final$group,result_final$pathwayname),]
  
  result_final$pathwayname <- gsub(" (","",result_final$pathwayname,fixed = T)
  result_final$pathwayname <- gsub(")","",result_final$pathwayname,fixed = T)
  result_final$pathwayname <- gsub("GO:[0-9]+| hsa[0-9]+","",result_final$pathwayname)
  result_final$pathwayname <- gsub(" Homo sapiens","",result_final$pathwayname)
  result_final$pathwayname <- firstup(result_final$pathwayname)
  result_final$pathwayname <- factor(result_final$pathwayname , levels = unique(result_final$pathwayname))
  
  newlevels <- unique(gsub("_Fruc|_HFHS","",gsub("oligodendrocytes precursor","Oligodendrocytes precursor",as.character(unique(result_final$celltype)))))
  result_final$celltype <- as.character(result_final$celltype)
  result_final$celltype <- gsub("_Fruc",":Fruc",result_final$celltype)
  result_final$celltype <- gsub("_HFHS",":HFHS",result_final$celltype)
  result_final$treatment <- sapply(strsplit(as.character(result_final$celltype),":"), `[[`, 2)
  result_final$celltype <- factor(sapply(strsplit(as.character(result_final$celltype),":"), `[[`, 1),levels = newlevels)
  p2 <- ggplot(result_final, aes(y = pathwayname,x = celltype)) +facet_grid(.~treatment, scale = "free",space = "free")+        ## global aes
    geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
    geom_point(aes(colour = fold, 
                   size = FDR))  +    ## geom_point for circle illusion
    scale_color_gradientn(limits = c(-3.1,3.1),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+   ## color of the corresponding aes          ## to tune the size of circles
    theme_bw()+theme(axis.text.y = element_text(size = 16),
                     axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                     strip.text.x = element_text(size = 20),
                     plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))+ guides(colour=guide_legend(title="-log(fold change)"),
                                                                         size =guide_legend(title="-log(FDR)"))
  
  
  
  p2

  
  assign( paste0("p2_",alltissue[Tissue]),p2)
  setwd("..")
  #pdf(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",alltissue[Tissue],"/pathway_dotplot.pdf"),width = 14,height = 8)
  #print(p2)
  #dev.off()
  
}
p2_SI
p2_Liver
library(gridExtra)
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/dotgeneplot.rda")
pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Combined_figure3_pathway2.pdf",width =15, height = 22)
grid.arrange(p2_SI,p2_SVF,p2_liver,p2_Hyp, heights = c(unit(1.1, "npc"),unit(1.45, "npc"),unit(1.55, "npc"),unit(1.1, "npc")),ncol = 1)
dev.off()

mat <- rbind(c(1, 1,	1,	2,	2,	2),c(1, 1,	1,	2,	2,	2),c(1, 1,	1,	2,	2,	2),c(1, 1,	1,	2,	2,	2),
                    c(3	,3	,4	,4	,4	,4),c(3	,3	,4	,4	,4	,4),
                    c(5	,5,	5	,6,	6	,6),c(5	,5,	5	,6,	6	,6),c(5	,5,	5	,6,	6	,6),
                    c(7	,7,	7	,7,	7	,7),c(7	,7,	7	,7,	7	,7))
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Dot_heatmap_selected_individual_pathwa.rda")

pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Figure4.pdf",width =28, height = 24)
grid.arrange(p2_liver,p2_SVF,p1_liver,p1_SVF,p2_SI,p2_Hyp,p1_hyp,layout_matrix = mat)
dev.off()



library(svglite)
finalgg <- grid.arrange(p2_Liver,p2_SVF,p1_liver,p1_SVF,p2_SI,p2_Hyp,p1_hyp,layout_matrix = mat)
ggsave(file="/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Figure4.svg", plot=finalgg,width =28, height = 24)

