#plot all top 5 genes
library(metap)
library(Seurat)
library(dplyr)
library(plyr)
library(heatmap3)
library(RColorBrewer)
library(gplots)
library(ggplot2)

allpath <- c("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined",
             "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined",
             "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined",
             "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER")
allSet <- c("Set2","Set1","Set1","Set1")
alltissue <- c("SI","SVF","Hyp","Liver")
resultinfo <- c()
for(Tissue in 1:4){
  setwd(allpath[Tissue])
  Set <- allSet[Tissue]
  currenttissue <- alltissue[Tissue ]
  if(Tissue == 4){
    load( "HumaninefoldFruc_monocle.rda")
    DEGlist_up_Fruc <- DEGlist_up
    DEGlist_down_Fruc <- DEGlist_down
    DEGlist_Fruc <- DEGlist
    load( "HumaninefoldHFHS_monocle.rda")
  }else{
    load( "Humaninefold.rda")
    DEGlist_up_Fruc <- DEGlist_up
    DEGlist_down_Fruc <- DEGlist_down
    DEGlist_Fruc <- DEGlist
    load( "HumaninefoldHFHS.rda")
  }
  allcelltypes <- union(names(DEGlist),names(DEGlist_Fruc))

  for(celltype in allcelltypes){
    finalcelltype <- gsub("Set1_|Set2_","",celltype)
    currentFruc <- DEGlist_Fruc[[celltype]]
    resultinfo <- c(resultinfo, paste(currenttissue,finalcelltype,"Fruc","combined",length(currentFruc), sep = ";"))
    currentFruc_up <- DEGlist_up_Fruc[[celltype]]
    resultinfo <- c(resultinfo, paste(currenttissue,finalcelltype,"Fruc","up",length(currentFruc_up), sep = ";"))
    currentFruc_down <- DEGlist_down_Fruc[[celltype]]
    resultinfo <- c(resultinfo, paste(currenttissue,finalcelltype,"Fruc","down",length(currentFruc_down), sep = ";"))
    
    currentHFHS <- DEGlist[[celltype]]
    resultinfo <- c(resultinfo, paste(currenttissue,finalcelltype,"HFHS","combined",length(currentHFHS), sep = ";"))
    currentHFHS_up <- DEGlist_up[[celltype]]
    resultinfo <- c(resultinfo, paste(currenttissue,finalcelltype,"HFHS","up",length(currentHFHS_up), sep = ";"))
    currentHFHS_down <- DEGlist_down[[celltype]]
    resultinfo <- c(resultinfo, paste(currenttissue,finalcelltype,"HFHS","down",length(currentHFHS_down), sep = ";"))
    
    resultinfo <- c(resultinfo, paste(currenttissue,finalcelltype,"shared","combined",length(intersect(currentHFHS,currentFruc)), sep = ";"))
    resultinfo <- c(resultinfo, paste(currenttissue,finalcelltype,"shared","up",length(intersect(currentHFHS_up,currentFruc_up)), sep = ";"))
    resultinfo <- c(resultinfo, paste(currenttissue,finalcelltype,"shared","down",length(intersect(currentHFHS_down,currentFruc_down)), sep = ";"))
  }
  
}
finalinfo <- strsplit(resultinfo,";")
finalframe <- data.frame(tissue = sapply(finalinfo, `[[`, 1),
                         celltype = sapply(finalinfo, `[[`, 2),
                         status = sapply(finalinfo, `[[`, 3),
                         direction = sapply(finalinfo, `[[`, 4),
                         count = as.numeric(sapply(finalinfo, `[[`, 5)))
finalframe <- finalframe[order(finalframe$count, decreasing = T),]
finalframe$celltype <- factor(finalframe$celltype, levels = unique(finalframe$celltype))
finalframe$tissue <- factor(finalframe$tissue, levels = c("SI","SVF","Liver","Hyp"))


combined <- ggplot(finalframe[finalframe$direction %in% "combined",], aes(celltype, count)) + facet_grid(.~tissue,scales = "free", space = "free")+   
  geom_bar(aes(fill = status), position = "dodge", stat="identity")+xlab("")+ylab("")+
  scale_fill_manual(values=c("blue", "red", "gray40"))+
  theme_bw()+theme(axis.text.y = element_text(size = 16),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
                   strip.text.x = element_text(size = 20),
                   plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

up <- ggplot(finalframe[finalframe$direction %in% "up",], aes(celltype, count)) + facet_grid(.~tissue,scales = "free", space = "free")+   
  geom_bar(aes(fill = status), position = "dodge", stat="identity")+xlab("")+ylab("")+
  scale_fill_manual(values=c("blue", "red", "gray40"))+
  theme_bw()+theme(axis.text.y = element_text(size = 16),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
                   strip.text.x = element_text(size = 20),
                   plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))
down <- ggplot(finalframe[finalframe$direction %in% "down",], aes(celltype, count)) + facet_grid(.~tissue,scales = "free", space = "free")+   
  geom_bar(aes(fill = status), position = "dodge", stat="identity")+xlab("")+ylab("")+
  scale_fill_manual(values=c("blue", "red", "gray40"))+
  theme_bw()+theme(axis.text.y = element_text(size = 16),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
                   strip.text.x = element_text(size = 20),
                   plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))

library(gridExtra)
#pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Fig4_Genecounts.pdf", width = 16, height = 16)
#grid.arrange(combined, up, down, ncol = 1)
#dev.off()

load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/EucdistanceFig.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/tSNEFig2.rda")


#Matrix_fig2 <- as.matrix(read.csv("~/Desktop/Box Sync/Fructose/Matrix_fig2.csv", header=FALSE))
Matrix_fig2 <- as.matrix(read.csv("~/Desktop/Box Sync/Fructose/Matrix_fig2_new.csv", header=FALSE))

pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Combined_figure3_mergeomics.pdf",width =32, height = 20)
grid.arrange(SI_Fig1,SVF_Fig1,Liver_Fig1,Hyp_Fig1,p5,p6,layout_matrix = Matrix_fig2)
dev.off()



pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/plot_for_resource_version/SupFig2_combined.pdf",width =14, height = 6)
print(combined)
dev.off()
