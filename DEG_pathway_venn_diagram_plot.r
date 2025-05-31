library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

intersectitems <- function(Set1, Set2){
  result <- list()
  result[["Set1"]] <- c(NA,NA)
  result[["Set1"]] <-setdiff(Set1,Set2)
  result[["inters"]] <- c(NA,NA)
  result[["inters"]] <-intersect(Set1,Set2)
  result[["Set2"]] <- c(NA,NA)
  result[["Set2"]] <-setdiff(Set2,Set1)
  return(result)
}

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
pathway_for_venn <- read.csv("~/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Fructose_shared_Yen-Wei/pathway_for_venn.csv")
pathway_for_venn$Tissue[pathway_for_venn$Tissue %in% "liver"] <- "Liver"
rm(overlappingsummary)
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
    currentHFHS <- DEGlist[[celltype]]
    if(length(currentFruc) < 10 |length(currentHFHS) < 10) next
    v <- venn.diagram(list(Fruc = currentFruc,HFHS = currentHFHS),resolution = 200,width = 1000, height = 1000,
                      fill = c("darkslategray2", "darkseagreen2"),main = paste0(currenttissue,"_",finalcelltype),
                      alpha = 0.5, cat.cex = 1.8, cex=1.8, main.cex = 2,hyper.test = T,total.population = 4000,lower.tail = F,filename = NULL)
    #filename=paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Venntesting/",currenttissue,"_",finalcelltype,"DEGvenn.png"))
    assign(paste0(currenttissue,"_",finalcelltype,"_DEG"),v)
    pathway_for_venn$celltype <- gsub("ADSC_","APC_",pathway_for_venn$celltype)
    pathway_for_venn$celltype <- gsub("Preadipocytes_","Mesothelial_",pathway_for_venn$celltype)
    
    currentpathway_Fruc <- pathway_for_venn$pathwayname[pathway_for_venn$Tissue %in% currenttissue & pathway_for_venn$celltype %in% paste0(finalcelltype,"_Fruc")]
    currentpathway_HFHS <- pathway_for_venn$pathwayname[pathway_for_venn$Tissue %in% currenttissue & pathway_for_venn$celltype %in% paste0(finalcelltype,"_HFHS")]
    
    if(length(currentpathway_Fruc) == 0 & length(currentpathway_HFHS) == 0){
      cat(paste0(currenttissue,"_",celltype," "))
      next
    }
    
    if(length(currentpathway_Fruc) == 0 & length(currentpathway_HFHS) != 0){
      #      v <- venn.diagram(list(Fruc = currentpathway_Fruc),resolution = 200,width = 1000, height = 1000,
      #                        fill = c("blue"),
      #                        alpha = 0.5, cat.cex = 1.5, cex=1.5,filename = NULL)
      #                        filename=paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Venntesting/",currenttissue,"_",finalcelltype,"pathwayvenn.png"))
      next     
    }else if(length(currentpathway_HFHS) == 0 & length(currentpathway_Fruc) != 0){
      #  v <- venn.diagram(list(HFHS = currentpathway_HFHS),resolution = 200,width = 1000, height = 1000,
      #                    fill = c("red"),
      #                    alpha = 0.5, cat.cex = 1.5, cex=1.5,filename = NULL)
      #                        filename=paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Venntesting/",currenttissue,"_",finalcelltype,"pathwayvenn.png"))
      next 
    }else{
      currentinters <- intersectitems(currentpathway_Fruc, currentpathway_HFHS)
      v <- venn.diagram(list(Fruc = currentpathway_Fruc,HFHS = currentpathway_HFHS),resolution = 200,width = 1000, height = 1000,
                        fill = c("darkslategray2", "darkseagreen2"),main = paste0(currenttissue,"_",finalcelltype,"_pathway"),
                        alpha = 0.5, cat.cex = 1.8, cex=1.8, main.cex = 2,hyper.test = T,total.population = 500,lower.tail = F,filename = NULL)
      #                     filename=paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Venntesting/",currenttissue,"_",finalcelltype,"pathwayvenn.png"))
      
    }
    frame <- data.frame(tissue  = currenttissue, celltype = finalcelltype, pathway = paste0(currentinters[["Set1"]], collapse = ","), condition = "Fruc specific")
    frame <- rbind.data.frame(frame,data.frame(tissue  = currenttissue, celltype = finalcelltype, pathway = paste0(currentinters[["inters"]], collapse = ","), condition = "shared in HFHS and Fruc"))
    frame <- rbind.data.frame(frame,data.frame(tissue  = currenttissue, celltype = finalcelltype, pathway = paste0(currentinters[["Set2"]], collapse = ","), condition = "HFHS specific"))
    if(!exists("overlappingsummary")){
      overlappingsummary <- frame
    }else{
      overlappingsummary <- rbind.data.frame(overlappingsummary, frame)
    }
    
    assign(paste0(currenttissue,"_",finalcelltype,"_Pathway"),v)
    assign(paste0(currenttissue,"_",finalcelltype,"_Pathway_items"),currentinters)
  }
  
}
overlappingsummary <- overlappingsummary[!overlappingsummary$pathway %in% "",]
write.csv(overlappingsummary, file = "/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/plot_for_resource_version/SupTable5_pathway_overlap.csv",row.names = F)
Hyp_Astrocytes_Pathway_items
`SI_Proximal enterocytes_Pathway_items`
library(gridExtra)
setwd("/Users/Tsai_Lab/Downloads/")
pdf("All Venn combined_SI.pdf",height = 12,width = 5)
grid.arrange(gTree(children = `SI_Proximal enterocytes_DEG`),gTree(children = SI_Goblet_DEG),gTree(children = `SI_Distal enterocytes_DEG`),gTree(children = `SI_T cells_DEG`), ncol = 1)
dev.off()

pdf("All Venn combined_Liver.pdf",height = 12,width = 10)
grid.arrange(gTree(children = `Liver_Pericentral hepatocytes_DEG`),gTree(children = `Liver_Periportal hepatocytes_DEG`),gTree(children = `Liver_Sinusoidal endothelial cells_DEG`),gTree(children = Liver_Kupffer_DEG),
             gTree(children = `Liver_NKT cells_DEG`),gTree(children = `Liver_Classical dendritic cells_DEG`),gTree(children = Liver_Cholangiocytes_DEG),gTree(children = `Liver_Dividing cells_DEG`),ncol = 2)
dev.off()

pdf("All Venn combined_SVF.pdf",height = 15,width = 5)
grid.arrange(gTree(children = SVF_APC_Agt_DEG),gTree(children = SVF_APC_Hsd11b1_DEG),gTree(children = SVF_APC_Pi16_DEG),gTree(children = `SVF_M2 Macrophage_DEG`), ncol = 1)
dev.off()

pdf("All Venn combined_Hyp.pdf",height = 12,width = 10)
grid.arrange(gTree(children = `Hyp_Glutamatergic neurons_DEG`),gTree(children = `Hyp_GABAergic neurons_DEG`),gTree(children = `Hyp_General neurons_DEG`),
             gTree(children = Hyp_Astrocytes_DEG),gTree(children = `Hyp_Endothelial cells_DEG`),gTree(children = `Hyp_Oligodendrocytes precursor_DEG`),gTree(children = `Hyp_Myelinating oligodendrocytes_DEG`),
             gTree(children = `Hyp_Ependymal cells_DEG`), ncol = 2)
dev.off()




setwd("/Users/Tsai_Lab/Downloads/")
pdf("All Venn combined_SI_Path.pdf",height = 12,width = 7)
grid.arrange(gTree(children = `SI_Proximal enterocytes_Pathway`),gTree(children = SI_Goblet_Pathway),gTree(children = `SI_T cells_Pathway`), ncol = 1)
dev.off()

pdf("All Venn combined_Liver_Path.pdf",height = 12,width = 16)
grid.arrange(gTree(children = `Liver_Pericentral hepatocytes_Pathway`),gTree(children = `Liver_Periportal hepatocytes_Pathway`),gTree(children = `Liver_Sinusoidal endothelial cells_Pathway`),gTree(children = Liver_Kupffer_Pathway),
            ncol = 2)
dev.off()

pdf("All Venn combined_SVF_Path.pdf",height = 15,width = 7)
grid.arrange(gTree(children = SVF_APC_Hsd11b1_Pathway),gTree(children = SVF_APC_Pi16_Pathway), ncol = 1)
dev.off()

pdf("All Venn combined_Hyp_Path.pdf",height = 12,width = 12)
grid.arrange(gTree(children = `Hyp_Glutamatergic neurons_Pathway`),gTree(children = `Hyp_GABAergic neurons_Pathway`),gTree(children = `Hyp_General neurons_Pathway`),
             gTree(children = Hyp_Astrocytes_Pathway),gTree(children = `Hyp_Oligodendrocytes precursor_Pathway`),
             ncol = 2)
dev.off()

test <- as.data.frame(table(pathway_for_venn$Tissue, pathway_for_venn$celltype))
test <- test[!test$Freq == 0,]
