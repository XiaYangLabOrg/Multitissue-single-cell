setwd("/Users/Tsai_Lab/Desktop/Box Sync/Mergeomics-master")
source("./Mergeomics.R")

allfiles <- list.files(path = "/Users/Tsai_Lab/Desktop/Box Sync/GWAS/Metabolic")
for(i in 10:42){
  job.ssea <- list()
  job.ssea$label <- allfiles[i]
  job.ssea$folder <- "./"
  job.ssea$genfile <- "/Users/Tsai_Lab/Desktop/Box Sync/Mergeomics-master/gene2loci.050kb.txt"		
  job.ssea$locfile <- paste0("/Users/Tsai_Lab/Desktop/Box Sync/GWAS/Metabolic/",allfiles[i])
  #job.ssea$modfile <- "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/DEGsforMSEA.txt"
  #job.ssea$inffile <- "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/ModdescrMSEA.txt"
  job.ssea$modfile <- "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/MSEA_combine_networktarget.txt"
  job.ssea$inffile <- "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/MSEA_combine_networktarget_descr.txt"
  
  job.ssea$permtype <- "gene"
  job.ssea$nperm <- 10000
  job.ssea <- ssea.start(job.ssea)
  job.ssea <- ssea.prepare(job.ssea)
  job.ssea <- ssea.control(job.ssea)
  job.ssea <- ssea.analyze(job.ssea,trim_start=0.005,trim_end=0.995)
  job.ssea <- ssea.finish(job.ssea)
  
}

#Loading up results from Hoffman2 server

library(dplyr)
setwd("/u/project/xyang123/xyang123-NOBACKUP/aleph999/Mergeomics-master/Part2/msea")
setwd("/u/project/xyang123/xyang123-NOBACKUP/aleph999/Mergeomics-master/msea")

allresultfiles <- list.files( pattern = "_.*results.txt")
rm(finalfile)
for(i in 1:length(allresultfiles)){
  currentfile = read.table(allresultfiles[i], sep = "\t",header = T,stringsAsFactors = F)
  currentfile = currentfile[!currentfile[,1] %in% c("_ctrlB","_ctrlA"),]
  currentfile = currentfile[,c(1,7)]
  currentfile$celltype <- gsub("_HFHS|_Fruc","",currentfile$MODULE)
  treatment <- strsplit(currentfile$MODULE,"_")
  currentfile$treatment <-  sapply(treatment,function(x) x[2])
  currentfile <- currentfile[,-1]
  colnames(currentfile)[1] <- gsub(".txt.results.txt","",allresultfiles[i])
  currentfile <- currentfile[,c(2,3,1)]
  currentfile[,3] <- -log10(currentfile[,3])
  if(!exists("finalfile")){
    finalfile <- currentfile
  }else{
    finalfile <- full_join(finalfile,currentfile)
  }
}
finalfile <- finalfile[order(finalfile$celltype,finalfile$treatment),]
write.csv(finalfile,"/u/project/xyang123/xyang123-NOBACKUP/aleph999/Mergeomics-master/PART2_50kb_MSEA_finalfile.csv",row.names = F)
write.csv(finalfile,"/u/project/xyang123/xyang123-NOBACKUP/aleph999/Mergeomics-master/PART1_50kb_MSEA_finalfile.csv",row.names = F)

write.csv(finalfile,"/u/project/xyang123/xyang123-NOBACKUP/aleph999/Mergeomics-master/50kb_MSEA_finalfile.csv",row.names = F)


library("gplots")
library("devtools")
library(readr)
#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

library(RColorBrewer)
cols <- c("white",colorRampPalette(brewer.pal(11,"OrRd"))(39))

#Set a working directory for output files

#finalfile <- read_csv("~/Desktop/Box Sync/Mergeomics-master/PART1_50kb_MSEA_finalfile.csv")
#finalfile <- read_csv("~/Desktop/Box Sync/Mergeomics-master/PART2_50kb_MSEA_finalfile.csv")


main_title="distance_tokn"
par(cex.main=1)
colnames(finalfile) <- gsub("combine_","",colnames(finalfile))
rownames(finalfile) <- paste0(finalfile$celltype,"_",finalfile$treatment)
finalmat <- as.matrix(finalfile[,-c(1,2)])
rownames(finalmat) <- finalfile$celltype
rowcol <- ifelse(finalfile$treatment == "Fruc","red","blue")
rowcol <- as.matrix(rowcol)
colnames(rowcol) = "treatment"


finalmat <- finalmat[,-which(colSums(finalmat > 1.3,na.rm = T) < 3)]
rownames(finalmat) <- gsub("_oligo","_Oligo",rownames(finalmat))
rownames(finalmat) <- gsub("Gabaergic","GABAergic",rownames(finalmat))
rownames(finalmat) <- gsub("_ADSC","_APC",rownames(finalmat))
rownames(finalmat) <- gsub("_Preadipocytes","_Mesothelial",rownames(finalmat))
rownames(finalmat) <- gsub("Macrophage_Pf4","M2 Macrophage",rownames(finalmat))
rownames(finalmat) <- gsub("Macrophage_Il1b","M1 Macrophage",rownames(finalmat))


#GWAStable <- read.csv("~/Desktop/Box Sync/Mergeomics-master/List of publicly available human GWAS_for_Fructose.csv")
GWAStable <- read.csv("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/List of publicly available human GWAS_for_Fructose.csv")
GWAStable <- GWAStable[!GWAStable$Remove %in% "*",]
GWAStable$name <- paste(GWAStable$Consortium,GWAStable$Trait,sep = "_")
setdiff(colnames(finalmat), GWAStable$name)
GWAStable <- GWAStable[GWAStable$name %in% colnames(finalmat),]
finalmat <- finalmat[,colnames(finalmat) %in% GWAStable$name]

test <- finalmat[which(rowcol %in% "red"),]
test <- as.data.frame(t(test))
colnames(test) <- paste0("Fruc_",colnames(test))
test$name <- rownames(test)
test_Fruc <- test

test <- finalmat[which(rowcol %in% "blue"),]
test <- as.data.frame(t(test))
colnames(test) <- paste0("HFHS_",colnames(test))
test$name <- rownames(test)
test <- cbind.data.frame(test_Fruc[,-ncol(test_Fruc)],test)
GWAStable <- inner_join(GWAStable,test)
GWAStable <- GWAStable[,!colnames(GWAStable) %in% "Remove"]

pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Mergeomics_heatmapFruc.pdf",height = 16, width = 16)
heatmap.3(finalmat[which(rowcol %in% "red"),], na.rm = TRUE, scale="none", dendrogram="col", margins=c(30,30),
                  Rowv=T, Colv=T,  symbreaks=FALSE, key=TRUE, symkey=FALSE, breaks = c(0,1.3,seq(1.3, max(finalmat,na.rm = T),length.out = 39)),
                  density.info="none", trace="none", main="Fruc", labCol=colnames(finalmat), cexRow=2.2,cexCol = 2.2, col=cols, KeyValueName="-log(FDR)")
dev.off()
