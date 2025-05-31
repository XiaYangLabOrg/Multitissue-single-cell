#Hoffman version
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  a = 1
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#Hoffman version

setwd("/u/project/xyang123/aleph999/Mergeomics-master")
source("./Mergeomics.R")


#merging eQTL files
library(dplyr)

#alltissues <- c("SI","SVF","HEP","NPC","Hyp")
#allloci <- c("/u/flashscratch/j/jading/resources/mapping/GTEx_v8_eQTL/Small_Intestine_Terminal_Ileum.txt",
#             "/u/flashscratch/j/jading/resources/mapping/GTEx_v8_eQTL/Adipose_Subcutaneous.txt",
#             "/u/flashscratch/j/jading/resources/mapping/GTEx_v8_eQTL/Liver.txt",
#             "/u/flashscratch/j/jading/resources/mapping/GTEx_v8_eQTL/Liver.txt",
#             "/u/flashscratch/j/jading/resources/mapping/GTEx_v8_eQTL/Brain_Hypothalamus.txt")
#test <- read.table("/u/project/xyang123-nobackup/aleph999/Mergeomics-master/gene2loci.050kb.txt", sep = "\t", header = T,stringsAsFactors = F) 
#for(i in 1:length(allloci)){
#  test2 <- read.table(allloci[i], sep = "\t", header = T,stringsAsFactors = F) 
#  test3 <- anti_join(test2,test)
#  test3 <- rbind.data.frame(test,test2)
#  #test3 <- test3[!duplicated(test3),]
#  test3 <- test3[order(test3$GENE),]  
#  write.table(test3, file = paste0("/u/scratch/a/aleph999/",alltissues[i],"_eQTL_dist.txt"),sep = "\t",quote = F, row.names = F)
#}


#MSEA stuffs
allfiles <- list.files(path = "/u/project/xyang123/shared/datasets/GWAS/GWAS_Cleaned/",pattern = ".txt", include.dirs = F)
#allfiles <- allfiles[-1]
#alltissues <- c("SI","SVF","HEP","NPC","Hyp","combine")
#allloci <- rep("/u/project/xyang123/xyang123-NOBACKUP/aleph999/Mergeomics-master/gene2loci.050kb.txt",6)
#alltissues <- c("SI","SVF","HEP","NPC","Hyp","combine")


#allloci <- c("/u/flashscratch/j/jading/resources/mapping/GTEx_v8_eQTL/Small_Intestine_Terminal_Ileum.txt",
#             "/u/flashscratch/j/jading/resources/mapping/GTEx_v8_eQTL/Adipose_Subcutaneous.txt",
#             "/u/flashscratch/j/jading/resources/mapping/GTEx_v8_eQTL/Liver.txt",
#             "/u/flashscratch/j/jading/resources/mapping/GTEx_v8_eQTL/Liver.txt",
#             "/u/flashscratch/j/jading/resources/mapping/GTEx_v8_eQTL/Brain_Hypothalamus.txt",
#             "/u/project/xyang123-nobackup/aleph999/Mergeomics-master/gene2loci.050kb.txt")
allloci <- c(             "/u/project/xyang123/aleph999/Mergeomics-master/gene2loci.050kb.txt")

alltissues <- c("combine")
#allloci <- paste0("/u/scratch/a/aleph999/",alltissues,"_eQTL_dist.txt")

for(i in 1:length(alltissues)){
  job.ssea <- list()
  job.ssea$label <- paste0(alltissues[i],"_",allfiles[a])
  job.ssea$folder <- "/u/home/a/aleph999/"
  job.ssea$genfile <- allloci
  job.ssea$locfile <- paste0("/u/project/xyang123/shared/datasets/GWAS/GWAS_Cleaned/",allfiles[a])
  job.ssea$modfile <- paste0("/u/project/xyang123/aleph999/Mergeomics-master/MSEA_",alltissues[i],"_DEG.txt")
  job.ssea$inffile <- paste0("/u/project/xyang123/aleph999/Mergeomics-master/MSEA_",alltissues[i],"_descr.txt")
  job.ssea$permtype <- "gene"
  job.ssea$nperm <- 10000
  job.ssea$maxgenes <- 700
  job.ssea <- ssea.start(job.ssea)
  job.ssea <- ssea.prepare(job.ssea)
  job.ssea <- ssea.control(job.ssea)
  job.ssea <- ssea.analyze(job.ssea,trim_start=0.005,trim_end=0.995)
  job.ssea <- ssea.finish(job.ssea)
  
}

