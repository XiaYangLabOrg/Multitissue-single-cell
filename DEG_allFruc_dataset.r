#conservativecompare_all Fruc
#conservative compare final
library(metap)
library(Seurat)
library(dplyr)
library(plyr)

rm(list=ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined")
#load("SI_Set1_DropEST.Object1.rda")
#load("SI_HFHS_demulti.rda") #added HFHS demultiplexed samples

#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "Proximal endothelium","1" = "Proximal endothelium","3" = "Proximal endothelium","2" = "Proximal endothelium", "4" = "T cells","7" = "Goblet", "5" = "Proximal endothelium",
#                                                    "6" = "Distal endothelium","8" = "Proximal endothelium","9" = "Distal endothelium_Pmp22","10" = "Proliferative cells","11" = "I cells","12" = "S cells","13" = "Dendritic cells","14" = "Tuft","15" = "B cell"))
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/SI_combined_Set2_final.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SI.rda")


minimump2 <- function(x){
  NAmean <- mean(is.na(x))
  if(NAmean > 0.25){return(NA)
  }else{
    return(minimump(na.omit(x))$p)
  }
}
maximump2 <- function(x){
  NAmean <- mean(is.na(x))
  if(NAmean > 0.25){return(NA)
  }else{
    return(maximump(na.omit(x))$p)
  }
}


dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")
allsamples <- unique(dropEST.combined.filtered@meta.data$orig.ident)
sampleinfo <- substr(allsamples,1,4)
treatmentinfo <- substr(allsamples,5,5)
treatmentinfo[treatmentinfo %in% c("1","2")] <- "Control"
treatmentinfo[treatmentinfo %in% c("3","4")] <- "Fruc"
treatmentinfo[treatmentinfo %in% c("5","6")] <- "HFHS"
treatmentinfotable <- data.frame(allsamples,sampleinfo,treatmentinfo, stringsAsFactors = F)

#For small intestine where Set3 is available
#if(sum(treatmentinfotable$sampleinfo %in% "Set3") > 0 ){
#  treatmentinfotable$treatmentinfo[grep("Set32",treatmentinfotable$allsamples)] <- "Fruc"
#}
#


DEGlist <- list()
DEGlist2 <- list()
topDEG <- list()
Humaninelist <- list()
Pu1list <- list()
Fcgrtlist <- list()

DEGlist_up <- list()
DEGlist_down <- list()
methodused <- list()
DEGlist3 <- list()
Khklist <- list()
Aldoblist <- list()
Glut5list <- list()
DEGfulltable_single <- list()
DEGfulltable_multiple <- list()
DEGfulltableind <- c()
rm(pvaluemerge_final,foldchangemerge_final)
for(i in 2:2){
  currenttreatment <- treatmentinfotable[treatmentinfotable$sampleinfo %in% paste0("Set",i),]
  Controlsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "Control"]
  Frucsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "Fruc"]
  #Frucsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "HFHS"]
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "orig.ident")
  Celltypesubset <- subset(dropEST.combined.filtered, idents = c(Controlsample,Frucsample))
  #Celltypesubset = SetIdent(Celltypesubset, value = "integrated_snn_res.0.5")
  Celltypesubset = SetIdent(Celltypesubset, value = "celltype")
  
  alllevels <- levels(Celltypesubset@active.ident)
  allleveltable <- as.data.frame(table(Celltypesubset@active.ident, Celltypesubset@meta.data$orig.ident))
  
  
  
  for(j in 1:length(alllevels)){
    currentlevel <- alllevels[j]
    
    
    print(j)
    currentlevel <- alllevels[j]
    Celltypesubset2 <- subset(Celltypesubset, idents = currentlevel)
    DefaultAssay(Celltypesubset2) <- "RNA"
    Celltypesubset2 <- SetIdent(Celltypesubset2, value = "orig.ident")
    #set minimum cell representation per sample
    currentsample <- as.character(allleveltable$Var2[allleveltable$Var1 %in% currentlevel & allleveltable$Freq >= 3])
    availcontrolsample <- currentsample[currentsample %in% Controlsample]
    availfrucsample <- currentsample[currentsample %in% Frucsample]
    if(length(availcontrolsample) < 1 | length(availfrucsample) < 1){
      next
    }
    
    
    
    ###
    #1-comparison status
    if(length(availfrucsample) == 1 & length(availcontrolsample) == 1){
      #if only one sample per group, no need to replicate
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample, ident.2 =availcontrolsample , logfc.threshold = 0.05)
      
      finalDEG <- rownames(controltable[controltable$p_val_adj < 0.05,])
      #finalDEG <- rownames(controltable[controltable$p_val < 0.01,])
      methodused[[paste0("Set",i,"_",currentlevel)]] <- "Seurat_default"
      if(length(finalDEG) > 0){
        #indrib <- grep("^Rps|^Rpl|^Gm[0-9]|^mt-Rnr[0-9]",finalDEG)
        #if(length(indrib) > 0){
        #  finalDEG <- finalDEG[-indrib]
        #}
        DEGfulltable_single[[paste0("Set",i,"_",currentlevel)]] <- controltable
        DEGfulltableind[paste0("Set",i,"_",currentlevel)] <- "single"
        DEGlist2[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
        Humaninelist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["mt-Rnr2",]$avg_logFC,":",controltable["mt-Rnr2",]$p_val_adj)
        Pu1list[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Spi1",]$avg_logFC,":",controltable["Spi1",]$p_val_adj)
        Fcgrtlist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Fcgrt",]$avg_logFC,":",controltable["Fcgrt",]$p_val_adj)
        
        Khklist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Khk",]$avg_logFC,":",controltable["Khk",]$p_val_adj)
        Aldoblist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Aldob",]$avg_logFC,":",controltable["Aldob",]$p_val_adj)
        Glut5list[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Slc2a5",]$avg_logFC,":",controltable["Slc2a5",]$p_val_adj)
        
        DEGlist[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
        topDEG[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable)[1:5]
        DEGlist_up[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable[controltable$p_val_adj < 0.05 & controltable$avg_logFC > 0,])
        DEGlist_down[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable[controltable$p_val_adj < 0.05 & controltable$avg_logFC < 0,])
      }
      controltable$genename <- rownames(controltable)
      foldchangemerge <- controltable[,c("genename","avg_log2FC")]
      pvaluemerge <- controltable[,c("genename","p_val_adj")]
      colnames(foldchangemerge)[2] <- paste0("Set",i,"_",currentlevel)
      colnames(pvaluemerge)[2] <- paste0("Set",i,"_",currentlevel)
      
      if(!exists("foldchangemerge_final")){
        foldchangemerge_final <- foldchangemerge
        pvaluemerge_final <- pvaluemerge
      }else{
        foldchangemerge_final <- full_join(foldchangemerge_final,foldchangemerge, by = "genename")
        pvaluemerge_final <- full_join(pvaluemerge_final,pvaluemerge, by = "genename")
      }
      next
    }
    ###
    
    
    methodused[[paste0("Set",i,"_",currentlevel)]] <- "Dougmethod"
    #DEGlist2 with adjusted p value only
    controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample, ident.2 = availcontrolsample)
    finalDEG <- rownames(controltable[controltable$p_val_adj < 0.05,])
    if(length(finalDEG) > 0){
      DEGlist2[[paste0("Set",i,"_",currentlevel)]] <- finalDEG 
    }
    
    #For whole DEG across batches
    rm(finalptable,finalfoldtable)
    for(k in 1:length(availcontrolsample)){
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample,ident.2 = availcontrolsample[k], logfc.threshold = 0.05)
      controltable <- controltable[,c(1,2,5)]
      colnames(controltable) <- paste0(colnames(controltable),availcontrolsample[k])
      controltable$genename <- rownames(controltable)
      if(!exists("finalptable")){
        #finalptable <- controltable[,c(1,4)] #using raw pvalue
        finalptable <- controltable[,c(3,4)] #using adj pvalue
        finalfoldtable <- controltable[,c(2,4)]
      }else{
        #no NAs allowed, only consider full samples
        #finalptable <- full_join(finalptable,controltable[,c(1,4)])
        finalptable <- full_join(finalptable,controltable[,c(3,4)])
        finalfoldtable <- full_join(finalfoldtable, controltable[,c(2,4)])
      }
    }
    
    for(k in 1:length(availfrucsample)){
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample[k] ,ident.2 = availcontrolsample, logfc.threshold = 0.05)
      controltable <- controltable[,c(1,2,5)]
      colnames(controltable) <- paste0(colnames(controltable),availfrucsample[k])
      controltable$genename <- rownames(controltable)
      #no NAs allowed, only consider full samples
      #finalptable <- full_join(finalptable,controltable[,c(1,4)]) #using raw pvalue
      finalptable <- full_join(finalptable,controltable[,c(3,4)])
      finalfoldtable <- full_join(finalfoldtable, controltable[,c(2,4)])
    }
    
    rownames(finalptable) <- finalptable$genename
    finalptable <- finalptable[,!colnames(finalptable) %in% "genename"]
    allcompare <- ncol(finalptable)
    rownames(finalfoldtable) <- finalfoldtable$genename
    finalfoldtable <- finalfoldtable[,!colnames(finalfoldtable) %in% "genename"]
    
    finalptable$minP = apply(as.matrix(finalptable[,1:allcompare]), 1, minimump2)
    finalptable$maxP = apply(as.matrix(finalptable[,1:allcompare]), 1, maximump2)
    #finalptable$minPadj = p.adjust(finalptable$minP,method = 'fdr',n = nrow(Celltypesubset2@assays$RNA@data))
    #finalptable$maxPadj = p.adjust(finalptable$maxP,method = 'fdr',n = nrow(Celltypesubset2@assays$RNA@data))
    finalptable$minPadj = finalptable$minP
    finalptable$maxPadj = finalptable$maxP
    
    #finalptable$finalp <- apply(finalptable[,c("minPadj","maxPadj")],1,max)
    finalptable$finalp <- apply(finalptable[,c("minPadj","maxPadj")],1,min)
    #finalptable$FDR <- p.adjust(finalptable$finalp, method = "fdr")
    finalptable <- finalptable[!is.na(finalptable$finalp),]
    finalptable$genename <- rownames(finalptable)
    
    DEGfulltableind[paste0("Set",i,"_",currentlevel)] <- "multiple"
    DEGfulltable_multiple[[paste0("foldptable_Set",i,"_",currentlevel)]] <- finalfoldtable
    DEGfulltable_multiple[[paste0("finalptable_Set",i,"_",currentlevel)]] <- finalptable
    
    Humaninelist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["mt-Rnr2",]), na.rm = T),":",finalptable["mt-Rnr2",]$finalp)
    Pu1list[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Spi1",]), na.rm = T),":",finalptable["Spi1",]$finalp)
    Fcgrtlist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Fcgrt",]), na.rm = T),":",finalptable["Fcgrt",]$finalp)
    
    Khklist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Khk",]), na.rm = T),":",finalptable["Khk",]$finalp)
    Aldoblist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Aldob",]), na.rm = T),":",finalptable["Aldob",]$finalp)
    Glut5list[[paste0("Set",i,"_",currentlevel)]] <-paste0(median(as.numeric(finalfoldtable["Slc2a5",]), na.rm = T),":",finalptable["Slc2a5",]$finalp)
    
    #matching direction
    finalfoldtable$directionless <- apply(finalfoldtable[,1:allcompare],1,function(x){mean(x < 0)})
    finalfoldtable$directionmore <- apply(finalfoldtable[,1:allcompare],1,function(x){mean(x > 0)})
    finalfoldtable$avgfold <- apply(finalfoldtable[,1:allcompare],1,median,na.rm = T)
    finalfoldtable$genename <- rownames(finalfoldtable)
    
    foldchangemerge <- finalfoldtable[,c("genename","avgfold")]
    pvaluemerge <- finalptable[,c("genename","finalp")]
    colnames(foldchangemerge)[2] <- paste0("Set",i,"_",currentlevel)
    colnames(pvaluemerge)[2] <- paste0("Set",i,"_",currentlevel)
    
    if(!exists("foldchangemerge_final")){
      foldchangemerge_final <- foldchangemerge
      pvaluemerge_final <- pvaluemerge
    }else{
      foldchangemerge_final <- full_join(foldchangemerge_final,foldchangemerge, by = "genename")
      pvaluemerge_final <- full_join(pvaluemerge_final,pvaluemerge, by = "genename")
    }
    
    
    matchedgenes <- rownames(finalfoldtable[which(finalfoldtable$directionless >= 0.75 | finalfoldtable$directionmore >= 0.75),])
    upgenes <- rownames(finalfoldtable[which(finalfoldtable$directionmore >= 0.75) ,])
    downgenes <- rownames(finalfoldtable[which(finalfoldtable$directionless >= 0.75) ,])
    
    finalDEG <- intersect(matchedgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    #finalDEG <- intersect(matchedgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    #finalptable <- finalptable[finalptable$finalp < 0.05,]
    finalptable <- finalptable[order(finalptable$finalp),]
    #save(finalptable,file = paste0("Set",i,"_",currentlevel,"_Fruc.rda"))
    #save(finalptable,file = paste0("Set",i,"_",currentlevel,"_HFHS.rda"))
    finalptable <- finalptable[rownames(finalptable) %in% matchedgenes,]
    finalDEG_up <- intersect(upgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    finalDEG_down <- intersect(downgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    
    #finalDEG_up <- intersect(upgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    #finalDEG_down <- intersect(downgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    
    if(length(finalDEG) > 0){
      topDEG[[paste0("Set",i,"_",currentlevel)]] <- rownames(finalptable)[1:5]
      #exlude ribosomal from DEG list
      #indrib <- grep("^Rps|^Rpl|^Gm[0-9]|^mt-Rnr[0-9]",finalDEG)
      #if(length(indrib) > 0){
      #  finalDEG <- finalDEG[-indrib]
      #}
      DEGlist[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
      DEGlist_up[[paste0("Set",i,"_",currentlevel)]] <- finalDEG_up
      DEGlist_down[[paste0("Set",i,"_",currentlevel)]] <- finalDEG_down
      
    }
    
    
  }
  
}
save(DEGlist,DEGlist_up,DEGlist_down, Humaninelist,Khklist,Aldoblist,Glut5list,Fcgrtlist,Pu1list,DEGfulltable_multiple,DEGfulltable_single,DEGfulltableind,foldchangemerge_final,pvaluemerge_final,topDEG,file = "Humaninefold.rda")

rm(list=ls())
load("Humaninefold.rda")
library(httr)
set_config(config(ssl_verifypeer = FALSE))
options(RCurlOptions = list(ssl_verifypeer = FALSE))
options(rsconnect.check.certificate = FALSE)
library(enrichR)
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

enrichrfunction <- function(siggenes){
  if(length(siggenes) == 0){return(NA)}
  res <- enrichr(siggenes,databases = c("KEGG_2016","GO_Biological_Process_2017"))
  for(j in 1:length(res)){
    currentframe <- res[[j]]
    currentframe <- currentframe[currentframe$Adjusted.P.value < 0.05,]
    if(nrow(currentframe) > 0) currentframe <- currentframe[sapply(strsplit(currentframe$Genes,";"),length) > 3,]
    if(nrow(currentframe)  == 0){next
    }else{
      if(!exists("drugframe")){
        drugframe <- currentframe
      }else{
        drugframe <- rbind.data.frame(drugframe,currentframe)
      }
    }
  }
  if(exists("drugframe")){
    return(drugframe)
  }else{
    return(NULL)
  }
}

allcomparisons <- names(DEGlist)
siggenes_set1 <- NULL
siggenes_set1num <- NULL
siggenes_set2 <- NULL
siggenes_set2num<- NULL
siggenes_shared <- NULL
siggenes_sharednum <- NULL
sigpathways_set1 <- NULL
sigpathways_set1num <- NULL
sigpathways_set2 <- NULL
sigpathways_set2num <- NULL
sigpathways_sharedpath <- NULL
sigpathways_sharedpathnum <- NULL
sigpathways_sharedgenes <- NULL
sigpathways_sharedgenesnum <- NULL

alllevels  <- unique(gsub("Set[0-9]_","",allcomparisons))
foldchangevector <- NULL
pvaluevector <- NULL
genevector <- NULL
foldchangemerge_final$genename  <- firstup(foldchangemerge_final$genename)

for(i in 1:length(alllevels)){
  currentlevel <- alllevels[i]
  #matchingset <- allcomparisons[grep(paste0("_",currentlevel,"$"), allcomparisons)]
  matchingset <- allcomparisons[grep(paste0("_",currentlevel,"$"), allcomparisons)]
  if(length(matchingset) > 1){
    siggenes_set1[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
    
    tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
    if(!is.null(tempresult)){
      for(r in 1:nrow(tempresult)){
        name <- tempresult$Term[r]
        genevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Genes[r]
        foldchangevector[paste(name,currentlevel,"Set1",sep = "__")] <- 
          median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set1_",currentlevel)], na.rm = T)
        pvaluevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Adjusted.P.value[r]
      }
    }
    pathwayset1 <- tempresult$Term
    
    sigpathways_set1[i] <- paste(pathwayset1,collapse = ";")
    siggenes_set1num[i] <- length(DEGlist[[matchingset[1]]])
    sigpathways_set1num[i] <- length(pathwayset1)
    
    siggenes_set2[i] <- paste(DEGlist[[matchingset[2]]], collapse = " ")
    
    tempresult <-  enrichrfunction(DEGlist[[matchingset[2]]])
    if(!is.null(tempresult)){
      for(r in 1:nrow(tempresult)){
        name <- tempresult$Term[r]
        genevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Genes[r]
        foldchangevector[paste(name,currentlevel,"Set2",sep = "__")] <- 
          median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set2_",currentlevel)], na.rm = T)
        pvaluevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Adjusted.P.value[r]
      }
    }
    pathwayset2 <- tempresult$Term
    sigpathways_set2[i] <- paste(pathwayset2, collapse = ";")
    siggenes_set2num[i] <- length(DEGlist[[matchingset[2]]])
    sigpathways_set2num[i] <- length(pathwayset2)
    
    sigpathways_sharedpath[i] <- paste(intersect(pathwayset1,pathwayset2),collapse = " ")
    sigpathways_sharedpathnum[i] <- length(intersect(pathwayset1,pathwayset2))
    shared <- intersect(DEGlist[[matchingset[1]]],DEGlist[[matchingset[2]]])
    siggenes_shared[i] <- paste(shared,collapse = " ")
    siggenes_sharednum[i] <- length(shared)
    
    if(length(shared) > 0){
      DEGlist[[paste0("Shared",i,"_",currentlevel)]] <- shared
      sharedpathways <- enrichrfunction(shared)$Term
      sigpathways_sharedgenes[i] <- paste(sharedpathways,collapse = " ")
      sigpathways_sharedgenesnum[i] <- length(sharedpathways)
    }else{
      sigpathways_sharedgenes[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
    }
  }else if(length(matchingset) == 1){
    setname <- substr(matchingset,1,4)
    if(setname %in% "Set1"){
      siggenes_set1[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
      
      tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
      if(!is.null(tempresult)){
        for(r in 1:nrow(tempresult)){
          name <- tempresult$Term[r]
          genevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Genes[r]
          foldchangevector[paste(name,currentlevel,"Set1",sep = "__")] <- 
            median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set1_",currentlevel)], na.rm = T)
          pvaluevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Adjusted.P.value[r]
        }
      }
      pathwayset1 <- tempresult$Term
      
      sigpathways_set1[i] <- paste(pathwayset1,collapse = ";")
      siggenes_set1num[i] <- length(DEGlist[[matchingset[1]]])
      sigpathways_set1num[i] <- length(pathwayset1)
      
      siggenes_set2[i] <- NA
      sigpathways_set2[i] <- NA
      sigpathways_sharedpath[i] <- NA
      siggenes_set2num[i] <- NA
      sigpathways_set2num[i] <- NA
      
      siggenes_shared[i] <- NA
      sigpathways_sharedgenes[i] <- NA
      siggenes_sharednum[i] <- NA
      sigpathways_sharedpathnum[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
    }else if(setname %in% "Set2"){
      siggenes_set1[i] <- NA
      sigpathways_set1[i] <- NA
      siggenes_set1num[i] <- NA
      sigpathways_set1num[i] <- NA
      
      siggenes_set2[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
      
      tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
      if(!is.null(tempresult)){
        for(r in 1:nrow(tempresult)){
          name <- tempresult$Term[r]
          genevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Genes[r]
          foldchangevector[paste(name,currentlevel,"Set2",sep = "__")] <- 
            median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set2_",currentlevel)], na.rm = T)
          pvaluevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Adjusted.P.value[r]
        }
      }
      
      
      pathwayset2 <- tempresult$Term
      sigpathways_set2[i] <- paste(pathwayset2, collapse = ";")
      siggenes_set2num[i] <- length(DEGlist[[matchingset[1]]])
      sigpathways_set2num[i] <- length(pathwayset2)
      
      sigpathways_sharedpath[i] <- NA 
      
      siggenes_shared[i] <- NA
      sigpathways_sharedgenes[i] <- NA
      siggenes_sharednum[i] <- NA
      sigpathways_sharedpathnum[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
      
    }
  }else{
    siggenes_set1[i] <- NA
    sigpathways_set1[i] <- NA
    siggenes_set1num[i] <- NA
    sigpathways_set1num[i] <- NA
    
    siggenes_set2[i] <- NA
    sigpathways_set2[i] <-NA
    siggenes_set2num[i] <- NA
    sigpathways_set2num[i] <- NA
    
    sigpathways_sharedpath[i] <- NA 
    siggenes_shared[i] <- NA
    sigpathways_sharedgenes[i] <- NA
    siggenes_sharednum[i] <- NA
    sigpathways_sharedpathnum[i] <- NA
    sigpathways_sharedgenesnum[i] <- NA
  }
  
  
}
topDEG2 = sapply(topDEG,paste,collapse = ", ")
topDEG2 <- topDEG2[paste0("Set1_",alllevels)]
topDEG2 <- topDEG2[paste0("Set2_",alllevels)]

sigpathways_set1 <- gsub("Homo sapiens_","",sigpathways_set1)
sigpathways_set2 <- gsub("Homo sapiens_","",sigpathways_set2)
sigpathways_sharedpath <- gsub("Homo sapiens_","",sigpathways_sharedpath)
sigpathways_sharedgenes <- gsub("Homo sapiens_","",sigpathways_sharedgenes)

finalframe <- data.frame(names = alllevels,alllevels,siggenes_set1, siggenes_set1num, sigpathways_set1, sigpathways_set1num,
                         siggenes_set2, siggenes_set2num, sigpathways_set2, sigpathways_set2num,
                         siggenes_shared,siggenes_sharednum,sigpathways_sharedpath,sigpathways_sharedpathnum,sigpathways_sharedgenes,sigpathways_sharedgenesnum)

identical( names(foldchangevector),names(pvaluevector))
result <- data.frame(name = names(foldchangevector),fold = foldchangevector,FDR = -log10(pvaluevector), genevector)

splitname <- strsplit(result$name,"__")
result$pathwayname <- sapply(splitname, `[[`, 1)
result$celltype <- sapply(splitname, `[[`, 2)
result$Set <- sapply(splitname, `[[`, 3)
result$set_celltype <- paste0(result$celltype,"__",result$Set)

result_final <- result %>% group_by(set_celltype) %>% top_n(n = 5, wt = FDR)
save(result,result_final, file= "pathwayinfo_Fruc.rda")
write.csv(finalframe,file = "cross_batch_DEG_FDR_Doug_Fruc.csv")



rm(list=ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined")
#load("HEP_Set1_DropEST.Object1.rda")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "ADSC_Pi16","1" = "ADSC_Hsd11b1","3" = "ADSC_Hsd11b1","5" = "ADSC_Agt", "8" = "Preadipocytes","7" = "Endtohelial cells", "10" = "NK",
#                                                    "11" = "Germ cells","9" = "Macrophage_Mmp12","13" = "B cells","12" = "Dividing cells","6" = "Macrophage_Il1b","4" = "Macrophage_Pf4",
#                                                    "2" = "Macrophage_Pf4" ,"14" = "Mast cells","15"="Unknown"))
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined/SVF_Final_Set1only.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SVF.rda")
minimump2 <- function(x){
  NAmean <- mean(is.na(x))
  if(NAmean > 0.25){return(NA)
  }else{
    return(minimump(na.omit(x))$p)
  }
}
maximump2 <- function(x){
  NAmean <- mean(is.na(x))
  if(NAmean > 0.25){return(NA)
  }else{
    return(maximump(na.omit(x))$p)
  }
}


dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")
allsamples <- unique(dropEST.combined.filtered@meta.data$orig.ident)
sampleinfo <- substr(allsamples,1,4)
treatmentinfo <- substr(allsamples,5,5)
treatmentinfo[treatmentinfo %in% c("1","2")] <- "Control"
treatmentinfo[treatmentinfo %in% c("3","4")] <- "Fruc"
treatmentinfo[treatmentinfo %in% c("5","6")] <- "HFHS"
treatmentinfotable <- data.frame(allsamples,sampleinfo,treatmentinfo, stringsAsFactors = F)

#For small intestine where Set3 is available
#if(sum(treatmentinfotable$sampleinfo %in% "Set3") > 0 ){
#  treatmentinfotable$treatmentinfo[grep("Set32",treatmentinfotable$allsamples)] <- "Fruc"
#}
#


DEGlist <- list()
DEGlist2 <- list()
topDEG <- list()
Humaninelist <- list()
Pu1list <- list()
Fcgrtlist <- list()

DEGlist_up <- list()
DEGlist_down <- list()
methodused <- list()
DEGlist3 <- list()
Khklist <- list()
Aldoblist <- list()
Glut5list <- list()
DEGfulltable_single <- list()
DEGfulltable_multiple <- list()
DEGfulltableind <- c()
rm(pvaluemerge_final,foldchangemerge_final)
for(i in 1:1){
  currenttreatment <- treatmentinfotable[treatmentinfotable$sampleinfo %in% paste0("Set",i),]
  Controlsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "Control"]
  Frucsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "Fruc"]
  #Frucsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "HFHS"]
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "orig.ident")
  Celltypesubset <- subset(dropEST.combined.filtered, idents = c(Controlsample,Frucsample))
  #Celltypesubset = SetIdent(Celltypesubset, value = "integrated_snn_res.0.5")
  Celltypesubset = SetIdent(Celltypesubset, value = "celltype")
  
  alllevels <- levels(Celltypesubset@active.ident)
  allleveltable <- as.data.frame(table(Celltypesubset@active.ident, Celltypesubset@meta.data$orig.ident))
  
  
  
  for(j in 1:length(alllevels)){
    currentlevel <- alllevels[j]
    
    
    print(j)
    currentlevel <- alllevels[j]
    Celltypesubset2 <- subset(Celltypesubset, idents = currentlevel)
    DefaultAssay(Celltypesubset2) <- "RNA"
    Celltypesubset2 <- SetIdent(Celltypesubset2, value = "orig.ident")
    #set minimum cell representation per sample
    currentsample <- as.character(allleveltable$Var2[allleveltable$Var1 %in% currentlevel & allleveltable$Freq >= 3])
    availcontrolsample <- currentsample[currentsample %in% Controlsample]
    availfrucsample <- currentsample[currentsample %in% Frucsample]
    if(length(availcontrolsample) < 1 | length(availfrucsample) < 1){
      next
    }
    
    
    
    ###
    #1-comparison status
    if(length(availfrucsample) == 1 & length(availcontrolsample) == 1){
      #if only one sample per group, no need to replicate
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample, ident.2 =availcontrolsample , logfc.threshold = 0.05)
      
      finalDEG <- rownames(controltable[controltable$p_val_adj < 0.05,])
      #finalDEG <- rownames(controltable[controltable$p_val < 0.01,])
      methodused[[paste0("Set",i,"_",currentlevel)]] <- "Seurat_default"
      if(length(finalDEG) > 0){
        #indrib <- grep("^Rps|^Rpl|^Gm[0-9]|^mt-Rnr[0-9]",finalDEG)
        #if(length(indrib) > 0){
        #  finalDEG <- finalDEG[-indrib]
        #}
        DEGfulltable_single[[paste0("Set",i,"_",currentlevel)]] <- controltable
        DEGfulltableind[paste0("Set",i,"_",currentlevel)] <- "single"
        DEGlist2[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
        Humaninelist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["mt-Rnr2",]$avg_logFC,":",controltable["mt-Rnr2",]$p_val_adj)
        Pu1list[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Spi1",]$avg_logFC,":",controltable["Spi1",]$p_val_adj)
        Fcgrtlist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Fcgrt",]$avg_logFC,":",controltable["Fcgrt",]$p_val_adj)
        
        Khklist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Khk",]$avg_logFC,":",controltable["Khk",]$p_val_adj)
        Aldoblist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Aldob",]$avg_logFC,":",controltable["Aldob",]$p_val_adj)
        Glut5list[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Slc2a5",]$avg_logFC,":",controltable["Slc2a5",]$p_val_adj)
        
        DEGlist[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
        topDEG[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable)[1:5]
        DEGlist_up[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable[controltable$p_val_adj < 0.05 & controltable$avg_logFC > 0,])
        DEGlist_down[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable[controltable$p_val_adj < 0.05 & controltable$avg_logFC < 0,])
      }
      controltable$genename <- rownames(controltable)
      foldchangemerge <- controltable[,c("genename","avg_log2FC")]
      pvaluemerge <- controltable[,c("genename","p_val_adj")]
      colnames(foldchangemerge)[2] <- paste0("Set",i,"_",currentlevel)
      colnames(pvaluemerge)[2] <- paste0("Set",i,"_",currentlevel)
      
      if(!exists("foldchangemerge_final")){
        foldchangemerge_final <- foldchangemerge
        pvaluemerge_final <- pvaluemerge
      }else{
        foldchangemerge_final <- full_join(foldchangemerge_final,foldchangemerge, by = "genename")
        pvaluemerge_final <- full_join(pvaluemerge_final,pvaluemerge, by = "genename")
      }
      next
    }
    ###
    
    
    methodused[[paste0("Set",i,"_",currentlevel)]] <- "Dougmethod"
    #DEGlist2 with adjusted p value only
    controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample, ident.2 = availcontrolsample)
    finalDEG <- rownames(controltable[controltable$p_val_adj < 0.05,])
    if(length(finalDEG) > 0){
      DEGlist2[[paste0("Set",i,"_",currentlevel)]] <- finalDEG 
    }
    
    #For whole DEG across batches
    rm(finalptable,finalfoldtable)
    for(k in 1:length(availcontrolsample)){
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample,ident.2 = availcontrolsample[k], logfc.threshold = 0.05)
      controltable <- controltable[,c(1,2,5)]
      colnames(controltable) <- paste0(colnames(controltable),availcontrolsample[k])
      controltable$genename <- rownames(controltable)
      if(!exists("finalptable")){
        #finalptable <- controltable[,c(1,4)] #using raw pvalue
        finalptable <- controltable[,c(3,4)] #using adj pvalue
        finalfoldtable <- controltable[,c(2,4)]
      }else{
        #no NAs allowed, only consider full samples
        #finalptable <- full_join(finalptable,controltable[,c(1,4)])
        finalptable <- full_join(finalptable,controltable[,c(3,4)])
        finalfoldtable <- full_join(finalfoldtable, controltable[,c(2,4)])
      }
    }
    
    for(k in 1:length(availfrucsample)){
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample[k] ,ident.2 = availcontrolsample, logfc.threshold = 0.05)
      controltable <- controltable[,c(1,2,5)]
      colnames(controltable) <- paste0(colnames(controltable),availfrucsample[k])
      controltable$genename <- rownames(controltable)
      #no NAs allowed, only consider full samples
      #finalptable <- full_join(finalptable,controltable[,c(1,4)]) #using raw pvalue
      finalptable <- full_join(finalptable,controltable[,c(3,4)])
      finalfoldtable <- full_join(finalfoldtable, controltable[,c(2,4)])
    }
    
    rownames(finalptable) <- finalptable$genename
    finalptable <- finalptable[,!colnames(finalptable) %in% "genename"]
    allcompare <- ncol(finalptable)
    rownames(finalfoldtable) <- finalfoldtable$genename
    finalfoldtable <- finalfoldtable[,!colnames(finalfoldtable) %in% "genename"]
    
    finalptable$minP = apply(as.matrix(finalptable[,1:allcompare]), 1, minimump2)
    finalptable$maxP = apply(as.matrix(finalptable[,1:allcompare]), 1, maximump2)
    #finalptable$minPadj = p.adjust(finalptable$minP,method = 'fdr',n = nrow(Celltypesubset2@assays$RNA@data))
    #finalptable$maxPadj = p.adjust(finalptable$maxP,method = 'fdr',n = nrow(Celltypesubset2@assays$RNA@data))
    finalptable$minPadj = finalptable$minP
    finalptable$maxPadj = finalptable$maxP
    
    #finalptable$finalp <- apply(finalptable[,c("minPadj","maxPadj")],1,max)
    finalptable$finalp <- apply(finalptable[,c("minPadj","maxPadj")],1,min)
    #finalptable$FDR <- p.adjust(finalptable$finalp, method = "fdr")
    finalptable <- finalptable[!is.na(finalptable$finalp),]
    finalptable$genename <- rownames(finalptable)
    
    DEGfulltableind[paste0("Set",i,"_",currentlevel)] <- "multiple"
    DEGfulltable_multiple[[paste0("foldptable_Set",i,"_",currentlevel)]] <- finalfoldtable
    DEGfulltable_multiple[[paste0("finalptable_Set",i,"_",currentlevel)]] <- finalptable
    
    Humaninelist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["mt-Rnr2",]), na.rm = T),":",finalptable["mt-Rnr2",]$finalp)
    Pu1list[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Spi1",]), na.rm = T),":",finalptable["Spi1",]$finalp)
    Fcgrtlist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Fcgrt",]), na.rm = T),":",finalptable["Fcgrt",]$finalp)
    
    Khklist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Khk",]), na.rm = T),":",finalptable["Khk",]$finalp)
    Aldoblist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Aldob",]), na.rm = T),":",finalptable["Aldob",]$finalp)
    Glut5list[[paste0("Set",i,"_",currentlevel)]] <-paste0(median(as.numeric(finalfoldtable["Slc2a5",]), na.rm = T),":",finalptable["Slc2a5",]$finalp)
    
    #matching direction
    finalfoldtable$directionless <- apply(finalfoldtable[,1:allcompare],1,function(x){mean(x < 0)})
    finalfoldtable$directionmore <- apply(finalfoldtable[,1:allcompare],1,function(x){mean(x > 0)})
    finalfoldtable$avgfold <- apply(finalfoldtable[,1:allcompare],1,median,na.rm = T)
    finalfoldtable$genename <- rownames(finalfoldtable)
    
    foldchangemerge <- finalfoldtable[,c("genename","avgfold")]
    pvaluemerge <- finalptable[,c("genename","finalp")]
    colnames(foldchangemerge)[2] <- paste0("Set",i,"_",currentlevel)
    colnames(pvaluemerge)[2] <- paste0("Set",i,"_",currentlevel)
    
    if(!exists("foldchangemerge_final")){
      foldchangemerge_final <- foldchangemerge
      pvaluemerge_final <- pvaluemerge
    }else{
      foldchangemerge_final <- full_join(foldchangemerge_final,foldchangemerge, by = "genename")
      pvaluemerge_final <- full_join(pvaluemerge_final,pvaluemerge, by = "genename")
    }
    
    
    matchedgenes <- rownames(finalfoldtable[which(finalfoldtable$directionless >= 0.75 | finalfoldtable$directionmore >= 0.75),])
    upgenes <- rownames(finalfoldtable[which(finalfoldtable$directionmore >= 0.75) ,])
    downgenes <- rownames(finalfoldtable[which(finalfoldtable$directionless >= 0.75) ,])
    
    finalDEG <- intersect(matchedgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    #finalDEG <- intersect(matchedgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    #finalptable <- finalptable[finalptable$finalp < 0.05,]
    finalptable <- finalptable[order(finalptable$finalp),]
    #save(finalptable,file = paste0("Set",i,"_",currentlevel,"_Fruc.rda"))
    #save(finalptable,file = paste0("Set",i,"_",currentlevel,"_HFHS.rda"))
    finalptable <- finalptable[rownames(finalptable) %in% matchedgenes,]
    finalDEG_up <- intersect(upgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    finalDEG_down <- intersect(downgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    
    #finalDEG_up <- intersect(upgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    #finalDEG_down <- intersect(downgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    
    if(length(finalDEG) > 0){
      topDEG[[paste0("Set",i,"_",currentlevel)]] <- rownames(finalptable)[1:5]
      #exlude ribosomal from DEG list
      #indrib <- grep("^Rps|^Rpl|^Gm[0-9]|^mt-Rnr[0-9]",finalDEG)
      #if(length(indrib) > 0){
      #  finalDEG <- finalDEG[-indrib]
      #}
      DEGlist[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
      DEGlist_up[[paste0("Set",i,"_",currentlevel)]] <- finalDEG_up
      DEGlist_down[[paste0("Set",i,"_",currentlevel)]] <- finalDEG_down
      
    }
    
    
  }
  
}
save(DEGlist,DEGlist_up,DEGlist_down, Humaninelist,Khklist,Aldoblist,Glut5list,Fcgrtlist,Pu1list,DEGfulltable_multiple,DEGfulltable_single,DEGfulltableind,foldchangemerge_final,pvaluemerge_final,topDEG,file = "Humaninefold.rda")

rm(list=ls())
load("Humaninefold.rda")
library(enrichR)
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

enrichrfunction <- function(siggenes){
  if(length(siggenes) == 0){return(NA)}
  res <- enrichr(siggenes,databases = c("KEGG_2016","GO_Biological_Process_2017"))
  for(j in 1:length(res)){
    currentframe <- res[[j]]
    currentframe <- currentframe[currentframe$Adjusted.P.value < 0.05,]
    if(nrow(currentframe) > 0) currentframe <- currentframe[sapply(strsplit(currentframe$Genes,";"),length) > 3,]
    if(nrow(currentframe)  == 0){next
    }else{
      if(!exists("drugframe")){
        drugframe <- currentframe
      }else{
        drugframe <- rbind.data.frame(drugframe,currentframe)
      }
    }
  }
  if(exists("drugframe")){
    return(drugframe)
  }else{
    return(NULL)
  }
}

allcomparisons <- names(DEGlist)
siggenes_set1 <- NULL
siggenes_set1num <- NULL
siggenes_set2 <- NULL
siggenes_set2num<- NULL
siggenes_shared <- NULL
siggenes_sharednum <- NULL
sigpathways_set1 <- NULL
sigpathways_set1num <- NULL
sigpathways_set2 <- NULL
sigpathways_set2num <- NULL
sigpathways_sharedpath <- NULL
sigpathways_sharedpathnum <- NULL
sigpathways_sharedgenes <- NULL
sigpathways_sharedgenesnum <- NULL

alllevels  <- unique(gsub("Set[0-9]_","",allcomparisons))
foldchangevector <- NULL
pvaluevector <- NULL
genevector <- NULL
foldchangemerge_final$genename  <- firstup(foldchangemerge_final$genename)

for(i in 1:length(alllevels)){
  currentlevel <- alllevels[i]
  #matchingset <- allcomparisons[grep(paste0("_",currentlevel,"$"), allcomparisons)]
  matchingset <- allcomparisons[grep(paste0("_",currentlevel,"$"), allcomparisons)]
  if(length(matchingset) > 1){
    siggenes_set1[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
    
    tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
    if(!is.null(tempresult)){
      for(r in 1:nrow(tempresult)){
        name <- tempresult$Term[r]
        genevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Genes[r]
        foldchangevector[paste(name,currentlevel,"Set1",sep = "__")] <- 
          median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set1_",currentlevel)], na.rm = T)
        pvaluevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Adjusted.P.value[r]
      }
    }
    pathwayset1 <- tempresult$Term
    
    sigpathways_set1[i] <- paste(pathwayset1,collapse = ";")
    siggenes_set1num[i] <- length(DEGlist[[matchingset[1]]])
    sigpathways_set1num[i] <- length(pathwayset1)
    
    siggenes_set2[i] <- paste(DEGlist[[matchingset[2]]], collapse = " ")
    
    tempresult <-  enrichrfunction(DEGlist[[matchingset[2]]])
    if(!is.null(tempresult)){
      for(r in 1:nrow(tempresult)){
        name <- tempresult$Term[r]
        genevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Genes[r]
        foldchangevector[paste(name,currentlevel,"Set2",sep = "__")] <- 
          median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set2_",currentlevel)], na.rm = T)
        pvaluevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Adjusted.P.value[r]
      }
    }
    pathwayset2 <- tempresult$Term
    sigpathways_set2[i] <- paste(pathwayset2, collapse = ";")
    siggenes_set2num[i] <- length(DEGlist[[matchingset[2]]])
    sigpathways_set2num[i] <- length(pathwayset2)
    
    sigpathways_sharedpath[i] <- paste(intersect(pathwayset1,pathwayset2),collapse = " ")
    sigpathways_sharedpathnum[i] <- length(intersect(pathwayset1,pathwayset2))
    shared <- intersect(DEGlist[[matchingset[1]]],DEGlist[[matchingset[2]]])
    siggenes_shared[i] <- paste(shared,collapse = " ")
    siggenes_sharednum[i] <- length(shared)
    
    if(length(shared) > 0){
      DEGlist[[paste0("Shared",i,"_",currentlevel)]] <- shared
      sharedpathways <- enrichrfunction(shared)$Term
      sigpathways_sharedgenes[i] <- paste(sharedpathways,collapse = " ")
      sigpathways_sharedgenesnum[i] <- length(sharedpathways)
    }else{
      sigpathways_sharedgenes[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
    }
  }else if(length(matchingset) == 1){
    setname <- substr(matchingset,1,4)
    if(setname %in% "Set1"){
      siggenes_set1[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
      
      tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
      if(!is.null(tempresult)){
        for(r in 1:nrow(tempresult)){
          name <- tempresult$Term[r]
          genevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Genes[r]
          foldchangevector[paste(name,currentlevel,"Set1",sep = "__")] <- 
            median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set1_",currentlevel)], na.rm = T)
          pvaluevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Adjusted.P.value[r]
        }
      }
      pathwayset1 <- tempresult$Term
      
      sigpathways_set1[i] <- paste(pathwayset1,collapse = ";")
      siggenes_set1num[i] <- length(DEGlist[[matchingset[1]]])
      sigpathways_set1num[i] <- length(pathwayset1)
      
      siggenes_set2[i] <- NA
      sigpathways_set2[i] <- NA
      sigpathways_sharedpath[i] <- NA
      siggenes_set2num[i] <- NA
      sigpathways_set2num[i] <- NA
      
      siggenes_shared[i] <- NA
      sigpathways_sharedgenes[i] <- NA
      siggenes_sharednum[i] <- NA
      sigpathways_sharedpathnum[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
    }else if(setname %in% "Set2"){
      siggenes_set1[i] <- NA
      sigpathways_set1[i] <- NA
      siggenes_set1num[i] <- NA
      sigpathways_set1num[i] <- NA
      
      siggenes_set2[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
      
      tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
      if(!is.null(tempresult)){
        for(r in 1:nrow(tempresult)){
          name <- tempresult$Term[r]
          genevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Genes[r]
          foldchangevector[paste(name,currentlevel,"Set2",sep = "__")] <- 
            median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set2_",currentlevel)], na.rm = T)
          pvaluevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Adjusted.P.value[r]
        }
      }
      
      
      pathwayset2 <- tempresult$Term
      sigpathways_set2[i] <- paste(pathwayset2, collapse = ";")
      siggenes_set2num[i] <- length(DEGlist[[matchingset[1]]])
      sigpathways_set2num[i] <- length(pathwayset2)
      
      sigpathways_sharedpath[i] <- NA 
      
      siggenes_shared[i] <- NA
      sigpathways_sharedgenes[i] <- NA
      siggenes_sharednum[i] <- NA
      sigpathways_sharedpathnum[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
      
    }
  }else{
    siggenes_set1[i] <- NA
    sigpathways_set1[i] <- NA
    siggenes_set1num[i] <- NA
    sigpathways_set1num[i] <- NA
    
    siggenes_set2[i] <- NA
    sigpathways_set2[i] <-NA
    siggenes_set2num[i] <- NA
    sigpathways_set2num[i] <- NA
    
    sigpathways_sharedpath[i] <- NA 
    siggenes_shared[i] <- NA
    sigpathways_sharedgenes[i] <- NA
    siggenes_sharednum[i] <- NA
    sigpathways_sharedpathnum[i] <- NA
    sigpathways_sharedgenesnum[i] <- NA
  }
  
  
}
topDEG2 = sapply(topDEG,paste,collapse = ", ")
topDEG2 <- topDEG2[paste0("Set1_",alllevels)]
topDEG2 <- topDEG2[paste0("Set2_",alllevels)]

sigpathways_set1 <- gsub("Homo sapiens_","",sigpathways_set1)
sigpathways_set2 <- gsub("Homo sapiens_","",sigpathways_set2)
sigpathways_sharedpath <- gsub("Homo sapiens_","",sigpathways_sharedpath)
sigpathways_sharedgenes <- gsub("Homo sapiens_","",sigpathways_sharedgenes)

finalframe <- data.frame(names = alllevels,alllevels,siggenes_set1, siggenes_set1num, sigpathways_set1, sigpathways_set1num,
                         siggenes_set2, siggenes_set2num, sigpathways_set2, sigpathways_set2num,
                         siggenes_shared,siggenes_sharednum,sigpathways_sharedpath,sigpathways_sharedpathnum,sigpathways_sharedgenes,sigpathways_sharedgenesnum)

identical( names(foldchangevector),names(pvaluevector))
result <- data.frame(name = names(foldchangevector),fold = foldchangevector,FDR = -log10(pvaluevector), genevector)

splitname <- strsplit(result$name,"__")
result$pathwayname <- sapply(splitname, `[[`, 1)
result$celltype <- sapply(splitname, `[[`, 2)
result$Set <- sapply(splitname, `[[`, 3)
result$set_celltype <- paste0(result$celltype,"__",result$Set)

result_final <- result %>% group_by(set_celltype) %>% top_n(n = 5, wt = FDR)
save(result,result_final, file= "pathwayinfo_Fruc.rda")
write.csv(finalframe,file = "cross_batch_DEG_FDR_Doug_Fruc.csv")


setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/NPC_liver.rda")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "PN_Hep","1" = "SEC","2" = "Kupffer","3" = "NKT", "4" = "CV_Hep","5" = "Macrophages", "6" = "PN_Hep",
#                                                    "7" = "Hepatic stellate","8" = "B cells","9" = "Dividing cells","10" = "SEC_Dnase1l3","11" = "Plasmacytoid dendritic cells","12" = "Cholangiocytes"))
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/liver.rda")
minimump2 <- function(x){
  NAmean <- mean(is.na(x))
  if(NAmean > 0.25){return(NA)
  }else{
    return(minimump(na.omit(x))$p)
  }
}
maximump2 <- function(x){
  NAmean <- mean(is.na(x))
  if(NAmean > 0.25){return(NA)
  }else{
    return(maximump(na.omit(x))$p)
  }
}


dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")
allsamples <- unique(dropEST.combined.filtered@meta.data$orig.ident)
sampleinfo <- substr(allsamples,1,4)
treatmentinfo <- substr(allsamples,5,5)
treatmentinfo[treatmentinfo %in% c("1","2")] <- "Control"
treatmentinfo[treatmentinfo %in% c("3","4")] <- "Fruc"
treatmentinfo[treatmentinfo %in% c("5","6")] <- "HFHS"
treatmentinfotable <- data.frame(allsamples,sampleinfo,treatmentinfo, stringsAsFactors = F)

#For small intestine where Set3 is available
#if(sum(treatmentinfotable$sampleinfo %in% "Set3") > 0 ){
#  treatmentinfotable$treatmentinfo[grep("Set32",treatmentinfotable$allsamples)] <- "Fruc"
#}
#


DEGlist <- list()
DEGlist2 <- list()
topDEG <- list()
Humaninelist <- list()
Pu1list <- list()
Fcgrtlist <- list()

DEGlist_up <- list()
DEGlist_down <- list()
methodused <- list()
DEGlist3 <- list()
Khklist <- list()
Aldoblist <- list()
Glut5list <- list()
DEGfulltable_single <- list()
DEGfulltable_multiple <- list()
DEGfulltableind <- c()
rm(pvaluemerge_final,foldchangemerge_final)
for(i in 1:2){
  currenttreatment <- treatmentinfotable[treatmentinfotable$sampleinfo %in% paste0("Set",i),]
  Controlsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "Control"]
  Frucsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "Fruc"]
  #Frucsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "HFHS"]
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "orig.ident")
  Celltypesubset <- subset(dropEST.combined.filtered, idents = c(Controlsample,Frucsample))
  #Celltypesubset = SetIdent(Celltypesubset, value = "integrated_snn_res.0.5")
  Celltypesubset = SetIdent(Celltypesubset, value = "celltype")
  
  alllevels <- levels(Celltypesubset@active.ident)
  allleveltable <- as.data.frame(table(Celltypesubset@active.ident, Celltypesubset@meta.data$orig.ident))
  
  
  
  for(j in 1:length(alllevels)){
    currentlevel <- alllevels[j]
    
    
    print(j)
    currentlevel <- alllevels[j]
    Celltypesubset2 <- subset(Celltypesubset, idents = currentlevel)
    DefaultAssay(Celltypesubset2) <- "RNA"
    Celltypesubset2 <- SetIdent(Celltypesubset2, value = "orig.ident")
    #set minimum cell representation per sample
    currentsample <- as.character(allleveltable$Var2[allleveltable$Var1 %in% currentlevel & allleveltable$Freq >= 3])
    availcontrolsample <- currentsample[currentsample %in% Controlsample]
    availfrucsample <- currentsample[currentsample %in% Frucsample]
    if(length(availcontrolsample) < 1 | length(availfrucsample) < 1){
      next
    }
    
    
    
    ###
    #1-comparison status
    if(length(availfrucsample) == 1 & length(availcontrolsample) == 1){
      #if only one sample per group, no need to replicate
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample, ident.2 =availcontrolsample , logfc.threshold = 0.05)
      
      finalDEG <- rownames(controltable[controltable$p_val_adj < 0.05,])
      #finalDEG <- rownames(controltable[controltable$p_val < 0.01,])
      methodused[[paste0("Set",i,"_",currentlevel)]] <- "Seurat_default"
      if(length(finalDEG) > 0){
        #indrib <- grep("^Rps|^Rpl|^Gm[0-9]|^mt-Rnr[0-9]",finalDEG)
        #if(length(indrib) > 0){
        #  finalDEG <- finalDEG[-indrib]
        #}
        DEGfulltable_single[[paste0("Set",i,"_",currentlevel)]] <- controltable
        DEGfulltableind[paste0("Set",i,"_",currentlevel)] <- "single"
        DEGlist2[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
        Humaninelist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["mt-Rnr2",]$avg_logFC,":",controltable["mt-Rnr2",]$p_val_adj)
        Pu1list[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Spi1",]$avg_logFC,":",controltable["Spi1",]$p_val_adj)
        Fcgrtlist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Fcgrt",]$avg_logFC,":",controltable["Fcgrt",]$p_val_adj)
        
        Khklist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Khk",]$avg_logFC,":",controltable["Khk",]$p_val_adj)
        Aldoblist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Aldob",]$avg_logFC,":",controltable["Aldob",]$p_val_adj)
        Glut5list[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Slc2a5",]$avg_logFC,":",controltable["Slc2a5",]$p_val_adj)
        
        DEGlist[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
        topDEG[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable)[1:5]
        DEGlist_up[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable[controltable$p_val_adj < 0.05 & controltable$avg_logFC > 0,])
        DEGlist_down[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable[controltable$p_val_adj < 0.05 & controltable$avg_logFC < 0,])
      }
      controltable$genename <- rownames(controltable)
      foldchangemerge <- controltable[,c("genename","avg_log2FC")]
      pvaluemerge <- controltable[,c("genename","p_val_adj")]
      colnames(foldchangemerge)[2] <- paste0("Set",i,"_",currentlevel)
      colnames(pvaluemerge)[2] <- paste0("Set",i,"_",currentlevel)
      
      if(!exists("foldchangemerge_final")){
        foldchangemerge_final <- foldchangemerge
        pvaluemerge_final <- pvaluemerge
      }else{
        foldchangemerge_final <- full_join(foldchangemerge_final,foldchangemerge, by = "genename")
        pvaluemerge_final <- full_join(pvaluemerge_final,pvaluemerge, by = "genename")
      }
      next
    }
    ###
    
    
    methodused[[paste0("Set",i,"_",currentlevel)]] <- "Dougmethod"
    #DEGlist2 with adjusted p value only
    controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample, ident.2 = availcontrolsample)
    finalDEG <- rownames(controltable[controltable$p_val_adj < 0.05,])
    if(length(finalDEG) > 0){
      DEGlist2[[paste0("Set",i,"_",currentlevel)]] <- finalDEG 
    }
    
    #For whole DEG across batches
    rm(finalptable,finalfoldtable)
    for(k in 1:length(availcontrolsample)){
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample,ident.2 = availcontrolsample[k], logfc.threshold = 0.05)
      controltable <- controltable[,c(1,2,5)]
      colnames(controltable) <- paste0(colnames(controltable),availcontrolsample[k])
      controltable$genename <- rownames(controltable)
      if(!exists("finalptable")){
        #finalptable <- controltable[,c(1,4)] #using raw pvalue
        finalptable <- controltable[,c(3,4)] #using adj pvalue
        finalfoldtable <- controltable[,c(2,4)]
      }else{
        #no NAs allowed, only consider full samples
        #finalptable <- full_join(finalptable,controltable[,c(1,4)])
        finalptable <- full_join(finalptable,controltable[,c(3,4)])
        finalfoldtable <- full_join(finalfoldtable, controltable[,c(2,4)])
      }
    }
    
    for(k in 1:length(availfrucsample)){
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample[k] ,ident.2 = availcontrolsample, logfc.threshold = 0.05)
      controltable <- controltable[,c(1,2,5)]
      colnames(controltable) <- paste0(colnames(controltable),availfrucsample[k])
      controltable$genename <- rownames(controltable)
      #no NAs allowed, only consider full samples
      #finalptable <- full_join(finalptable,controltable[,c(1,4)]) #using raw pvalue
      finalptable <- full_join(finalptable,controltable[,c(3,4)])
      finalfoldtable <- full_join(finalfoldtable, controltable[,c(2,4)])
    }
    
    rownames(finalptable) <- finalptable$genename
    finalptable <- finalptable[,!colnames(finalptable) %in% "genename"]
    allcompare <- ncol(finalptable)
    rownames(finalfoldtable) <- finalfoldtable$genename
    finalfoldtable <- finalfoldtable[,!colnames(finalfoldtable) %in% "genename"]
    
    finalptable$minP = apply(as.matrix(finalptable[,1:allcompare]), 1, minimump2)
    finalptable$maxP = apply(as.matrix(finalptable[,1:allcompare]), 1, maximump2)
    #finalptable$minPadj = p.adjust(finalptable$minP,method = 'fdr',n = nrow(Celltypesubset2@assays$RNA@data))
    #finalptable$maxPadj = p.adjust(finalptable$maxP,method = 'fdr',n = nrow(Celltypesubset2@assays$RNA@data))
    finalptable$minPadj = finalptable$minP
    finalptable$maxPadj = finalptable$maxP
    
    #finalptable$finalp <- apply(finalptable[,c("minPadj","maxPadj")],1,max)
    finalptable$finalp <- apply(finalptable[,c("minPadj","maxPadj")],1,min)
    #finalptable$FDR <- p.adjust(finalptable$finalp, method = "fdr")
    finalptable <- finalptable[!is.na(finalptable$finalp),]
    finalptable$genename <- rownames(finalptable)
    
    DEGfulltableind[paste0("Set",i,"_",currentlevel)] <- "multiple"
    DEGfulltable_multiple[[paste0("foldptable_Set",i,"_",currentlevel)]] <- finalfoldtable
    DEGfulltable_multiple[[paste0("finalptable_Set",i,"_",currentlevel)]] <- finalptable
    
    Humaninelist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["mt-Rnr2",]), na.rm = T),":",finalptable["mt-Rnr2",]$finalp)
    Pu1list[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Spi1",]), na.rm = T),":",finalptable["Spi1",]$finalp)
    Fcgrtlist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Fcgrt",]), na.rm = T),":",finalptable["Fcgrt",]$finalp)
    
    Khklist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Khk",]), na.rm = T),":",finalptable["Khk",]$finalp)
    Aldoblist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Aldob",]), na.rm = T),":",finalptable["Aldob",]$finalp)
    Glut5list[[paste0("Set",i,"_",currentlevel)]] <-paste0(median(as.numeric(finalfoldtable["Slc2a5",]), na.rm = T),":",finalptable["Slc2a5",]$finalp)
    
    #matching direction
    finalfoldtable$directionless <- apply(finalfoldtable[,1:allcompare],1,function(x){mean(x < 0)})
    finalfoldtable$directionmore <- apply(finalfoldtable[,1:allcompare],1,function(x){mean(x > 0)})
    finalfoldtable$avgfold <- apply(finalfoldtable[,1:allcompare],1,median,na.rm = T)
    finalfoldtable$genename <- rownames(finalfoldtable)
    
    foldchangemerge <- finalfoldtable[,c("genename","avgfold")]
    pvaluemerge <- finalptable[,c("genename","finalp")]
    colnames(foldchangemerge)[2] <- paste0("Set",i,"_",currentlevel)
    colnames(pvaluemerge)[2] <- paste0("Set",i,"_",currentlevel)
    
    if(!exists("foldchangemerge_final")){
      foldchangemerge_final <- foldchangemerge
      pvaluemerge_final <- pvaluemerge
    }else{
      foldchangemerge_final <- full_join(foldchangemerge_final,foldchangemerge, by = "genename")
      pvaluemerge_final <- full_join(pvaluemerge_final,pvaluemerge, by = "genename")
    }
    
    
    matchedgenes <- rownames(finalfoldtable[which(finalfoldtable$directionless >= 0.75 | finalfoldtable$directionmore >= 0.75),])
    upgenes <- rownames(finalfoldtable[which(finalfoldtable$directionmore >= 0.75) ,])
    downgenes <- rownames(finalfoldtable[which(finalfoldtable$directionless >= 0.75) ,])
    
    finalDEG <- intersect(matchedgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    #finalDEG <- intersect(matchedgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    #finalptable <- finalptable[finalptable$finalp < 0.05,]
    finalptable <- finalptable[order(finalptable$finalp),]
    #save(finalptable,file = paste0("Set",i,"_",currentlevel,"_Fruc.rda"))
    #save(finalptable,file = paste0("Set",i,"_",currentlevel,"_HFHS.rda"))
    finalptable <- finalptable[rownames(finalptable) %in% matchedgenes,]
    finalDEG_up <- intersect(upgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    finalDEG_down <- intersect(downgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    
    #finalDEG_up <- intersect(upgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    #finalDEG_down <- intersect(downgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    
    if(length(finalDEG) > 0){
      topDEG[[paste0("Set",i,"_",currentlevel)]] <- rownames(finalptable)[1:5]
      #exlude ribosomal from DEG list
      #indrib <- grep("^Rps|^Rpl|^Gm[0-9]|^mt-Rnr[0-9]",finalDEG)
      #if(length(indrib) > 0){
      #  finalDEG <- finalDEG[-indrib]
      #}
      DEGlist[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
      DEGlist_up[[paste0("Set",i,"_",currentlevel)]] <- finalDEG_up
      DEGlist_down[[paste0("Set",i,"_",currentlevel)]] <- finalDEG_down
      
    }
    
    
  }
  
}
save(DEGlist,DEGlist_up,DEGlist_down, Humaninelist,Khklist,Aldoblist,Glut5list,Fcgrtlist,Pu1list,DEGfulltable_multiple,DEGfulltable_single,DEGfulltableind,foldchangemerge_final,pvaluemerge_final,topDEG,file = "Humaninefold.rda")

rm(list=ls())
load("Humaninefold.rda")
library(enrichR)
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
enrichrfunction <- function(siggenes){
  if(length(siggenes) == 0){return(NA)}
  res <- enrichr(siggenes,databases = c("KEGG_2016","GO_Biological_Process_2017"))
  for(j in 1:length(res)){
    currentframe <- res[[j]]
    currentframe <- currentframe[currentframe$Adjusted.P.value < 0.05,]
    if(nrow(currentframe) > 0) currentframe <- currentframe[sapply(strsplit(currentframe$Genes,";"),length) > 3,]
    if(nrow(currentframe)  == 0){next
    }else{
      if(!exists("drugframe")){
        drugframe <- currentframe
      }else{
        drugframe <- rbind.data.frame(drugframe,currentframe)
      }
    }
  }
  if(exists("drugframe")){
    return(drugframe)
  }else{
    return(NULL)
  }
}

allcomparisons <- names(DEGlist)
siggenes_set1 <- NULL
siggenes_set1num <- NULL
siggenes_set2 <- NULL
siggenes_set2num<- NULL
siggenes_shared <- NULL
siggenes_sharednum <- NULL
sigpathways_set1 <- NULL
sigpathways_set1num <- NULL
sigpathways_set2 <- NULL
sigpathways_set2num <- NULL
sigpathways_sharedpath <- NULL
sigpathways_sharedpathnum <- NULL
sigpathways_sharedgenes <- NULL
sigpathways_sharedgenesnum <- NULL

alllevels  <- unique(gsub("Set[0-9]_","",allcomparisons))
foldchangevector <- NULL
pvaluevector <- NULL
genevector <- NULL
foldchangemerge_final$genename  <- firstup(foldchangemerge_final$genename)

for(i in 1:length(alllevels)){
  currentlevel <- alllevels[i]
  #matchingset <- allcomparisons[grep(paste0("_",currentlevel,"$"), allcomparisons)]
  matchingset <- allcomparisons[grep(paste0("_",currentlevel,"$"), allcomparisons)]
  if(length(matchingset) > 1){
    siggenes_set1[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
    
    tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
    if(!is.null(tempresult)){
      for(r in 1:nrow(tempresult)){
        name <- tempresult$Term[r]
        genevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Genes[r]
        foldchangevector[paste(name,currentlevel,"Set1",sep = "__")] <- 
          median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set1_",currentlevel)], na.rm = T)
        pvaluevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Adjusted.P.value[r]
      }
    }
    pathwayset1 <- tempresult$Term
    
    sigpathways_set1[i] <- paste(pathwayset1,collapse = ";")
    siggenes_set1num[i] <- length(DEGlist[[matchingset[1]]])
    sigpathways_set1num[i] <- length(pathwayset1)
    
    siggenes_set2[i] <- paste(DEGlist[[matchingset[2]]], collapse = " ")
    
    tempresult <-  enrichrfunction(DEGlist[[matchingset[2]]])
    if(!is.null(tempresult)){
      for(r in 1:nrow(tempresult)){
        name <- tempresult$Term[r]
        genevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Genes[r]
        foldchangevector[paste(name,currentlevel,"Set2",sep = "__")] <- 
          median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set2_",currentlevel)], na.rm = T)
        pvaluevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Adjusted.P.value[r]
      }
    }
    pathwayset2 <- tempresult$Term
    sigpathways_set2[i] <- paste(pathwayset2, collapse = ";")
    siggenes_set2num[i] <- length(DEGlist[[matchingset[2]]])
    sigpathways_set2num[i] <- length(pathwayset2)
    
    sigpathways_sharedpath[i] <- paste(intersect(pathwayset1,pathwayset2),collapse = " ")
    sigpathways_sharedpathnum[i] <- length(intersect(pathwayset1,pathwayset2))
    shared <- intersect(DEGlist[[matchingset[1]]],DEGlist[[matchingset[2]]])
    siggenes_shared[i] <- paste(shared,collapse = " ")
    siggenes_sharednum[i] <- length(shared)
    
    if(length(shared) > 0){
      DEGlist[[paste0("Shared",i,"_",currentlevel)]] <- shared
      sharedpathways <- enrichrfunction(shared)$Term
      sigpathways_sharedgenes[i] <- paste(sharedpathways,collapse = " ")
      sigpathways_sharedgenesnum[i] <- length(sharedpathways)
    }else{
      sigpathways_sharedgenes[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
    }
  }else if(length(matchingset) == 1){
    setname <- substr(matchingset,1,4)
    if(setname %in% "Set1"){
      siggenes_set1[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
      
      tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
      if(!is.null(tempresult)){
        for(r in 1:nrow(tempresult)){
          name <- tempresult$Term[r]
          genevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Genes[r]
          foldchangevector[paste(name,currentlevel,"Set1",sep = "__")] <- 
            median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set1_",currentlevel)], na.rm = T)
          pvaluevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Adjusted.P.value[r]
        }
      }
      pathwayset1 <- tempresult$Term
      
      sigpathways_set1[i] <- paste(pathwayset1,collapse = ";")
      siggenes_set1num[i] <- length(DEGlist[[matchingset[1]]])
      sigpathways_set1num[i] <- length(pathwayset1)
      
      siggenes_set2[i] <- NA
      sigpathways_set2[i] <- NA
      sigpathways_sharedpath[i] <- NA
      siggenes_set2num[i] <- NA
      sigpathways_set2num[i] <- NA
      
      siggenes_shared[i] <- NA
      sigpathways_sharedgenes[i] <- NA
      siggenes_sharednum[i] <- NA
      sigpathways_sharedpathnum[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
    }else if(setname %in% "Set2"){
      siggenes_set1[i] <- NA
      sigpathways_set1[i] <- NA
      siggenes_set1num[i] <- NA
      sigpathways_set1num[i] <- NA
      
      siggenes_set2[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
      
      tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
      if(!is.null(tempresult)){
        for(r in 1:nrow(tempresult)){
          name <- tempresult$Term[r]
          genevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Genes[r]
          foldchangevector[paste(name,currentlevel,"Set2",sep = "__")] <- 
            median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set2_",currentlevel)], na.rm = T)
          pvaluevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Adjusted.P.value[r]
        }
      }
      
      
      pathwayset2 <- tempresult$Term
      sigpathways_set2[i] <- paste(pathwayset2, collapse = ";")
      siggenes_set2num[i] <- length(DEGlist[[matchingset[1]]])
      sigpathways_set2num[i] <- length(pathwayset2)
      
      sigpathways_sharedpath[i] <- NA 
      
      siggenes_shared[i] <- NA
      sigpathways_sharedgenes[i] <- NA
      siggenes_sharednum[i] <- NA
      sigpathways_sharedpathnum[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
      
    }
  }else{
    siggenes_set1[i] <- NA
    sigpathways_set1[i] <- NA
    siggenes_set1num[i] <- NA
    sigpathways_set1num[i] <- NA
    
    siggenes_set2[i] <- NA
    sigpathways_set2[i] <-NA
    siggenes_set2num[i] <- NA
    sigpathways_set2num[i] <- NA
    
    sigpathways_sharedpath[i] <- NA 
    siggenes_shared[i] <- NA
    sigpathways_sharedgenes[i] <- NA
    siggenes_sharednum[i] <- NA
    sigpathways_sharedpathnum[i] <- NA
    sigpathways_sharedgenesnum[i] <- NA
  }
  
  
}
topDEG2 = sapply(topDEG,paste,collapse = ", ")
topDEG2 <- topDEG2[paste0("Set1_",alllevels)]
topDEG2 <- topDEG2[paste0("Set2_",alllevels)]

sigpathways_set1 <- gsub("Homo sapiens_","",sigpathways_set1)
sigpathways_set2 <- gsub("Homo sapiens_","",sigpathways_set2)
sigpathways_sharedpath <- gsub("Homo sapiens_","",sigpathways_sharedpath)
sigpathways_sharedgenes <- gsub("Homo sapiens_","",sigpathways_sharedgenes)

finalframe <- data.frame(names = alllevels,alllevels,siggenes_set1, siggenes_set1num, sigpathways_set1, sigpathways_set1num,
                         siggenes_set2, siggenes_set2num, sigpathways_set2, sigpathways_set2num,
                         siggenes_shared,siggenes_sharednum,sigpathways_sharedpath,sigpathways_sharedpathnum,sigpathways_sharedgenes,sigpathways_sharedgenesnum)

identical( names(foldchangevector),names(pvaluevector))
result <- data.frame(name = names(foldchangevector),fold = foldchangevector,FDR = -log10(pvaluevector), genevector)

splitname <- strsplit(result$name,"__")
result$pathwayname <- sapply(splitname, `[[`, 1)
result$celltype <- sapply(splitname, `[[`, 2)
result$Set <- sapply(splitname, `[[`, 3)
result$set_celltype <- paste0(result$celltype,"__",result$Set)

result_final <- result %>% group_by(set_celltype) %>% top_n(n = 5, wt = FDR)
save(result,result_final, file= "pathwayinfo_Fruc.rda")
write.csv(finalframe,file = "cross_batch_DEG_FDR_Doug_Fruc.csv")




rm(list=ls())
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
#load("SI_Set1_DropEST.Object1.rda")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "Astrocytes","1" = "Myelinating oligodendrocyte","2" = "Gabaergic neurons","3" = "Endothelial cells", "4" = "Oligodendrocyte precursor","5" = "General neurons", "6" = "Glutamatergic neurons",
#                                                    "7" = "Ependymal cells","8" = "Gabaergic neurons","9" = "Microglial cells","10" = "Microglial cells_Nfkbiz","11" = "Pericytes","12" = "Myelinating oligodendrocyte",
#                                                    "13" = "Tanycytes","14" = "Muralpericytes_Acta2","15"="Glutamatergic neurons","16"="Glutamatergic neurons","17" = "Pericyes_Kcnj8","18" = "Vascular Leptomenigeal cells", "19" = "Gabaergic neurons",
#                                                    "20" = "Pars tuberalis","21" = "Macrophages","22" = "Endothelial cells","23" = "Newly formed oligodendrocyte","24" = "Myelinating oligodendrocyte"))
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Hyp_combined_Set1_final.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/Hyp.rda")
minimump2 <- function(x){
  NAmean <- mean(is.na(x))
  if(NAmean > 0.25){return(NA)
  }else{
    return(minimump(na.omit(x))$p)
  }
}
maximump2 <- function(x){
  NAmean <- mean(is.na(x))
  if(NAmean > 0.25){return(NA)
  }else{
    return(maximump(na.omit(x))$p)
  }
}


dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")
allsamples <- unique(dropEST.combined.filtered@meta.data$orig.ident)
sampleinfo <- substr(allsamples,1,4)
treatmentinfo <- substr(allsamples,5,5)
treatmentinfo[treatmentinfo %in% c("1","2")] <- "Control"
treatmentinfo[treatmentinfo %in% c("3","4")] <- "Fruc"
treatmentinfo[treatmentinfo %in% c("5","6")] <- "HFHS"
treatmentinfotable <- data.frame(allsamples,sampleinfo,treatmentinfo, stringsAsFactors = F)

#For small intestine where Set3 is available
#if(sum(treatmentinfotable$sampleinfo %in% "Set3") > 0 ){
#  treatmentinfotable$treatmentinfo[grep("Set32",treatmentinfotable$allsamples)] <- "Fruc"
#}
#


DEGlist <- list()
DEGlist2 <- list()
topDEG <- list()
Humaninelist <- list()
Pu1list <- list()
Fcgrtlist <- list()

DEGlist_up <- list()
DEGlist_down <- list()
methodused <- list()
DEGlist3 <- list()
Khklist <- list()
Aldoblist <- list()
Glut5list <- list()
DEGfulltable_single <- list()
DEGfulltable_multiple <- list()
DEGfulltableind <- c()
rm(pvaluemerge_final,foldchangemerge_final)
for(i in 1:1){
  currenttreatment <- treatmentinfotable[treatmentinfotable$sampleinfo %in% paste0("Set",i),]
  Controlsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "Control"]
  Frucsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "Fruc"]
  #Frucsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "HFHS"]
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "orig.ident")
  Celltypesubset <- subset(dropEST.combined.filtered, idents = c(Controlsample,Frucsample))
  #Celltypesubset = SetIdent(Celltypesubset, value = "integrated_snn_res.0.5")
  Celltypesubset = SetIdent(Celltypesubset, value = "celltype")
  
  alllevels <- levels(Celltypesubset@active.ident)
  allleveltable <- as.data.frame(table(Celltypesubset@active.ident, Celltypesubset@meta.data$orig.ident))
  
  
  
  for(j in 1:length(alllevels)){
    currentlevel <- alllevels[j]
    
    
    print(j)
    currentlevel <- alllevels[j]
    Celltypesubset2 <- subset(Celltypesubset, idents = currentlevel)
    DefaultAssay(Celltypesubset2) <- "RNA"
    Celltypesubset2 <- SetIdent(Celltypesubset2, value = "orig.ident")
    #set minimum cell representation per sample
    currentsample <- as.character(allleveltable$Var2[allleveltable$Var1 %in% currentlevel & allleveltable$Freq >= 3])
    availcontrolsample <- currentsample[currentsample %in% Controlsample]
    availfrucsample <- currentsample[currentsample %in% Frucsample]
    if(length(availcontrolsample) < 1 | length(availfrucsample) < 1){
      next
    }
    
    
    
    ###
    #1-comparison status
    if(length(availfrucsample) == 1 & length(availcontrolsample) == 1){
      #if only one sample per group, no need to replicate
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample, ident.2 =availcontrolsample , logfc.threshold = 0.05)
      
      finalDEG <- rownames(controltable[controltable$p_val_adj < 0.05,])
      #finalDEG <- rownames(controltable[controltable$p_val < 0.01,])
      methodused[[paste0("Set",i,"_",currentlevel)]] <- "Seurat_default"
      if(length(finalDEG) > 0){
        #indrib <- grep("^Rps|^Rpl|^Gm[0-9]|^mt-Rnr[0-9]",finalDEG)
        #if(length(indrib) > 0){
        #  finalDEG <- finalDEG[-indrib]
        #}
        DEGfulltable_single[[paste0("Set",i,"_",currentlevel)]] <- controltable
        DEGfulltableind[paste0("Set",i,"_",currentlevel)] <- "single"
        DEGlist2[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
        Humaninelist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["mt-Rnr2",]$avg_logFC,":",controltable["mt-Rnr2",]$p_val_adj)
        Pu1list[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Spi1",]$avg_logFC,":",controltable["Spi1",]$p_val_adj)
        Fcgrtlist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Fcgrt",]$avg_logFC,":",controltable["Fcgrt",]$p_val_adj)
        
        Khklist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Khk",]$avg_logFC,":",controltable["Khk",]$p_val_adj)
        Aldoblist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Aldob",]$avg_logFC,":",controltable["Aldob",]$p_val_adj)
        Glut5list[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Slc2a5",]$avg_logFC,":",controltable["Slc2a5",]$p_val_adj)
        
        DEGlist[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
        topDEG[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable)[1:5]
        DEGlist_up[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable[controltable$p_val_adj < 0.05 & controltable$avg_logFC > 0,])
        DEGlist_down[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable[controltable$p_val_adj < 0.05 & controltable$avg_logFC < 0,])
      }
      controltable$genename <- rownames(controltable)
      foldchangemerge <- controltable[,c("genename","avg_log2FC")]
      pvaluemerge <- controltable[,c("genename","p_val_adj")]
      colnames(foldchangemerge)[2] <- paste0("Set",i,"_",currentlevel)
      colnames(pvaluemerge)[2] <- paste0("Set",i,"_",currentlevel)
      
      if(!exists("foldchangemerge_final")){
        foldchangemerge_final <- foldchangemerge
        pvaluemerge_final <- pvaluemerge
      }else{
        foldchangemerge_final <- full_join(foldchangemerge_final,foldchangemerge, by = "genename")
        pvaluemerge_final <- full_join(pvaluemerge_final,pvaluemerge, by = "genename")
      }
      next
    }
    ###
    
    
    methodused[[paste0("Set",i,"_",currentlevel)]] <- "Dougmethod"
    #DEGlist2 with adjusted p value only
    controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample, ident.2 = availcontrolsample)
    finalDEG <- rownames(controltable[controltable$p_val_adj < 0.05,])
    if(length(finalDEG) > 0){
      DEGlist2[[paste0("Set",i,"_",currentlevel)]] <- finalDEG 
    }
    
    #For whole DEG across batches
    rm(finalptable,finalfoldtable)
    for(k in 1:length(availcontrolsample)){
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample,ident.2 = availcontrolsample[k], logfc.threshold = 0.05)
      controltable <- controltable[,c(1,2,5)]
      colnames(controltable) <- paste0(colnames(controltable),availcontrolsample[k])
      controltable$genename <- rownames(controltable)
      if(!exists("finalptable")){
        #finalptable <- controltable[,c(1,4)] #using raw pvalue
        finalptable <- controltable[,c(3,4)] #using adj pvalue
        finalfoldtable <- controltable[,c(2,4)]
      }else{
        #no NAs allowed, only consider full samples
        #finalptable <- full_join(finalptable,controltable[,c(1,4)])
        finalptable <- full_join(finalptable,controltable[,c(3,4)])
        finalfoldtable <- full_join(finalfoldtable, controltable[,c(2,4)])
      }
    }
    
    for(k in 1:length(availfrucsample)){
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample[k] ,ident.2 = availcontrolsample, logfc.threshold = 0.05)
      controltable <- controltable[,c(1,2,5)]
      colnames(controltable) <- paste0(colnames(controltable),availfrucsample[k])
      controltable$genename <- rownames(controltable)
      #no NAs allowed, only consider full samples
      #finalptable <- full_join(finalptable,controltable[,c(1,4)]) #using raw pvalue
      finalptable <- full_join(finalptable,controltable[,c(3,4)])
      finalfoldtable <- full_join(finalfoldtable, controltable[,c(2,4)])
    }
    
    rownames(finalptable) <- finalptable$genename
    finalptable <- finalptable[,!colnames(finalptable) %in% "genename"]
    allcompare <- ncol(finalptable)
    rownames(finalfoldtable) <- finalfoldtable$genename
    finalfoldtable <- finalfoldtable[,!colnames(finalfoldtable) %in% "genename"]
    
    finalptable$minP = apply(as.matrix(finalptable[,1:allcompare]), 1, minimump2)
    finalptable$maxP = apply(as.matrix(finalptable[,1:allcompare]), 1, maximump2)
    #finalptable$minPadj = p.adjust(finalptable$minP,method = 'fdr',n = nrow(Celltypesubset2@assays$RNA@data))
    #finalptable$maxPadj = p.adjust(finalptable$maxP,method = 'fdr',n = nrow(Celltypesubset2@assays$RNA@data))
    finalptable$minPadj = finalptable$minP
    finalptable$maxPadj = finalptable$maxP
    
    #finalptable$finalp <- apply(finalptable[,c("minPadj","maxPadj")],1,max)
    finalptable$finalp <- apply(finalptable[,c("minPadj","maxPadj")],1,min)
    #finalptable$FDR <- p.adjust(finalptable$finalp, method = "fdr")
    finalptable <- finalptable[!is.na(finalptable$finalp),]
    finalptable$genename <- rownames(finalptable)
    
    DEGfulltableind[paste0("Set",i,"_",currentlevel)] <- "multiple"
    DEGfulltable_multiple[[paste0("foldptable_Set",i,"_",currentlevel)]] <- finalfoldtable
    DEGfulltable_multiple[[paste0("finalptable_Set",i,"_",currentlevel)]] <- finalptable
    
    Humaninelist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["mt-Rnr2",]), na.rm = T),":",finalptable["mt-Rnr2",]$finalp)
    Pu1list[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Spi1",]), na.rm = T),":",finalptable["Spi1",]$finalp)
    Fcgrtlist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Fcgrt",]), na.rm = T),":",finalptable["Fcgrt",]$finalp)
    
    Khklist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Khk",]), na.rm = T),":",finalptable["Khk",]$finalp)
    Aldoblist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Aldob",]), na.rm = T),":",finalptable["Aldob",]$finalp)
    Glut5list[[paste0("Set",i,"_",currentlevel)]] <-paste0(median(as.numeric(finalfoldtable["Slc2a5",]), na.rm = T),":",finalptable["Slc2a5",]$finalp)
    
    #matching direction
    finalfoldtable$directionless <- apply(finalfoldtable[,1:allcompare],1,function(x){mean(x < 0)})
    finalfoldtable$directionmore <- apply(finalfoldtable[,1:allcompare],1,function(x){mean(x > 0)})
    finalfoldtable$avgfold <- apply(finalfoldtable[,1:allcompare],1,median,na.rm = T)
    finalfoldtable$genename <- rownames(finalfoldtable)
    
    foldchangemerge <- finalfoldtable[,c("genename","avgfold")]
    pvaluemerge <- finalptable[,c("genename","finalp")]
    colnames(foldchangemerge)[2] <- paste0("Set",i,"_",currentlevel)
    colnames(pvaluemerge)[2] <- paste0("Set",i,"_",currentlevel)
    
    if(!exists("foldchangemerge_final")){
      foldchangemerge_final <- foldchangemerge
      pvaluemerge_final <- pvaluemerge
    }else{
      foldchangemerge_final <- full_join(foldchangemerge_final,foldchangemerge, by = "genename")
      pvaluemerge_final <- full_join(pvaluemerge_final,pvaluemerge, by = "genename")
    }
    
    
    matchedgenes <- rownames(finalfoldtable[which(finalfoldtable$directionless >= 0.75 | finalfoldtable$directionmore >= 0.75),])
    upgenes <- rownames(finalfoldtable[which(finalfoldtable$directionmore >= 0.75) ,])
    downgenes <- rownames(finalfoldtable[which(finalfoldtable$directionless >= 0.75) ,])
    
    finalDEG <- intersect(matchedgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    #finalDEG <- intersect(matchedgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    #finalptable <- finalptable[finalptable$finalp < 0.05,]
    finalptable <- finalptable[order(finalptable$finalp),]
    #save(finalptable,file = paste0("Set",i,"_",currentlevel,"_Fruc.rda"))
    #save(finalptable,file = paste0("Set",i,"_",currentlevel,"_HFHS.rda"))
    finalptable <- finalptable[rownames(finalptable) %in% matchedgenes,]
    finalDEG_up <- intersect(upgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    finalDEG_down <- intersect(downgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    
    #finalDEG_up <- intersect(upgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    #finalDEG_down <- intersect(downgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    
    if(length(finalDEG) > 0){
      topDEG[[paste0("Set",i,"_",currentlevel)]] <- rownames(finalptable)[1:5]
      #exlude ribosomal from DEG list
      #indrib <- grep("^Rps|^Rpl|^Gm[0-9]|^mt-Rnr[0-9]",finalDEG)
      #if(length(indrib) > 0){
      #  finalDEG <- finalDEG[-indrib]
      #}
      DEGlist[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
      DEGlist_up[[paste0("Set",i,"_",currentlevel)]] <- finalDEG_up
      DEGlist_down[[paste0("Set",i,"_",currentlevel)]] <- finalDEG_down
      
    }
    
    
  }
  
}
save(DEGlist,DEGlist_up,DEGlist_down, Humaninelist,Khklist,Aldoblist,Glut5list,Fcgrtlist,Pu1list,DEGfulltable_multiple,DEGfulltable_single,DEGfulltableind,foldchangemerge_final,pvaluemerge_final,topDEG,file = "Humaninefold.rda")

rm(list=ls())
load("Humaninefold.rda")
library(enrichR)
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
enrichrfunction <- function(siggenes){
  if(length(siggenes) == 0){return(NA)}
  res <- enrichr(siggenes,databases = c("KEGG_2016","GO_Biological_Process_2017"))
  for(j in 1:length(res)){
    currentframe <- res[[j]]
    currentframe <- currentframe[currentframe$Adjusted.P.value < 0.05,]
    if(nrow(currentframe) > 0) currentframe <- currentframe[sapply(strsplit(currentframe$Genes,";"),length) > 3,]
    if(nrow(currentframe)  == 0){next
    }else{
      if(!exists("drugframe")){
        drugframe <- currentframe
      }else{
        drugframe <- rbind.data.frame(drugframe,currentframe)
      }
    }
  }
  if(exists("drugframe")){
    return(drugframe)
  }else{
    return(NULL)
  }
}

allcomparisons <- names(DEGlist)
siggenes_set1 <- NULL
siggenes_set1num <- NULL
siggenes_set2 <- NULL
siggenes_set2num<- NULL
siggenes_shared <- NULL
siggenes_sharednum <- NULL
sigpathways_set1 <- NULL
sigpathways_set1num <- NULL
sigpathways_set2 <- NULL
sigpathways_set2num <- NULL
sigpathways_sharedpath <- NULL
sigpathways_sharedpathnum <- NULL
sigpathways_sharedgenes <- NULL
sigpathways_sharedgenesnum <- NULL

alllevels  <- unique(gsub("Set[0-9]_","",allcomparisons))
foldchangevector <- NULL
pvaluevector <- NULL
genevector <- NULL
foldchangemerge_final$genename  <- firstup(foldchangemerge_final$genename)

for(i in 1:length(alllevels)){
  currentlevel <- alllevels[i]
  #matchingset <- allcomparisons[grep(paste0("_",currentlevel,"$"), allcomparisons)]
  matchingset <- allcomparisons[grep(paste0("_",currentlevel,"$"), allcomparisons)]
  if(length(matchingset) > 1){
    siggenes_set1[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
    
    tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
    if(!is.null(tempresult)){
      for(r in 1:nrow(tempresult)){
        name <- tempresult$Term[r]
        genevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Genes[r]
        foldchangevector[paste(name,currentlevel,"Set1",sep = "__")] <- 
          median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set1_",currentlevel)], na.rm = T)
        pvaluevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Adjusted.P.value[r]
      }
    }
    pathwayset1 <- tempresult$Term
    
    sigpathways_set1[i] <- paste(pathwayset1,collapse = ";")
    siggenes_set1num[i] <- length(DEGlist[[matchingset[1]]])
    sigpathways_set1num[i] <- length(pathwayset1)
    
    siggenes_set2[i] <- paste(DEGlist[[matchingset[2]]], collapse = " ")
    
    tempresult <-  enrichrfunction(DEGlist[[matchingset[2]]])
    if(!is.null(tempresult)){
      for(r in 1:nrow(tempresult)){
        name <- tempresult$Term[r]
        genevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Genes[r]
        foldchangevector[paste(name,currentlevel,"Set2",sep = "__")] <- 
          median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set2_",currentlevel)], na.rm = T)
        pvaluevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Adjusted.P.value[r]
      }
    }
    pathwayset2 <- tempresult$Term
    sigpathways_set2[i] <- paste(pathwayset2, collapse = ";")
    siggenes_set2num[i] <- length(DEGlist[[matchingset[2]]])
    sigpathways_set2num[i] <- length(pathwayset2)
    
    sigpathways_sharedpath[i] <- paste(intersect(pathwayset1,pathwayset2),collapse = " ")
    sigpathways_sharedpathnum[i] <- length(intersect(pathwayset1,pathwayset2))
    shared <- intersect(DEGlist[[matchingset[1]]],DEGlist[[matchingset[2]]])
    siggenes_shared[i] <- paste(shared,collapse = " ")
    siggenes_sharednum[i] <- length(shared)
    
    if(length(shared) > 0){
      DEGlist[[paste0("Shared",i,"_",currentlevel)]] <- shared
      sharedpathways <- enrichrfunction(shared)$Term
      sigpathways_sharedgenes[i] <- paste(sharedpathways,collapse = " ")
      sigpathways_sharedgenesnum[i] <- length(sharedpathways)
    }else{
      sigpathways_sharedgenes[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
    }
  }else if(length(matchingset) == 1){
    setname <- substr(matchingset,1,4)
    if(setname %in% "Set1"){
      siggenes_set1[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
      
      tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
      if(!is.null(tempresult)){
        for(r in 1:nrow(tempresult)){
          name <- tempresult$Term[r]
          genevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Genes[r]
          foldchangevector[paste(name,currentlevel,"Set1",sep = "__")] <- 
            median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set1_",currentlevel)], na.rm = T)
          pvaluevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Adjusted.P.value[r]
        }
      }
      pathwayset1 <- tempresult$Term
      
      sigpathways_set1[i] <- paste(pathwayset1,collapse = ";")
      siggenes_set1num[i] <- length(DEGlist[[matchingset[1]]])
      sigpathways_set1num[i] <- length(pathwayset1)
      
      siggenes_set2[i] <- NA
      sigpathways_set2[i] <- NA
      sigpathways_sharedpath[i] <- NA
      siggenes_set2num[i] <- NA
      sigpathways_set2num[i] <- NA
      
      siggenes_shared[i] <- NA
      sigpathways_sharedgenes[i] <- NA
      siggenes_sharednum[i] <- NA
      sigpathways_sharedpathnum[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
    }else if(setname %in% "Set2"){
      siggenes_set1[i] <- NA
      sigpathways_set1[i] <- NA
      siggenes_set1num[i] <- NA
      sigpathways_set1num[i] <- NA
      
      siggenes_set2[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
      
      tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
      if(!is.null(tempresult)){
        for(r in 1:nrow(tempresult)){
          name <- tempresult$Term[r]
          genevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Genes[r]
          foldchangevector[paste(name,currentlevel,"Set2",sep = "__")] <- 
            median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set2_",currentlevel)], na.rm = T)
          pvaluevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Adjusted.P.value[r]
        }
      }
      
      
      pathwayset2 <- tempresult$Term
      sigpathways_set2[i] <- paste(pathwayset2, collapse = ";")
      siggenes_set2num[i] <- length(DEGlist[[matchingset[1]]])
      sigpathways_set2num[i] <- length(pathwayset2)
      
      sigpathways_sharedpath[i] <- NA 
      
      siggenes_shared[i] <- NA
      sigpathways_sharedgenes[i] <- NA
      siggenes_sharednum[i] <- NA
      sigpathways_sharedpathnum[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
      
    }
  }else{
    siggenes_set1[i] <- NA
    sigpathways_set1[i] <- NA
    siggenes_set1num[i] <- NA
    sigpathways_set1num[i] <- NA
    
    siggenes_set2[i] <- NA
    sigpathways_set2[i] <-NA
    siggenes_set2num[i] <- NA
    sigpathways_set2num[i] <- NA
    
    sigpathways_sharedpath[i] <- NA 
    siggenes_shared[i] <- NA
    sigpathways_sharedgenes[i] <- NA
    siggenes_sharednum[i] <- NA
    sigpathways_sharedpathnum[i] <- NA
    sigpathways_sharedgenesnum[i] <- NA
  }
  
  
}
topDEG2 = sapply(topDEG,paste,collapse = ", ")
topDEG2 <- topDEG2[paste0("Set1_",alllevels)]
topDEG2 <- topDEG2[paste0("Set2_",alllevels)]

sigpathways_set1 <- gsub("Homo sapiens_","",sigpathways_set1)
sigpathways_set2 <- gsub("Homo sapiens_","",sigpathways_set2)
sigpathways_sharedpath <- gsub("Homo sapiens_","",sigpathways_sharedpath)
sigpathways_sharedgenes <- gsub("Homo sapiens_","",sigpathways_sharedgenes)

finalframe <- data.frame(names = alllevels,alllevels,siggenes_set1, siggenes_set1num, sigpathways_set1, sigpathways_set1num,
                         siggenes_set2, siggenes_set2num, sigpathways_set2, sigpathways_set2num,
                         siggenes_shared,siggenes_sharednum,sigpathways_sharedpath,sigpathways_sharedpathnum,sigpathways_sharedgenes,sigpathways_sharedgenesnum)

identical( names(foldchangevector),names(pvaluevector))
result <- data.frame(name = names(foldchangevector),fold = foldchangevector,FDR = -log10(pvaluevector), genevector)

splitname <- strsplit(result$name,"__")
result$pathwayname <- sapply(splitname, `[[`, 1)
result$celltype <- sapply(splitname, `[[`, 2)
result$Set <- sapply(splitname, `[[`, 3)
result$set_celltype <- paste0(result$celltype,"__",result$Set)

result_final <- result %>% group_by(set_celltype) %>% top_n(n = 5, wt = FDR)
save(result,result_final, file= "pathwayinfo_Fruc.rda")
write.csv(finalframe,file = "cross_batch_DEG_FDR_Doug_Fruc.csv")



#Hyp neurons
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Neuroncelltype.rda")
#dropEST.combined.filtered <- CellTypeSubset
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.1.5,
#                                                  c("0" = "General neurons","1" = "General neurons","2" = "General neurons","3" = "Trh_GABAergic", 
#                                                    "4" = "General neurons","5" = "Grp_Glutamatergic","6" = "Prok2_GABAergic","7" = "Agrp_GABAergic",
#                                                    "8" = "Pomc_Glutamatergic","9" = "Sst_Bcl11b_GABAergic","10" = "Prlr_Glutamatergic","11" = "Sim1_Glutamatergic","12" = "Ghrh_GABAergic",
#                                                    "13" = "Crabp1_GABAergic","14" = "Kl_GABAergic","15" = "Avp_Glutamatergic","16" = "Lhx1_GABAergic",
#                                                    "17" = "Lhx6_GABAergic","18" = "Prph_Histaminergic","19" = "Cx3cr1_GABAergic","20" = "Ndnf_Glutamatergic"))
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/Hyp_neurons.rda")
DefaultAssay(dropEST.combined.filtered) <- "RNA"
dropEST.combined.filtered = SetIdent(dropEST.combined.filtered,value="batch")
dropEST.combined.filtered = subset(dropEST.combined.filtered,idents="Set1")

minimump2 <- function(x){
  NAmean <- mean(is.na(x))
  if(NAmean > 0.25){return(NA)
  }else{
    return(minimump(na.omit(x))$p)
  }
}
maximump2 <- function(x){
  NAmean <- mean(is.na(x))
  if(NAmean > 0.25){return(NA)
  }else{
    return(maximump(na.omit(x))$p)
  }
}


dropEST.combined.filtered = SetIdent(dropEST.combined.filtered, value = "integrated_snn_res.0.5")
allsamples <- unique(dropEST.combined.filtered@meta.data$orig.ident)
sampleinfo <- substr(allsamples,1,4)
treatmentinfo <- substr(allsamples,5,5)
treatmentinfo[treatmentinfo %in% c("1","2")] <- "Control"
treatmentinfo[treatmentinfo %in% c("3","4")] <- "Fruc"
treatmentinfo[treatmentinfo %in% c("5","6")] <- "HFHS"
treatmentinfotable <- data.frame(allsamples,sampleinfo,treatmentinfo, stringsAsFactors = F)

#For small intestine where Set3 is available
#if(sum(treatmentinfotable$sampleinfo %in% "Set3") > 0 ){
#  treatmentinfotable$treatmentinfo[grep("Set32",treatmentinfotable$allsamples)] <- "Fruc"
#}
#


DEGlist <- list()
DEGlist2 <- list()
topDEG <- list()
Humaninelist <- list()
Pu1list <- list()
Fcgrtlist <- list()

DEGlist_up <- list()
DEGlist_down <- list()
methodused <- list()
DEGlist3 <- list()
Khklist <- list()
Aldoblist <- list()
Glut5list <- list()
DEGfulltable_single <- list()
DEGfulltable_multiple <- list()
DEGfulltableind <- c()
rm(pvaluemerge_final,foldchangemerge_final)
for(i in 1:1){
  currenttreatment <- treatmentinfotable[treatmentinfotable$sampleinfo %in% paste0("Set",i),]
  Controlsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "Control"]
  Frucsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "Fruc"]
  #Frucsample <- currenttreatment$allsamples[currenttreatment$treatmentinfo %in% "HFHS"]
  dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "orig.ident")
  Celltypesubset <- subset(dropEST.combined.filtered, idents = c(Controlsample,Frucsample))
  #Celltypesubset = SetIdent(Celltypesubset, value = "integrated_snn_res.0.5")
  Celltypesubset = SetIdent(Celltypesubset, value = "celltype")
  
  alllevels <- levels(Celltypesubset@active.ident)
  allleveltable <- as.data.frame(table(Celltypesubset@active.ident, Celltypesubset@meta.data$orig.ident))
  
  
  
  for(j in 1:length(alllevels)){
    currentlevel <- alllevels[j]
    
    
    print(j)
    currentlevel <- alllevels[j]
    Celltypesubset2 <- subset(Celltypesubset, idents = currentlevel)
    DefaultAssay(Celltypesubset2) <- "RNA"
    Celltypesubset2 <- SetIdent(Celltypesubset2, value = "orig.ident")
    #set minimum cell representation per sample
    currentsample <- as.character(allleveltable$Var2[allleveltable$Var1 %in% currentlevel & allleveltable$Freq >= 3])
    availcontrolsample <- currentsample[currentsample %in% Controlsample]
    availfrucsample <- currentsample[currentsample %in% Frucsample]
    if(length(availcontrolsample) < 1 | length(availfrucsample) < 1){
      next
    }
    
    
    
    ###
    #1-comparison status
    if(length(availfrucsample) == 1 & length(availcontrolsample) == 1){
      #if only one sample per group, no need to replicate
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample, ident.2 =availcontrolsample , logfc.threshold = 0.05)
      
      finalDEG <- rownames(controltable[controltable$p_val_adj < 0.05,])
      #finalDEG <- rownames(controltable[controltable$p_val < 0.01,])
      methodused[[paste0("Set",i,"_",currentlevel)]] <- "Seurat_default"
      if(length(finalDEG) > 0){
        #indrib <- grep("^Rps|^Rpl|^Gm[0-9]|^mt-Rnr[0-9]",finalDEG)
        #if(length(indrib) > 0){
        #  finalDEG <- finalDEG[-indrib]
        #}
        DEGfulltable_single[[paste0("Set",i,"_",currentlevel)]] <- controltable
        DEGfulltableind[paste0("Set",i,"_",currentlevel)] <- "single"
        DEGlist2[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
        Humaninelist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["mt-Rnr2",]$avg_logFC,":",controltable["mt-Rnr2",]$p_val_adj)
        Pu1list[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Spi1",]$avg_logFC,":",controltable["Spi1",]$p_val_adj)
        Fcgrtlist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Fcgrt",]$avg_logFC,":",controltable["Fcgrt",]$p_val_adj)
        
        Khklist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Khk",]$avg_logFC,":",controltable["Khk",]$p_val_adj)
        Aldoblist[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Aldob",]$avg_logFC,":",controltable["Aldob",]$p_val_adj)
        Glut5list[[paste0("Set",i,"_",currentlevel)]] <- paste0(controltable["Slc2a5",]$avg_logFC,":",controltable["Slc2a5",]$p_val_adj)
        
        DEGlist[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
        topDEG[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable)[1:5]
        DEGlist_up[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable[controltable$p_val_adj < 0.05 & controltable$avg_logFC > 0,])
        DEGlist_down[[paste0("Set",i,"_",currentlevel)]] <- rownames(controltable[controltable$p_val_adj < 0.05 & controltable$avg_logFC < 0,])
      }
      controltable$genename <- rownames(controltable)
      foldchangemerge <- controltable[,c("genename","avg_log2FC")]
      pvaluemerge <- controltable[,c("genename","p_val_adj")]
      colnames(foldchangemerge)[2] <- paste0("Set",i,"_",currentlevel)
      colnames(pvaluemerge)[2] <- paste0("Set",i,"_",currentlevel)
      
      if(!exists("foldchangemerge_final")){
        foldchangemerge_final <- foldchangemerge
        pvaluemerge_final <- pvaluemerge
      }else{
        foldchangemerge_final <- full_join(foldchangemerge_final,foldchangemerge, by = "genename")
        pvaluemerge_final <- full_join(pvaluemerge_final,pvaluemerge, by = "genename")
      }
      next
    }
    ###
    
    
    methodused[[paste0("Set",i,"_",currentlevel)]] <- "Dougmethod"
    #DEGlist2 with adjusted p value only
    controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample, ident.2 = availcontrolsample)
    finalDEG <- rownames(controltable[controltable$p_val_adj < 0.05,])
    if(length(finalDEG) > 0){
      DEGlist2[[paste0("Set",i,"_",currentlevel)]] <- finalDEG 
    }
    
    #For whole DEG across batches
    rm(finalptable,finalfoldtable)
    for(k in 1:length(availcontrolsample)){
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample,ident.2 = availcontrolsample[k], logfc.threshold = 0.05)
      controltable <- controltable[,c(1,2,5)]
      colnames(controltable) <- paste0(colnames(controltable),availcontrolsample[k])
      controltable$genename <- rownames(controltable)
      if(!exists("finalptable")){
        #finalptable <- controltable[,c(1,4)] #using raw pvalue
        finalptable <- controltable[,c(3,4)] #using adj pvalue
        finalfoldtable <- controltable[,c(2,4)]
      }else{
        #no NAs allowed, only consider full samples
        #finalptable <- full_join(finalptable,controltable[,c(1,4)])
        finalptable <- full_join(finalptable,controltable[,c(3,4)])
        finalfoldtable <- full_join(finalfoldtable, controltable[,c(2,4)])
      }
    }
    
    for(k in 1:length(availfrucsample)){
      controltable <- FindMarkers(Celltypesubset2, ident.1 = availfrucsample[k] ,ident.2 = availcontrolsample, logfc.threshold = 0.05)
      controltable <- controltable[,c(1,2,5)]
      colnames(controltable) <- paste0(colnames(controltable),availfrucsample[k])
      controltable$genename <- rownames(controltable)
      #no NAs allowed, only consider full samples
      #finalptable <- full_join(finalptable,controltable[,c(1,4)]) #using raw pvalue
      finalptable <- full_join(finalptable,controltable[,c(3,4)])
      finalfoldtable <- full_join(finalfoldtable, controltable[,c(2,4)])
    }
    
    rownames(finalptable) <- finalptable$genename
    finalptable <- finalptable[,!colnames(finalptable) %in% "genename"]
    allcompare <- ncol(finalptable)
    rownames(finalfoldtable) <- finalfoldtable$genename
    finalfoldtable <- finalfoldtable[,!colnames(finalfoldtable) %in% "genename"]
    
    finalptable$minP = apply(as.matrix(finalptable[,1:allcompare]), 1, minimump2)
    finalptable$maxP = apply(as.matrix(finalptable[,1:allcompare]), 1, maximump2)
    #finalptable$minPadj = p.adjust(finalptable$minP,method = 'fdr',n = nrow(Celltypesubset2@assays$RNA@data))
    #finalptable$maxPadj = p.adjust(finalptable$maxP,method = 'fdr',n = nrow(Celltypesubset2@assays$RNA@data))
    finalptable$minPadj = finalptable$minP
    finalptable$maxPadj = finalptable$maxP
    
    #finalptable$finalp <- apply(finalptable[,c("minPadj","maxPadj")],1,max)
    finalptable$finalp <- apply(finalptable[,c("minPadj","maxPadj")],1,min)
    #finalptable$FDR <- p.adjust(finalptable$finalp, method = "fdr")
    finalptable <- finalptable[!is.na(finalptable$finalp),]
    finalptable$genename <- rownames(finalptable)
    
    DEGfulltableind[paste0("Set",i,"_",currentlevel)] <- "multiple"
    DEGfulltable_multiple[[paste0("foldptable_Set",i,"_",currentlevel)]] <- finalfoldtable
    DEGfulltable_multiple[[paste0("finalptable_Set",i,"_",currentlevel)]] <- finalptable
    
    Humaninelist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["mt-Rnr2",]), na.rm = T),":",finalptable["mt-Rnr2",]$finalp)
    Pu1list[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Spi1",]), na.rm = T),":",finalptable["Spi1",]$finalp)
    Fcgrtlist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Fcgrt",]), na.rm = T),":",finalptable["Fcgrt",]$finalp)
    
    Khklist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Khk",]), na.rm = T),":",finalptable["Khk",]$finalp)
    Aldoblist[[paste0("Set",i,"_",currentlevel)]] <- paste0(median(as.numeric(finalfoldtable["Aldob",]), na.rm = T),":",finalptable["Aldob",]$finalp)
    Glut5list[[paste0("Set",i,"_",currentlevel)]] <-paste0(median(as.numeric(finalfoldtable["Slc2a5",]), na.rm = T),":",finalptable["Slc2a5",]$finalp)
    
    #matching direction
    finalfoldtable$directionless <- apply(finalfoldtable[,1:allcompare],1,function(x){mean(x < 0)})
    finalfoldtable$directionmore <- apply(finalfoldtable[,1:allcompare],1,function(x){mean(x > 0)})
    finalfoldtable$avgfold <- apply(finalfoldtable[,1:allcompare],1,median,na.rm = T)
    finalfoldtable$genename <- rownames(finalfoldtable)
    
    foldchangemerge <- finalfoldtable[,c("genename","avgfold")]
    pvaluemerge <- finalptable[,c("genename","finalp")]
    colnames(foldchangemerge)[2] <- paste0("Set",i,"_",currentlevel)
    colnames(pvaluemerge)[2] <- paste0("Set",i,"_",currentlevel)
    
    if(!exists("foldchangemerge_final")){
      foldchangemerge_final <- foldchangemerge
      pvaluemerge_final <- pvaluemerge
    }else{
      foldchangemerge_final <- full_join(foldchangemerge_final,foldchangemerge, by = "genename")
      pvaluemerge_final <- full_join(pvaluemerge_final,pvaluemerge, by = "genename")
    }
    
    
    matchedgenes <- rownames(finalfoldtable[which(finalfoldtable$directionless >= 0.75 | finalfoldtable$directionmore >= 0.75),])
    upgenes <- rownames(finalfoldtable[which(finalfoldtable$directionmore >= 0.75) ,])
    downgenes <- rownames(finalfoldtable[which(finalfoldtable$directionless >= 0.75) ,])
    
    finalDEG <- intersect(matchedgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    #finalDEG <- intersect(matchedgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    #finalptable <- finalptable[finalptable$finalp < 0.05,]
    finalptable <- finalptable[order(finalptable$finalp),]
    #save(finalptable,file = paste0("Set",i,"_",currentlevel,"_Fruc.rda"))
    #save(finalptable,file = paste0("Set",i,"_",currentlevel,"_HFHS.rda"))
    finalptable <- finalptable[rownames(finalptable) %in% matchedgenes,]
    finalDEG_up <- intersect(upgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    finalDEG_down <- intersect(downgenes,rownames(finalptable[which(finalptable$finalp < 0.05),]))
    
    #finalDEG_up <- intersect(upgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    #finalDEG_down <- intersect(downgenes,rownames(finalptable[finalptable$finalp < 0.01,]))
    
    if(length(finalDEG) > 0){
      topDEG[[paste0("Set",i,"_",currentlevel)]] <- rownames(finalptable)[1:5]
      #exlude ribosomal from DEG list
      #indrib <- grep("^Rps|^Rpl|^Gm[0-9]|^mt-Rnr[0-9]",finalDEG)
      #if(length(indrib) > 0){
      #  finalDEG <- finalDEG[-indrib]
      #}
      DEGlist[[paste0("Set",i,"_",currentlevel)]] <- finalDEG
      DEGlist_up[[paste0("Set",i,"_",currentlevel)]] <- finalDEG_up
      DEGlist_down[[paste0("Set",i,"_",currentlevel)]] <- finalDEG_down
      
    }
    
    
  }
  
}
save(DEGlist,DEGlist_up,DEGlist_down, Humaninelist,Khklist,Aldoblist,Glut5list,Fcgrtlist,Pu1list,DEGfulltable_multiple,DEGfulltable_single,DEGfulltableind,foldchangemerge_final,pvaluemerge_final,topDEG,file = "Humaninefold_subset.rda")

rm(list=ls())
load("Humaninefold_subset.rda")
library(enrichR)
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
enrichrfunction <- function(siggenes){
  if(length(siggenes) == 0){return(NA)}
  res <- enrichr(siggenes,databases = c("KEGG_2016","GO_Biological_Process_2017"))
  for(j in 1:length(res)){
    currentframe <- res[[j]]
    currentframe <- currentframe[currentframe$Adjusted.P.value < 0.05,]
    if(nrow(currentframe) > 0) currentframe <- currentframe[sapply(strsplit(currentframe$Genes,";"),length) > 3,]
    if(nrow(currentframe)  == 0){next
    }else{
      if(!exists("drugframe")){
        drugframe <- currentframe
      }else{
        drugframe <- rbind.data.frame(drugframe,currentframe)
      }
    }
  }
  if(exists("drugframe")){
    return(drugframe)
  }else{
    return(NULL)
  }
}

allcomparisons <- names(DEGlist)
siggenes_set1 <- NULL
siggenes_set1num <- NULL
siggenes_set2 <- NULL
siggenes_set2num<- NULL
siggenes_shared <- NULL
siggenes_sharednum <- NULL
sigpathways_set1 <- NULL
sigpathways_set1num <- NULL
sigpathways_set2 <- NULL
sigpathways_set2num <- NULL
sigpathways_sharedpath <- NULL
sigpathways_sharedpathnum <- NULL
sigpathways_sharedgenes <- NULL
sigpathways_sharedgenesnum <- NULL

alllevels  <- unique(gsub("Set[0-9]_","",allcomparisons))
foldchangevector <- NULL
pvaluevector <- NULL
genevector <- NULL
foldchangemerge_final$genename  <- firstup(foldchangemerge_final$genename)

for(i in 1:length(alllevels)){
  currentlevel <- alllevels[i]
  #matchingset <- allcomparisons[grep(paste0("_",currentlevel,"$"), allcomparisons)]
  matchingset <- allcomparisons[grep(paste0("_",currentlevel,"$"), allcomparisons)]
  if(length(matchingset) > 1){
    siggenes_set1[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
    
    tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
    if(!is.null(tempresult)){
      for(r in 1:nrow(tempresult)){
        name <- tempresult$Term[r]
        genevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Genes[r]
        foldchangevector[paste(name,currentlevel,"Set1",sep = "__")] <- 
          median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set1_",currentlevel)], na.rm = T)
        pvaluevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Adjusted.P.value[r]
      }
    }
    pathwayset1 <- tempresult$Term
    
    sigpathways_set1[i] <- paste(pathwayset1,collapse = ";")
    siggenes_set1num[i] <- length(DEGlist[[matchingset[1]]])
    sigpathways_set1num[i] <- length(pathwayset1)
    
    siggenes_set2[i] <- paste(DEGlist[[matchingset[2]]], collapse = " ")
    
    tempresult <-  enrichrfunction(DEGlist[[matchingset[2]]])
    if(!is.null(tempresult)){
      for(r in 1:nrow(tempresult)){
        name <- tempresult$Term[r]
        genevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Genes[r]
        foldchangevector[paste(name,currentlevel,"Set2",sep = "__")] <- 
          median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set2_",currentlevel)], na.rm = T)
        pvaluevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Adjusted.P.value[r]
      }
    }
    pathwayset2 <- tempresult$Term
    sigpathways_set2[i] <- paste(pathwayset2, collapse = ";")
    siggenes_set2num[i] <- length(DEGlist[[matchingset[2]]])
    sigpathways_set2num[i] <- length(pathwayset2)
    
    sigpathways_sharedpath[i] <- paste(intersect(pathwayset1,pathwayset2),collapse = " ")
    sigpathways_sharedpathnum[i] <- length(intersect(pathwayset1,pathwayset2))
    shared <- intersect(DEGlist[[matchingset[1]]],DEGlist[[matchingset[2]]])
    siggenes_shared[i] <- paste(shared,collapse = " ")
    siggenes_sharednum[i] <- length(shared)
    
    if(length(shared) > 0){
      DEGlist[[paste0("Shared",i,"_",currentlevel)]] <- shared
      sharedpathways <- enrichrfunction(shared)$Term
      sigpathways_sharedgenes[i] <- paste(sharedpathways,collapse = " ")
      sigpathways_sharedgenesnum[i] <- length(sharedpathways)
    }else{
      sigpathways_sharedgenes[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
    }
  }else if(length(matchingset) == 1){
    setname <- substr(matchingset,1,4)
    if(setname %in% "Set1"){
      siggenes_set1[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
      
      tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
      if(!is.null(tempresult)){
        for(r in 1:nrow(tempresult)){
          name <- tempresult$Term[r]
          genevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Genes[r]
          foldchangevector[paste(name,currentlevel,"Set1",sep = "__")] <- 
            median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set1_",currentlevel)], na.rm = T)
          pvaluevector[paste(name,currentlevel,"Set1",sep = "__")] <- tempresult$Adjusted.P.value[r]
        }
      }
      pathwayset1 <- tempresult$Term
      
      sigpathways_set1[i] <- paste(pathwayset1,collapse = ";")
      siggenes_set1num[i] <- length(DEGlist[[matchingset[1]]])
      sigpathways_set1num[i] <- length(pathwayset1)
      
      siggenes_set2[i] <- NA
      sigpathways_set2[i] <- NA
      sigpathways_sharedpath[i] <- NA
      siggenes_set2num[i] <- NA
      sigpathways_set2num[i] <- NA
      
      siggenes_shared[i] <- NA
      sigpathways_sharedgenes[i] <- NA
      siggenes_sharednum[i] <- NA
      sigpathways_sharedpathnum[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
    }else if(setname %in% "Set2"){
      siggenes_set1[i] <- NA
      sigpathways_set1[i] <- NA
      siggenes_set1num[i] <- NA
      sigpathways_set1num[i] <- NA
      
      siggenes_set2[i] <- paste(DEGlist[[matchingset[1]]], collapse = " ")
      
      tempresult <-  enrichrfunction(DEGlist[[matchingset[1]]])
      if(!is.null(tempresult)){
        for(r in 1:nrow(tempresult)){
          name <- tempresult$Term[r]
          genevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Genes[r]
          foldchangevector[paste(name,currentlevel,"Set2",sep = "__")] <- 
            median(foldchangemerge_final[foldchangemerge_final$genename %in% firstup(strsplit(tempresult$Genes[r],";")[[1]]),paste0("Set2_",currentlevel)], na.rm = T)
          pvaluevector[paste(name,currentlevel,"Set2",sep = "__")] <- tempresult$Adjusted.P.value[r]
        }
      }
      
      
      pathwayset2 <- tempresult$Term
      sigpathways_set2[i] <- paste(pathwayset2, collapse = ";")
      siggenes_set2num[i] <- length(DEGlist[[matchingset[1]]])
      sigpathways_set2num[i] <- length(pathwayset2)
      
      sigpathways_sharedpath[i] <- NA 
      
      siggenes_shared[i] <- NA
      sigpathways_sharedgenes[i] <- NA
      siggenes_sharednum[i] <- NA
      sigpathways_sharedpathnum[i] <- NA
      sigpathways_sharedgenesnum[i] <- NA
      
    }
  }else{
    siggenes_set1[i] <- NA
    sigpathways_set1[i] <- NA
    siggenes_set1num[i] <- NA
    sigpathways_set1num[i] <- NA
    
    siggenes_set2[i] <- NA
    sigpathways_set2[i] <-NA
    siggenes_set2num[i] <- NA
    sigpathways_set2num[i] <- NA
    
    sigpathways_sharedpath[i] <- NA 
    siggenes_shared[i] <- NA
    sigpathways_sharedgenes[i] <- NA
    siggenes_sharednum[i] <- NA
    sigpathways_sharedpathnum[i] <- NA
    sigpathways_sharedgenesnum[i] <- NA
  }
  
  
}
topDEG2 = sapply(topDEG,paste,collapse = ", ")
topDEG2 <- topDEG2[paste0("Set1_",alllevels)]
topDEG2 <- topDEG2[paste0("Set2_",alllevels)]

sigpathways_set1 <- gsub("Homo sapiens_","",sigpathways_set1)
sigpathways_set2 <- gsub("Homo sapiens_","",sigpathways_set2)
sigpathways_sharedpath <- gsub("Homo sapiens_","",sigpathways_sharedpath)
sigpathways_sharedgenes <- gsub("Homo sapiens_","",sigpathways_sharedgenes)

finalframe <- data.frame(names = alllevels,alllevels,siggenes_set1, siggenes_set1num, sigpathways_set1, sigpathways_set1num,
                         siggenes_set2, siggenes_set2num, sigpathways_set2, sigpathways_set2num,
                         siggenes_shared,siggenes_sharednum,sigpathways_sharedpath,sigpathways_sharedpathnum,sigpathways_sharedgenes,sigpathways_sharedgenesnum)

identical( names(foldchangevector),names(pvaluevector))
result <- data.frame(name = names(foldchangevector),fold = foldchangevector,FDR = -log10(pvaluevector), genevector)

splitname <- strsplit(result$name,"__")
result$pathwayname <- sapply(splitname, `[[`, 1)
result$celltype <- sapply(splitname, `[[`, 2)
result$Set <- sapply(splitname, `[[`, 3)
result$set_celltype <- paste0(result$celltype,"__",result$Set)

result_final <- result %>% group_by(set_celltype) %>% top_n(n = 5, wt = FDR)
save(result,result_final, file= "pathwayinfo_Fruc_subset.rda")
write.csv(finalframe,file = "cross_batch_DEG_FDR_Doug_Fruc_subset.csv")

