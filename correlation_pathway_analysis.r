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
  if(length(siggenes) < 20){return(NA)}
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

library(Seurat)
library(enrichR)
#Seurat correlation in hypothalamus neuron
#load("/Users/Tsai_Lab/Desktop/shared_slides/filtered_object/Hyp_neurons.rda")
load("/Users/Tsai_Lab/Desktop/shared_slides/filtered_object/Hyp.rda")

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered,value = "data.diseasestatus")
allvalues <- unique(dropEST.combined.filtered$data.diseasestatus)
rm("final")
for(trt in allvalues){
  currenttrt <- subset(dropEST.combined.filtered,idents = trt)
  currenttrt <- SetIdent(currenttrt,value = "celltype")
  allcelltypes <- as.character(unique(currenttrt$celltype))
  rm("mergedframe")
  for(celltype in allcelltypes){
    currentcell <- subset(currenttrt, idents = celltype)
    seurat.data = currentcell[["RNA"]]@data
    validind <- rowMeans(seurat.data == 0)
    seurat.data <- seurat.data[-which(validind >= 0.9),,drop = F]
    if(ncol(seurat.data) < 3){next}
    allgenes <- rownames(seurat.data)
    pval <- NULL
    corval <- NULL
    generes <- NULL
    rm("frame")
    for(gene in allgenes){
      if(gene %in% "mt-Rnr2"){next}
      ind <- which( allgenes %in% gene)
      current <- t(as.matrix(seurat.data[c("mt-Rnr2",gene),]))
      #if(mean(current[,2] == 0) >= 0.9){next}
      res <- cor.test(current[,1],current[,2] , method = "spearman")
      pval[ind] <- res$p.value
      corval[ind] <-  res$estimate
      generes[ind] <- gene
    }
    frame <- data.frame(generes,pval, corval, celltype)
    frame$pval <- p.adjust(frame$pval, method = "BH")
    if(!exists("mergedframe")){
      mergedframe <- frame
    }else{
      mergedframe <- rbind.data.frame(mergedframe, frame)
    }
    #seurat.data.cor.spearman = cor(t(as.matrix(seurat.data)), 
    #                               method = "spearman")
  }
  mergedframe$trt <- trt
  #mergedframe$pval <- p.adjust(mergedframe$pval, method = "BH")
  if(!exists("final")){
    final <- mergedframe
  }else{
    final <- rbind.data.frame(final, mergedframe)
  }
}
#Doing correlation across all genes with Humanin in each cell type
final_sig <- final[final$pval < 0.05,,]
final_sig <- final_sig[!final_sig$generes %in% "Malat1",]
final_sig <- final_sig[!is.na(final_sig$generes),]
#Significant ones we do pathway enrichment analysis
alltrt <- unique(final$trt)
rm("finalpath")
for(trt in alltrt){
  currentframe <- final[final$trt %in% trt,]
  currentframe <- currentframe[currentframe$pval < 0.05,]
  currentframe <- currentframe[!is.na(currentframe$generes),]
  allcelltypes <- unique(currentframe$celltype)
  rm("mergedpath")
  for(cell in allcelltypes){
    currentgene <- currentframe[currentframe$celltype %in% cell,]
    #if(nrow(currentgene) <= 5){next}
    currentpath <- enrichrfunction(currentgene$generes)
    if(is.null(currentpath)){next}
    if(is.na(currentpath)){next}
    fold <- NULL
    for(i in 1:nrow(currentpath)){
      availgene <- unlist(strsplit(currentpath$Genes[i],";"))
      fold[i] <- median(currentgene$corval[currentgene$generes %in% firstup(availgene)],na.rm = T)
    }
    currentpath$fold <- fold
    currentpath$celltype <- cell
    if(!exists("mergedpath")){
      mergedpath <- currentpath
    }else{
      mergedpath <- rbind.data.frame(mergedpath, currentpath)
    }
  }
  if(!exists("mergedpath")){next}
  mergedpath$trt <- trt
  if(!exists("finalpath")){
    finalpath <- mergedpath
  }else{
    finalpath <- rbind.data.frame(finalpath, mergedpath)
  }
}
finalpath <- finalpath[!is.na(finalpath$fold),]
#dotplot, color is median correlation and size is FDR
selectedpath <- c("chemical synaptic transmission (GO:0007268)","Dopaminergic synapse Homo sapiens hsa04728","mitochondrial electron transport, NADH to ubiquinone (GO:0006120)","Oxytocin signaling pathway Homo sapiens hsa04921",
                  "rRNA processing (GO:0006364)","negative regulation of cell proliferation (GO:0008285)")
result_final <- finalpath[finalpath$Term %in% selectedpath,]
result_final$Term <- gsub(" (","",result_final$Term,fixed = T)
result_final$Term <- gsub(")","",result_final$Term,fixed = T)
result_final$Term <- gsub("GO:[0-9]+| hsa[0-9]+","",result_final$Term)
result_final$Term <- gsub(" Homo sapiens","",result_final$Term)
result_final$Term <- firstup(result_final$Term)
result_final$FDR <- -log10(result_final$Adjusted.P.value)
result_final$FDR[result_final$FDR > 3 ] <- 3
write.csv(final_sig, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/significant_correlation_genes_hyp.csv")
write.csv(result_final, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/significant_correlation_pathways_hyp.csv")
write.csv(finalpath, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/significant_correlation_all_pathways_hyp.csv")

final_sig$celltype_tissue <- paste0(final_sig$celltype,"_",final_sig$trt,"_Hyp")
final_sig_Hyp <- final_sig
upregtable_Hyp <- final_sig[final_sig$corval > 0,]
downregtable_Hyp <- final_sig[final_sig$corval < 0,]

library(RColorBrewer)
library(ggplot2)
p_Hyp <- ggplot(result_final, aes(y = Term,
                              x = celltype)) +        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = fold, 
                 size = FDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-0.35,0.35),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+  facet_grid(.~trt)+
  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                   strip.text.x = element_text(size = 20))+ guides(colour=guide_legend(title="median(Corr)"),
                                                                   size =guide_legend(title="-log10(FDR)"))

pdf("/Users/Tsai_Lab/Downloads/Stuffs_folder/Hyp.pdf", width= 16, height = 10)
print(p_Hyp)
dev.off()


#Seurat correlation in SI
#load("/Users/Tsai_Lab/Desktop/shared_slides/filtered_object/Hyp_neurons.rda")
load("/Users/Tsai_Lab/Desktop/shared_slides/filtered_object/SI.rda")

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered,value = "data.diseasestatus")
allvalues <- unique(dropEST.combined.filtered$data.diseasestatus)
rm("final")
for(trt in allvalues){
  currenttrt <- subset(dropEST.combined.filtered,idents = trt)
  currenttrt <- SetIdent(currenttrt,value = "celltype")
  allcelltypes <- as.character(unique(currenttrt$celltype))
  rm("mergedframe")
  for(celltype in allcelltypes){
    currentcell <- subset(currenttrt, idents = celltype)
    seurat.data = currentcell[["RNA"]]@data
    validind <- rowMeans(seurat.data == 0)
    seurat.data <- seurat.data[-which(validind >= 0.9),,drop = F]
    if(ncol(seurat.data) < 3){next}
    allgenes <- rownames(seurat.data)
    pval <- NULL
    corval <- NULL
    generes <- NULL
    rm("frame")
    for(gene in allgenes){
      if(gene %in% "mt-Rnr2"){next}
      ind <- which( allgenes %in% gene)
      current <- t(as.matrix(seurat.data[c("mt-Rnr2",gene),]))
      #if(mean(current[,2] == 0) >= 0.9){next}
      res <- cor.test(current[,1],current[,2] , method = "spearman")
      pval[ind] <- res$p.value
      corval[ind] <-  res$estimate
      generes[ind] <- gene
    }
    frame <- data.frame(generes,pval, corval, celltype)
    frame$pval <- p.adjust(frame$pval, method = "BH")
    if(!exists("mergedframe")){
      mergedframe <- frame
    }else{
      mergedframe <- rbind.data.frame(mergedframe, frame)
    }
    #seurat.data.cor.spearman = cor(t(as.matrix(seurat.data)), 
    #                               method = "spearman")
  }
  mergedframe$trt <- trt
  #mergedframe$pval <- p.adjust(mergedframe$pval, method = "BH")
  if(!exists("final")){
    final <- mergedframe
  }else{
    final <- rbind.data.frame(final, mergedframe)
  }
}
#Doing correlation across all genes with Humanin in each cell type
final_sig <- final[final$pval < 0.05,,]
final_sig <- final_sig[!final_sig$generes %in% "Malat1",]
final_sig <- final_sig[!is.na(final_sig$generes),]
#Significant ones we do pathway enrichment analysis
alltrt <- unique(final$trt)
rm("finalpath")
for(trt in alltrt){
  currentframe <- final[final$trt %in% trt,]
  currentframe <- currentframe[currentframe$pval < 0.05,]
  currentframe <- currentframe[!is.na(currentframe$generes),]
  allcelltypes <- unique(currentframe$celltype)
  rm("mergedpath")
  for(cell in allcelltypes){
    currentgene <- currentframe[currentframe$celltype %in% cell,]
    #if(nrow(currentgene) <= 5){next}
    currentpath <- enrichrfunction(currentgene$generes)
    if(is.null(currentpath)){next}
    if(is.na(currentpath)){next}
    fold <- NULL
    for(i in 1:nrow(currentpath)){
      availgene <- unlist(strsplit(currentpath$Genes[i],";"))
      fold[i] <- median(currentgene$corval[currentgene$generes %in% firstup(availgene)],na.rm = T)
    }
    currentpath$fold <- fold
    currentpath$celltype <- cell
    if(!exists("mergedpath")){
      mergedpath <- currentpath
    }else{
      mergedpath <- rbind.data.frame(mergedpath, currentpath)
    }
  }
  if(!exists("mergedpath")){next}
  mergedpath$trt <- trt
  if(!exists("finalpath")){
    finalpath <- mergedpath
  }else{
    finalpath <- rbind.data.frame(finalpath, mergedpath)
  }
}
finalpath <- finalpath[!is.na(finalpath$fold),]
#dotplot, color is median correlation and size is FDR
selectedpath <- c("Antigen processing and presentation Homo sapiens hsa04612","aerobic respiration (GO:0009060)","Carbon metabolism Homo sapiens hsa01200",
                  "cholesterol efflux (GO:0033344)","very-low-density lipoprotein particle assembly (GO:0034379)","Fructose and mannose metabolism Homo sapiens hsa00051","Peroxisome Homo sapiens hsa04146",
                  "Pyruvate metabolism Homo sapiens hsa00620")
result_final <- finalpath[finalpath$Term %in% selectedpath,]
result_final$Term <- gsub(" (","",result_final$Term,fixed = T)
result_final$Term <- gsub(")","",result_final$Term,fixed = T)
result_final$Term <- gsub("GO:[0-9]+| hsa[0-9]+","",result_final$Term)
result_final$Term <- gsub(" Homo sapiens","",result_final$Term)
result_final$Term <- firstup(result_final$Term)
result_final$FDR <- -log10(result_final$Adjusted.P.value)
result_final$FDR[result_final$FDR > 3 ] <- 3
write.csv(final_sig, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/significant_correlation_genes_SI.csv")
write.csv(result_final, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/significant_correlation_pathways_SI.csv")
write.csv(finalpath, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/significant_correlation_all_pathways_SI.csv")

final_sig$celltype_tissue <- paste0(final_sig$celltype,"_",final_sig$trt,"_SI")
final_sig_SI <- final_sig
upregtable_SI <- final_sig[final_sig$corval > 0,]
downregtable_SI <- final_sig[final_sig$corval < 0,]

library(RColorBrewer)
library(ggplot2)
p_SI <- ggplot(result_final, aes(y = Term,
                              x = celltype)) +        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = fold, 
                 size = FDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-0.40,0.40),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+  facet_grid(.~trt)+
  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                   strip.text.x = element_text(size = 20))+ guides(colour=guide_legend(title="median(Corr)"),
                                                                   size =guide_legend(title="-log10(FDR)"))

pdf("/Users/Tsai_Lab/Downloads/Stuffs_folder/SI.pdf", width = 10, height = 10)
print(p_SI)
dev.off()


#
#Seurat correlation in Liver
#load("/Users/Tsai_Lab/Desktop/shared_slides/filtered_object/Hyp_neurons.rda")
load("/Users/Tsai_Lab/Desktop/shared_slides/filtered_object/Liver.rda")

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered,value = "data.diseasestatus")
allvalues <- unique(dropEST.combined.filtered$data.diseasestatus)
rm("final")
for(trt in allvalues){
  currenttrt <- subset(dropEST.combined.filtered,idents = trt)
  currenttrt <- SetIdent(currenttrt,value = "celltype")
  allcelltypes <- as.character(unique(currenttrt$celltype))
  rm("mergedframe")
  for(celltype in allcelltypes){
    currentcell <- subset(currenttrt, idents = celltype)
    seurat.data = currentcell[["RNA"]]@data
    validind <- rowMeans(seurat.data == 0)
    seurat.data <- seurat.data[-which(validind >= 0.9),,drop = F]
    if(ncol(seurat.data) < 3){next}
    allgenes <- rownames(seurat.data)
    pval <- NULL
    corval <- NULL
    generes <- NULL
    rm("frame")
    for(gene in allgenes){
      if(gene %in% "mt-Rnr2"){next}
      ind <- which( allgenes %in% gene)
      current <- t(as.matrix(seurat.data[c("mt-Rnr2",gene),]))
      #if(mean(current[,2] == 0) >= 0.9){next}
      res <- cor.test(current[,1],current[,2] , method = "spearman")
      pval[ind] <- res$p.value
      corval[ind] <-  res$estimate
      generes[ind] <- gene
    }
    frame <- data.frame(generes,pval, corval, celltype)
    frame$pval <- p.adjust(frame$pval, method = "BH")
    if(!exists("mergedframe")){
      mergedframe <- frame
    }else{
      mergedframe <- rbind.data.frame(mergedframe, frame)
    }
    #seurat.data.cor.spearman = cor(t(as.matrix(seurat.data)), 
    #                               method = "spearman")
  }
  mergedframe$trt <- trt
  #mergedframe$pval <- p.adjust(mergedframe$pval, method = "BH")
  if(!exists("final")){
    final <- mergedframe
  }else{
    final <- rbind.data.frame(final, mergedframe)
  }
}
#Doing correlation across all genes with Humanin in each cell type
final_sig <- final[final$pval < 0.05,,]
final_sig <- final_sig[!final_sig$generes %in% "Malat1",]
final_sig <- final_sig[!is.na(final_sig$generes),]
#significant ones we do pathway enrichment analyLivers
alltrt <- unique(final$trt)
rm("finalpath")
for(trt in alltrt){
  currentframe <- final[final$trt %in% trt,]
  currentframe <- currentframe[currentframe$pval < 0.05,]
  currentframe <- currentframe[!is.na(currentframe$generes),]
  allcelltypes <- unique(currentframe$celltype)
  rm("mergedpath")
  for(cell in allcelltypes){
    currentgene <- currentframe[currentframe$celltype %in% cell,]
    #if(nrow(currentgene) <= 5){next}
    currentpath <- enrichrfunction(currentgene$generes)
    if(is.null(currentpath)){next}
    if(is.na(currentpath)){next}
    fold <- NULL
    for(i in 1:nrow(currentpath)){
      availgene <- unlist(strsplit(currentpath$Genes[i],";"))
      fold[i] <- median(currentgene$corval[currentgene$generes %in% firstup(availgene)],na.rm = T)
    }
    currentpath$fold <- fold
    currentpath$celltype <- cell
    if(!exists("mergedpath")){
      mergedpath <- currentpath
    }else{
      mergedpath <- rbind.data.frame(mergedpath, currentpath)
    }
  }
  if(!exists("mergedpath")){next}
  mergedpath$trt <- trt
  if(!exists("finalpath")){
    finalpath <- mergedpath
  }else{
    finalpath <- rbind.data.frame(finalpath, mergedpath)
  }
}
finalpath <- finalpath[!is.na(finalpath$fold),]
#dotplot, color is median correlation and Liverze is FDR
selectedpath <- c("Apoptosis Homo sapiens hsa04210","Biosynthesis of amino acids Homo sapiens hsa01230",
                  "Carbon metabolism Homo sapiens hsa01200","cholesterol efflux (GO:0033344)","cellular response to oxidative stress (GO:0034599)",
                  "extracellular matrix organization (GO:0030198)","neutrophil degranulation (GO:0043312)","Ribosome Homo sapiens hsa03010",
                  "very-low-density lipoprotein particle assembly (GO:0034379)")
result_final <- finalpath[finalpath$Term %in% selectedpath,]
result_final$Term <- gsub(" (","",result_final$Term,fixed = T)
result_final$Term <- gsub(")","",result_final$Term,fixed = T)
result_final$Term <- gsub("GO:[0-9]+| hsa[0-9]+","",result_final$Term)
result_final$Term <- gsub(" Homo sapiens","",result_final$Term)
result_final$Term <- firstup(result_final$Term)
result_final$FDR <- -log10(result_final$Adjusted.P.value)
result_final$FDR[result_final$FDR > 3 ] <- 3
write.csv(final_sig, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/Livergnificant_correlation_genes_Liver.csv")
write.csv(result_final, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/Livergnificant_correlation_pathways_Liver.csv")
write.csv(finalpath, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/significant_correlation_all_pathways_Liver.csv")

final_sig$celltype_tissue <- paste0(final_sig$celltype,"_",final_sig$trt,"_Liver")
final_sig_Liver <- final_sig
upregtable_Liver <- final_sig[final_sig$corval > 0,]
downregtable_Liver <- final_sig[final_sig$corval < 0,]


library(RColorBrewer)
library(ggplot2)
p_Liver <- ggplot(result_final, aes(y = Term,
                              x = celltype)) +        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = fold, 
                 size = FDR))  +    ## geom_point for circle illuLiveron
  scale_color_gradientn(limits = c(-0.5,0.5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+  facet_grid(.~trt)+
  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                   strip.text.x = element_text(size = 20))+ guides(colour=guide_legend(title="median(Corr)"),
                                                                   size =guide_legend(title="-log10(FDR)"))

pdf("/Users/Tsai_Lab/Downloads/Stuffs_folder/Liver.pdf", width = 16, height = 10)
print(p_Liver)
dev.off()


#SVF


#Seurat correlation in SVF
#load("/Users/Tsai_Lab/Desktop/shared_slides/filtered_object/Hyp_neurons.rda")
load("/Users/Tsai_Lab/Desktop/shared_slides/filtered_object/SVF.rda")

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered,value = "data.diseasestatus")
allvalues <- unique(dropEST.combined.filtered$data.diseasestatus)
rm("final")
for(trt in allvalues){
  currenttrt <- subset(dropEST.combined.filtered,idents = trt)
  currenttrt <- SetIdent(currenttrt,value = "celltype")
  allcelltypes <- as.character(unique(currenttrt$celltype))
  rm("mergedframe")
  for(celltype in allcelltypes){
    currentcell <- subset(currenttrt, idents = celltype)
    seurat.data = currentcell[["RNA"]]@data
    validind <- rowMeans(seurat.data == 0)
    seurat.data <- seurat.data[-which(validind >= 0.9),,drop = F]
    if(ncol(seurat.data) < 3){next}
    allgenes <- rownames(seurat.data)
    pval <- NULL
    corval <- NULL
    generes <- NULL
    rm("frame")
    for(gene in allgenes){
      if(gene %in% "mt-Rnr2"){next}
      ind <- which( allgenes %in% gene)
      current <- t(as.matrix(seurat.data[c("mt-Rnr2",gene),]))
      #if(mean(current[,2] == 0) >= 0.9){next}
      res <- cor.test(current[,1],current[,2] , method = "spearman")
      pval[ind] <- res$p.value
      corval[ind] <-  res$estimate
      generes[ind] <- gene
    }
    frame <- data.frame(generes,pval, corval, celltype)
    frame$pval <- p.adjust(frame$pval, method = "BH")
    if(!exists("mergedframe")){
      mergedframe <- frame
    }else{
      mergedframe <- rbind.data.frame(mergedframe, frame)
    }
    #seurat.data.cor.spearman = cor(t(as.matrix(seurat.data)), 
    #                               method = "spearman")
  }
  mergedframe$trt <- trt
  #mergedframe$pval <- p.adjust(mergedframe$pval, method = "BH")
  if(!exists("final")){
    final <- mergedframe
  }else{
    final <- rbind.data.frame(final, mergedframe)
  }
}
#Doing correlation across all genes with Humanin in each cell type
final_sig <- final[final$pval < 0.05,,]
final_sig <- final_sig[!final_sig$generes %in% "Malat1",]
final_sig <- final_sig[!is.na(final_sig$generes),]
#Significant ones we do pathway enrichment analysis
alltrt <- unique(final$trt)
rm("finalpath")
for(trt in alltrt){
  currentframe <- final[final$trt %in% trt,]
  currentframe <- currentframe[currentframe$pval < 0.05,]
  currentframe <- currentframe[!is.na(currentframe$generes),]
  allcelltypes <- unique(currentframe$celltype)
  rm("mergedpath")
  for(cell in allcelltypes){
    currentgene <- currentframe[currentframe$celltype %in% cell,]
    #if(nrow(currentgene) <= 5){next}
    currentpath <- enrichrfunction(currentgene$generes)
    if(is.null(currentpath)){next}
    if(is.na(currentpath)){next}
    fold <- NULL
    for(i in 1:nrow(currentpath)){
      availgene <- unlist(strsplit(currentpath$Genes[i],";"))
      fold[i] <- median(currentgene$corval[currentgene$generes %in% firstup(availgene)],na.rm = T)
    }
    currentpath$fold <- fold
    currentpath$celltype <- cell
    if(!exists("mergedpath")){
      mergedpath <- currentpath
    }else{
      mergedpath <- rbind.data.frame(mergedpath, currentpath)
    }
  }
  if(!exists("mergedpath")){next}
  mergedpath$trt <- trt
  if(!exists("finalpath")){
    finalpath <- mergedpath
  }else{
    finalpath <- rbind.data.frame(finalpath, mergedpath)
  }
}
finalpath <- finalpath[!is.na(finalpath$fold),]
#dotplot, color is median correlation and size is FDR
selectedpath <- c("Antigen processing and presentation Homo sapiens hsa04612","cellular protein metabolic process (GO:0044267)",
                  "collagen catabolic process (GO:0030574)","ECM-receptor interaction Homo sapiens hsa04512",
                  "extracellular matrix disassembly (GO:0022617)","Lysosome Homo sapiens hsa04142","low-density lipoprotein particle clearance (GO:0034383)",
                  "neutrophil degranulation (GO:0043312)","Ribosome Homo sapiens hsa03010")
result_final <- finalpath[finalpath$Term %in% selectedpath,]
result_final$Term <- gsub(" (","",result_final$Term,fixed = T)
result_final$Term <- gsub(")","",result_final$Term,fixed = T)
result_final$Term <- gsub("GO:[0-9]+| hsa[0-9]+","",result_final$Term)
result_final$Term <- gsub(" Homo sapiens","",result_final$Term)
result_final$Term <- firstup(result_final$Term)
result_final$FDR <- -log10(result_final$Adjusted.P.value)
result_final$FDR[result_final$FDR > 3 ] <- 3
write.csv(final_sig, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/significant_correlation_genes_SVF.csv")
write.csv(result_final, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/significant_correlation_pathways_SVF.csv")
write.csv(finalpath, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/significant_correlation_all_pathways_SVF.csv")

final_sig$celltype_tissue <- paste0(final_sig$celltype,"_",final_sig$trt,"_SVF")
final_sig_SVF <- final_sig
upregtable_SVF <- final_sig[final_sig$corval > 0,]
downregtable_SVF <- final_sig[final_sig$corval < 0,]


library(RColorBrewer)
library(ggplot2)
p_SVF <- ggplot(result_final, aes(y = Term,
                              x = celltype)) +        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = fold, 
                 size = FDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-0.35,0.35),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+  facet_grid(.~trt)+
  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                   strip.text.x = element_text(size = 20))+ guides(colour=guide_legend(title="median(Corr)"),
                                                                   size =guide_legend(title="-log10(FDR)"))

pdf("/Users/Tsai_Lab/Downloads/Stuffs_folder/SVF.pdf",width = 10, height = 10)
print(p_SVF)
dev.off()

library(gridExtra)
pdf("/Users/Tsai_Lab/Downloads/Stuffs_folder/corr_combined.pdf",width = 28, height = 16)
grid.arrange(p_SI,p_Liver,p_SVF,p_Hyp, nrow = 2)
dev.off()

#Global analysis

upregtable_all <- rbind.data.frame(upregtable_Hyp,upregtable_SI,upregtable_Liver,upregtable_SVF)
result_up <- as.data.frame(table(upregtable_all$generes))
celltype_SI <- NULL
celltype_Liver <- NULL
celltype_SVF <- NULL
celltype_Hyp <- NULL
for(i in 1:nrow(result_up)){
  currentcelltype <- upregtable_all$celltype_tissue[upregtable_all$generes %in% result_up$Var1[i]]
  celltype_SI[i] <- paste0(gsub("_SI","",currentcelltype[grep("_SI",currentcelltype)]), collapse = ";")
  celltype_SVF[i] <- paste0(gsub("_SVF","",currentcelltype[grep("_SVF",currentcelltype)]), collapse = ";")
  celltype_Liver[i] <- paste0(gsub("_Liver","",currentcelltype[grep("_Liver",currentcelltype)]), collapse = ";")
  celltype_Hyp[i] <- paste0(gsub("_Hyp","",currentcelltype[grep("_Hyp",currentcelltype)]), collapse = ";")
}
result_up$celltype_SI <- celltype_SI
result_up$celltype_SVF <- celltype_SVF
result_up$celltype_Liver <- celltype_Liver
result_up$celltype_Hyp <- celltype_Hyp

downregtable_all <- rbind.data.frame(downregtable_Hyp,downregtable_SI,downregtable_Liver,downregtable_SVF)
result_down <- as.data.frame(table(downregtable_all$generes))
celltype_SI <- NULL
celltype_Liver <- NULL
celltype_SVF <- NULL
celltype_Hyp <- NULL
for(i in 1:nrow(result_down)){
  currentcelltype <- downregtable_all$celltype_tissue[downregtable_all$generes %in% result_down$Var1[i]]
  celltype_SI[i] <- paste0(gsub("_SI","",currentcelltype[grep("_SI",currentcelltype)]), collapse = ";")
  celltype_SVF[i] <- paste0(gsub("_SVF","",currentcelltype[grep("_SVF",currentcelltype)]), collapse = ";")
  celltype_Liver[i] <- paste0(gsub("_Liver","",currentcelltype[grep("_Liver",currentcelltype)]), collapse = ";")
  celltype_Hyp[i] <- paste0(gsub("_Hyp","",currentcelltype[grep("_Hyp",currentcelltype)]), collapse = ";")
}
result_down$celltype_SI <- celltype_SI
result_down$celltype_SVF <- celltype_SVF
result_down$celltype_Liver <- celltype_Liver
result_down$celltype_Hyp <- celltype_Hyp

#Now visualize some genes
result_up <- result_up[order(result_up$Freq, decreasing = T),]
result_down <- result_down[order(result_down$Freq,decreasing = T),]
write.csv(result_up, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/up_corr_gene_summary.csv")
write.csv(result_down, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/down_corr_gene_summary.csv")

selectedgenes <- unique(c(as.character(result_up$Var1[1:10]),as.character(result_down$Var1[1:10])))
plotsig <- rbind.data.frame(final_sig_Liver[final_sig_Liver$generes %in% selectedgenes,],
                            final_sig_SVF[final_sig_SVF$generes %in% selectedgenes,],
                            final_sig_SI[final_sig_SI$generes %in% selectedgenes,],
                            final_sig_Hyp[final_sig_Hyp$generes %in% selectedgenes,])
plotsig$celltype_tissue2 <- gsub("_Control|_Fruc|_HFHS","_",plotsig$celltype_tissue)
plotsig$celltype_tissue2 <- factor(plotsig$celltype_tissue2, levels = unique(plotsig$celltype_tissue2))
plotsig$corval[plotsig$corval > 0.5] <- 0.5
plotsig$corval[plotsig$corval < -0.5] <- -0.5
plotsig$logFDR <- log10(plotsig$pval)
plotsig$logFDR[plotsig$logFDR < -10] <- -10
p_all_gene <- ggplot(plotsig, aes(y = generes,
                         x = celltype_tissue2)) +        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = corval, 
                 size = -logFDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-0.5,0.5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+  facet_grid(.~trt)+
  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                   strip.text.x = element_text(size = 20))+ guides(colour=guide_legend(title="median(Corr)"),
                                                                   size =guide_legend(title="-log10(FDR)"))

#In each tissue
result_up <- as.data.frame(table(upregtable_Hyp$generes))
celltype_Hyp <- NULL
for(i in 1:nrow(result_up)){
  currentcelltype <- upregtable_all$celltype_tissue[upregtable_all$generes %in% result_up$Var1[i]]
  celltype_Hyp[i] <- paste0(gsub("_Hyp","",currentcelltype[grep("_Hyp",currentcelltype)]), collapse = ";")
}
result_up$celltype <- celltype_Hyp

result_down <- as.data.frame(table(downregtable_Hyp$generes))
celltype_Hyp <- NULL
for(i in 1:nrow(result_down)){
  currentcelltype <- downregtable_all$celltype_tissue[downregtable_all$generes %in% result_down$Var1[i]]
  celltype_Hyp[i] <- paste0(gsub("_Hyp","",currentcelltype[grep("_Hyp",currentcelltype)]), collapse = ";")
}
result_down$celltype <- celltype_Hyp

result_up <- result_up[order(result_up$Freq, decreasing = T),]
result_down <- result_down[order(result_down$Freq,decreasing = T),]
write.csv(result_up, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/up_corr_gene_summary_Hyp.csv")
write.csv(result_down, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/down_corr_gene_summary_Hyp.csv")

selectedgenes <- unique(c(as.character(result_up$Var1[1:10]),as.character(result_down$Var1[1:10])))
plotsig <- final_sig_Hyp[final_sig_Hyp$generes %in% selectedgenes,]
#plotsig$celltype_tissue2 <- gsub("_Control|_Fruc|_HFHS","_",plotsig$celltype_tissue)
plotsig$corval[plotsig$corval > 0.5] <- 0.5
plotsig$corval[plotsig$corval < -0.5] <- -0.5
plotsig$logFDR <- log10(plotsig$pval)
plotsig$logFDR[plotsig$logFDR < -10] <- -10
p_Hyp_gene <- ggplot(plotsig, aes(y = generes,
                                  x = celltype)) +        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = corval, 
                 size = -logFDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-0.5,0.5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+  facet_grid(.~trt)+
  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                   strip.text.x = element_text(size = 20))+ guides(colour=guide_legend(title="median(Corr)"),
                                                                   size =guide_legend(title="-log10(FDR)"))


result_up <- as.data.frame(table(upregtable_Liver$generes))
celltype_Liver <- NULL
for(i in 1:nrow(result_up)){
  currentcelltype <- upregtable_all$celltype_tissue[upregtable_all$generes %in% result_up$Var1[i]]
  celltype_Liver[i] <- paste0(gsub("_Liver","",currentcelltype[grep("_Liver",currentcelltype)]), collapse = ";")
}
result_up$celltype <- celltype_Liver

result_down <- as.data.frame(table(downregtable_Liver$generes))
celltype_Liver <- NULL
for(i in 1:nrow(result_down)){
  currentcelltype <- downregtable_all$celltype_tissue[downregtable_all$generes %in% result_down$Var1[i]]
  celltype_Liver[i] <- paste0(gsub("_Liver","",currentcelltype[grep("_Liver",currentcelltype)]), collapse = ";")
}
result_down$celltype <- celltype_Liver

result_up <- result_up[order(result_up$Freq, decreasing = T),]
result_down <- result_down[order(result_down$Freq,decreasing = T),]

write.csv(result_up, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/up_corr_gene_summary_Liver.csv")
write.csv(result_down, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/down_corr_gene_summary_Liver.csv")

selectedgenes <- unique(c(as.character(result_up$Var1[1:10]),as.character(result_down$Var1[1:10])))
plotsig <- final_sig_Liver[final_sig_Liver$generes %in% selectedgenes,]
#plotsig$celltype_tissue2 <- gsub("_Control|_Fruc|_HFHS","_",plotsig$celltype_tissue)
plotsig$corval[plotsig$corval > 0.5] <- 0.5
plotsig$corval[plotsig$corval < -0.5] <- -0.5
plotsig$logFDR <- log10(plotsig$pval)
plotsig$logFDR[plotsig$logFDR < -10] <- -10
p_Liver_gene <- ggplot(plotsig, aes(y = generes,
                                  x = celltype)) +        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = corval, 
                 size = -logFDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-0.5,0.5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+  facet_grid(.~trt)+
  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                   strip.text.x = element_text(size = 20))+ guides(colour=guide_legend(title="median(Corr)"),
                                                                   size =guide_legend(title="-log10(FDR)"))


result_up <- as.data.frame(table(upregtable_SI$generes))
celltype_SI <- NULL
for(i in 1:nrow(result_up)){
  currentcelltype <- upregtable_all$celltype_tissue[upregtable_all$generes %in% result_up$Var1[i]]
  celltype_SI[i] <- paste0(gsub("_SI","",currentcelltype[grep("_SI",currentcelltype)]), collapse = ";")
}
result_up$celltype <- celltype_SI

result_down <- as.data.frame(table(downregtable_SI$generes))
celltype_SI <- NULL
for(i in 1:nrow(result_down)){
  currentcelltype <- downregtable_all$celltype_tissue[downregtable_all$generes %in% result_down$Var1[i]]
  celltype_SI[i] <- paste0(gsub("_SI","",currentcelltype[grep("_SI",currentcelltype)]), collapse = ";")
}
result_down$celltype <- celltype_SI

result_up <- result_up[order(result_up$Freq, decreasing = T),]
result_down <- result_down[order(result_down$Freq,decreasing = T),]

write.csv(result_up, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/up_corr_gene_summary_SI.csv")
write.csv(result_down, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/down_corr_gene_summary_SI.csv")

selectedgenes <- unique(c(as.character(result_up$Var1[1:10]),as.character(result_down$Var1[1:10])))
plotsig <- final_sig_SI[final_sig_SI$generes %in% selectedgenes,]
#plotsig$celltype_tissue2 <- gsub("_Control|_Fruc|_HFHS","_",plotsig$celltype_tissue)
plotsig$corval[plotsig$corval > 0.5] <- 0.5
plotsig$corval[plotsig$corval < -0.5] <- -0.5
plotsig$logFDR <- log10(plotsig$pval)
plotsig$logFDR[plotsig$logFDR < -10] <- -10
p_SI_gene <- ggplot(plotsig, aes(y = generes,
                                  x = celltype)) +        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = corval, 
                 size = -logFDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-0.5,0.5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+  facet_grid(.~trt)+
  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                   strip.text.x = element_text(size = 20))+ guides(colour=guide_legend(title="median(Corr)"),
                                                                   size =guide_legend(title="-log10(FDR)"))



result_up <- as.data.frame(table(upregtable_SVF$generes))
celltype_SVF <- NULL
for(i in 1:nrow(result_up)){
  currentcelltype <- upregtable_all$celltype_tissue[upregtable_all$generes %in% result_up$Var1[i]]
  celltype_SVF[i] <- paste0(gsub("_SVF","",currentcelltype[grep("_SVF",currentcelltype)]), collapse = ";")
}
result_up$celltype <- celltype_SVF

result_down <- as.data.frame(table(downregtable_SVF$generes))
celltype_SVF <- NULL
for(i in 1:nrow(result_down)){
  currentcelltype <- downregtable_all$celltype_tissue[downregtable_all$generes %in% result_down$Var1[i]]
  celltype_SVF[i] <- paste0(gsub("_SVF","",currentcelltype[grep("_SVF",currentcelltype)]), collapse = ";")
}
result_down$celltype <- celltype_SVF

result_up <- result_up[order(result_up$Freq, decreasing = T),]
result_down <- result_down[order(result_down$Freq,decreasing = T),]

write.csv(result_up, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/up_corr_gene_summary_SVF.csv")
write.csv(result_down, file = "/Users/Tsai_Lab/Downloads/Stuffs_folder/down_corr_gene_summary_SVF.csv")

selectedgenes <- unique(c(as.character(result_up$Var1[1:10]),as.character(result_down$Var1[1:10])))
plotsig <- final_sig_SVF[final_sig_SVF$generes %in% selectedgenes,]
#plotsig$celltype_tissue2 <- gsub("_Control|_Fruc|_HFHS","_",plotsig$celltype_tissue)
plotsig$corval[plotsig$corval > 0.5] <- 0.5
plotsig$corval[plotsig$corval < -0.5] <- -0.5
plotsig$logFDR <- log10(plotsig$pval)
plotsig$logFDR[plotsig$logFDR < -10] <- -10
p_SVF_gene <- ggplot(plotsig, aes(y = generes,
                                  x = celltype)) +        ## global aes
  geom_tile(fill = "white") +  xlab("")+ylab("")+       ## to get the rect filled
  geom_point(aes(colour = corval, 
                 size = -logFDR))  +    ## geom_point for circle illusion
  scale_color_gradientn(limits = c(-0.5,0.5),colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40))+  facet_grid(.~trt)+
  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
                   strip.text.x = element_text(size = 20))+ guides(colour=guide_legend(title="median(Corr)"),
                                                                   size =guide_legend(title="-log10(FDR)"))



library(gridExtra)
pdf("/Users/Tsai_Lab/Downloads/Stuffs_folder/corr_gene_combined.pdf",width = 36, height = 36)
grid.arrange(p_all_gene,p_SI_gene,p_Liver_gene,p_SVF_gene,p_Hyp_gene, layout_matrix = matrix(c(1,2,4,1,3,5),ncol = 2))
dev.off()



