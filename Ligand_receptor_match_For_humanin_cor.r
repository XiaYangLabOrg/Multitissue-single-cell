setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/")
alltissues <- c("LR_hyp_","Liver__","SI__","SVF__")
orignote <- c("Orig_hyp","Orig_Liver","Orig_SI","Orig_SVF")
self = c("Su_Hyp_LR_result_Orig_hyp.rda","Su_Liver_result_Orig_Liver.rda","Su_SI_result_Orig_SI.rda",
         "Su_SVF_result_Orig_SVF.rda")
Tissue_receptor_collector <- c("Hyp","Liver","SI","SVF")
Tissue_ligand_collector <- c("HYP","Liver","SI","SVF")

HYP_rec_HFHS <- NULL
HYP_rec_intra_HFHS <- NULL
Liver_rec_HFHS <- NULL
Liver_rec_intra_HFHS <- NULL
SI_rec_HFHS <- NULL
SI_rec_intra_HFHS <- NULL
SVF_rec_HFHS <- NULL
SVF_rec_intra_HFHS <- NULL
HYP_rec_Fruc <- NULL
HYP_rec_intra_Fruc <- NULL
Liver_rec_Fruc <- NULL
Liver_rec_intra_Fruc <- NULL
SI_rec_Fruc <- NULL
SI_rec_intra_Fruc <- NULL
SVF_rec_Fruc <- NULL
SVF_rec_intra_Fruc <- NULL

HYP_lig_HFHS <- NULL
HYP_lig_intra_HFHS <- NULL
Liver_lig_HFHS <- NULL
Liver_lig_intra_HFHS <- NULL
SI_lig_HFHS <- NULL
SI_lig_intra_HFHS <- NULL
SVF_lig_HFHS <- NULL
SVF_lig_intra_HFHS <- NULL
HYP_lig_Fruc <- NULL
HYP_lig_intra_Fruc <- NULL
Liver_lig_Fruc <- NULL
Liver_lig_intra_Fruc <- NULL
SI_lig_Fruc <- NULL
SI_lig_intra_Fruc <- NULL
SVF_lig_Fruc <- NULL
SVF_lig_intra_Fruc <- NULL

for(i in 1:length(alltissues)){
  tissue = alltissues[i]
  currentnode <- read.csv(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",tissue,"Fruc_nodes.csv"))
  currentedge <- read.csv(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",tissue,"Fruc_edges.csv"))
  currentligand <- currentnode$shared.name[currentnode$Tissue %in% "P"]
  assign(paste0(Tissue_ligand_collector[i],"_lig_Fruc"),unique(c(get(paste0(Tissue_ligand_collector[i],"_lig_Fruc")), currentligand)))
  
  allorigs <- list.files(pattern = orignote[i])
  allorigs <- allorigs[!allorigs %in% self[i]]
  #select all related receptors
  
  for(file in allorigs){
    load(file)
    currentreceptor <- unlist(strsplit(file, "_"))[2]
    if(currentreceptor %in% "Hyp"){currentreceptor <- "HYP"}
    currentligand <- currentligand[currentligand %in% currentedge$from[currentedge$tissue %in% currentreceptor]]
    
    for(j in 1:length(result_Fruc_list)){
      currentmat <- result_Fruc_list[[j]][[5]]
      currentreceptors <- currentmat[currentligand[currentligand %in% rownames(currentmat)],,drop = F]
      currentreceptors <- colSums(currentreceptors)
      currentreceptors <- names(currentreceptors[currentreceptors > 0])
      if(length(currentreceptors) == 0){next}
      assign(paste0(currentreceptor,"_rec_Fruc"),unique(c(get(paste0(currentreceptor,"_rec_Fruc")), currentreceptors)))
    }
  }
  
  currentnode <- read.csv(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",tissue,"HFHS_nodes.csv"))
  currentedge <- read.csv(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",tissue,"HFHS_edges.csv"))
  
  currentligand <- currentnode$shared.name[currentnode$Tissue %in% "P"]
  assign(paste0(Tissue_ligand_collector[i],"_lig_HFHS"),unique(c(get(paste0(Tissue_ligand_collector[i],"_lig_HFHS")), currentligand)))
  
  for(file in allorigs){
    load(file)
    currentreceptor <- unlist(strsplit(file, "_"))[2]
    if(currentreceptor %in% "Hyp"){currentreceptor <- "HYP"}
    currentligand <- currentligand[currentligand %in% currentedge$from[currentedge$tissue %in% currentreceptor]]
    
    for(j in 1:length(result_HFHS_list)){
      currentmat <- result_HFHS_list[[j]][[5]]
      currentreceptors <- currentmat[currentligand[currentligand %in% rownames(currentmat)],,drop = F]
      currentreceptors <- colSums(currentreceptors)
      currentreceptors <- names(currentreceptors[currentreceptors > 0])
      if(length(currentreceptors) == 0){next}
      assign(paste0(currentreceptor,"_rec_HFHS"),unique(c(get(paste0(currentreceptor,"_rec_HFHS")), currentreceptors)))
    }
  }
 
  
  
  
  
  #inter_net only
  currentnode <- read.csv(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",tissue,"Fruc_nodes_intra.csv"))
  currentedge <- read.csv(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",tissue,"Fruc_edges_intra.csv"))
  currentligand <- currentnode$shared.name[currentnode$Tissue %in% "P"]
  assign(paste0(Tissue_ligand_collector[i],"_lig_intra_Fruc"),unique(c(get(paste0(Tissue_ligand_collector[i],"_lig_intra_Fruc")), currentligand)))
  
  allorigs = self[i]
  for(file in allorigs){
    load(file)
    currentreceptor <- unlist(strsplit(file, "_"))[2]
    if(currentreceptor %in% "Hyp"){currentreceptor <- "HYP"}
    #currentligand <- currentligand[currentligand %in% currentedge$from[currentedge$tissue %in% currentreceptor]]
    
    for(j in 1:length(result_Fruc_list)){
      currentmat <- result_Fruc_list[[j]][[5]]
      currentreceptors <- currentmat[currentligand[currentligand %in% rownames(currentmat)],,drop = F]
      currentreceptors <- colSums(currentreceptors)
      currentreceptors <- names(currentreceptors[currentreceptors > 0])
      if(length(currentreceptors) == 0){next}
      assign(paste0(currentreceptor,"_rec_intra_Fruc"),unique(c(get(paste0(currentreceptor,"_rec_intra_Fruc")), currentreceptors)))
    }
  }
  
  currentnode <- read.csv(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",tissue,"HFHS_nodes_intra.csv"))
  currentedge <- read.csv(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",tissue,"HFHS_edges_intra.csv"))
  currentligand <- currentnode$shared.name[currentnode$Tissue %in% "P"]
  assign(paste0(Tissue_ligand_collector[i],"_lig_intra_HFHS"),unique(c(get(paste0(Tissue_ligand_collector[i],"_lig_intra_HFHS")), currentligand)))
  
  
  for(file in allorigs){
    load(file)
    currentreceptor <- unlist(strsplit(file, "_"))[2]
    if(currentreceptor %in% "Hyp"){currentreceptor <- "HYP"}
    #currentligand <- currentligand[currentligand %in% currentedge$from[currentedge$tissue %in% currentreceptor]]
    
    for(j in 1:length(result_HFHS_list)){
      currentmat <- result_HFHS_list[[j]][[5]]
      currentreceptors <- currentmat[currentligand[currentligand %in% rownames(currentmat)],,drop = F]
      currentreceptors <- colSums(currentreceptors)
      currentreceptors <- names(currentreceptors[currentreceptors > 0])
      if(length(currentreceptors) == 0){next}
      assign(paste0(currentreceptor,"_rec_intra_HFHS"),unique(c(get(paste0(currentreceptor,"_rec_intra_HFHS")), currentreceptors)))
    }
  }
}


load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/liver.rda")
universe <- rownames(dropEST.combined.filtered@assays$RNA@data)
library(clusterProfiler)
library(GeneOverlap)
setwd("/Users/Tsai_Lab/Downloads/Stuffs_folder/")
alltissue <- c("hyp","SVF","SI","Liver")
treatment <- c("HFHS","Fruc")
rm(final)
receptorsize <- NULL
gsetsize <- NULL
overlapP <- NULL
intersectionsize <- NULL
intersection <- NULL
finaltissue <- NULL
finaltrt <- NULL
finalnetwork <- NULL
type <- NULL
rm(merged_ligand)
rm(merged_ligand_intra)
for(tissue in alltissue){
  significant_correlation_genes <- read.csv(paste0("~/Downloads/Stuffs_folder/significant_correlation_genes_",tissue,".csv"))
  if(tissue %in% "hyp"){tissue = "HYP"}
  for(trt in treatment){
    currentgset <- unique(significant_correlation_genes$generes[significant_correlation_genes$trt %in% trt])
    frame1 <- get(paste0(tissue,"_rec_",trt))
    go.obj <- newGeneOverlap(frame1,
                             currentgset,
                             genome.size=23000)
    go.obj <- testGeneOverlap(go.obj)
    
    receptorsize <- c( receptorsize, length(frame1))
    gsetsize <- c(gsetsize, length(currentgset))
    overlapP <- c(overlapP, go.obj@pval)
    intersectionsize <- c(intersectionsize, length(go.obj@intersection))
    intersection <- c(intersection , paste0(go.obj@intersection, collapse = ";"))
    finaltissue  <- c(finaltissue , tissue)
    finaltrt <- c(finaltrt, trt)
    finalnetwork <- c(finalnetwork, "inter")
    type <- c(type, "Receptor")
    
    
    frame2 <- get(paste0(tissue,"_rec_intra_",trt))
    go.obj <- newGeneOverlap(frame2,
                             currentgset,
                             genome.size=23000)
    go.obj <- testGeneOverlap(go.obj)
    

    
    receptorsize <- c( receptorsize, length(frame2))
    gsetsize <- c(gsetsize, length(currentgset))
    overlapP <- c(overlapP, go.obj@pval)
    intersectionsize <- c(intersectionsize, length(go.obj@intersection))
    intersection <- c(intersection , paste0(go.obj@intersection, collapse = ";"))
    finaltissue  <- c(finaltissue , tissue)
    finaltrt <- c(finaltrt, trt)    
    finalnetwork <- c(finalnetwork, "intra")
    type <- c(type, "Receptor")
    
    
    
    
    ligand <-  get(paste0(tissue,"_lig_",trt))
    ligand_intra <-  get(paste0(tissue,"_lig_intra_",trt))
    
    go.obj <- newGeneOverlap(ligand,
                             currentgset,
                             genome.size=23000)
    go.obj <- testGeneOverlap(go.obj)
    
    receptorsize <- c( receptorsize, length(ligand))
    gsetsize <- c(gsetsize, length(currentgset))
    overlapP <- c(overlapP, go.obj@pval)
    intersectionsize <- c(intersectionsize, length(go.obj@intersection))
    intersection <- c(intersection , paste0(go.obj@intersection, collapse = ";"))
    finaltissue  <- c(finaltissue , tissue)
    finaltrt <- c(finaltrt, trt)
    finalnetwork <- c(finalnetwork, "inter")
    type <- c(type, "ligand")
    
    
    go.obj <- newGeneOverlap(ligand_intra,
                             currentgset,
                             genome.size=23000)
    go.obj <- testGeneOverlap(go.obj)
    
    receptorsize <- c( receptorsize, length(ligand_intra))
    gsetsize <- c(gsetsize, length(currentgset))
    overlapP <- c(overlapP, go.obj@pval)
    intersectionsize <- c(intersectionsize, length(go.obj@intersection))
    intersection <- c(intersection , paste0(go.obj@intersection, collapse = ";"))
    finaltissue  <- c(finaltissue , tissue)
    finaltrt <- c(finaltrt, trt)
    finalnetwork <- c(finalnetwork, "intra")
    type <- c(type, "ligand")
    
    significant_ligands <- significant_correlation_genes[significant_correlation_genes$trt %in% trt,]
    significant_ligands <- significant_ligands[significant_ligands$generes %in% ligand,]
    if(nrow(significant_ligands) > 0){
      significant_ligands$tissue <- tissue
      significant_ligands$trt <- trt
      significant_ligands$network <- "inter"
      if(!exists("merged_ligand")){
        merged_ligand <- significant_ligands
      }else{
        rbind.data.frame(merged_ligand, significant_ligands)
      }
      
    }
   
    
    significant_ligands_intra <- significant_correlation_genes[significant_correlation_genes$trt %in% trt,]
    significant_ligands_intra <- significant_ligands_intra[significant_ligands_intra$generes %in% ligand_intra,]
    if(nrow(significant_ligands_intra) > 0){
      significant_ligands_intra$tissue <- tissue
      significant_ligands_intra$trt <- trt
      significant_ligands_intra$network <- "inter"
      if(!exists("merged_ligand_intra")){
        merged_ligand_intra <- significant_ligands_intra
      }else{
        merged_ligand_intra <- rbind.data.frame(merged_ligand_intra, significant_ligands_intra)
      }
    }
  
  }
}

enrichframe <- data.frame(receptorsize ,
                          gsetsize,
                          overlapP,
                          intersectionsize,
                          intersection ,
                          finaltissue,
                          finaltrt ,finalnetwork,type)
write.csv(enrichframe[enrichframe$type %in% "Receptor",], file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Humanin_network_receptor_entichment.csv")

write.csv(enrichframe[enrichframe$type %in% "ligand",], file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Humanin_network_ligand_entichment.csv")
Receptorframe <- enrichframe[enrichframe$type %in% "Receptor",]
Receptorframe <- Receptorframe[,c(1,3,4,5,2,6:8)]
colnames(Receptorframe)[1:4] <- paste0("Receptor_",colnames(Receptorframe)[1:4])


ligandframe <- enrichframe[enrichframe$type %in% "ligand",]
ligandframe <- ligandframe[,c(1,3,4,5,2,6:8)]
colnames(ligandframe)[1:4] <- paste0("ligand_",colnames(ligandframe)[1:4])

final <- full_join(Receptorframe,ligandframe)
write.csv(final[order(final$Receptor_overlapP),], file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/final_network_humanin_enrich.csv")
##Matching ligand