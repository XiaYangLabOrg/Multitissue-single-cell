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

library(readr)
Mouse_symbols <- read_delim("/Users/Tsai_Lab/Desktop/Box Sync/PharmOmics_code_data/human_mouse_hcop_fifteen_column.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
Mouse_symbols  <- Mouse_symbols[grep("Ensembl|HGNC",Mouse_symbols$support),]
Mouse_symbols2 <- as.data.frame(Mouse_symbols[,c(5,12)])

alltissues <- c("SVF_","SI_","Liver_","LR_hyp")
alltreatment <- c("HFHS","Fruc")

rm("allnodesresult")
for(currenttissue in alltissues){
  for(treatment in alltreatment){
    current <- read.csv(paste0("~/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/combined_",treatment,"_",currenttissue,"network_final_ver2.csv"))
    allnodes <- unique(c(current$from,current$to))
    Tissues = NULL
    Label.Location = NULL
    for(i in 1:length(allnodes)){
      if(allnodes[i] %in% Mouse_symbols2$mouse_symbol){Tissues[i] = "P"
      Label.Location[i] <- "C"
      }else if(allnodes[i] %in% current$from){
        Tissues[i] <- currenttissue
        Label.Location[i] <- "R"
      }else{
        Tissues[i] <- unique(current$tissue[which(current$to %in% allnodes[i])])
        Label.Location[i] <- "L"
      }
    }
    newnodetable <- data.frame("shared name" = allnodes,name = allnodes, Tissue = Tissues, Label.Location  = Label.Location ,name2 = allnodes)
    colnames(newnodetable)[1] <- "shared name"
    #then filter for the secretome
    proteinatlas <- read.delim("~/Downloads/proteinatlas.tsv")
    #detected <- proteinatlas[!is.na(proteinatlas$Blood.concentration...Conc..blood.IM..pg.L.) | !is.na(proteinatlas$Blood.concentration...Conc..blood.MS..pg.L.),]
    detected <- proteinatlas[proteinatlas$Secretome.location %in% "Secreted to blood",]
    allproteins <- newnodetable$name[newnodetable$Tissue %in% "P"]
    final <- toupper(allproteins)
    finalligands <- NULL
    to_remove <- NULL
    for(ligand in final){
      grepmatch <- c(paste0("^",ligand,"| ",ligand))
      
      if(ligand %in% detected$Gene | length(grep(grepmatch,detected$Gene.synonym)) > 0){
        finalligands <- c(finalligands, ligand)
      }else{
        to_remove <- c(to_remove, ligand)
      }
    }

    newnodetable$Tissue[newnodetable$Tissue %in% "LR_hyp"] <- "Hyp"
    newnodetable$Tissue[newnodetable$Tissue %in% "HYP"] <- "Hyp"
    newnodetable$Tissue <- gsub("_","",newnodetable$Tissue)
    oldnodetable <- newnodetable
    currenttissue2 <- currenttissue
    currenttissue2[currenttissue2 %in% c("LR_hyp")] <- "HYP"
    currenttissue2 <- gsub("_","",currenttissue2)

    oldcurrent <- current[current$tissue %in% currenttissue2,]
    oldnodes <- unique(c(oldcurrent$from,oldcurrent$to))
    oldnodetable <- oldnodetable[oldnodetable$`shared name` %in% oldnodes,]
    oldcurrent$tissue[oldcurrent$from %in% oldnodetable$`shared name`[oldnodetable$Tissue %in% "P"]] <- "P"
    
    #further segregating input and output level
    ligands <- oldnodetable$name[oldnodetable$Tissue %in% "P"]
    oldcurrent$to[!oldcurrent$to %in% ligands] <- paste0(oldcurrent$to[!oldcurrent$to %in% ligands],"_1")
    nodenames_to_add <- unique(oldcurrent$to[!oldcurrent$to %in% ligands])
    tissuename <- unique(oldnodetable$Tissue)
    tissuename <- tissuename[!tissuename %in% "P"]
    frame_to_add <-  data.frame("shared name" = nodenames_to_add, name = nodenames_to_add, Tissue = tissuename, Label.Location = "L",name2 = gsub("_1","",nodenames_to_add))
    colnames(frame_to_add)[1] <- "shared name"
    oldnodetable <- rbind.data.frame(oldnodetable,frame_to_add)
    write.csv(oldcurrent, file = paste0("~/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",currenttissue,"_",treatment,"_edges_intra.csv"),quote = F)
    write.csv(oldnodetable, file = paste0("~/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",currenttissue,"_",treatment,"_nodes_intra.csv"), quote = F)
    
    
    
    newnodetable <- newnodetable[!newnodetable$`shared name` %in% firstup(tolower(to_remove)),]
    current <- current[-union(c(which(current$from %in% firstup(tolower(to_remove)))) , 
                              c(which(current$to %in% firstup(tolower(to_remove))))),]
    current <- current[!duplicated(paste0(current$from,"_",current$to)),]
    write.csv(current, file = paste0("~/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",currenttissue,"_",treatment,"_edges.csv"),quote = F)
    
    write.csv(newnodetable, file = paste0("~/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",currenttissue,"_",treatment,"_nodes.csv"), quote = F)
    cat(currenttissue," with node",nrow(newnodetable)," with edge ",nrow(current))
    cat("\n")
  }
}

#Prepare final ligand list per intra and inter


rm("allnodesresult_intra")
rm("allnodesresult")
for(currenttissue in alltissues){
  for(treatment in alltreatment){
    current <- read.csv(paste0("~/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/combined_",treatment,"_",currenttissue,"network_final_ver2.csv"))
    allnodes <- unique(c(current$from,current$to))
    Tissues = NULL
    Label.Location = NULL
    for(i in 1:length(allnodes)){
      if(allnodes[i] %in% Mouse_symbols2$mouse_symbol){Tissues[i] = "P"
      Label.Location[i] <- "C"
      }else if(allnodes[i] %in% current$from){
        Tissues[i] <- currenttissue
        Label.Location[i] <- "R"
      }else{
        Tissues[i] <- unique(current$tissue[which(current$to %in% allnodes[i])])
        Label.Location[i] <- "L"
      }
    }
    newnodetable <- data.frame("shared name" = allnodes,name = allnodes, Tissue = Tissues, Label.Location  = Label.Location ,name2 = allnodes)
    colnames(newnodetable)[1] <- "shared name"
    #then filter for the secretome
    proteinatlas <- read.delim("~/Downloads/proteinatlas.tsv")
    #detected <- proteinatlas[!is.na(proteinatlas$Blood.concentration...Conc..blood.IM..pg.L.) | !is.na(proteinatlas$Blood.concentration...Conc..blood.MS..pg.L.),]
    detected <- proteinatlas[proteinatlas$Secretome.location %in% "Secreted to blood",]
    allproteins <- newnodetable$name[newnodetable$Tissue %in% "P"]
    final <- toupper(allproteins)
    finalligands <- NULL
    to_remove <- NULL
    for(ligand in final){
      grepmatch <- c(paste0("^",ligand,"| ",ligand))
      
      if(ligand %in% detected$Gene | length(grep(grepmatch,detected$Gene.synonym)) > 0){
        finalligands <- c(finalligands, ligand)
      }else{
        to_remove <- c(to_remove, ligand)
      }
    }
    
    newnodetable$Tissue[newnodetable$Tissue %in% "LR_hyp"] <- "Hyp"
    newnodetable$Tissue[newnodetable$Tissue %in% "HYP"] <- "Hyp"
    newnodetable$Tissue <- gsub("_","",newnodetable$Tissue)
    oldnodetable <- newnodetable
    currenttissue2 <- currenttissue
    currenttissue2[currenttissue2 %in% c("LR_hyp")] <- "HYP"
    currenttissue2 <- gsub("_","",currenttissue2)
    
    oldcurrent <- current[current$tissue %in% currenttissue2,]
    oldnodes <- unique(c(oldcurrent$from,oldcurrent$to))
    oldnodetable <- oldnodetable[oldnodetable$`shared name` %in% oldnodes,]
    oldcurrent$tissue[oldcurrent$from %in% oldnodetable$`shared name`[oldnodetable$Tissue %in% "P"]] <- "P"
    
    #further segregating input and output level
    ligands <- oldnodetable$name[oldnodetable$Tissue %in% "P"]
    oldcurrent$to[!oldcurrent$to %in% ligands] <- paste0(oldcurrent$to[!oldcurrent$to %in% ligands],"_1")
    nodenames_to_add <- unique(oldcurrent$to[!oldcurrent$to %in% ligands])
    tissuename <- unique(oldnodetable$Tissue)
    tissuename <- tissuename[!tissuename %in% "P"]
    frame_to_add <-  data.frame("shared name" = nodenames_to_add, name = nodenames_to_add, Tissue = tissuename, Label.Location = "L",name2 = gsub("_1","",nodenames_to_add))
    colnames(frame_to_add)[1] <- "shared name"
    oldnodetable <- rbind.data.frame(oldnodetable,frame_to_add)
    oldnodetable <- oldnodetable[oldnodetable$Tissue %in% "P",c(2,3)]
    colnames(oldnodetable ) <- c("Protein name","tissue")
    oldnodetable$treatment <- treatment
    oldnodetable$tissue <- currenttissue
    
    if(!exists("allnodesresult_intra")){
      allnodesresult_intra <- oldnodetable
    }else{
      allnodesresult_intra <- rbind.data.frame(allnodesresult_intra, oldnodetable)
    }
    #write.csv(oldcurrent, file = paste0("~/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",currenttissue,"_",treatment,"_edges_intra.csv"),quote = F)
    #write.csv(oldnodetable, file = paste0("~/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/",currenttissue,"_",treatment,"_nodes_intra.csv"), quote = F)
    
    
    
    newnodetable <- newnodetable[!newnodetable$`shared name` %in% firstup(tolower(to_remove)),]
    
    newnodetable <- newnodetable[newnodetable$Tissue %in% "P",c(2,3)]
    colnames(newnodetable ) <- c("Protein name","tissue")
    newnodetable$treatment <- treatment
    newnodetable$tissue <- currenttissue
    
    if(!exists("allnodesresult")){
      allnodesresult <- newnodetable
    }else{
      allnodesresult <- rbind.data.frame(allnodesresult, newnodetable)
    }
    
    }
}

allnodesresult$tissue <- gsub("_","",allnodesresult$tissue)
allnodesresult$tissue[allnodesresult$tissue %in% "LR_hyp"] <- "Hyp"
allnodesresult$treatment[allnodesresult$treatment %in% "Fruc"] <- "Fructose"
allnodesresult$`Protein name` <- toupper(allnodesresult$`Protein name`)
write.csv(allnodesresult, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Final_fructose_ligand_new.csv", row.names = F)

allnodesresult_intra$tissue <- gsub("_","",allnodesresult_intra$tissue)
allnodesresult_intra$tissue[allnodesresult_intra$tissue %in% "LR_hyp"] <- "Hyp"
allnodesresult_intra$treatment[allnodesresult_intra$treatment %in% "Fruc"] <- "Fructose"
allnodesresult_intra$`Protein name` <- toupper(allnodesresult_intra$`Protein name`)

write.csv(allnodesresult_intra, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Final_fructose_ligand_new_intra.csv",row.names = F)
