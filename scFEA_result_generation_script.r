#allfiles <- list.files("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/scFEA_input")
allfiles <- c("Hyp_All.csv" , "Liver_All.csv"  ,  "Liver_Set1_All.csv", 
              "Liver_Set2_All.csv","SI_All.csv"  ,  "SVF_All.csv" )
allfiles <- gsub(".csv","",allfiles,fixed = T)
finalstring <- NULL

#Non-impute version
string <- "python src/scFEA.py --data_dir data --input_dir input --test_file TO_BE_CHANGED.csv --moduleGene_file module_gene_complete_mouse_m168.csv --stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv --cName_file cName_complete_mouse_c70_m168.csv --sc_imputation False --output_flux_file TO_BE_CHANGED_flux.csv --output_balance_file TO_BE_CHANGED_balance.csv"

frontstring <- readLines("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/scFEA_test_template.sh")

for(i in 1:length(allfiles)){
  finalstring[i] <- gsub("TO_BE_CHANGED",allfiles[i],string)
}
#for(i in 1:9){
#  file.remove(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/scFEA_test_",i,".sh"))
#  frontstring2 <- gsub("TO_BE_CHANGED",paste0("scFEA_",i),frontstring)
#  writeLines(c(frontstring2,finalstring[seq(((i-1)*2+1),i*2,1)]), con = paste0("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/scFEA_test_",i,".sh"))
#}
for(i in 1:6){
  file.remove(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/scFEA_test_",i,".sh"))
  frontstring2 <- gsub("TO_BE_CHANGED",paste0("scFEA_",i),frontstring)
  writeLines(c(frontstring2,finalstring[i]), con = paste0("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/scFEA_test_",i,".sh"))
}

#impute version
string <- "python src/scFEA.py --data_dir data --input_dir input --test_file TO_BE_CHANGED.csv --moduleGene_file module_gene_complete_mouse_m168.csv --stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv --cName_file cName_complete_mouse_c70_m168.csv --sc_imputation True --output_flux_file TO_BE_CHANGED_flux_impute.csv --output_balance_file TO_BE_CHANGED_balance_impute.csv"
string <- "python src/scFEA.py --data_dir data --input_dir input --test_file TO_BE_CHANGED.csv --moduleGene_file module_gene_complete_mouse_m168_Fructose.csv --stoichiometry_matrix cmMat_complete_mouse_c70_m168_Fructose.csv --cName_file cName_complete_mouse_c70_m168_Fructose.csv --sc_imputation True --output_flux_file TO_BE_CHANGED_flux_Fructose_impute.csv --output_balance_file TO_BE_CHANGED_balance_Fructose_impute.csv"

for(i in 1:length(allfiles)){
  finalstring[i] <- gsub("TO_BE_CHANGED",allfiles[i],string)
}
#for(i in 1:9){
#  file.remove(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/scFEA_test_",i,"_impute.sh"))
#  frontstring2 <- gsub("TO_BE_CHANGED",paste0("scFEA_impute",i),frontstring)
#  writeLines(c(frontstring2,finalstring[seq(((i-1)*2+1),i*2,1)]), con = paste0("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/scFEA_test_",i,"_impute.sh"))
#}
for(i in 1:6){
  file.remove(paste0("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/scFEA_test_",i,"_impute.sh"))
  frontstring2 <- gsub("TO_BE_CHANGED",paste0("scFEA_impute",i),frontstring)
  writeLines(c(frontstring2,finalstring[i]), con = paste0("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/scFEA_test_",i,"_impute.sh"))
}
#in Hoffman2 R
setwd("/u/home/a/aleph999/scFEA")
allfiles <- list.files(pattern = ".sh")
for(i in 1:length(allfiles)){
  system(paste0("qsub ",allfiles[i]))
}
setwd("/u/home/a/aleph999/scFEA")
allfiles <- list.files(pattern = "impute.sh")
for(i in 1:length(allfiles)){
  system(paste0("qsub ",allfiles[i]))
}
cp -a /u/home/a/aleph999/scFEA_old/input/. /u/home/a/aleph999/scFEA/input
#Obtaining and analyzing results
library(tidyr)
library(plyr)
library(ggplot2)
library(ggridges)
library(Seurat)
library(dplyr)
library(ggpubr)
options(stringsAsFactors = F)
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/scFEA_result")
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/scFEA_result/old_version")
allfiles <- list.files()



#SI
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/SI_combined_Set2_final.rda")
load("/Users/Tsai_Lab/Desktop/shared_slides/filtered_object/SI.rda")
#module_info <- read.csv("./module_info.csv")
module_info <- read.csv("./module_info_Fructose.csv")
colnames(module_info)[1] <- "module_name"
module_info$module_final <- paste0(module_info$Compound_IN_name,"->",module_info$Compound_OUT_name)
currentmeta <- dropEST.combined.filtered@meta.data
currentmeta$name <- rownames(currentmeta)

#result <- read.csv("./SI_All_flux_impute.csv", row.names=1)
result <- read.csv("./SI_All_flux_Fructose_impute.csv", row.names=1)
result2 <- result[,module_info$module_name[module_info$selected %in% "*"]]
result2$name <- rownames(result2)
result2 <- gather(result2, module_name, flux,-name)
result3 <- inner_join(result2, module_info)
result3 <- inner_join(result3, currentmeta)
finalcombined <- result3
allmodules <- unique(finalcombined$module_name)

finalmodules <- NULL
celltype <- NULL
FC <- NULL
pvalue <- NULL
Treatment <- NULL
for(i in 1:length(allmodules)){
  currentframe <- finalcombined[finalcombined$module_name %in% allmodules[i],]
  allcelltypes <- unique(currentframe$celltype)
  for(j in 1:length(allcelltypes)){
    currentframe2 <- currentframe[currentframe$celltype %in% allcelltypes[j],]
  
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , module_info$module_final[module_info$module_name %in%  allmodules[i]])
      FC <- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"Fruc")
    }
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , module_info$module_final[module_info$module_name %in%  allmodules[i]])
      FC <- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"HFHS")
    }
    
  }
  
  #aa <- ggplot(finalcombined[finalcombined$module_name %in% allmodules[i],], aes(y = flux, x = celltype, fill = data.diseasestatus)) + 
  #  geom_boxplot(outlier.size = 0.1) +
  #  stat_compare_means(method = "anova",label="p.signif")+
  #  theme_classic() + ggtitle(module_info$module_final[module_info$module_name %in%  allmodules[i]])+
  #  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))
  #plot(aa)
}
resultvalue_flux <- data.frame(finalmodules,celltype,FC,pvalue,Treatment)
resultvalue_flux_selected <- resultvalue_flux[resultvalue_flux$celltype %in% c("Distal enterocytes","Proximal enterocytes","Goblet","Macrophages"),]
resultvalue_flux_selected$FC[resultvalue_flux_selected$FC > 5] <- 5
resultvalue_flux_selected$FC[resultvalue_flux_selected$FC < -5] <- -5
resultvalue_flux_selected_SI <- resultvalue_flux_selected
ggplot(resultvalue_flux_selected, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_grid(.~celltype)+
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))



selected <- c("Glucose","G6P","Fatty.Acid","Fructose")
#result <- read.csv("./SI_All_balance_impute.csv", row.names=1)
result <- read.csv("./SI_All_balance_Fructose_impute.csv", row.names=1)

result2 <- result[,selected]
result2$name <- rownames(result2)
result2 <- gather(result2, module_name, flux,-name)
result3 <- inner_join(result2, currentmeta)
finalcombined <- result3
allmodules <- unique(finalcombined$module_name)

finalmodules <- NULL
celltype <- NULL
FC <- NULL
pvalue <- NULL
Treatment <- NULL
for(i in 1:length(allmodules)){
  currentframe <- finalcombined[finalcombined$module_name %in% allmodules[i],]
  allcelltypes <- unique(currentframe$celltype)
  for(j in 1:length(allcelltypes)){
    currentframe2 <- currentframe[currentframe$celltype %in% allcelltypes[j],]
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
    celltype <- c(celltype, allcelltypes[j])
    finalmodules  <- c(finalmodules , allmodules[i])
    FC<- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
    pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
    Treatment <- c(Treatment,"Fruc")}
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
    celltype <- c(celltype, allcelltypes[j])
    finalmodules  <- c(finalmodules , allmodules[i])
    FC <- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
    pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
    Treatment <- c(Treatment,"HFHS")
    }
  }
  
  #aa <- ggplot(finalcombined[finalcombined$module_name %in% allmodules[i],], aes(y = flux, x = celltype, fill = data.diseasestatus)) + 
  #  geom_boxplot(outlier.size = 0.1) +
  #  ggtitle(allmodules[i])+ylab("metabolite balance")+theme_ridges()+
  #  stat_compare_means(method = "anova",label="p.signif")+
  #  theme_classic()+
  #  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))
  #plot(aa)
}
resultvalue_balance <- data.frame(finalmodules,celltype,FC,pvalue,Treatment)
resultvalue_balance_selected <- resultvalue_balance[resultvalue_balance$celltype %in% c("Distal enterocytes","Proximal enterocytes","Goblet","Macrophages"),]
resultvalue_balance_selected_SI <- resultvalue_balance_selected 
ggplot(resultvalue_balance_selected, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_grid(.~celltype)+
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))


#Liver
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Combined_allliver.rda")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "Sinusoidal endothelial cells","1" = "Periportal hepatocytes","2" = "Periportal hepatocytes","3" = "Kupffer", "4" = "Pericentral hepatocytes","5" = "Periportal hepatocytes", "6" = "NKT cells",
#                                                    "7" = "Periportal hepatocytes","8" = "Classical dendritic cells","9" = "B cells","10" = "Hepatic stellate cells","11" = "Cholangiocytes","12" = "Dividing cells","13" = "Plasmacytoid dendritic cells"))
#colnames(dropEST.combined.filtered) <- gsub("^_","",colnames(dropEST.combined.filtered))
#dropEST.combined.filtered <- RenameCells(dropEST.combined.filtered, new.names =  gsub("^_","",colnames(dropEST.combined.filtered)))
load("/Users/Tsai_Lab/Desktop/shared_slides/filtered_object/liver.rda")
module_info <- read.csv("./module_info_Fructose.csv")
colnames(module_info)[1] <- "module_name"
module_info$module_final <- paste0(module_info$Compound_IN_name,"->",module_info$Compound_OUT_name)

currentmeta <- dropEST.combined.filtered@meta.data
currentmeta$name <- rownames(currentmeta)
currentmeta$name <- gsub("^_","",currentmeta$name)
currentmeta$celltype <- as.character(currentmeta$celltype)

result <- read.csv("./Liver_All_flux_Fructose_impute.csv", row.names=1)
result2 <- result[,module_info$module_name[module_info$selected %in% "*"]]
result2$name <- rownames(result2)
result2 <- gather(result2, module_name, flux,-name)
result3 <- inner_join(result2, module_info)
result3$name <- gsub("X_","",result3$name)
result3 <- inner_join(result3, currentmeta)
finalcombined <- result3
allmodules <- unique(finalcombined$module_name)
finalmodules <- NULL
celltype <- NULL
FC <- NULL
pvalue <- NULL
Treatment <- NULL
for(i in 1:length(allmodules)){
  currentframe <- finalcombined[finalcombined$module_name %in% allmodules[i],]
  allcelltypes <- unique(currentframe$celltype)
  for(j in 1:length(allcelltypes)){
    currentframe2 <- currentframe[currentframe$celltype %in% allcelltypes[j],]
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , module_info$module_final[module_info$module_name %in%  allmodules[i]])
      FC <- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"Fruc")
    }
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , module_info$module_final[module_info$module_name %in%  allmodules[i]])
      FC <- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"HFHS")
    }
    
  }
  
  #aa <- ggplot(finalcombined[finalcombined$module_name %in% allmodules[i],], aes(y = flux, x = celltype, fill = data.diseasestatus)) + 
  #  geom_boxplot(outlier.size = 0.1) +
  #  stat_compare_means(method = "anova",label="p.signif")+
  #  theme_classic() + ggtitle(module_info$module_final[module_info$module_name %in%  allmodules[i]])+
  #  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))
  #plot(aa)
}
resultvalue_flux <- data.frame(finalmodules,celltype,FC,pvalue,Treatment)
resultvalue_flux_selected <- resultvalue_flux[resultvalue_flux$celltype %in% c("Pericentral hepatocytes","Periportal hepatocytes","Kupffer","Sinusoidal endothelial cells"),]
resultvalue_flux_selected_Liver <- resultvalue_flux_selected
ggplot(resultvalue_flux_selected, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_grid(.~celltype)+
  geom_bar(stat = "identity",position = "dodge") +geom_hline(yintercept = 1)+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))


selected <- c("Glucose","G6P","Fatty.Acid","Fructose")
result <- read.csv("./Liver_All_balance_Fructose_impute.csv", row.names=1)
result2 <- result[,selected]
result2$name <- rownames(result2)
result2 <- gather(result2, module_name, flux,-name)
result2$name <- gsub("X_","",result2$name)
result3 <- inner_join(result2, currentmeta)
finalcombined <- result3
allmodules <- unique(finalcombined$module_name)
finalmodules <- NULL
celltype <- NULL
FC <- NULL
pvalue <- NULL
Treatment <- NULL
for(i in 1:length(allmodules)){
  currentframe <- finalcombined[finalcombined$module_name %in% allmodules[i],]
  allcelltypes <- unique(currentframe$celltype)
  for(j in 1:length(allcelltypes)){
    currentframe2 <- currentframe[currentframe$celltype %in% allcelltypes[j],]
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , allmodules[i])
      FC<- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"Fruc")}
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , allmodules[i])
      FC <- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"HFHS")
    }
  }
  
  #aa <- ggplot(finalcombined[finalcombined$module_name %in% allmodules[i],], aes(y = flux, x = celltype, fill = data.diseasestatus)) + 
  #  geom_boxplot(outlier.size = 0.1) +
  #  ggtitle(allmodules[i])+ylab("metabolite balance")+theme_ridges()+
  #  stat_compare_means(method = "anova",label="p.signif")+
  #  theme_classic()+
  #  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))
  #plot(aa)
}
resultvalue_balance <- data.frame(finalmodules,celltype,FC,pvalue,Treatment)
resultvalue_balance$FC[resultvalue_balance$FC > 5] <- 5
resultvalue_balance_selected <- resultvalue_balance[resultvalue_balance$celltype %in% c("Pericentral hepatocytes","Periportal hepatocytes","Kupffer","Sinusoidal endothelial cells"),]
resultvalue_balance_selected_Liver <- resultvalue_balance_selected
ggplot(resultvalue_balance_selected, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_grid(.~celltype)+
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))



#SVF
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined/SVF_Final_Set1only.rda")
load("/Users/Tsai_Lab/Desktop/shared_slides/filtered_object/SVF.rda")
module_info <- read.csv("./module_info_Fructose.csv")
colnames(module_info)[1] <- "module_name"
module_info$module_final <- paste0(module_info$Compound_IN_name,"->",module_info$Compound_OUT_name)

currentmeta <- dropEST.combined.filtered@meta.data
currentmeta$name <- rownames(currentmeta)

result <- read.csv("./SVF_All_flux_Fructose_impute.csv", row.names=1)
result2 <- result[,module_info$module_name[module_info$selected %in% "*"]]
result2$name <- rownames(result2)
result2 <- gather(result2, module_name, flux,-name)
result3 <- inner_join(result2, module_info)
result3 <- inner_join(result3, currentmeta)
finalcombined <- result3
allmodules <- unique(finalcombined$module_name)

finalmodules <- NULL
celltype <- NULL
FC <- NULL
pvalue <- NULL
Treatment <- NULL
for(i in 1:length(allmodules)){
  currentframe <- finalcombined[finalcombined$module_name %in% allmodules[i],]
  allcelltypes <- unique(currentframe$celltype)
  for(j in 1:length(allcelltypes)){
    currentframe2 <- currentframe[currentframe$celltype %in% allcelltypes[j],]
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , module_info$module_final[module_info$module_name %in%  allmodules[i]])
      FC <- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"Fruc")
    }
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , module_info$module_final[module_info$module_name %in%  allmodules[i]])
      FC <- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"HFHS")
    }
    
  }
  
  #aa <- ggplot(finalcombined[finalcombined$module_name %in% allmodules[i],], aes(y = flux, x = celltype, fill = data.diseasestatus)) + 
  #  geom_boxplot(outlier.size = 0.1) +
  #  stat_compare_means(method = "anova",label="p.signif")+
  #  theme_classic() + ggtitle(module_info$module_final[module_info$module_name %in%  allmodules[i]])+
  #  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))
  #plot(aa)
}
resultvalue_flux <- data.frame(finalmodules,celltype,FC,pvalue,Treatment)
resultvalue_flux_selected <- resultvalue_flux[resultvalue_flux$celltype %in% c("APC_Hsd11b1","APC_Pi16","Goblet","APC_Agt", "Endothelial cells","M2 Macrophage"),]
resultvalue_flux_selected_SVF <- resultvalue_flux_selected
ggplot(resultvalue_flux_selected, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_grid(.~celltype)+
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))


selected <- c("Glucose","G6P","Fatty.Acid","Fructose")
result <- read.csv("./SVF_All_balance_Fructose_impute.csv", row.names=1)
result2 <- result[,selected]
result2$name <- rownames(result2)
result2 <- gather(result2, module_name, flux,-name)
result3 <- inner_join(result2, currentmeta)
finalcombined <- result3
allmodules <- unique(finalcombined$module_name)
finalmodules <- NULL
celltype <- NULL
FC <- NULL
pvalue <- NULL
Treatment <- NULL
for(i in 1:length(allmodules)){
  currentframe <- finalcombined[finalcombined$module_name %in% allmodules[i],]
  allcelltypes <- unique(currentframe$celltype)
  for(j in 1:length(allcelltypes)){
    currentframe2 <- currentframe[currentframe$celltype %in% allcelltypes[j],]
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , allmodules[i])
      FC<- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"Fruc")}
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , allmodules[i])
      FC <- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"HFHS")
    }
  }
  
  #aa <- ggplot(finalcombined[finalcombined$module_name %in% allmodules[i],], aes(y = flux, x = celltype, fill = data.diseasestatus)) + 
  #  geom_boxplot(outlier.size = 0.1) +
  #  ggtitle(allmodules[i])+ylab("metabolite balance")+theme_ridges()+
  #  stat_compare_means(method = "anova",label="p.signif")+
  #  theme_classic()+
  #  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))
  #plot(aa)
}
resultvalue_balance <- data.frame(finalmodules,celltype,FC,pvalue,Treatment)
resultvalue_balance_selected <- resultvalue_balance[resultvalue_balance$celltype %in% c("APC_Hsd11b1","APC_Pi16","Goblet","APC_Agt", "Endothelial cells","M2 Macrophage"),]
resultvalue_balance_selected_SVF <- resultvalue_balance_selected
ggplot(resultvalue_balance_selected, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_grid(.~celltype)+
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))




#Hyp
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Hyp_combined_Set1_final.rda")
load("/Users/Tsai_Lab/Desktop/shared_slides/filtered_object/Hyp.rda")
#module_info <- read.csv("./module_info_Fructose_Hyp.csv")
module_info <- read.csv("./module_info_Fructose_Hyp.csv")
#module_info <- read.csv("./module_info_Fructose.csv")

colnames(module_info)[1] <- "module_name"
module_info$module_final <- paste0(module_info$Compound_IN_name,"->",module_info$Compound_OUT_name)

currentmeta <- dropEST.combined.filtered@meta.data
currentmeta$name <- rownames(currentmeta)

result <- read.csv("./Hyp_All_flux_Fructose_impute.csv", row.names=1)
result2 <- result[,module_info$module_name[module_info$selected %in% "*"]]
result2$name <- rownames(result2)
result2 <- gather(result2, module_name, flux,-name)
result3 <- inner_join(result2, module_info)
result3 <- inner_join(result3, currentmeta)
finalcombined <- result3
allmodules <- unique(finalcombined$module_name)
finalmodules <- NULL
celltype <- NULL
FC <- NULL
pvalue <- NULL
Treatment <- NULL
for(i in 1:length(allmodules)){
  currentframe <- finalcombined[finalcombined$module_name %in% allmodules[i],]
  allcelltypes <- unique(currentframe$celltype)
  for(j in 1:length(allcelltypes)){
    currentframe2 <- currentframe[currentframe$celltype %in% allcelltypes[j],]
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , module_info$module_final[module_info$module_name %in%  allmodules[i]])
      FC <- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"Fruc")
    }
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , module_info$module_final[module_info$module_name %in%  allmodules[i]])
      FC <- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"HFHS")
    }
    
  }
  
  #aa <- ggplot(finalcombined[finalcombined$module_name %in% allmodules[i],], aes(y = flux, x = celltype, fill = data.diseasestatus)) + 
  #  geom_boxplot(outlier.size = 0.1) +
  #  stat_compare_means(method = "anova",label="p.signif")+
  #  theme_classic() + ggtitle(module_info$module_final[module_info$module_name %in%  allmodules[i]])+
  #  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))
  #plot(aa)
}
resultvalue_flux <- data.frame(finalmodules,celltype,FC,pvalue,Treatment)
resultvalue_flux$FC[resultvalue_flux$FC > 5] <- 5
resultvalue_flux$FC[resultvalue_flux$FC < -5] <- -5
resultvalue_flux_selected <- resultvalue_flux[resultvalue_flux$celltype %in% c("GABAergic neurons","Glutamatergic neurons","Astrocytes","Endothelial cells"),]
resultvalue_flux_selected_Hyp <- resultvalue_flux_selected
ggplot(resultvalue_flux_selected_Hyp, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_grid(.~celltype)+
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))


selected <- c("Glucose","G6P","Fatty.Acid","GABA","Fructose")
result <- read.csv("./Hyp_All_balance_Fructose_impute.csv", row.names=1)
result2 <- result[,selected]
result2$name <- rownames(result2)
result2 <- gather(result2, module_name, flux,-name)
result3 <- inner_join(result2, currentmeta)
finalcombined <- result3
allmodules <- unique(finalcombined$module_name)
finalmodules <- NULL
celltype <- NULL
FC <- NULL
pvalue <- NULL
Treatment <- NULL
for(i in 1:length(allmodules)){
  currentframe <- finalcombined[finalcombined$module_name %in% allmodules[i],]
  allcelltypes <- unique(currentframe$celltype)
  for(j in 1:length(allcelltypes)){
    currentframe2 <- currentframe[currentframe$celltype %in% allcelltypes[j],]
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , allmodules[i])
      FC<- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "Fruc"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"Fruc")}
    
    if(length(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"]) > 10 & length(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]) > 10){
      celltype <- c(celltype, allcelltypes[j])
      finalmodules  <- c(finalmodules , allmodules[i])
      FC <- c(FC,(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"])-median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"]))/abs(median(currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])))
      pvalue <- c(pvalue, wilcox.test(currentframe2$flux[currentframe2$data.diseasestatus %in% "HFHS"],currentframe2$flux[currentframe2$data.diseasestatus %in% "Control"])$p.value)
      Treatment <- c(Treatment,"HFHS")
    }
  }
  
  #aa <- ggplot(finalcombined[finalcombined$module_name %in% allmodules[i],], aes(y = flux, x = celltype, fill = data.diseasestatus)) + 
  #  geom_boxplot(outlier.size = 0.1) +
  #  ggtitle(allmodules[i])+ylab("metabolite balance")+theme_ridges()+
  #  stat_compare_means(method = "anova",label="p.signif")+
  #  theme_classic()+
  #  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))
  #plot(aa)
}
resultvalue_balance <- data.frame(finalmodules,celltype,FC,pvalue,Treatment)
resultvalue_balance$FC[resultvalue_balance$FC < -5] <- -5
resultvalue_balance$FC[resultvalue_balance$FC > 5] <- 5

resultvalue_balance_selected <- resultvalue_balance[resultvalue_balance$celltype %in% c("GABAergic neurons","Glutamatergic neurons","Astrocytes","Endothelial cells"),]
resultvalue_balance_selected_Hyp <- resultvalue_balance_selected
ggplot(resultvalue_balance_selected, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_grid(.~celltype)+
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))




resultvalue_balance_selected_SI$tissue <- "SI"
resultvalue_balance_selected_SVF$tissue <- "SVF"
resultvalue_balance_selected_Liver$tissue <- "Liver"
resultvalue_balance_selected_Hyp$tissue <- "Hyp"
resultvalue_balance_selected <- rbind.data.frame(resultvalue_balance_selected_SI,resultvalue_balance_selected_SVF,
                                                 resultvalue_balance_selected_Liver,resultvalue_balance_selected_Hyp)
resultvalue_balance_selected$FC[resultvalue_balance_selected $FC < -5] <- -5
ggplot(resultvalue_balance_selected , aes(y = FC, x = finalmodules, fill = Treatment)) + facet_wrap(tissue~celltype , nrow = 2, scales = "free_x")+
  geom_bar(stat = "identity",position = "dodge") +geom_hline(yintercept = c(-1,1))+
  theme_classic() +ylab("Fold increment")+xlab("")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1))


resultvalue_flux_selected_SI$tissue <- "SI"
resultvalue_flux_selected_SVF$tissue <- "SVF"
resultvalue_flux_selected_Liver$tissue <- "Liver"
resultvalue_flux_selected_Hyp$tissue <- "Hyp"

resultvalue_flux_selected <- rbind.data.frame(resultvalue_flux_selected_SI,resultvalue_flux_selected_SVF,
                                                 resultvalue_flux_selected_Liver,resultvalue_flux_selected_Hyp)
resultvalue_flux_selected$finalmodules[resultvalue_flux_selected$finalmodules %in% "Pyruvate->Acetyl-Coa"] <- "Pyruvate->Acetyl-CoA"
resultvalue_flux_selected_metabolic <- resultvalue_flux_selected[-grep("_in",resultvalue_flux_selected$finalmodules),]
resultvalue_flux_selected <- resultvalue_flux_selected[grep("_in",resultvalue_flux_selected$finalmodules),]
resultvalue_flux_selected$finalmodules <- sapply(strsplit(resultvalue_flux_selected$finalmodules,"->",fixed = T), `[[`, 2)
resultvalue_flux_selected$FC[resultvalue_flux_selected$FC > 5] <- 5
resultvalue_flux_selected$category <- "Intake"
p2 <- ggplot(resultvalue_flux_selected, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_wrap(.~tissue+celltype, scales = "free_x", nrow = 2)+
  geom_bar(stat = "identity",position = "dodge") +geom_hline(yintercept = c(-1,1))+
  theme_classic() +ylab("Fold increment")+xlab("")+
  theme(legend.text=element_text(size=30),legend.title =element_text(size=30),plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1),
        strip.text.x = element_text(size = 15),axis.title=element_text(size=24))

resultvalue_flux_selected_metabolic$finalmodules <- factor(resultvalue_flux_selected_metabolic$finalmodules,levels = unique(resultvalue_flux_selected_metabolic$finalmodules ))
resultvalue_flux_selected_metabolic$FC[resultvalue_flux_selected_metabolic$FC > 5] <- 5
resultvalue_flux_selected_metabolic$category <- "Metabolic pathway"
p3 <- ggplot(resultvalue_flux_selected_metabolic, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_wrap(.~tissue+celltype, scales = "free_x", nrow = 2)+
  geom_bar(stat = "identity",position = "dodge") +geom_hline(yintercept = c(-1,1))+
  theme_classic() +ylab("Fold increment\n\n")+xlab("")+
  theme(legend.text=element_text(size=30),legend.title =element_text(size=30),plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18,angle = 45, hjust=1),
        plot.margin = margin(20, 20, 20, 50),strip.text.x = element_text(size = 15),axis.title=element_text(size=24))

SI_part <- rbind.data.frame(resultvalue_flux_selected[resultvalue_flux_selected$tissue %in% "SI",]
  ,resultvalue_flux_selected_metabolic[resultvalue_flux_selected_metabolic$tissue %in% "SI",])
SI_part$finalmodules <- factor(SI_part$finalmodules, levels = unique(SI_part$finalmodules))
SI_plot <- ggplot(SI_part, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_grid(celltype~category, scales = "free_x",space = "free_x", labeller = label_wrap_gen(width=10))+
  geom_bar(stat = "identity",position = "dodge") +geom_hline(yintercept = c(-1,1))+
  theme_classic() +ylab("Fold increment")+xlab("")+ggtitle("SI")+
  theme(plot.title = element_text(size=24,hjust = 0.5),legend.text=element_text(size=30),legend.title =element_text(size=30),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 20,angle = 45, hjust=1),
        strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15),axis.title=element_text(size=24),legend.position = "none") + 
  scale_fill_manual(values = c("Fruc" = "#00BFC4","HFHS" = "#F8766D"))

Liver_part <- rbind.data.frame(resultvalue_flux_selected[resultvalue_flux_selected$tissue %in% "Liver",]
                            ,resultvalue_flux_selected_metabolic[resultvalue_flux_selected_metabolic$tissue %in% "Liver",])
Liver_part$FC[Liver_part$FC < -5] <- -5
Liver_part$finalmodules <- factor(Liver_part$finalmodules, levels = unique(Liver_part$finalmodules))
Liver_plot <- ggplot(Liver_part, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_grid(celltype~category, scales = "free_x",space = "free_x", labeller = label_wrap_gen(width=10))+
  geom_bar(stat = "identity",position = "dodge") +geom_hline(yintercept = c(-1,1))+
  theme_classic() +ylab("")+xlab("")+ggtitle("Liver")+
  theme(plot.title = element_text(size=24,hjust = 0.5),legend.text=element_text(size=30),legend.title =element_text(size=30),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 20,angle = 45, hjust=1),
        strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15),axis.title=element_text(size=24),legend.position = "none") + 
  scale_fill_manual(values = c("Fruc" = "#00BFC4","HFHS" = "#F8766D"))

SVF_part <- rbind.data.frame(resultvalue_flux_selected[resultvalue_flux_selected$tissue %in% "SVF",]
                            ,resultvalue_flux_selected_metabolic[resultvalue_flux_selected_metabolic$tissue %in% "SVF",])
SVF_part$finalmodules <- factor(SVF_part$finalmodules, levels = unique(SVF_part$finalmodules))
SVF_plot <- ggplot(SVF_part, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_grid(celltype~category, scales = "free_x",space = "free_x", labeller = label_wrap_gen(width=10))+
  geom_bar(stat = "identity",position = "dodge") +geom_hline(yintercept = c(-1,1))+
  theme_classic() +ylab("")+xlab("")+ggtitle("SVF")+
  theme(plot.title = element_text(size=24,hjust = 0.5),legend.text=element_text(size=30),legend.title =element_text(size=30),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 20,angle = 45, hjust=1),
        strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15),axis.title=element_text(size=24),legend.position = "none") + 
  scale_fill_manual(values = c("Fruc" = "#00BFC4","HFHS" = "#F8766D"))
Hyp_part <- rbind.data.frame(resultvalue_flux_selected[resultvalue_flux_selected$tissue %in% "Hyp",]
                            ,resultvalue_flux_selected_metabolic[resultvalue_flux_selected_metabolic$tissue %in% "Hyp",])
Hyp_part$finalmodules <- factor(Hyp_part$finalmodules, levels = unique(Hyp_part$finalmodules))
Hyp_plot <- ggplot(Hyp_part, aes(y = FC, x = finalmodules, fill = Treatment)) + facet_grid(celltype~category, scales = "free_x",space = "free_x", labeller = label_wrap_gen(width=10))+
  geom_bar(stat = "identity",position = "dodge") +geom_hline(yintercept = c(-1,1))+
  theme_classic() +ylab("")+xlab("")+ggtitle("Hyp")+
  theme(plot.title = element_text(size=24,hjust = 0.5),legend.text=element_text(size=30),legend.title =element_text(size=30),axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 20,angle = 45, hjust=1),
        strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15),axis.title=element_text(size=24)) + 
  scale_fill_manual(values = c("Fruc" = "#00BFC4","HFHS" = "#F8766D"))
library(gridExtra)
#pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/combined_scFEA_result_new.pdf",width = 18, height = 24)
pdf("/Users/Tsai_Lab/Downloads/old_combined2.pdf",width = 18, height = 24)
grid.arrange(SI_plot,Liver_plot,SVF_plot,Hyp_plot,nrow = 2)
dev.off()
library(svglite)
finalgg <- grid.arrange(SI_plot,Liver_plot,SVF_plot,Hyp_plot,nrow = 2,heights = c(1,1.4), widths = c(1,1.4))
ggsave(file="/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/combined_scFEA_result.svg", plot=finalgg,width =18, height = 24)

#p3
pdf("flux_intake.pdf",width = 20,height = 10)
print(p2)
dev.off()
pdf("flux_path.pdf",width = 28,height = 15)
print(p3)
dev.off()
#Old version, the Treatments were processed seperately
module_info <- read.csv("./module_info_Fructose.csv")
colnames(module_info)[1] <- "module_name"

rm("finalcombined")
for(i in 1:length(allTreatment)){
  result <- read.csv(paste0("Hyp_",allTreatment[i],"_flux.csv"), row.names=1)
  currentmeta <- dropEST.combined.filtered@meta.data
  currentmeta <- currentmeta[currentmeta$data.diseasestatus %in% allTreatment[i],]
  currentmeta <- currentmeta[,c("data.diseasestatus","celltype")]
  currentmeta$name <- rownames(currentmeta)
  result2 <- result[,module_info$module_name[module_info$selected %in% "*"]]
  result2$name <- rownames(result2)
  result2 <- gather(result2, module_name, flux,-name)

  result3 <- inner_join(result2, module_info)
  result3 <- inner_join(result3, currentmeta)
  
  result3$Treatment <- allTreatment[i]
  if(!exists("finalcombined")){
    finalcombined <- result3
  }else{
    finalcombined <- rbind.data.frame(finalcombined, result3)
  }
}

allmodules <- unique(finalcombined$module_name)


rm("finalcombined")
selected <- c("Fatty.Acid","G6P","Glucose","Lactate","Cholesterol","G3P","Glucose.1.phosphate","Citrate")
for(i in 1:length(allTreatment)){
  result <- read.csv(paste0("Hyp_",allTreatment[i],"_balance.csv"), row.names=1)
  currentmeta <- dropEST.combined.filtered@meta.data
  currentmeta <- currentmeta[currentmeta$data.diseasestatus %in% allTreatment[i],]
  currentmeta <- currentmeta[,c("data.diseasestatus","celltype")]
  currentmeta$name <- rownames(currentmeta)
  result2 <- result[,selected]
  result2$name <- rownames(result2)
  result2 <- gather(result2, module_name, balance,-name)
  result3 <- inner_join(result2, currentmeta)
  result3$Treatment <- allTreatment[i]
  if(!exists("finalcombined")){
    finalcombined <- result3
  }else{
    finalcombined <- rbind.data.frame(finalcombined, result3)
  }
}

allmodules <- unique(finalcombined$module_name)
for(i in 1:length(allmodules)){
  aa <- ggplot(finalcombined[finalcombined$module_name %in% allmodules[i],], aes(x = balance, y = celltype, fill= Treatment)) + 
    geom_density_ridges(alpha = 0.2) + xlim(c(-0.1,0.1))+
    theme_ridges() + ggtitle(allmodules[i])+
    theme(plot.title = element_text(hjust = 0.5))
  plot(aa)
}
