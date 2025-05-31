#plotting treatment effect
library(ggplot2)
library(ggrepel)
library(dplyr)
library(grid)
library(gtable)
library("ggplotify")
#wdlist <- c("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined",
#            "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined",
#            "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_NPC",
#            "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1",
#            "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
#tissuename <- c("SI","SVF","NPC","HEP","HYP","HYP_Neuron")

wdlist <- c("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/",
            "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined/",
              "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Set1_",
            "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Set2_",
            "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/")
tissuename <- c("SI","SVF","Liver_Set1","Liver_Set2","HYP")

Fructoseloc = paste0(wdlist,"Fructosedistanceresultzscore.rda")
#Fructoseloc = c(Fructoseloc,"/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Fructosedistanceresult_subsetzscore.rda")
rm(Fructosetoplot, finalorigdistane, combinealldistance)

for(i in 1:length(Fructoseloc )){
  load(Fructoseloc[i])
  toPlot$tissue = tissuename[i]
  toPlot$celltype <- rownames(toPlot)
  tissuedistance = melted.combinedPermutationDistances[,c(2,4)]
  tissuedistance <- tissuedistance[!duplicated(tissuedistance),]
  tissuedistance$tissue <- tissuename[i]
  
  if(!exists("Fructosetoplot")){
    Fructosetoplot = toPlot
    finalorigdistane <- tissuedistance
    combinealldistance = unlist(cellTypePermutationDistances)
    
  }else{
    Fructosetoplot = rbind.data.frame(Fructosetoplot, toPlot)
    finalorigdistane = rbind.data.frame(finalorigdistane, tissuedistance)
    combinealldistance <- c( combinealldistance, unlist(cellTypePermutationDistances))
  }
}

Liverpart <- Fructosetoplot[grep("Liver",Fructosetoplot$tissue),]
#Liverpart$celltype <- gsub("1|2","",rownames(Liverpart))
Liverpart <- Liverpart %>% 
  dplyr::group_by(celltype) %>% 
  dplyr::summarise(logFC = mean(logFC),pval = mean(pval))
Liverpart$tissue = "Liver"
Liverpart <- as.data.frame(Liverpart)
#rownames(Liverpart) <- Liverpart$celltype

Fructosetoplot <- Fructosetoplot[-grep("Liver",Fructosetoplot$tissue),] 
Fructosetoplot <- rbind.data.frame(Fructosetoplot, Liverpart)
#FDR <- NULL
#for(i in 1:nrow(finalorigdistane)){
#  FDR[i] <- (sum((finalorigdistane$distances[i] < combinealldistance))+1)/length(combinealldistance)
#}
#finalorigdistane$FDR <- FDR
Fructosetoplot$L1 <- Fructosetoplot$celltype
#Fructosetoplot <- inner_join(Fructosetoplot, finalorigdistane)

Fructosetoplot <- Fructosetoplot[-grep("UMI|Unk",Fructosetoplot$L1),]

#Fructosetoplot$FDR <- -log10(Fructosetoplot$FDR)

Fructosetoplot$FDR2 <- -log10(p.adjust(10^(-Fructosetoplot$pval), "bonferroni"))


p2 <- ggplot(Fructosetoplot, aes(x = pval, y = logFC, color = tissue)) +
  geom_point() +
  theme_classic(base_size = 25) +
  xlab("-log10 FDR") + 
  geom_text_repel(aes(label = Fructosetoplot$celltype), size = 3.5) 
print(p2)


HFHSloc = paste0(wdlist,"HFHSdistanceresultscore.rda")
#HFHSloc = c(HFHSloc,"/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/HFHSdistanceresult_subsetzscore.rda")
rm(HFHStoplot, finalorigdistane, combinealldistance)

for(i in 1:length(HFHSloc )){
  load(HFHSloc[i])
  toPlot$tissue = tissuename[i]
  toPlot$celltype <- rownames(toPlot)
  
  tissuedistance = melted.combinedPermutationDistances[,c(2,4)]
  tissuedistance <- tissuedistance[!duplicated(tissuedistance),]
  tissuedistance$tissue <- tissuename[i]
  
  if(!exists("HFHStoplot")){
    HFHStoplot = toPlot
    finalorigdistane <- tissuedistance
    combinealldistance = unlist(cellTypePermutationDistances)
    
  }else{
    HFHStoplot = rbind.data.frame(HFHStoplot, toPlot)
    finalorigdistane = rbind.data.frame(finalorigdistane, tissuedistance)
    combinealldistance <- c( combinealldistance, unlist(cellTypePermutationDistances))
  }
}

Liverpart <- HFHStoplot[grep("Liver",HFHStoplot$tissue),]
#Liverpart$celltype <- gsub("1|2","",rownames(Liverpart))
Liverpart <- Liverpart %>% 
  dplyr::group_by(celltype) %>% 
  dplyr::summarise(logFC = mean(logFC),pval = mean(pval))
Liverpart$tissue = "Liver"
Liverpart <- as.data.frame(Liverpart)
#rownames(Liverpart) <- Liverpart$celltype
#Liverpart <- Liverpart[,c("logFC","pval","tissue")]

HFHStoplot <- HFHStoplot[-grep("Liver",HFHStoplot$tissue),] 
HFHStoplot <- rbind.data.frame(HFHStoplot, Liverpart)

#FDR <- NULL
#for(i in 1:nrow(finalorigdistane)){
#  FDR[i] <- (sum((finalorigdistane$distances[i] < combinealldistance))+1)/length(combinealldistance)
#}
#finalorigdistane$FDR <- FDR
HFHStoplot$L1 <- HFHStoplot$celltype
#HFHStoplot <- inner_join(HFHStoplot, finalorigdistane)
HFHStoplot$FDR2 <- -log10(p.adjust(10^(-HFHStoplot$pval), "bonferroni"))

HFHStoplot <- HFHStoplot[-grep("UMI|Unk",HFHStoplot$L1),]
#HFHStoplot$celltype <- HFHStoplot$L1
#HFHStoplot$FDR <- -log10(HFHStoplot$FDR)


p3 <- ggplot(HFHStoplot, aes(x = pval, y = logFC, color = tissue)) +
  geom_point() +
  theme_classic(base_size = 25) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = HFHStoplot$celltype), size = 3.5) 
print(p3)







Fructosetoplot$treatment <- "Fructose"
Fructosetoplot <- Fructosetoplot[!Fructosetoplot$tissue %in% "HYP_Neuron",]
HFHStoplot$treatment <- "HFHS"
HFHStoplot <- HFHStoplot[!HFHStoplot$tissue %in% "HYP_Neuron",]
#FrucvHFHStoplot$treatment <- "FrucvHFHS"
#FrucvHFHStoplot <- FrucvHFHStoplot[!FrucvHFHStoplot$tissue %in% "HYP_Neuron",]
final <- rbind.data.frame(Fructosetoplot,HFHStoplot)
#final <- rbind.data.frame(Fructosetoplot,HFHStoplot,FrucvHFHStoplot)

#final$tissue[final$tissue %in% c("NPC","HEP")] <- "Liver"
#final <- final[!final$tissue %in% "HYP_Neuron",]
final$tissue[final$tissue %in% c("HYP","HYP_Neuron")] <- "Hypothalamus"
final$tissue[final$tissue %in% c("SVF")] <- "Adipose"

#final <- final[!final$L1 %in% c("hepatocyte"),]
#final$L1[final$L1 %in% c("Neural-cells")] <- "Total Neural cells"
#final$L1[final$L1 %in% c("Centralvein")] <- "Hepatocyte_Centralvein"
#final$L1[final$L1 %in% c("Perinodal")] <- "Hepatocyte_Perinodal"
final <- final[order(final$logFC,decreasing = T),]


final$L1[final$L1 %in% "oligodendrocytes precursor"] <- "Oligodendrocytes precursor"
final$L1[final$L1 %in% "Gabaergic neurons"] <- "GABAergic neurons"
allcelltypes <- unique(final$L1)
final$L1 <- factor(final$L1, levels = allcelltypes)
final$logFC <- log2(10^(final$logFC))
final$treatment[final$FDR2 < 1.3] <- paste0(final$treatment[final$FDR2 < 1.3],"_non-significant")

#p4 <- ggplot(final, aes(x = FDR2, y = logFC, color = tissue)) +
#  geom_point() +
#  theme_classic(base_size = 25) +
#  xlab("-log10 FDR") + 
#  geom_text_repel(aes(label = final$celltype), size = 4) + facet_grid(. ~ treatment)

test <- final[!final$tissue %in% "Hypothalamus",]

p5 <- ggplot(final[!final$tissue %in% "Hypothalamus",], aes(L1, logFC, color = treatment)) +
  geom_point(aes(size=pval)) +
  #geom_point(data = nonSignificant, aes(size=FDR2), color="grey") +
  # geom_point(data = nonSignificant, aes(size=FDR2, color=significance)) +
  facet_grid(.~tissue,scales="free", space="free") +
  #scale_color_manual(values=c("blue", "red", "gray"))+
  scale_color_manual(values=c("blue","cyan", "red", "bisque"))+
  #scale_x_discrete(labels= function(x) str_replace(x, "_.+$", "")) +
  theme_bw(base_size = 20) +
  xlab("") +
  guides(color=guide_legend(title="Treatment",override.aes = list(size=6)),
         size=guide_legend(title="-log10(p value)",override.aes = list(color="black"))) +
  theme(  axis.text=element_text(size=18),plot.margin = margin(0.25, 0.25, 0, 2, "cm"),
    axis.text.x = element_text( angle = 45, hjust = 1, vjust = 1), legend.position = "none")
p5
p6 <- ggplot(final[final$tissue %in% "Hypothalamus",], aes(L1, logFC, color = treatment)) +
  geom_point(aes(size=pval)) +
  #geom_point(data = nonSignificant, aes(size=FDR2), color="grey") +
  # geom_point(data = nonSignificant, aes(size=FDR2, color=significance)) +
  facet_grid(.~tissue,scales="free", space="free") +
  #scale_color_manual(values=c("blue", "red", "gray"))+
  scale_color_manual(values=c("blue","cyan", "red", "bisque"))+
  #scale_x_discrete(labels= function(x) str_replace(x, "_.+$", "")) +
  theme_bw(base_size = 20) +
  xlab("") +
  guides(color=guide_legend(title="Treatment",override.aes = list(size=6)),
         size=guide_legend(title="-log10(p value)",override.aes = list(color="black"))) +
  theme(  axis.text=element_text(size=18),
          axis.text.x = element_text( angle = 45, hjust = 1, vjust = 1))
library(gridExtra)


#save(p5,p6, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/EucdistanceFig.rda")
pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Euclideanplot_zscore_newcluster.pdf",width = 24, height = 12)
print(grid.arrange(p5,p6,layout_matrix = rbind(c(1,1,1,1,1),c(2,2,2,2,NA))))
dev.off()


sort_frame <- function(dataframe,tissue){
  dataframe <- dataframe[dataframe$tissue %in% tissue,]
  dataframe <- dataframe[order(dataframe$logFC, decreasing = T),]
  dataframe$celltype <- as.character(dataframe$celltype)
  dataframe$celltype <- factor(dataframe$celltype, levels = unique(dataframe$celltype))
  p_current <- ggplot(dataframe, aes(celltype, logFC, color = treatment)) +
    geom_point(aes(size=pval)) +
    #geom_point(data = nonSignificant, aes(size=FDR2), color="grey") +
    # geom_point(data = nonSignificant, aes(size=FDR2, color=significance)) +
    facet_grid(.~tissue,scales="free", space="free") +
    #scale_color_manual(values=c("blue", "red", "gray"))+
    scale_color_manual(values=c("blue","cyan", "red", "bisque"))+
    #scale_x_discrete(labels= function(x) str_replace(x, "_.+$", "")) +
    theme_bw(base_size = 20) +
    xlab("") +
    guides(color=guide_legend(title="Treatment",override.aes = list(size=6)),
           size=guide_legend(title="-log10(p value)",override.aes = list(color="black"))) +
    theme(  axis.text=element_text(size=18),plot.margin = margin(0.25, 0.25, 1, 0.25, "cm"),
            axis.text.x = element_text( angle = 45, hjust = 1, vjust = 1), legend.position = "none")
  
  
  return(p_current)
}
sort_frame_Adipose <- function(dataframe,tissue){
  dataframe <- dataframe[dataframe$tissue %in% tissue,]
  dataframe <- dataframe[order(dataframe$logFC, decreasing = T),]
  dataframe$celltype <- as.character(dataframe$celltype)
  dataframe$celltype <- factor(dataframe$celltype, levels = unique(dataframe$celltype))
  p_current <- ggplot(dataframe, aes(celltype, logFC, color = treatment)) +
    geom_point(aes(size=pval)) +
    #geom_point(data = nonSignificant, aes(size=FDR2), color="grey") +
    # geom_point(data = nonSignificant, aes(size=FDR2, color=significance)) +
    facet_grid(.~tissue,scales="free", space="free") +
    #scale_color_manual(values=c("blue", "red", "gray"))+
    scale_color_manual(values=c("blue", "red", "bisque"))+
    #scale_x_discrete(labels= function(x) str_replace(x, "_.+$", "")) +
    theme_bw(base_size = 20) +
    xlab("") +
    guides(color=guide_legend(title="Treatment",override.aes = list(size=6)),
           size=guide_legend(title="-log10(p value)",override.aes = list(color="black"))) +
    theme(  axis.text=element_text(size=18),plot.margin = margin(0.25, 0.25, 1, 2, "cm"),
            axis.text.x = element_text( angle = 45, hjust = 1, vjust = 1), legend.position = "none")
  
  
  return(p_current)
}
sort_frame_Hyp <- function(dataframe,tissue){
  dataframe <- dataframe[dataframe$tissue %in% tissue,]
  dataframe <- dataframe[order(dataframe$logFC, decreasing = T),]
  dataframe$celltype <- as.character(dataframe$celltype)
  dataframe$celltype <- factor(dataframe$celltype, levels = unique(dataframe$celltype))
  p_current <- ggplot(dataframe, aes(celltype, logFC, color = treatment)) +
    geom_point(aes(size=pval)) +
    #geom_point(data = nonSignificant, aes(size=FDR2), color="grey") +
    # geom_point(data = nonSignificant, aes(size=FDR2, color=significance)) +
    facet_grid(.~tissue,scales="free", space="free") +
    #scale_color_manual(values=c("blue", "red", "gray"))+
    scale_color_manual(values=c("blue","cyan", "red", "bisque"))+
    #scale_x_discrete(labels= function(x) str_replace(x, "_.+$", "")) +
    theme_bw(base_size = 20) +
    xlab("") +
    guides(color=guide_legend(title="Treatment",override.aes = list(size=6)),
           size=guide_legend(title="-log10(p value)",override.aes = list(color="black"))) +
    theme(  axis.text=element_text(size=18),plot.margin = margin(0.25, 0.25, 0, 2, "cm"),
            axis.text.x = element_text( angle = 45, hjust = 1, vjust = 1))
  
  return(p_current)
}
p6 <- sort_frame_Hyp(final,"Hypothalamus")
final_Adipose <- ggplotGrob(sort_frame_Adipose(final,"Adipose"))
final_Liver <- ggplotGrob(sort_frame(final,"Liver"))
final_SI <- ggplotGrob(sort_frame(final,"SI"))

test <- gtable:::cbind.gtable(final_Adipose,final_Liver,final_SI, size = "first")
p5 <- as.ggplot(test)


save(p5,p6, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/EucdistanceFig.rda")

#pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Euclideanplot_zscore_newcluster2.pdf",width = 12, height = 10)
#print(grid.arrange(p5,p6,layout_matrix = rbind(c(1,1,1),c(2,2,2))))
#dev.off()
pdf("/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/Euclideanplot_zscore_newcluster3.pdf",width = 16, height = 10)
print(grid.arrange(p5,p6,layout_matrix = rbind(c(1,1,1),c(2,2,2))))
dev.off()



p = arrangeGrob(ggarrange(list(final_Adipose,final_Liver,final_SI,final_Hyp), ncol=length(plist)),
                textGrob("Item"), heights=c(20,1))
grid.newpage()
grid.draw(p)


sigFruc <- Fructosetoplot$L1[Fructosetoplot$logFC > 0.5 & Fructosetoplot$pval > 2]
sigHFHS <- HFHStoplot$L1[HFHStoplot$logFC > 0.5 & HFHStoplot$pval > 2]
intersect(sigFruc,sigHFHS)
setdiff(sigFruc,sigHFHS)
setdiff(sigHFHS,sigFruc)



#Fruc vs HFHS part

#Frucvs HFHS part
FrucvHFHSloc = paste0(wdlist,"/FrucvHFHSdistanceresultzscore.rda")
FrucvHFHSloc = c(FrucvHFHSloc,"/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/FrucvHFHSdistanceresult_subsetzscore.rda")
rm(FrucvHFHStoplot, finalorigdistane, combinealldistance)

for(i in 1:length(FrucvHFHSloc )){
  load(FrucvHFHSloc[i])
  toPlot$tissue = tissuename[i]
  
  tissuedistance = melted.combinedPermutationDistances[,c(2,4)]
  tissuedistance <- tissuedistance[!duplicated(tissuedistance),]
  tissuedistance$tissue <- tissuename[i]
  
  if(!exists("FrucvHFHStoplot")){
    FrucvHFHStoplot = toPlot
    finalorigdistane <- tissuedistance
    combinealldistance = unlist(cellTypePermutationDistances)
    
  }else{
    FrucvHFHStoplot = rbind.data.frame(FrucvHFHStoplot, toPlot)
    finalorigdistane = rbind.data.frame(finalorigdistane, tissuedistance)
    combinealldistance <- c( combinealldistance, unlist(cellTypePermutationDistances))
  }
}

FDR <- NULL
for(i in 1:nrow(finalorigdistane)){
  FDR[i] <- (sum((finalorigdistane$distances[i] < combinealldistance))+1)/length(combinealldistance)
}
finalorigdistane$FDR <- FDR
FrucvHFHStoplot$L1 <- rownames(FrucvHFHStoplot)
FrucvHFHStoplot <- inner_join(FrucvHFHStoplot, finalorigdistane)
FrucvHFHStoplot$FDR2 <- -log10(p.adjust(10^(-FrucvHFHStoplot$pval), "bonferroni"))

FrucvHFHStoplot <- FrucvHFHStoplot[-grep("UMI|Unk",FrucvHFHStoplot$L1),]
FrucvHFHStoplot$celltype <- gsub('1|2',"",FrucvHFHStoplot$L1)
FrucvHFHStoplot$FDR <- -log10(FrucvHFHStoplot$FDR)


p7 <- ggplot(FrucvHFHStoplot, aes(x = pval, y = logFC, color = tissue)) +
  geom_point() +
  theme_classic(base_size = 25) +
  xlab("-log10 pval") + 
  geom_text_repel(aes(label = FrucvHFHStoplot$celltype), size = 3.5) 
print(p7)

