library(dplyr)
#ligand_matching
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/")
#Final_fructose_ligand <- read.csv("./Final_fructose_ligand.csv")
Final_fructose_ligand <- read.csv("./Final_fructose_ligand_new.csv")

colnames(Final_fructose_ligand)[1] <- "Gene"
Table_3a <- read.csv("./42255_2021_478_MOESM4_ESM_3a.csv")
Table_3a <- Table_3a[!duplicated(Table_3a$Entrez.ID),]
unique(Table_3a$PGS)
#Table_3a <- Table_3a[!Table_3a$PGS %in% c("CKD","AF"),]
#Table_3a <- Table_3a[Table_3a$PGS %in% "T2D",]
#Table_3a <- Table_3a[Table_3a$PGS %in% "CAD",]
#Table_3a <- Table_3a[Table_3a$PGS %in% "IS",]
#Table_3a <- Table_3a[Table_3a$PGS %in% "CKD",]

length(unique(Table_3a$Gene))
combined <- left_join(Final_fructose_ligand,Table_3a)
unique(combined$Gene)
unique(Final_fructose_ligand$Gene)[!unique(Final_fructose_ligand$Gene) %in% unique(combined$Gene[!is.na(combined$PGS)])]
wilcox.test(combined$P.value,Table_3a$P.value, na.rm = T)
wilcox.test(combined$P.value[combined$treatment %in% "Fructose"],Table_3a$P.value, na.rm = T)
wilcox.test(combined$P.value[combined$treatment %in% "HFHS"],Table_3a$P.value, na.rm = T)
combined2 <- combined[combined$P.value < 0.05,]
mean(combined$P.value < 0.05,na.rm = T)
mean(Table_3a$P.value < 0.05)
testmat <- matrix(c(sum(combined$P.value < 0.05,na.rm = T), sum(combined$P.value > 0.05,na.rm = T), 
                  sum(Table_3a$P.value < 0.05),sum(Table_3a$P.value > 0.05)),ncol = 2)
fisher.test(testmat)


combined4 <- combined[combined$treatment %in% "Fruc",]
testmat <- matrix(c(sum(combined4$P.value < 0.05,na.rm = T), sum(combined4$P.value > 0.05,na.rm = T), 
                    sum(Table_3a$P.value < 0.05),sum(Table_3a$P.value > 0.05)),ncol = 2)
fisher.test(testmat)

combined4 <- combined[combined$treatment %in% "HFHS",]
testmat <- matrix(c(sum(combined$P.value < 0.05,na.rm = T), sum(combined$P.value > 0.05,na.rm = T), 
                    sum(Table_3a$P.value < 0.05),sum(Table_3a$P.value > 0.05)),ncol = 2)
fisher.test(testmat)
#11.2% vs 7.5% pvalue< 0.05

setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/")
alldisease <- c("T2D","CAD","AF","IS","CKD")
alltreatment <- list(c("Fructose","HFHS"),"HFHS","Fructose")
allwilcoxPcomp <- NULL
allpropcomp <- NULL
#allpvals <- NULL
sigpropsc <- NULL
sigpropall <- NULL
finaldisease <- NULL
finaltreatment <- NULL

for(a in alldisease){
  Final_fructose_ligand <- read.csv("./Final_fructose_ligand_new.csv")
  colnames(Final_fructose_ligand)[1] <- "Gene"
  Table_3a <- read.csv("./42255_2021_478_MOESM4_ESM_3a.csv")
  #Table_3a <- Table_3a[!duplicated(Table_3a$Entrez.ID),]
  unique(Table_3a$PGS)
  Table_3a <- Table_3a[Table_3a$PGS %in% a,]
  Table_3a <- Table_3a[!duplicated(Table_3a$Entrez.ID),]
  
  for(j in 1:length(alltreatment)){
    combined <- left_join(Final_fructose_ligand,Table_3a)
    finaldisease <- c(finaldisease, a)
    currenttreatment <- alltreatment[[j]]
    finaltreatment <- c(finaltreatment , paste0(currenttreatment,collapse = "|"))
    combined <- combined[combined$treatment %in% currenttreatment,]
    allwilcoxPcomp <- c(allwilcoxPcomp, wilcox.test(combined$P.value,Table_3a$P.value, na.rm = T)$p.value)
    sigpropsc <- c(sigpropsc, mean(combined$P.value < 0.05,na.rm = T))
    sigpropall <- c(sigpropall,  mean(Table_3a$P.value < 0.05))
    
    testmat <- matrix(c(sum(combined$P.value < 0.05,na.rm = T), sum(combined$P.value > 0.05,na.rm = T), 
                        sum(Table_3a$P.value < 0.05),sum(Table_3a$P.value > 0.05)),ncol = 2)
    allpropcomp  <- c(allpropcomp , fisher.test(testmat)$p.value)
    
  }
}
allresult <- data.frame(finaltreatment,finaldisease,
                        sigpropsc, sigpropall,
                        allpropcomp,allwilcoxPcomp )

#Going through all tissue
alldisease <- c("T2D","CAD","AF","IS","CKD")
alltreatment <- list(c("Fructose","HFHS"),"HFHS","Fructose")
allwilcoxPcomp <- NULL
allpropcomp <- NULL
#allpvals <- NULL
sigpropsc <- NULL
sigpropall <- NULL
finaldisease <- NULL
finaltreatment <- NULL
finaltissues <- NULL
for(a in alldisease){
  Final_fructose_ligand <- read.csv("./Final_fructose_ligand.csv")
  colnames(Final_fructose_ligand)[1] <- "Gene"
  Table_3a <- read.csv("./42255_2021_478_MOESM4_ESM_3a.csv")
  unique(Table_3a$PGS)
  Table_3a <- Table_3a[Table_3a$PGS %in% a,]
  
  
  for(j in 1:length(alltreatment)){
    combined <- left_join(Final_fructose_ligand,Table_3a)
    finaldisease <- c(finaldisease, a)
    currenttreatment <- alltreatment[[j]]
    finaltreatment <- c(finaltreatment , paste0(currenttreatment,collapse = "|"))
    combined <- combined[combined$treatment %in% currenttreatment,]
    alltissues <- unique(combined$tissue)
    combined_orig <- combined
    for(t in alltissues){
      finaltissues <- c(finaltissues,t)
      combined <- combined_orig[combined_orig$tissue %in% t,]
      allwilcoxPcomp <- c(allwilcoxPcomp, wilcox.test(combined$P.value,Table_3a$P.value, na.rm = T)$p.value)
      sigpropsc <- c(sigpropsc, mean(combined$P.value < 0.05,na.rm = T))
      sigpropall <- c(sigpropall,  mean(Table_3a$P.value < 0.05))
      
      testmat <- matrix(c(sum(combined$P.value < 0.05,na.rm = T), sum(combined$P.value > 0.05,na.rm = T), 
                          sum(Table_3a$P.value < 0.05),sum(Table_3a$P.value > 0.05)),ncol = 2)
      allpropcomp  <- c(allpropcomp , fisher.test(testmat)$p.value)

    }
   }
}
allresult_organspec <- data.frame(finaltreatment,finaldisease,
                        sigpropsc, sigpropall,
                        allpropcomp,allwilcoxPcomp,finaltissues )


#11.2% vs 7.5% pvalue< 0.05




alldisease <- unique(Table_3a$PGS)
allres <- NULL
alldisease2 <- NULL
percentsc <- NULL
percentack <- NULL
wilxocP <- NULL
for(a in alldisease){
  alldisease2  <- c(alldisease2  , a)
  currentcombined <- combined[combined$PGS %in% a,]
  currenttable <- Table_3a[Table_3a$PGS %in% a,]
  combined3 <- currentcombined %>% dplyr::group_by(Gene) %>% dplyr::slice(which.min(P.value))
  Table_3a2 <- currenttable  %>% dplyr::group_by(Gene) %>% dplyr::slice(which.min(P.value))
  testmat <- matrix(c(sum(combined3$P.value < 0.05,na.rm = T), sum(combined3$P.value > 0.05,na.rm = T), 
                      sum(Table_3a2$P.value < 0.05),sum(Table_3a2$P.value > 0.05)),ncol = 2)
  
  wilxocP <- c(wilxocP, wilcox.test(combined3$P.value,Table_3a2$P.value, na.rm = T)$p.value)
  allres <- c(allres,fisher.test(testmat)$p.value)
  percentsc <- c(percentsc , mean(combined3$P.value < 0.05,na.rm = T))
  percentack <- c(percentack, mean(Table_3a2$P.value < 0.05))
}
final_all <- data.frame(alldisease2,wilxocP,Prop005sigP = allres,percentsc,percentack)




alldisease <- unique(Table_3a$PGS)
allres <- NULL
alldisease2 <- NULL
percentsc <- NULL
percentack <- NULL
wilxocP <- NULL
for(a in alldisease){
  alldisease2  <- c(alldisease2  , a)
  currentcombined <- combined[combined$PGS %in% a & combined$treatment %in% "Fructose",]
  currenttable <- Table_3a[Table_3a$PGS %in% a,]
  combined3 <- currentcombined %>% dplyr::group_by(Gene) %>% dplyr::slice(which.min(P.value))
  Table_3a2 <- currenttable  %>% dplyr::group_by(Gene) %>% dplyr::slice(which.min(P.value))
  testmat <- matrix(c(sum(combined3$P.value < 0.05,na.rm = T), sum(combined3$P.value > 0.05,na.rm = T), 
                      sum(Table_3a2$P.value < 0.05),sum(Table_3a2$P.value > 0.05)),ncol = 2)
 
  wilxocP <- c(wilxocP, wilcox.test(combined3$P.value,Table_3a2$P.value, na.rm = T)$p.value)
  allres <- c(allres,fisher.test(testmat)$p.value)
  percentsc <- c(percentsc , mean(combined3$P.value < 0.05,na.rm = T))
  percentack <- c(percentack, mean(Table_3a2$P.value < 0.05))
}
final_Fruc <- data.frame(alldisease2,wilxocP,Prop005sigP = allres,percentsc,percentack)


allres <- NULL
alldisease2 <- NULL
percentsc <- NULL
percentack <- NULL
wilxocP <- NULL
for(a in alldisease){
  alldisease2  <- c(alldisease2  , a)
  currentcombined <- combined[combined$PGS %in% a & combined$treatment %in% "HFHS",]
  currenttable <- Table_3a[Table_3a$PGS %in% a,]
  combined3 <- currentcombined %>% dplyr::group_by(Gene) %>% dplyr::slice(which.min(P.value))
  Table_3a2 <- currenttable  %>% dplyr::group_by(Gene) %>% dplyr::slice(which.min(P.value))
  testmat <- matrix(c(sum(combined3$P.value < 0.05,na.rm = T), sum(combined3$P.value > 0.05,na.rm = T), 
                      sum(Table_3a2$P.value < 0.05),sum(Table_3a2$P.value > 0.05)),ncol = 2)
  
  wilxocP <- c(wilxocP, wilcox.test(combined3$P.value,Table_3a2$P.value, na.rm = T)$p.value)
  allres <- c(allres,fisher.test(testmat)$p.value)
  percentsc <- c(percentsc , mean(combined3$P.value < 0.05,na.rm = T))
  percentack <- c(percentack, mean(Table_3a2$P.value < 0.05))
}
final_HFHS <- data.frame(alldisease2,wilxocP,Prop005sigP = allres,percentsc,percentack)

#combined3 <- combined %>% dplyr::group_by(Gene) %>% dplyr::slice(which.min(P.value))
#Table_3a2 <- Table_3a %>% dplyr::group_by(Gene) %>% dplyr::slice(which.min(P.value))
#testmat <- matrix(c(sum(combined3$P.value < 0.05,na.rm = T), sum(combined3$P.value > 0.05,na.rm = T), 
#                    sum(Table_3a2$P.value < 0.05),sum(Table_3a2$P.value > 0.05)),ncol = 2)
#fisher.test(testmat)


Table_3c <- read.csv("./42255_2021_478_MOESM4_ESM_3c.csv")
unique(Table_3c$Protein)
colnames(Final_fructose_ligand)[1] <- "Protein"
combined2 <- left_join(Final_fructose_ligand,Table_3c)
