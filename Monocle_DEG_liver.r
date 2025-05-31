library(monocle)
#library(monocle3)
library(Seurat)
library(plyr)
library(dplyr)
library(enrichR)
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER")
#Monocle new analysis with regards in correcting batch effect
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
differentialGeneTest_extract_coef2 <- function (cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)", 
                                                reducedModelFormulaStr = "~1", relative_expr = TRUE, cores = 1, 
                                                verbose = FALSE) 
{
  status <- NA
  if (class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  all_vars <- c(all.vars(formula(fullModelFormulaStr)), all.vars(formula(reducedModelFormulaStr)))
  pd <- pData(cds)
  for (i in all_vars) {
    x <- pd[, i]
    if (any((c(Inf, NaN, NA) %in% x))) {
      stop("Error: Inf, NaN, or NA values were located in pData of cds in columns mentioned in model terms")
    }
  }
  if (relative_expr && cds@expressionFamily@vfamily %in% c("negbinomial", 
                                                           "negbinomial.size")) {
    if (is.null(sizeFactors(cds)) || sum(is.na(sizeFactors(cds)))) {
      stop("Error: to call this function with relative_expr==TRUE, you must first call estimateSizeFactors() on the CellDataSet.")
    }
  }
  if (cores > 1) {
    diff_test_res <- mcesApply(cds, 1, diff_test_helper_extract_coef2, 
                               c("BiocGenerics", "VGAM", "Matrix"), cores = cores, 
                               fullModelFormulaStr = fullModelFormulaStr, reducedModelFormulaStr = reducedModelFormulaStr, 
                               expressionFamily = cds@expressionFamily, relative_expr = relative_expr, 
                               disp_func = cds@dispFitInfo[["blind"]]$disp_func, 
                               verbose = verbose)
    diff_test_res
  }
  else {
    diff_test_res <- monocle:::smartEsApply(cds, 1, diff_test_helper_extract_coef2, 
                                            convert_to_dense = TRUE, fullModelFormulaStr = fullModelFormulaStr, 
                                            reducedModelFormulaStr = reducedModelFormulaStr, 
                                            expressionFamily = cds@expressionFamily, relative_expr = relative_expr, 
                                            disp_func = cds@dispFitInfo[["blind"]]$disp_func, 
                                            verbose = verbose)
    diff_test_res
  }
  diff_test_res <- do.call(rbind.data.frame, diff_test_res)
  diff_test_res$qval <- 1
  diff_test_res$qval[which(diff_test_res$status == "OK")] <- p.adjust(subset(diff_test_res, 
                                                                             status == "OK")[, "pval"], method = "BH")
  diff_test_res <- merge(diff_test_res, fData(cds), by = "row.names")
  row.names(diff_test_res) <- diff_test_res[, 1]
  diff_test_res[, 1] <- NULL
  diff_test_res[row.names(cds), ]
}

diff_test_helper_extract_coef2 <- function (x, fullModelFormulaStr, reducedModelFormulaStr, expressionFamily, 
                                            relative_expr, weights, disp_func = NULL, verbose = FALSE) 
{
  reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, 
                                  sep = "")
  fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, 
                               sep = "")
  x_orig <- x
  disp_guess <- 0
  if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    if (relative_expr == TRUE) {
      x <- x/Size_Factor
    }
    f_expression <- round(x)
    if (is.null(disp_func) == FALSE) {
      disp_guess <- monocle:::calculate_NB_dispersion_hint(disp_func, 
                                                           round(x_orig))
      if (is.null(disp_guess) == FALSE && disp_guess > 
          0 && is.na(disp_guess) == FALSE) {
        if (expressionFamily@vfamily == "negbinomial") 
          expressionFamily <- negbinomial(isize = 1/disp_guess)
        else expressionFamily <- negbinomial.size(size = 1/disp_guess)
      }
    }
  }
  else if (expressionFamily@vfamily %in% c("uninormal")) {
    f_expression <- x
  }
  else if (expressionFamily@vfamily %in% c("binomialff")) {
    f_expression <- x
  }
  else {
    f_expression <- log10(x)
  }
  test_res <- tryCatch({
    if (expressionFamily@vfamily %in% c("binomialff")) {
      if (verbose) {
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), 
                                     epsilon = 0.1, family = expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), 
                                        epsilon = 0.1, family = expressionFamily)
      }
      else {
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), 
                                                      epsilon = 0.1, family = expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), 
                                                         epsilon = 0.1, family = expressionFamily))
      }
    }
    else {
      if (verbose) {
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), 
                                     epsilon = 0.1, family = expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), 
                                        epsilon = 0.1, family = expressionFamily)
      }
      else {
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), 
                                                      epsilon = 0.1, family = expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), 
                                                         epsilon = 0.1, family = expressionFamily))
      }
    }
    #print(coef(full_model_fit))
    cbind.data.frame(monocle::compareModels(list(full_model_fit), list(reduced_model_fit)),data.frame(coef = coef(full_model_fit)["data.diseasestatusFruc"], 
                                                                                                      FC =(coef(full_model_fit)["data.diseasestatusFruc"]+coef(full_model_fit)["(Intercept)"])/coef(full_model_fit)["(Intercept)"] ))
  }, error = function(e) {
    if (verbose) 
      print(e)
    data.frame(status = "FAIL", family = expressionFamily@vfamily, 
               pval = 1, qval = 1)
  })
  test_res
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


#Monocle liver only test
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Combined_allliver.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/liver.rda")
test <- dropEST.combined.filtered@meta.data

#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "Sinusoidal endothelial cells","1" = "Periportal hepatocytes","2" = "Periportal hepatocytes","3" = "Kupffer", "4" = "Pericentral hepatocytes","5" = "Periportal hepatocytes", "6" = "NKT cell",
#                                                    "7" = "Periportal hepatocytes","8" = "Classical dendritic cells","9" = "B cells","10" = "Hepatic stellate cells","11" = "Cholangiocytes","12" = "Dividing cells","13" = "Plasmacytoid dendritic cells"))

dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents = c("Control","Fruc"))
DefaultAssay(dropEST.combined.filtered) <- "RNA"
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "celltype")
alllevels <- as.character(unique(dropEST.combined.filtered$celltype))

alltissues <- NULL
allDEGs <- NULL
DEGnumber <- NULL
allpathways <- NULL
pathwaynumbers <- NULL
rm(finalptable,finalfoldtable,result_final)
Humaninelist <- list()
Pu1list <- list()
Fcgrtlist <- list()
DEGlist3 <- list()
Khklist <- list()
Aldoblist <- list()
Glut5list <- list()
top5pathways <- NULL
top5genes <- NULL
DEGlist <- list()
DEGlist_up <- list()
DEGlist_down <- list()

currentsubset <- subset(dropEST.combined.filtered, idents = "Pericentral hepatocytes")
currentsubset <- SetIdent(currentsubset, value= "data.diseasestatus")
VlnPlot(currentsubset, features = c("Apoa1","Fga","C3","Igf1"))
VlnPlot(currentsubset, features = c("Cyp3a11"))

currentsubset <- subset(dropEST.combined.filtered, idents = "Periportal hepatocytes")
currentsubset <- SetIdent(currentsubset, value= "data.diseasestatus")
VlnPlot(currentsubset, features = c("Apoa1","Fga","C3","Igf1"))

for(i in alllevels){
  currentsubset <- subset(dropEST.combined.filtered, idents = i)
  if(ncol(currentsubset) < 50){
    alltissues <- c(alltissues,i)
    allDEGs <- c(allDEGs,NA)
    DEGnumber <- c(DEGnumber,NA)
    allpathways <- c(allpathways,NA)
    pathwaynumbers <- c(pathwaynumbers, NA)
    top5genes <- c(top5genes, NA)
    top5pathways <- c(top5pathways, NA)
    next
  }
  
  cds <- newCellDataSet(as.matrix(currentsubset@assays$RNA@counts),
                        phenoData = new("AnnotatedDataFrame",
                                        data = currentsubset@meta.data ),
                        expressionFamily = negbinomial.size())
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- detectGenes(cds)
  print(head(fData(cds)))
  expressed_genes <- row.names(subset(fData(cds),
                                      num_cells_expressed >=  ncol(currentsubset) *0.25))
  
  diff_test_res <- differentialGeneTest_extract_coef2(cds[expressed_genes,], fullModelFormulaStr = "~batch+data.diseasestatus+nCount_RNA+nFeature_RNA",reducedModelFormulaStr = "~batch+nCount_RNA+nFeature_RNA",cores = 1)
  #For liver currently
  #diff_test_res <- differentialGeneTest_extract_coef2(cds[expressed_genes,], fullModelFormulaStr = "~data.diseasestatus+nCount_RNA+nFeature_RNA",reducedModelFormulaStr = "~nCount_RNA+nFeature_RNA",cores = 1)
  #diff_test_res <- differentialGeneTest(cds[expressed_genes[1:5],], fullModelFormulaStr = "~data.diseasestatus+nCount_RNA+nFeature_RNA",reducedModelFormulaStr = "~nCount_RNA+nFeature_RNA",cores = 1)
  
  diff_test_res$genename <- rownames(diff_test_res)
  diff_test_res_original <- diff_test_res
  if(diff_test_res$status[1] %in% "FAIL"){
    stop()
    alltissues <- c(alltissues,i)
    allDEGs <- c(allDEGs,NA)
    DEGnumber <- c(DEGnumber,0)
    allpathways <- c(allpathways,NA)
    pathwaynumbers <- c(pathwaynumbers, NA)
    top5genes <- c(top5genes, NA)
    top5pathways <- c(top5pathways, NA)
    next
  }
  
  if(!exists("finalptable")){
    #finalptable <- controltable[,c(1,4)] #using raw pvalue
    finalptable <- diff_test_res[,c(4,8)] #using adj pvalue
    finalfoldtable <- diff_test_res[c(5,8)]
    colnames(finalptable)[1] <- i
    colnames(finalfoldtable)[1] <- i
  }else{
    #no NAs allowed, only consider full samples
    #finalptable <- full_join(finalptable,controltable[,c(1,4)])
    ptable <- diff_test_res[,c(4,8)] #using adj pvalue
    foldtable <- diff_test_res[,c(5,8)]
    colnames(ptable)[1] <- i
    colnames(foldtable)[1] <- i
    finalptable <- full_join(finalptable,ptable)
    finalfoldtable <- full_join(finalfoldtable, foldtable)
  }
  
  Humaninelist[[i]] <- paste0(diff_test_res["mt-Rnr2",]$coef,":",diff_test_res["mt-Rnr2",]$qval)
  Pu1list[[i]] <- paste0(diff_test_res["Spi1",]$coef,":",diff_test_res["Spi1",]$qval)
  Fcgrtlist[[i]] <- paste0(diff_test_res["Fcgrt",]$coef,":",diff_test_res["Fcgrt",]$qval)
  
  Khklist[[i]] <- paste0(diff_test_res["Khk",]$coef,":",diff_test_res["Khk",]$qval)
  Aldoblist[[i]] <- paste0(diff_test_res["Aldob",]$coef,":",diff_test_res["Aldob",]$qval)
  Glut5list[[i]] <- paste0(diff_test_res["Slc2a5",]$coef,":",diff_test_res["Slc2a5",]$qval)
  
  #Requiring at least 10% of total cells expressing for DEG analysis
  diff_test_res <- diff_test_res[diff_test_res$num_cells_expressed >= ncol(currentsubset) *0.25,]
  diff_test_res <- diff_test_res[order(diff_test_res$qval),]
  
  
  
  DEGs <- rownames(diff_test_res[diff_test_res$qval < 0.05,])
  DEGs_up <- rownames(diff_test_res[diff_test_res$qval < 0.05 & diff_test_res$coef > 0,])
  DEGs_down <- rownames(diff_test_res[diff_test_res$qval < 0.05 & diff_test_res$coef < 0,])
  
  if(length(DEGs) > 0){
    DEGlist[[i]] <- DEGs
    if(length(DEGs_up) > 0) {DEGlist_up[[i]] <- DEGs_up
    }#else{DEGlist_up[[i]] <- NA}
    if(length(DEGs_down) > 0) {DEGlist_down[[i]] <- DEGs_down
    }#else{DEGlist_down[[i]] <- NA}
    
    pathways <- enrichrfunction(DEGs)
    if(is.null(pathways) | is.na(pathways)){
      alltissues <- c(alltissues,i)
      allDEGs <- c(allDEGs,paste(DEGs, collapse = ", "))
      DEGnumber <- c(DEGnumber,length(DEGs))
      top5genes <- c(top5genes, paste(DEGs[1:min(5,length(DEGs))], collapse = ", "))
      pathwaynumbers <- c(pathwaynumbers, 0)
      top5pathways <- c(top5pathways, NA)
      allpathways <- c(allpathways,NA)
    }else{
      pathways <- pathways[order(pathways$Adjusted.P.value),]
      alltissues <- c(alltissues,i)
      allDEGs <- c(allDEGs,paste(DEGs, collapse = ", "))
      DEGnumber <- c(DEGnumber,length(DEGs))
      top5genes <- c(top5genes, paste(DEGs[1:min(5,length(DEGs))], collapse = ", "))
      pathwaynumbers <- c(pathwaynumbers, nrow(pathways))
      top5pathways <- c(top5pathways, paste(pathways$Term[1:(min(5,nrow(pathways)))], collapse = ", "))
      allpathways <- c(allpathways,paste0(pathways$Term,collapse = ";"))
      if(nrow(pathways) > 0){
        pathways$celltype <- i
        pathways$FDR <- -log10(pathways$Adjusted.P.value)
        pathways$pathwayname <- pathways$Term
        fold <- NULL
        for(r in 1:nrow(pathways)){
          fold[r] <- median(diff_test_res_original[firstup(unlist(strsplit(pathways$Genes[r],";"))),]$coef,na.rm = T)
        }
        pathways$fold <- fold
        if(!exists("result_final")){
          result_final <- pathways
        }else{
          result_final <- rbind.data.frame(result_final,pathways)
        }
      }
    }
  }else{
    alltissues <- c(alltissues,i)
    allDEGs <- c(allDEGs,NA)
    DEGnumber <- c(DEGnumber,0)
    allpathways <- c(allpathways,NA)
    pathwaynumbers <- c(pathwaynumbers, NA)
    top5genes <- c(top5genes, NA)
    top5pathways <- c(top5pathways, NA)
    
  }
  
}

foldchangemerge_final <- finalfoldtable
pvaluemerge_final <- finalptable
finalframe <- data.frame(celltypes = alltissues,DEGnumber, allDEGs, pathwaynumbers, allpathways)
top5frame <- data.frame(celltypes = alltissues,DEGnumber, top5genes, pathwaynumbers, top5pathways)
save(result_final,Humaninelist,DEGlist,Khklist,Aldoblist,Glut5list,Fcgrtlist,foldchangemerge_final,pvaluemerge_final,DEGlist_up,DEGlist_down,file = "HumaninefoldFruc_monocle.rda")
write.csv(finalframe,file = "monocle_combined_Fruc.csv")
write.csv(top5frame,file = "monocle_combined_top_Fruc.csv")






#Liver HFHS
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
differentialGeneTest_extract_coef2 <- function (cds, fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)", 
                                                reducedModelFormulaStr = "~1", relative_expr = TRUE, cores = 1, 
                                                verbose = FALSE) 
{
  status <- NA
  if (class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  all_vars <- c(all.vars(formula(fullModelFormulaStr)), all.vars(formula(reducedModelFormulaStr)))
  pd <- pData(cds)
  for (i in all_vars) {
    x <- pd[, i]
    if (any((c(Inf, NaN, NA) %in% x))) {
      stop("Error: Inf, NaN, or NA values were located in pData of cds in columns mentioned in model terms")
    }
  }
  if (relative_expr && cds@expressionFamily@vfamily %in% c("negbinomial", 
                                                           "negbinomial.size")) {
    if (is.null(sizeFactors(cds)) || sum(is.na(sizeFactors(cds)))) {
      stop("Error: to call this function with relative_expr==TRUE, you must first call estimateSizeFactors() on the CellDataSet.")
    }
  }
  if (cores > 1) {
    diff_test_res <- mcesApply(cds, 1, diff_test_helper_extract_coef2, 
                               c("BiocGenerics", "VGAM", "Matrix"), cores = cores, 
                               fullModelFormulaStr = fullModelFormulaStr, reducedModelFormulaStr = reducedModelFormulaStr, 
                               expressionFamily = cds@expressionFamily, relative_expr = relative_expr, 
                               disp_func = cds@dispFitInfo[["blind"]]$disp_func, 
                               verbose = verbose)
    diff_test_res
  }
  else {
    diff_test_res <- monocle:::smartEsApply(cds, 1, diff_test_helper_extract_coef2, 
                                            convert_to_dense = TRUE, fullModelFormulaStr = fullModelFormulaStr, 
                                            reducedModelFormulaStr = reducedModelFormulaStr, 
                                            expressionFamily = cds@expressionFamily, relative_expr = relative_expr, 
                                            disp_func = cds@dispFitInfo[["blind"]]$disp_func, 
                                            verbose = verbose)
    diff_test_res
  }
  diff_test_res <- do.call(rbind.data.frame, diff_test_res)
  diff_test_res$qval <- 1
  diff_test_res$qval[which(diff_test_res$status == "OK")] <- p.adjust(subset(diff_test_res, 
                                                                             status == "OK")[, "pval"], method = "BH")
  diff_test_res <- merge(diff_test_res, fData(cds), by = "row.names")
  row.names(diff_test_res) <- diff_test_res[, 1]
  diff_test_res[, 1] <- NULL
  diff_test_res[row.names(cds), ]
}

diff_test_helper_extract_coef2 <- function (x, fullModelFormulaStr, reducedModelFormulaStr, expressionFamily, 
                                            relative_expr, weights, disp_func = NULL, verbose = FALSE) 
{
  reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, 
                                  sep = "")
  fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, 
                               sep = "")
  x_orig <- x
  disp_guess <- 0
  if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    if (relative_expr == TRUE) {
      x <- x/Size_Factor
    }
    f_expression <- round(x)
    if (is.null(disp_func) == FALSE) {
      disp_guess <- monocle:::calculate_NB_dispersion_hint(disp_func, 
                                                           round(x_orig))
      if (is.null(disp_guess) == FALSE && disp_guess > 
          0 && is.na(disp_guess) == FALSE) {
        if (expressionFamily@vfamily == "negbinomial") 
          expressionFamily <- negbinomial(isize = 1/disp_guess)
        else expressionFamily <- negbinomial.size(size = 1/disp_guess)
      }
    }
  }
  else if (expressionFamily@vfamily %in% c("uninormal")) {
    f_expression <- x
  }
  else if (expressionFamily@vfamily %in% c("binomialff")) {
    f_expression <- x
  }
  else {
    f_expression <- log10(x)
  }
  test_res <- tryCatch({
    if (expressionFamily@vfamily %in% c("binomialff")) {
      if (verbose) {
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), 
                                     epsilon = 0.1, family = expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), 
                                        epsilon = 0.1, family = expressionFamily)
      }
      else {
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), 
                                                      epsilon = 0.1, family = expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), 
                                                         epsilon = 0.1, family = expressionFamily))
      }
    }
    else {
      if (verbose) {
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), 
                                     epsilon = 0.1, family = expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), 
                                        epsilon = 0.1, family = expressionFamily)
      }
      else {
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), 
                                                      epsilon = 0.1, family = expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), 
                                                         epsilon = 0.1, family = expressionFamily))
      }
    }
    #print(coef(full_model_fit))
    cbind.data.frame(monocle::compareModels(list(full_model_fit), list(reduced_model_fit)),data.frame(coef = coef(full_model_fit)["data.diseasestatusHFHS"], 
                                                                                                      FC =(coef(full_model_fit)["data.diseasestatusHFHS"]+coef(full_model_fit)["(Intercept)"])/coef(full_model_fit)["(Intercept)"] ))
  }, error = function(e) {
    if (verbose) 
      print(e)
    data.frame(status = "FAIL", family = expressionFamily@vfamily, 
               pval = 1, qval = 1)
  })
  test_res
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


#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Combined_allliver.rda")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "Sinusoidal endothelial cells","1" = "Periportal hepatocytes","2" = "Periportal hepatocytes","3" = "Kupffer", "4" = "Pericentral hepatocytes","5" = "Periportal hepatocytes", "6" = "NKT cell",
#                                                    "7" = "Periportal hepatocytes","8" = "Classical dendritic cells","9" = "B cells","10" = "Hepatic stellate cells","11" = "Cholangiocytes","12" = "Dividing cells","13" = "Plasmacytoid dendritic cells"))
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/liver.rda")
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "data.diseasestatus")
dropEST.combined.filtered <- subset(dropEST.combined.filtered, idents = c("Control","HFHS"))
DefaultAssay(dropEST.combined.filtered) <- "RNA"
dropEST.combined.filtered <- SetIdent(dropEST.combined.filtered, value = "celltype")
alllevels <- as.character(unique(dropEST.combined.filtered$celltype))

alltissues <- NULL
allDEGs <- NULL
DEGnumber <- NULL
allpathways <- NULL
pathwaynumbers <- NULL
rm(finalptable,finalfoldtable,result_final)
Humaninelist <- list()
Pu1list <- list()
Fcgrtlist <- list()
DEGlist3 <- list()
Khklist <- list()
Aldoblist <- list()
Glut5list <- list()
DEGlist <- list()
top5pathways <- NULL
top5genes <- NULL
DEGlist_up <- list()
DEGlist_down <- list()
for(i in alllevels){
  currentsubset <- subset(dropEST.combined.filtered, idents = i)
  if(ncol(currentsubset) < 50){
    alltissues <- c(alltissues,i)
    allDEGs <- c(allDEGs,NA)
    DEGnumber <- c(DEGnumber,NA)
    allpathways <- c(allpathways,NA)
    pathwaynumbers <- c(pathwaynumbers, NA)
    top5genes <- c(top5genes, NA)
    top5pathways <- c(top5pathways, NA)
    next
  }
  
  cds <- newCellDataSet(as.matrix(currentsubset@assays$RNA@counts),
                        phenoData = new("AnnotatedDataFrame",
                                        data = currentsubset@meta.data ),
                        expressionFamily = negbinomial.size())
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- detectGenes(cds)
  print(head(fData(cds)))
  expressed_genes <- row.names(subset(fData(cds),
                                      num_cells_expressed >=  ncol(currentsubset) *0.25))
  
  #diff_test_res <- differentialGeneTest_extract_coef2(cds[expressed_genes,], fullModelFormulaStr = "~batch+data.diseasestatus+nCount_RNA+nFeature_RNA",reducedModelFormulaStr = "~batch+nCount_RNA+nFeature_RNA",cores = 1)
  #For liver currently
  diff_test_res <- differentialGeneTest_extract_coef2(cds[expressed_genes,], fullModelFormulaStr = "~data.diseasestatus+nCount_RNA+nFeature_RNA",reducedModelFormulaStr = "~nCount_RNA+nFeature_RNA",cores = 1)
  
  diff_test_res$genename <- rownames(diff_test_res)
  diff_test_res_original <- diff_test_res
  if(diff_test_res$status[1] %in% "FAIL"){
    alltissues <- c(alltissues,i)
    allDEGs <- c(allDEGs,NA)
    DEGnumber <- c(DEGnumber,0)
    allpathways <- c(allpathways,NA)
    pathwaynumbers <- c(pathwaynumbers, NA)
    top5genes <- c(top5genes, NA)
    top5pathways <- c(top5pathways, NA)
    next
  }
  
  if(!exists("finalptable")){
    #finalptable <- controltable[,c(1,4)] #using raw pvalue
    finalptable <- diff_test_res[,c(4,8)] #using adj pvalue
    finalfoldtable <- diff_test_res[c(5,8)]
    colnames(finalptable)[1] <- i
    colnames(finalfoldtable)[1] <- i
  }else{
    #no NAs allowed, only consider full samples
    #finalptable <- full_join(finalptable,controltable[,c(1,4)])
    ptable <- diff_test_res[,c(4,8)] #using adj pvalue
    foldtable <- diff_test_res[,c(5,8)]
    colnames(ptable)[1] <- i
    colnames(foldtable)[1] <- i
    finalptable <- full_join(finalptable,ptable)
    finalfoldtable <- full_join(finalfoldtable, foldtable)
  }
  
  Humaninelist[[i]] <- paste0(diff_test_res["mt-Rnr2",]$coef,":",diff_test_res["mt-Rnr2",]$qval)
  Pu1list[[i]] <- paste0(diff_test_res["Spi1",]$coef,":",diff_test_res["Spi1",]$qval)
  Fcgrtlist[[i]] <- paste0(diff_test_res["Fcgrt",]$coef,":",diff_test_res["Fcgrt",]$qval)
  
  Khklist[[i]] <- paste0(diff_test_res["Khk",]$coef,":",diff_test_res["Khk",]$qval)
  Aldoblist[[i]] <- paste0(diff_test_res["Aldob",]$coef,":",diff_test_res["Aldob",]$qval)
  Glut5list[[i]] <- paste0(diff_test_res["Slc2a5",]$coef,":",diff_test_res["Slc2a5",]$qval)
  
  #Requiring at least 10% of total cells expressing for DEG analysis
  diff_test_res <- diff_test_res[diff_test_res$num_cells_expressed >= ncol(currentsubset) *0.25,]
  diff_test_res <- diff_test_res[order(diff_test_res$qval),]
  
  
  
  DEGs <- rownames(diff_test_res[diff_test_res$qval < 0.05,])
  DEGs_up <- rownames(diff_test_res[diff_test_res$qval < 0.05 & diff_test_res$coef > 0,])
  DEGs_down <- rownames(diff_test_res[diff_test_res$qval < 0.05 & diff_test_res$coef < 0,])
  
  if(length(DEGs) > 0){
    DEGlist[[i]] <- DEGs
    if(length(DEGs_up) > 0) {DEGlist_up[[i]] <- DEGs_up
    }#else{DEGlist_up[[i]] <- NA}
    if(length(DEGs_down) > 0) {DEGlist_down[[i]] <- DEGs_down
    }#else{DEGlist_down[[i]] <- NA}
    
    pathways <- enrichrfunction(DEGs)
    if(is.na(pathways)){
      alltissues <- c(alltissues,i)
      allDEGs <- c(allDEGs,paste(DEGs, collapse = ", "))
      DEGnumber <- c(DEGnumber,length(DEGs))
      top5genes <- c(top5genes, paste(DEGs[1:min(5,length(DEGs))], collapse = ", "))
      pathwaynumbers <- c(pathwaynumbers, 0)
      top5pathways <- c(top5pathways, NA)
      allpathways <- c(allpathways,NA)
    }else{
      pathways <- pathways[order(pathways$Adjusted.P.value),]
      alltissues <- c(alltissues,i)
      allDEGs <- c(allDEGs,paste(DEGs, collapse = ", "))
      DEGnumber <- c(DEGnumber,length(DEGs))
      top5genes <- c(top5genes, paste(DEGs[1:min(5,length(DEGs))], collapse = ", "))
      pathwaynumbers <- c(pathwaynumbers, nrow(pathways))
      top5pathways <- c(top5pathways, paste(pathways$Term[1:(min(5,nrow(pathways)))], collapse = ", "))
      allpathways <- c(allpathways,paste0(pathways$Term,collapse = ";"))
      if(nrow(pathways) > 0){
        pathways$celltype <- i
        pathways$FDR <- -log10(pathways$Adjusted.P.value)
        pathways$pathwayname <- pathways$Term
        fold <- NULL
        for(r in 1:nrow(pathways)){
          fold[r] <- median(diff_test_res_original[firstup(unlist(strsplit(pathways$Genes[r],";"))),]$coef,na.rm = T)
        }
        pathways$fold <- fold
        if(!exists("result_final")){
          result_final <- pathways
        }else{
          result_final <- rbind.data.frame(result_final,pathways)
        }
      }
    }
  }else{
    alltissues <- c(alltissues,i)
    allDEGs <- c(allDEGs,NA)
    DEGnumber <- c(DEGnumber,0)
    allpathways <- c(allpathways,NA)
    pathwaynumbers <- c(pathwaynumbers, NA)
    top5genes <- c(top5genes, NA)
    top5pathways <- c(top5pathways, NA)
    
  }
}

foldchangemerge_final <- finalfoldtable
pvaluemerge_final <- finalptable
finalframe <- data.frame(celltypes = alltissues, DEGnumber, allDEGs, pathwaynumbers, allpathways)
top5frame <- data.frame(celltypes = alltissues, DEGnumber, top5genes, pathwaynumbers, top5pathways)
save(result_final,Humaninelist,DEGlist,Khklist,Aldoblist,Glut5list,Fcgrtlist,foldchangemerge_final,pvaluemerge_final,DEGlist_up,DEGlist_down,file = "HumaninefoldHFHS_monocle.rda")
write.csv(finalframe,file = "monocle_combined_HFHS.csv")
write.csv(top5frame,file = "monocle_top_HFHS.csv")



#Comparing with One batch wilcoxon results
library(VennDiagram)
#results <- read.csv("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/cross_batch_DEG_FDR_Doug_HFHS.csv")
#results2 <- read.csv("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/monocle_combined_HFHS.csv")
#Pericentral
#test <- list(Wilcoxin_Set1 = unlist(strsplit(results$siggenes_set1[5]," ")),Monocle_combined = unlist(strsplit(results2$allDEGs[1],", ")))
#venn.diagram(test, filename = "Venn_CV_HFHS.png")
#test2 <- list(Wilcoxin_Set1 = unlist(strsplit(results$sigpathways_set1[5],";")),Monocle_combined = unlist(strsplit(results2$allpathways[1],";")))
#venn.diagram(test2, filename = "Venn_CV_path_HFHS.png")

#Perinodal
#test <- list(Wilcoxin_Set1 = unlist(strsplit(results$siggenes_set1[1]," ")),Monocle_combined = unlist(strsplit(results2$allDEGs[2],", ")))
#venn.diagram(test, filename = "Venn_PN_HFHS.png")
#test2 <- list(Wilcoxin_Set1 = unlist(strsplit(results$sigpathways_set1[1],";")),Monocle_combined = unlist(strsplit(results2$allpathways[2],";")))
#venn.diagram(test2, filename = "Venn_PN_Path_HFHS.png")

#SEC
#test <- list(Wilcoxin_Set1 = unlist(strsplit(results$siggenes_set1[2]," ")),Monocle_combined = unlist(strsplit(results2$allDEGs[5],", ")))
#venn.diagram(test, filename = "Venn_SEC_HFHS.png")
#test2 <- list(Wilcoxin_Set1 = unlist(strsplit(results$sigpathways_set1[2],";")),Monocle_combined = unlist(strsplit(results2$allpathways[5],";")))
#venn.diagram(test2, filename = "Venn_SEC_Path_HFHS.png")

#Kupffer
#test <- list(Wilcoxin_Set1 = unlist(strsplit(results$siggenes_set1[3]," ")),Monocle_combined = unlist(strsplit(results2$allDEGs[3],", ")))
#venn.diagram(test, filename = "Venn_Kupffer_HFHS.png")
#test2 <- list(Wilcoxin_Set1 = unlist(strsplit(results$sigpathways_set1[3],";")),Monocle_combined = unlist(strsplit(results2$allpathways[3],";")))
#venn.diagram(test2, filename = "Venn_Kupffer_Path_HFHS.png")



#Fructose
#results <- read.csv("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/cross_batch_DEG_FDR_Doug_Fruc.csv")
#results2 <- read.csv("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/monocle_combined_Fruc.csv")
#Pericentral
#test <- list(Wilcoxin_Set1 = unlist(strsplit(results$siggenes_set1[5]," ")),Monocle_combined = unlist(strsplit(results2$allDEGs[1],", ")))
#venn.diagram(test, filename = "Venn_CV_Fruc.png")
#test2 <- list(Wilcoxin_Set1 = unlist(strsplit(results$sigpathways_set1[5],";")),Monocle_combined = unlist(strsplit(results2$allpathways[1],";")))
#venn.diagram(test2, filename = "Venn_CV_path_Fruc.png")

#Perinodal
#test <- list(Wilcoxin_Set1 = unlist(strsplit(results$siggenes_set1[1]," ")),Monocle_combined = unlist(strsplit(results2$allDEGs[2],", ")))
#venn.diagram(test, filename = "Venn_PN_Fruc.png")
#test2 <- list(Wilcoxin_Set1 = unlist(strsplit(results$sigpathways_set1[1],";")),Monocle_combined = unlist(strsplit(results2$allpathways[2],";")))
#venn.diagram(test2, filename = "Venn_PN_Path_Fruc.png")

#SEC
#test <- list(Wilcoxin_Set1 = unlist(strsplit(results$siggenes_set1[2]," ")),Monocle_combined = unlist(strsplit(results2$allDEGs[5],", ")))
#venn.diagram(test, filename = "Venn_SEC_Fruc.png")
#test2 <- list(Wilcoxin_Set1 = unlist(strsplit(results$sigpathways_set1[2],";")),Monocle_combined = unlist(strsplit(results2$allpathways[5],";")))
#venn.diagram(test2, filename = "Venn_SEC_Path_Fruc.png")

#Kupffer
#test <- list(Wilcoxin_Set1 = unlist(strsplit(results$siggenes_set1[3]," ")),Monocle_combined = unlist(strsplit(results2$allDEGs[3],", ")))
#venn.diagram(test, filename = "Venn_Kupffer_Fruc.png")
#test2 <- list(Wilcoxin_Set1 = unlist(strsplit(results$sigpathways_set1[3],";")),Monocle_combined = unlist(strsplit(results2$allpathways[3],";")))
#venn.diagram(test2, filename = "Venn_Kupffer_Path_Fruc.png")
