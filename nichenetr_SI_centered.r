library(plyr)
library(nichenetr)
library(Seurat)
library(tidyverse)
buildnetwork <- function(frame, matchingresult,currentset,status){
  from <- NULL
  to <- NULL
  matchingresult <- matchingresult[grep(currentset,names(matchingresult))]
  for(i in 1:length(matchingresult)){
    currentname <- names(matchingresult)[i]
    currentname <- gsub(paste0(currentset,"_"),"",currentname)
    currentgenes <- matchingresult[[i]]
    if(length(currentgenes) >0){
      for(r in 1:length(currentgenes)){
        from <- c(from,currentname)
        to <- c(to,currentgenes[r])
      }
    }
  }
  allDEGs <- unique(to)
  
  receiver <- strsplit(as.character(frame[,status]),",")
  for(i in 1:length(receiver)){
    currentname <- as.character(frame$receiver_cells[i])
    currentgenes <- receiver[[i]]
    currentgenes <- currentgenes[currentgenes %in% allDEGs]
    if(length(currentgenes) >0){
      for(r in 1:length(currentgenes)){
        from <- c(from,currentgenes[r])
        to <- c(to,currentname)
      }
    }
  }
  return(data.frame(from,to))
}

testnichenet <- function(seuratObj, currentset,seuratObj_reciever,currentset_receiver,
                         Treatment,receiver,sender_celltypes,
                         ligand_target_matrix,lr_network,weighted_networks,DEGlist,foldchangemerge_final){
  
  seuratObj_reciever <- SetIdent(seuratObj_reciever, value = "batch")
  seuratObj_reciever <- subset(seuratObj_reciever, idents  = currentset_receiver)
  seuratObj_reciever <- SetIdent(seuratObj_reciever, value = "data.diseasestatus")
  seuratObj_reciever <- subset(seuratObj_reciever, idents = c(Treatment,"Control"))
  seuratObj_reciever <- SetIdent(seuratObj_reciever, value = "celltype")
  DefaultAssay(seuratObj_reciever) <- "RNA"
  
  
  seuratObj <- SetIdent(seuratObj, value = "batch")
  seuratObj <- subset(seuratObj, idents  = currentset)
  seuratObj<- SetIdent(seuratObj, value = "data.diseasestatus")
  seuratObj <- subset(seuratObj, idents = c(Treatment,"Control"))
  seuratObj <- SetIdent(seuratObj, value = "celltype")
  DefaultAssay(seuratObj) <- "RNA"
  
  
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
  colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
  ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
  weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
  
  
  
  #Main analysis part
  
  expressed_genes_receiver = get_expressed_genes(receiver, seuratObj_reciever, pct = 0.05,assay_oi = "RNA")
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  
  list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05,assay_oi = "RNA") # lapply to get the expressed genes of every sender cell type separately here
  expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
  
  #Setting condition
  seurat_obj_receiver= subset(seuratObj_reciever, idents = receiver)
  seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = "data.diseasestatus")
  
  condition_oi = Treatment
  condition_reference = "Control" 
  
  #DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")
  
  #geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= 0.25) %>% pull(gene)
  #geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  geneset_oi <- unique(unlist(DEGlist[paste0(currentset_receiver,"_",receiver)]))
  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  
  #4)
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
  p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
    geom_histogram(color="black", fill="darkorange")  + 
    # geom_density(alpha=.1, fill="orange") +
    geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
    labs(x="ligand activity (PCC)", y = "# ligands") +
    theme_classic()
  #p_hist_lig_activity
  
  #5)
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
  
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
  #p_ligand_target_network
  
  
  #Receptors of top-ranked ligands
  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% dplyr:: select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
  #p_ligand_receptor_network
  
  #Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
  #lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  lr_network_strict = lr_network
  ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
  receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
  
  lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
  lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
  
  lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
  lr_network_top_matrix_strict = lr_network_top_df_strict %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
  
  dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
  
  vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
  
  p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
  #p_ligand_receptor_network_strict
  
  
  return(list(NA,p_ligand_receptor_network_strict,
              colnames(vis_ligand_receptor_network_strict),
              rownames(vis_ligand_receptor_network_strict),
              vis_ligand_target ))
}

#Liver version, specific for monocle
testnichenet_monocle <- function(seuratObj, currentset,seuratObj_reciever,currentset_receiver,
                                 Treatment,receiver,sender_celltypes,
                                 ligand_target_matrix,lr_network,weighted_networks,DEGlist,foldchangemerge_final){
  
  seuratObj_reciever <- SetIdent(seuratObj_reciever, value = "batch")
  seuratObj_reciever <- subset(seuratObj_reciever, idents  = currentset_receiver)
  seuratObj_reciever <- SetIdent(seuratObj_reciever, value = "data.diseasestatus")
  seuratObj_reciever <- subset(seuratObj_reciever, idents = c(Treatment,"Control"))
  seuratObj_reciever <- SetIdent(seuratObj_reciever, value = "celltype")
  DefaultAssay(seuratObj_reciever) <- "RNA"
  
  
  seuratObj <- SetIdent(seuratObj, value = "batch")
  seuratObj <- subset(seuratObj, idents  = currentset)
  seuratObj<- SetIdent(seuratObj, value = "data.diseasestatus")
  seuratObj <- subset(seuratObj, idents = c(Treatment,"Control"))
  seuratObj <- SetIdent(seuratObj, value = "celltype")
  DefaultAssay(seuratObj) <- "RNA"
  
  
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
  colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
  ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
  weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
  
  
  
  #Main analysis part
  
  expressed_genes_receiver = get_expressed_genes(receiver, seuratObj_reciever, pct = 0.05,assay_oi = "RNA")
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  
  list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05,assay_oi = "RNA") # lapply to get the expressed genes of every sender cell type separately here
  expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
  
  #Setting condition
  seurat_obj_receiver= subset(seuratObj_reciever, idents = receiver)
  seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = "data.diseasestatus")
  
  condition_oi = Treatment
  condition_reference = "Control" 
  
  #DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")
  
  #geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= 0.25) %>% pull(gene)
  #geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  geneset_oi <- unique(unlist(DEGlist[paste0(receiver)]))
  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  
  #4)
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
  p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
    geom_histogram(color="black", fill="darkorange")  + 
    # geom_density(alpha=.1, fill="orange") +
    geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
    labs(x="ligand activity (PCC)", y = "# ligands") +
    theme_classic()
  #p_hist_lig_activity
  
  #5)
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
  
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
  
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
  #p_ligand_target_network
  
  
  #Receptors of top-ranked ligands
  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% dplyr:: select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
  #p_ligand_receptor_network
  
  #Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
  #lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  lr_network_strict = lr_network
  ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
  receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
  
  lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
  lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
  
  lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
  lr_network_top_matrix_strict = lr_network_top_df_strict %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
  
  dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
  
  vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
  
  p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
  #p_ligand_receptor_network_strict
  
  
  return(list(NA,p_ligand_receptor_network_strict,colnames(vis_ligand_receptor_network_strict),
              rownames(vis_ligand_receptor_network_strict),vis_ligand_target ))
}




setwd("/Users/Tsai_Lab/Downloads/")
ligand_target_matrix = readRDS("./ligand_target_matrix.rds")
lr_network = readRDS("./lr_network.rds")
weighted_networks = readRDS("./weighted_networks.rds")

setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/SI_combined_Set2_final.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SI.rda")
currentset <- "Set2"
sender_celltypes = dropEST.combined.filtered$celltype
seuratObj <- dropEST.combined.filtered



#Liver
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/NPC_liver.rda")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "PN_Hep","1" = "SEC","2" = "Kupffer","3" = "NKT", "4" = "CV_Hep","5" = "Macrophages", "6" = "PN_Hep",
#                                                    "7" = "Hepatic stellate","8" = "B cells","9" = "Dividing cells","10" = "SEC_Dnase1l3","11" = "Plasmacytoid dendritic cells","12" = "Cholangiocytes"))
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Combined_allliver.rda")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "Sinusoidal endothelial cells","1" = "Periportal hepatocytes","2" = "Periportal hepatocytes","3" = "Kupffer", "4" = "Pericentral hepatocytes","5" = "Periportal hepatocytes", "6" = "NKT cells",
#                                                    "7" = "Periportal hepatocytes","8" = "Classical dendritic cells","9" = "B cells","10" = "Hepatic stellate cells","11" = "Cholangiocytes","12" = "Dividing cells","13" = "Plasmacytoid dendritic cells"))
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/liver.rda")

currentset_receiver <- c("Set1","Set2")
#sender_celltypes =  as.character(unique(dropEST.combined.filtered$celltype))
receiver_cells <-  as.character(unique(dropEST.combined.filtered$celltype))

shared <- NULL
Fructose_specific <- NULL
HFHS_specific <- NULL
Fructose <- NULL
HFHS <- NULL
Fructose_count <- NULL
HFHS_count <- NULL
ligand_target_frame_Fruc <- data.frame(ligand = character(0),target = character(0))
ligand_target_frame_HFHS <- data.frame(ligand = character(0),target = character(0), celltype= character(0))
result_HFHS_list <- list()
result_Fruc_list <- list()
for(i in 1:length(receiver_cells)){
  Treatment <- "Fruc"
  #Old version, combining from two
  #load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/Humaninefold.rda")
  #DEGlist2 <- DEGlist
  #foldchangemerge_final2 <- foldchangemerge_final
  #load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_NPC/Humaninefold (a23451352@yahoo.com.tw).rda")
  #DEGlist <- c(DEGlist,DEGlist2)
  #foldchangemerge_final <- full_join(foldchangemerge_final, foldchangemerge_final2)
  #load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/Humaninefold.rda")
  
  #geneset_oi <- unique(unlist(DEGlist[paste0(currentset,"_",receiver_cells[[i]])]))
  load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/HumaninefoldFruc_monocle.rda")
  geneset_oi <- unique(unlist(DEGlist[receiver_cells[[i]]]))
  
  result_Fruc <- list(NA,NA,character(0),NA)
  if(length(geneset_oi) > 20){
    try(result_Fruc <-  testnichenet_monocle(seuratObj, currentset,dropEST.combined.filtered,currentset_receiver,
                                             Treatment,receiver_cells[i],sender_celltypes,
                                             ligand_target_matrix,lr_network,weighted_networks,DEGlist,foldchangemerge_final))
  }else(result_Fruc <- list(NA,NA,character(0),NA,NA))
  
  if(!is.na(result_Fruc[[5]])){
    currentframe <-  as.data.frame.table(`row.names<-`(as.matrix(result_Fruc[[5]]),
                                                       rownames(result_Fruc[[5]])))
    
    currentframe <- currentframe[!currentframe$Freq == 0,]
    currentframe <- currentframe[currentframe$Var1 %in% result_Fruc[[3]],]
    currentframe[,1] <- as.character(currentframe[,1] )
    currentframe[,2] <- as.character(currentframe[,2] )
    currentframe$celltype <- receiver_cells[[i]]
    ligand_target_frame_Fruc <- rbind.data.frame(ligand_target_frame_Fruc , currentframe[,c(1,2,4)])
    
    result_Fruc_list[[receiver_cells[i]]] <- result_Fruc
  }
  
  
  Treatment <- "HFHS"
  #load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/HumaninefoldHFHS.rda")
  #DEGlist2 <- DEGlist
  #foldchangemerge_final2 <- foldchangemerge_final
  #load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_NPC/HumaninefoldHFHS.rda")
  #DEGlist <- c(DEGlist,DEGlist2)
  #foldchangemerge_final <- full_join(foldchangemerge_final, foldchangemerge_final2)
  #load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/HumaninefoldHFHS.rda")
  
  #geneset_oi <- unique(unlist(DEGlist[paste0(currentset,"_",receiver_cells[[i]])]))
  load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HEP_set1/NPC_LIVER/HumaninefoldHFHS_monocle.rda")
  geneset_oi <- unique(unlist(DEGlist[receiver_cells[[i]]]))
  
  result_HFHS <- list(NA,NA,character(0),NA)
  if(length(geneset_oi) > 20){
    try(result_HFHS <-  testnichenet_monocle(seuratObj, currentset,dropEST.combined.filtered,currentset_receiver,
                                             Treatment,receiver_cells[i],sender_celltypes,
                                             ligand_target_matrix,lr_network,weighted_networks,DEGlist,foldchangemerge_final))
  }else{result_HFHS <- list(NA,NA,character(0),NA,NA)}
  
  
  if(!is.na(result_HFHS[[5]])){
    currentframe <-  as.data.frame.table(`row.names<-`(as.matrix(result_HFHS[[5]]),
                                                       rownames(result_HFHS[[5]])))
    
    currentframe <- currentframe[!currentframe$Freq == 0,]
    currentframe <- currentframe[currentframe$Var1 %in% result_HFHS[[3]],]
    currentframe[,1] <- as.character(currentframe[,1] )
    currentframe[,2] <- as.character(currentframe[,2] )
    currentframe$celltype <- receiver_cells[[i]]
    ligand_target_frame_HFHS <- rbind.data.frame(ligand_target_frame_HFHS , currentframe[,c(1,2,4)])
    result_HFHS_list[[receiver_cells[i]]] <- result_HFHS
    }
  
  shared[i] <- paste(intersect(result_Fruc[[3]],result_HFHS[[3]]),collapse = ",")
  Fructose_specific[i] <- paste(setdiff(result_Fruc[[3]],result_HFHS[[3]]),collapse = ",")
  Fructose[i] <- paste(result_Fruc[[3]],collapse = ",")
  HFHS_specific[i] <- paste(setdiff(result_HFHS[[3]],result_Fruc[[3]]),collapse = ",")
  HFHS[i] <- paste(result_HFHS[[3]],collapse = ",")
  Fructose_count <- c(Fructose_count,result_Fruc[[3]] )
  HFHS_count <- c(HFHS_count,result_HFHS[[3]])
}
save(result_HFHS_list, result_Fruc_list, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Su_Liver_result_Orig_SI.rda")

frame_Liver <- data.frame(receiver_cells,shared,Fructose,Fructose_specific,HFHS,HFHS_specific)
sort(table(Fructose_count))
sort(table(HFHS_count))

load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/Humaninefold.rda")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Humaninefold_subset (a23451352@yahoo.com.tw).rda")
Fructose_match_Liver<- lapply(DEGlist, function(y) y[y %in% unique(Fructose_count)])

load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/HumaninefoldHFHS.rda")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/HumaninefoldHFHS_subset (a23451352@yahoo.com.tw).rda")

HFHS_match_Liver <- lapply(DEGlist, function(y) y[y %in% unique(HFHS_count)])


network_Fruc_Liver <- buildnetwork(frame_Liver,Fructose_match_Liver,currentset,"Fructose")
network_Fruc_Liver$treatment <- "Fruc"
network_HFHS_Liver <- buildnetwork(frame_Liver,HFHS_match_Liver,currentset,"HFHS")
network_HFHS_Liver$treatment <- "HFHS"
write.csv(network_Fruc_Liver,file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/network_Liver_test_SI_Fruc.csv",quote = F, row.names = F)
write.csv(network_HFHS_Liver,file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/network_Liver_test_SI_HFHS.csv",quote = F, row.names = F)


network_Liver <- rbind.data.frame(network_Fruc_Liver,network_HFHS_Liver)
network_Liver$treatment[which(duplicated(network_Liver[,c(1,2)]))] <- "shared"
network_Liver <- network_Liver[-which(duplicated(network_Liver[,c(1,2)],fromLast = T)),]
#write.csv(network_Liver,file = "network_Liver_test.csv",quote = F, row.names = F)
write.csv(network_Liver,file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/network_Liver_test_SI.csv",quote = F, row.names = F)
save(ligand_target_frame_Fruc,ligand_target_frame_HFHS, file = "/Users/Tsai_Lab/Desktop/Box Sync/Mergeomics-master/Network_SI_LR_Liver.rda")



#SVF
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined")
#load("HEP_Set1_DropEST.Object1.rda")
#dropEST.combined.filtered[["celltype"]] = revalue(dropEST.combined.filtered@meta.data$integrated_snn_res.0.5,
#                                                  c("0" = "ADSC_Pi16","1" = "ADSC_Hsd11b1","3" = "ADSC_Hsd11b1","5" = "ADSC_Agt", "8" = "Preadipocytes","7" = "Endothelial cells", "10" = "NK",
#                                                    "11" = "Germ cells","9" = "Macrophage_Mmp12","13" = "B cells","12" = "Dividing cells","6" = "Macrophage_Il1b","4" = "Macrophage_Pf4",
#                                                    "2" = "Macrophage_Pf4" ,"14" = "Mast cells","15"="Unknown"))
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SVF_combined/SVF_Final_Set1only.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SVF.rda")
#DimPlot(dropEST.combined.filtered, group.by = "celltype",reduction="tsne",label = T)
currentset_receiver <- "Set1"

#sender_celltypes =  as.character(unique(dropEST.combined.filtered$celltype))
receiver_cells <-  as.character(unique(dropEST.combined.filtered$celltype))

shared <- NULL
Fructose_specific <- NULL
HFHS_specific <- NULL
Fructose <- NULL
HFHS <- NULL
Fructose_count <- NULL
HFHS_count <- NULL
ligand_target_frame_Fruc <- data.frame(ligand = character(0),target = character(0))
ligand_target_frame_HFHS <- data.frame(ligand = character(0),target = character(0), celltype= character(0))
result_HFHS_list <- list()
result_Fruc_list <- list()
for(i in 1:length(receiver_cells)){
  Treatment <- "Fruc"
  info <- ifelse(Treatment %in% "Fruc","Humaninefold.rda","HumaninefoldHFHS.rda")
  load(info)
  geneset_oi <- unique(unlist(DEGlist[paste0(currentset_receiver,"_",receiver_cells[[i]])]))
  if(length(geneset_oi) > 20){
    try(result_Fruc <-  testnichenet(seuratObj, currentset,dropEST.combined.filtered,currentset_receiver,
                                     Treatment,receiver_cells[i],sender_celltypes,
                                     ligand_target_matrix,lr_network,weighted_networks,DEGlist,foldchangemerge_final))
  }else(result_Fruc <- list(NA,NA,character(0),NA,NA))
  
  if(!is.na(result_Fruc[[5]])){
    currentframe <-  as.data.frame.table(`row.names<-`(as.matrix(result_Fruc[[5]]),
                                                       rownames(result_Fruc[[5]])))
    
    currentframe <- currentframe[!currentframe$Freq == 0,]
    currentframe <- currentframe[currentframe$Var1 %in% result_Fruc[[3]],]
    currentframe[,1] <- as.character(currentframe[,1] )
    currentframe[,2] <- as.character(currentframe[,2] )
    currentframe$celltype <- receiver_cells[[i]]
    ligand_target_frame_Fruc <- rbind.data.frame(ligand_target_frame_Fruc , currentframe[,c(1,2,4)])
    result_Fruc_list[[receiver_cells[i]]] <- result_Fruc
    }
  
  
  
  Treatment <- "HFHS"
  info <- ifelse(Treatment %in% "Fruc","Humaninefold.rda","HumaninefoldHFHS.rda")
  load(info)
  geneset_oi <- unique(unlist(DEGlist[paste0(currentset_receiver,"_",receiver_cells[[i]])]))
  if(length(geneset_oi) > 20){
    try(result_HFHS <-  testnichenet(seuratObj, currentset,dropEST.combined.filtered,currentset_receiver,
                                     Treatment,receiver_cells[i],sender_celltypes,
                                     ligand_target_matrix,lr_network,weighted_networks,DEGlist,foldchangemerge_final))
  }else{result_HFHS <- list(NA,NA,character(0),NA,NA)}
  
  
  if(!is.na(result_HFHS[[5]])){
    currentframe <-  as.data.frame.table(`row.names<-`(as.matrix(result_HFHS[[5]]),
                                                       rownames(result_HFHS[[5]])))
    
    currentframe <- currentframe[!currentframe$Freq == 0,]
    currentframe <- currentframe[currentframe$Var1 %in% result_HFHS[[3]],]
    currentframe[,1] <- as.character(currentframe[,1] )
    currentframe[,2] <- as.character(currentframe[,2] )
    currentframe$celltype <- receiver_cells[[i]]
    ligand_target_frame_HFHS <- rbind.data.frame(ligand_target_frame_HFHS , currentframe[,c(1,2,4)])
    result_HFHS_list[[receiver_cells[i]]] <- result_HFHS
    }
  
  shared[i] <- paste(intersect(result_Fruc[[3]],result_HFHS[[3]]),collapse = ",")
  Fructose_specific[i] <- paste(setdiff(result_Fruc[[3]],result_HFHS[[3]]),collapse = ",")
  Fructose[i] <- paste(result_Fruc[[3]],collapse = ",")
  HFHS_specific[i] <- paste(setdiff(result_HFHS[[3]],result_Fruc[[3]]),collapse = ",")
  HFHS[i] <- paste(result_HFHS[[3]],collapse = ",")
  Fructose_count <- c(Fructose_count,result_Fruc[[3]] )
  HFHS_count <- c(HFHS_count,result_HFHS[[3]])
}
save(result_HFHS_list, result_Fruc_list, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Su_SVF_result_Orig_SI.rda")

frame_SVF <- data.frame(receiver_cells,shared,Fructose,Fructose_specific,HFHS,HFHS_specific)
sort(table(Fructose_count))
sort(table(HFHS_count))

load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/Humaninefold.rda")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Humaninefold_subset (a23451352@yahoo.com.tw).rda")
Fructose_match_SVF <- lapply(DEGlist, function(y) y[y %in% unique(Fructose_count)])

load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/HumaninefoldHFHS.rda")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/HumaninefoldHFHS_subset (a23451352@yahoo.com.tw).rda")
HFHS_match_SVF <- lapply(DEGlist, function(y) y[y %in% unique(HFHS_count)])


network_Fruc_SVF <- buildnetwork(frame_SVF,Fructose_match_SVF,currentset,"Fructose")
network_Fruc_SVF$treatment <- "Fruc"
network_HFHS_SVF <- buildnetwork(frame_SVF,HFHS_match_SVF,currentset,"HFHS")
network_HFHS_SVF$treatment <- "HFHS"
write.csv(network_Fruc_SVF ,file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/network_SVF_test_SI_Fruc.csv",quote = F, row.names = F)
write.csv(network_HFHS_SVF,file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/network_SVF_test_SI_HFHS.csv",quote = F, row.names = F)


network_SVF <- rbind.data.frame(network_Fruc_SVF,network_HFHS_SVF)
network_SVF$treatment[which(duplicated(network_SVF[,c(1,2)]))] <- "shared"
#network_SVF <- network_SVF[-which(duplicated(network_SVF[,c(1,2)],fromLast = T)),]
#write.csv(network_SVF,file = "network_SVF_test.csv",quote = F, row.names = F)
write.csv(network_SVF,file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/network_SVF_test_SI.csv",quote = F, row.names = F)
save(ligand_target_frame_Fruc,ligand_target_frame_HFHS, file = "/Users/Tsai_Lab/Desktop/Box Sync/Mergeomics-master/Network_SI_LR_SVF.rda")


#hyp
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Hyp_combined_Set1_final.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/Hyp.rda")
#DimPlot(dropEST.combined.filtered, group.by = "celltype",reduction="tsne",label = T)
currentset_receiver <- "Set1"

#sender_celltypes =  as.character(unique(dropEST.combined.filtered$celltype))
receiver_cells <-  as.character(unique(dropEST.combined.filtered$celltype))

shared <- NULL
Fructose_specific <- NULL
HFHS_specific <- NULL
Fructose <- NULL
HFHS <- NULL
Fructose_count <- NULL
HFHS_count <- NULL
ligand_target_frame_Fruc <- data.frame(ligand = character(0),target = character(0))
ligand_target_frame_HFHS <- data.frame(ligand = character(0),target = character(0), celltype= character(0))
result_HFHS_list <- list()
result_Fruc_list <- list()
for(i in 1:length(receiver_cells)){
  Treatment <- "Fruc"
  info <- ifelse(Treatment %in% "Fruc","Humaninefold.rda","HumaninefoldHFHS.rda")
  load(info)
  geneset_oi <- unique(unlist(DEGlist[paste0(currentset_receiver,"_",receiver_cells[[i]])]))
  if(length(geneset_oi) > 20){
    try(result_Fruc <-  testnichenet(seuratObj, currentset,dropEST.combined.filtered,currentset_receiver,
                                     Treatment,receiver_cells[i],sender_celltypes,
                                     ligand_target_matrix,lr_network,weighted_networks,DEGlist,foldchangemerge_final))
  }else(result_Fruc <- list(NA,NA,character(0),NA,NA))
  
  if(!is.na(result_Fruc[[5]])){
    currentframe <-  as.data.frame.table(`row.names<-`(as.matrix(result_Fruc[[5]]),
                                                       rownames(result_Fruc[[5]])))
    
    currentframe <- currentframe[!currentframe$Freq == 0,]
    currentframe <- currentframe[currentframe$Var1 %in% result_Fruc[[3]],]
    currentframe[,1] <- as.character(currentframe[,1] )
    currentframe[,2] <- as.character(currentframe[,2] )
    currentframe$celltype <- receiver_cells[[i]]
    ligand_target_frame_Fruc <- rbind.data.frame(ligand_target_frame_Fruc , currentframe[,c(1,2,4)])
    result_Fruc_list[[receiver_cells[i]]] <- result_Fruc
    }

  Treatment <- "HFHS"
  info <- ifelse(Treatment %in% "Fruc","Humaninefold.rda","HumaninefoldHFHS.rda")
  load(info)
  geneset_oi <- unique(unlist(DEGlist[paste0(currentset_receiver,"_",receiver_cells[[i]])]))
  if(length(geneset_oi) > 20){
    try(result_HFHS <-  testnichenet(seuratObj, currentset,dropEST.combined.filtered,currentset_receiver,
                                     Treatment,receiver_cells[i],sender_celltypes,
                                     ligand_target_matrix,lr_network,weighted_networks,DEGlist,foldchangemerge_final))
  }else{result_HFHS <- list(NA,NA,character(0),NA,NA)}
  
  
  if(!is.na(result_HFHS[[5]])){
    currentframe <-  as.data.frame.table(`row.names<-`(as.matrix(result_HFHS[[5]]),
                                                       rownames(result_HFHS[[5]])))
    
    currentframe <- currentframe[!currentframe$Freq == 0,]
    currentframe <- currentframe[currentframe$Var1 %in% result_HFHS[[3]],]
    currentframe[,1] <- as.character(currentframe[,1] )
    currentframe[,2] <- as.character(currentframe[,2] )
    currentframe$celltype <- receiver_cells[[i]]
    ligand_target_frame_HFHS <- rbind.data.frame(ligand_target_frame_HFHS , currentframe[,c(1,2,4)])
    result_HFHS_list[[receiver_cells[i]]] <- result_HFHS
    }

  shared[i] <- paste(intersect(result_Fruc[[3]],result_HFHS[[3]]),collapse = ",")
  Fructose_specific[i] <- paste(setdiff(result_Fruc[[3]],result_HFHS[[3]]),collapse = ",")
  Fructose[i] <- paste(result_Fruc[[3]],collapse = ",")
  HFHS_specific[i] <- paste(setdiff(result_HFHS[[3]],result_Fruc[[3]]),collapse = ",")
  HFHS[i] <- paste(result_HFHS[[3]],collapse = ",")
  Fructose_count <- c(Fructose_count,result_Fruc[[3]] )
  HFHS_count <- c(HFHS_count,result_HFHS[[3]])
}
save(result_HFHS_list, result_Fruc_list, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Su_Hyp_result_Orig_SI.rda")

frame_HYP <- data.frame(receiver_cells,shared,Fructose,Fructose_specific,HFHS,HFHS_specific)
sort(table(Fructose_count))
sort(table(HFHS_count))

load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/Humaninefold.rda")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Humaninefold_subset (a23451352@yahoo.com.tw).rda")
Fructose_match_HYP <- lapply(DEGlist, function(y) y[y %in% unique(Fructose_count)])

load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/HumaninefoldHFHS.rda")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/HumaninefoldHFHS_subset (a23451352@yahoo.com.tw).rda")
HFHS_match_HYP <- lapply(DEGlist, function(y) y[y %in% unique(HFHS_count)])


network_Fruc_HYP <- buildnetwork(frame_HYP,Fructose_match_HYP,currentset,"Fructose")
network_Fruc_HYP$treatment <- "Fruc"
network_HFHS_HYP <- buildnetwork(frame_HYP,HFHS_match_HYP,currentset,"HFHS")
network_HFHS_HYP$treatment <- "HFHS"
write.csv(network_Fruc_HYP ,file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/network_HYP_test_SI_Fruc.csv",quote = F, row.names = F)
write.csv(network_HFHS_HYP,file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/network_HYP_test_SI_HFHS.csv",quote = F, row.names = F)


network_HYP <- rbind.data.frame(network_Fruc_HYP,network_HFHS_HYP)
network_HYP$treatment[which(duplicated(network_HYP[,c(1,2)]))] <- "shared"
#network_HYP <- network_HYP[-which(duplicated(network_HYP[,c(1,2)],fromLast = T)),]
#write.csv(network_HYP,file = "network_HYP_test.csv",quote = F, row.names = F)
write.csv(network_HYP,file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/network_HYP_test_SI.csv",quote = F, row.names = F)

save(ligand_target_frame_Fruc,ligand_target_frame_HFHS, file = "/Users/Tsai_Lab/Desktop/Box Sync/Mergeomics-master/Network_SI_LR_HYP.rda")



#SI
setwd("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/SI_combined_Set2_final.rda")
load("/Users/Tsai_Lab/Desktop/Box Sync/shared_slides/filtered_object/SI.rda")
#DimPlot(dropEST.combined.filtered, group.by = "celltype",reduction="tsne",label = T)
currentset_receiver <- "Set2"
#sender_celltypes =  as.character(unique(dropEST.combined.filtered$celltype))
receiver_cells <-  as.character(unique(dropEST.combined.filtered$celltype))

shared <- NULL
Fructose_specific <- NULL
HFHS_specific <- NULL
Fructose <- NULL
HFHS <- NULL
Fructose_count <- NULL
HFHS_count <- NULL
ligand_target_frame_Fruc <- data.frame(ligand = character(0),target = character(0))
ligand_target_frame_HFHS <- data.frame(ligand = character(0),target = character(0), celltype= character(0))
result_HFHS_list <- list()
result_Fruc_list <- list()
for(i in 1:length(receiver_cells)){
  Treatment <- "Fruc"
  info <- ifelse(Treatment %in% "Fruc","Humaninefold.rda","HumaninefoldHFHS.rda")
  load(info)
  geneset_oi <- unique(unlist(DEGlist[paste0(currentset_receiver,"_",receiver_cells[[i]])]))
  if(length(geneset_oi) > 20){
    try(result_Fruc <-  testnichenet(seuratObj, currentset,dropEST.combined.filtered,currentset_receiver,
                                     Treatment,receiver_cells[i],sender_celltypes,
                                     ligand_target_matrix,lr_network,weighted_networks,DEGlist,foldchangemerge_final))
  }else(result_Fruc <- list(NA,NA,character(0),NA,NA))
  
  if(!is.na(result_Fruc[[5]])){
    currentframe <-  as.data.frame.table(`row.names<-`(as.matrix(result_Fruc[[5]]),
                                                       rownames(result_Fruc[[5]])))
    
    currentframe <- currentframe[!currentframe$Freq == 0,]
    currentframe <- currentframe[currentframe$Var1 %in% result_Fruc[[3]],]
    currentframe[,1] <- as.character(currentframe[,1] )
    currentframe[,2] <- as.character(currentframe[,2] )
    currentframe$celltype <- receiver_cells[[i]]
    ligand_target_frame_Fruc <- rbind.data.frame(ligand_target_frame_Fruc , currentframe[,c(1,2,4)])
    result_Fruc_list[[receiver_cells[i]]] <- result_Fruc
    }
  
  Treatment <- "HFHS"
  info <- ifelse(Treatment %in% "Fruc","Humaninefold.rda","HumaninefoldHFHS.rda")
  load(info)
  geneset_oi <- unique(unlist(DEGlist[paste0(currentset_receiver,"_",receiver_cells[[i]])]))
  if(length(geneset_oi) > 20){
    try(result_HFHS <-  testnichenet(seuratObj, currentset,dropEST.combined.filtered,currentset_receiver,
                                     Treatment,receiver_cells[i],sender_celltypes,
                                     ligand_target_matrix,lr_network,weighted_networks,DEGlist,foldchangemerge_final))
  }else{result_HFHS <- list(NA,NA,character(0),NA,NA)}
  
  
  if(!is.na(result_HFHS[[5]])){
    currentframe <-  as.data.frame.table(`row.names<-`(as.matrix(result_HFHS[[5]]),
                                                       rownames(result_HFHS[[5]])))
    
    currentframe <- currentframe[!currentframe$Freq == 0,]
    currentframe <- currentframe[currentframe$Var1 %in% result_HFHS[[3]],]
    currentframe[,1] <- as.character(currentframe[,1] )
    currentframe[,2] <- as.character(currentframe[,2] )
    currentframe$celltype <- receiver_cells[[i]]
    ligand_target_frame_HFHS <- rbind.data.frame(ligand_target_frame_HFHS , currentframe[,c(1,2,4)])
    result_HFHS_list[[receiver_cells[i]]] <- result_HFHS
    }
  
  shared[i] <- paste(intersect(result_Fruc[[3]],result_HFHS[[3]]),collapse = ",")
  Fructose_specific[i] <- paste(setdiff(result_Fruc[[3]],result_HFHS[[3]]),collapse = ",")
  Fructose[i] <- paste(result_Fruc[[3]],collapse = ",")
  HFHS_specific[i] <- paste(setdiff(result_HFHS[[3]],result_Fruc[[3]]),collapse = ",")
  HFHS[i] <- paste(result_HFHS[[3]],collapse = ",")
  Fructose_count <- c(Fructose_count,result_Fruc[[3]] )
  HFHS_count <- c(HFHS_count,result_HFHS[[3]])
}
save(result_HFHS_list, result_Fruc_list, file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/Su_SI_result_Orig_SI.rda")
frame_SI <- data.frame(receiver_cells,shared,Fructose,Fructose_specific,HFHS,HFHS_specific)
sort(table(Fructose_count))
sort(table(HFHS_count))

load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/Humaninefold.rda")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/Humaninefold_subset (a23451352@yahoo.com.tw).rda")
Fructose_match_SI <- lapply(DEGlist, function(y) y[y %in% unique(Fructose_count)])


load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/SI_combined/HumaninefoldHFHS.rda")
#load("/Users/Tsai_Lab/Desktop/Box Sync/Fructose/HYP_combined/HumaninefoldHFHS_subset (a23451352@yahoo.com.tw).rda")
HFHS_match_SI <- lapply(DEGlist, function(y) y[y %in% unique(HFHS_count)])


network_Fruc_SI <- buildnetwork(frame_SI,Fructose_match_SI,currentset,"Fructose")
network_Fruc_SI$treatment <- "Fruc"
network_HFHS_SI <- buildnetwork(frame_SI,HFHS_match_SI,currentset,"HFHS")
network_HFHS_SI$treatment <- "HFHS"
write.csv(network_Fruc_SI,file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/network_SI_test_SI_Fruc.csv",quote = F, row.names = F)
write.csv(network_HFHS_SI,file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/network_SI_test_SI_HFHS.csv",quote = F, row.names = F)

network_SI <- rbind.data.frame(network_Fruc_SI,network_HFHS_SI)
network_SI$treatment[which(duplicated(network_SI[,c(1,2)]))] <- "shared"
network_SI <- network_SI[-which(duplicated(network_SI[,c(1,2)],fromLast = T)),]
write.csv(network_SI,file = "/Users/Tsai_Lab/Desktop/Box Sync/Fructose/network_SI_test_SI.csv",quote = F, row.names = F)
save(ligand_target_frame_Fruc,ligand_target_frame_HFHS, file = "/Users/Tsai_Lab/Desktop/Box Sync/Mergeomics-master/Network_SI_LR_SI.rda")


network_Fruc_Liver$tissue <- "Liver"
network_HFHS_Liver$tissue <- "Liver"
network_Fruc_HYP$tissue <- "HYP"
network_HFHS_HYP$tissue <- "HYP"
network_Fruc_SVF$tissue <- "SVF"
network_HFHS_SVF$tissue <- "SVF"
network_Fruc_SI$tissue <- "SI"
network_HFHS_SI$tissue <- "SI"
combined_network_fruc <- rbind.data.frame(network_Fruc_Liver,network_Fruc_HYP,network_Fruc_SVF,network_Fruc_SI)
combined_network_HFHS <- rbind.data.frame(network_HFHS_Liver,network_HFHS_HYP,network_HFHS_SVF,network_HFHS_SI)
write.csv(combined_network_fruc,file = "/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/combined_Fruc_SI_network_final_ver2.csv",quote = F)
write.csv(combined_network_HFHS,file = "/Users/Tsai_Lab/Desktop/Box Sync/Yang_Lab_doc_share/projects/Fructose_single_cell/combined_HFHS_SI_network_final_ver2.csv", quote = F)

