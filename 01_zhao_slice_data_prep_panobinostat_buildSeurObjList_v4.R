# data cleaning

library(ggplot2)
library(cowplot)
library(GEOquery)
library(Seurat)
library(Hmisc)
library(dplyr)
library(pbapply)

################################################################################
# ID relevant metadata

GPL18573_series <- getGEO(file = "GSE148842-GPL18573_series_matrix.txt")
GPL18573_pheno <- GPL18573_series@phenoData@data

GPL24676_series <- getGEO(file = "GSE148842-GPL24676_series_matrix.txt")
GPL24676_pheno <- GPL24676_series@phenoData@data

all_pheno <- plyr::rbind.fill(GPL18573_pheno, GPL24676_pheno)

for (i in 1:length(rownames(all_pheno))){
  if (all_pheno$characteristics_ch1.5[i] == ""){
    all_pheno$characteristics_ch1.5[i] <- all_pheno$characteristics_ch1.4[i]
    all_pheno$characteristics_ch1.4[i] <- all_pheno$characteristics_ch1.3[i]
    all_pheno$characteristics_ch1.3[i] <- all_pheno$characteristics_ch1.2[i]
    all_pheno$characteristics_ch1.2[i] <- all_pheno$characteristics_ch1.1[i]
    all_pheno$characteristics_ch1.1[i] <- all_pheno$characteristics_ch1[i]
    
  }
  all_pheno$patientID[i] <- unlist(strsplit(all_pheno$title[i], split = "-", fixed = T))[1]
}

samplesOI <- subset(all_pheno, all_pheno$`treatment:ch1` %in% c("0.2 uM panobinostat"))
samplesOI2 <- subset(all_pheno, all_pheno$patientID %in% samplesOI$patientID)
samplesOI2 <- subset(samplesOI2, samplesOI2$`treatment:ch1` %in% c("vehicle (DMSO)", "0.2 uM panobinostat"))

# # individual ctc.txt files correspond to samples in GEO series matrices
# 
# cts_files <- dir(pattern = "*.cts.txt")
# 
# cts_df <- as.data.frame(cts_files)
# for (i in 1:length(rownames(cts_df))){
#   cts_df$geo_accession[i] <- unlist(strsplit(cts_df$cts_files[i], split = "_", fixed = T))[1]  
# }
# 
# cts_keep <- subset(cts_df, cts_df$geo_accession %in% unique(samplesOI2$geo_accession))$cts_files
# cts_keep
# 
# ################################################################################
# 
# ################################################################################
# #
# ################################################################################
# 
# ##----
# # ensembl to gene aggregation on individual counts mats
# 
# # what cts files are already done?
# finished <- dir(path = "count_geneNameAgg/")
# # need to remove end of filenames in finished to match those in source files. 
# 
# toRun <- cts_keep[which(cts_keep %nin% gsub(pattern = "_count_agg.RDS", x = finished, replacement = ""))]
# 
# 
# 
# # testing
# # for (file in toRun[1]){
# #   counts <- read.table(file = file, header = T, row.names = 1)
# #   duplicated <- counts$gene[which(duplicated(counts$gene) == T)]
# #   counts2 <- counts[counts$gene %in% duplicated,]
# #   nondup <- counts[counts$gene %nin% duplicated,]
# #   rownames(nondup) <- nondup$gene
# #   nondup$gene <- NULL
# #   
# #   # aggregation of duplicated genes
# #   cols <- colnames(counts2)
# #   agg_data <- pbmapply(
# #     split(counts2[, cols, drop = FALSE], counts2$gene),
# #     FUN = function(sub_data){
# #       mean_values <- as.data.frame(colMeans(sub_data[, -1, drop = FALSE]))
# #       colnames(mean_values) <- sub_data$gene[1]
# #       mean_values <- as.data.frame(t(mean_values))
# #       mean_values$gene_name <- sub_data$gene[1]
# #       return(mean_values)
# #     }
# #   )
# #   
# #   # prepare for rbind to final gene aggregate expression matrix
# #   agg_data <- t(agg_data)
# #   agg_df <- as.data.frame(agg_data)
# #   agg_df$gene_name <- NULL
# #   dim(nondup)
# #   dim(agg_df)
# #   final_agg <- rbind(nondup, agg_df)
# # }
# 
# ################################################################################
# #
# ################################################################################
# 
# pbo <- pboptions()
# count.list <- list()
# for (file in toRun){
#   counts <- read.table(file = file, header = T, row.names = 1)
#   if (colnames(counts)[1] == "gene"){
#     duplicated <- counts$gene[which(duplicated(counts$gene) == T)]
#     dupcounts <- counts[counts$gene %in% duplicated,]
#     nondup <- counts[counts$gene %nin% duplicated,]
#     rm(duplicated)
#     rownames(nondup) <- nondup$gene
#     nondup$gene <- NULL
#     
#     # aggregation of duplicated genes
#     cols <- colnames(dupcounts)
#     agg_data <- pbmapply(
#       split(dupcounts[, cols, drop = FALSE], dupcounts$gene),
#       FUN = function(sub_data){
#         mean_values <- as.data.frame(colMeans(sub_data[, -1, drop = FALSE]))
#         colnames(mean_values) <- sub_data$gene[1]
#         mean_values <- as.data.frame(t(mean_values))
#         mean_values$gene_name <- sub_data$gene[1]
#         return(mean_values)
#       }
#     )
#     rm(dupcounts)
#     
#     # prepare for rbind to final gene aggregate expression matrix
#     agg_df <- as.data.frame(t(agg_data))
#     rm(agg_data)
#     agg_df$gene_name <- NULL
#     final_agg <- rbind(nondup, agg_df)
#     saveRDS(final_agg, file = paste0("count_geneNameAgg/", file, "_count_agg.RDS"))
#     
#   } else { # "Gene"
#     duplicated <- counts$gene[which(duplicated(counts$Gene) == T)]
#     dupcounts <- counts[counts$Gene %in% duplicated,]
#     nondup <- counts[counts$Gene %nin% duplicated,]
#     rm(duplicated)
#     rownames(nondup) <- nondup$Gene
#     nondup$Gene <- NULL
#     
#     # aggregation of duplicated genes
#     cols <- colnames(dupcounts)
#     agg_data <- pbmapply(
#       split(dupcounts[, cols, drop = FALSE], dupcounts$Gene),
#       FUN = function(sub_data){
#         mean_values <- as.data.frame(colMeans(sub_data[, -1, drop = FALSE]))
#         colnames(mean_values) <- sub_data$Gene[1]
#         mean_values <- as.data.frame(t(mean_values))
#         mean_values$gene_name <- sub_data$Gene[1]
#         return(mean_values)
#       }
#     )
#     rm(dupcounts)
#     
#     # prepare for rbind to final gene aggregate expression matrix
#     agg_df <- as.data.frame(t(agg_data))
#     rm(agg_data)
#     agg_df$gene_name <- NULL
#     final_agg <- rbind(nondup, agg_df)
#     saveRDS(final_agg, file = paste0("count_geneNameAgg/", file, "_count_agg.RDS"))
#   }
#   # count.list[[file]] <- agg_data
#   rm(counts)
#   rm(agg_data)
#   gc()
# }
# 
# ################################################################################
# #
# ################################################################################
# 
# ##----
# # create seurat object list
# 
# cts_files <- dir(path = "count_geneNameAgg/", pattern = "*.RDS")
# 
# obj.list <- list()
# for (file in cts_files){
#   file <- cts_files[5]
#   counts <- readRDS(file = paste0("count_geneNameAgg/", file))
#   round(counts)
#   counts <- as.matrix(counts)
#   
#   # if (colnames(counts)[1] == "gene"){
#   #   cols <- colnames(counts)[which(colnames(counts) %nin% c("gene"))] 
#   #   gene_agg <- counts %>% 
#   #     group_by(gene) %>% # because of this, need to nest inside of if-else 
#   #     summarise(across(all_of(cols), mean), .groups = 'drop') %>% 
#   #     distinct() %>% 
#   #     as.data.frame()
#   #   rownames(gene_agg) <- gene_agg$gene
#   #   gene_agg$gene <- NULL
#   # } else {
#   #   cols <- colnames(counts)[which(colnames(counts) %nin% c("Gene"))] 
#   #   gene_agg <- counts %>% 
#   #     group_by(Gene) %>% # because of this, need to nest inside of if-else 
#   #     summarise(across(all_of(cols), mean), .groups = 'drop') %>% 
#   #     distinct() %>% 
#   #     as.data.frame()
#   #   rownames(gene_agg) <- gene_agg$Gene
#   #   gene_agg$Gene <- NULL
#   # }
#   # if (colnames(counts)[1] == "gene"){
#   #   ens_mapping <- as.data.frame(counts[,"gene"])
#   #   colnames(ens_mapping) <- c("gene")
#   #   rownames(ens_mapping) <- rownames(counts)
#   #   counts$gene <- NULL
#   # } else {
#   #   ens_mapping <- as.data.frame(counts[,"Gene"])
#   #   colnames(ens_mapping) <- c("Gene")
#   #   rownames(ens_mapping) <- rownames(counts)
#   #   counts$Gene <- NULL
#   # }
#   # rm(counts)
#   # gc()
#   obj <- CreateSeuratObject(counts = round(counts))
#   # rm(gene_agg)
#   obj$cts_file <- file
#   obj <- NormalizeData(object = obj)
#   obj <- FindVariableFeatures(object = obj)
#   obj <- ScaleData(object = obj)
#   obj <- RunPCA(object = obj)
#   obj <- FindNeighbors(object = obj)
#   obj <- FindClusters(object = obj)
#   obj <- RunTSNE(object = obj)
#   # DimPlot(object = obj, reduction = "tsne")
#   obj.list[[file]] <- obj
#   print("List Gen Progress! -------------------")
#   print(paste0("List Gen Progress! --------------   ", length(obj.list) / length(cts_files) * 100, " % Complete" ))
#   rm(obj)
#   gc()
# }

seurObj_filesList <- dir(path = "individualSeuratObjs/", pattern = "*.RDS")
obj.list <- list()
for (file in seurObj_filesList){
  obj.list[[file]] <- readRDS(file = paste0("individualSeuratObjs/", file))
}
saveRDS(obj.list, file = "panobinostat_obj.list.RDS")

plotlist <- list()
for (i in 1:length(obj.list)){
  plotlist[[i]] <- DimPlot(obj.list[[i]]) + NoLegend() + coord_fixed() + ggtitle(names(obj.list)[i])
}
pdf(file = "ind_obj_tsne.pdf", height = 70, width = 21)
cowplot::plot_grid(plotlist = plotlist, labels = "AUTO", ncol = 3)
dev.off()

##----
# need to add relevant metadata to each object prior to integration. 

for (i in 1:length(obj.list)){
  print(names(obj.list)[i])
  obj.list[[i]]$cts_file <- names(obj.list)[i]
  obj.list[[i]]$gsm_id <- unlist(strsplit(obj.list[[i]]$cts_file[1], split = "_", fixed = T))[1]
  obj.list[[i]]$pw_id <- unlist(strsplit(obj.list[[i]]$cts_file[1], split = "_", fixed = T))[2]
}

saveRDS(obj.list, file = "zhao_individual_object_list_preFilt.RDS")

##----
# QC filtering?

# generate pre-filter qc metric plots
# need to use gene name aggregates instead of ensembl-ids to do %^MT- filtering
prefilt_qc <- list()
postfilt_qc <- list()
for (i in 1:length(obj.list)){
  # obj.list[[i]]$percent_mt <- PercentageFeatureSet(obj.list[[i]], pattern = "^MT-")
  prefilt_qc[[i]] <- FeatureScatter(obj.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(names(obj.list)[i])
  cutoff <- quantile(obj.list[[i]]$nFeature_RNA)[4]
  obj.list[[i]] <- subset(obj.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < cutoff)
  postfilt_qc[[i]] <- FeatureScatter(obj.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(paste(names(obj.list)[i], "nFeature cut-off value = ", cutoff))
}

pdf(file = "prefilt_qc.pdf", height = 70, width = 21)
cowplot::plot_grid(plotlist = prefilt_qc, labels = "AUTO", ncol = 3)
dev.off()

pdf(file = "postfilt_qc.pdf", height = 70, width = 21)
cowplot::plot_grid(plotlist = postfilt_qc, labels = "AUTO", ncol = 3)
dev.off()

saveRDS(obj.list, file = "zhao_individual_object_list.RDS")
