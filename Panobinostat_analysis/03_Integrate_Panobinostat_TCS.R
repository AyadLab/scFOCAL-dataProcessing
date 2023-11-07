
library(ggplot2)
library(cowplot)
library(GEOquery)
library(Seurat)
library(Hmisc)
library(dplyr)
library(scrabble)

################################################################################
# ID relevant metadata

GPL18573_series <- getGEO(file = "GSE148842-GPL18573_series_matrix.txt")
GPL18573_pheno <- GPL18573_series@phenoData@data

GPL24676_series <- getGEO(file = "GSE148842-GPL24676_series_matrix.txt")
GPL24676_pheno <- GPL24676_series@phenoData@data

################################################################################

obj <- readRDS(file = "zhao_integrated_varFeatures.RDS")

meta_df <- obj@meta.data
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
# all_pheno$characteristics_ch1.5
meta_df$cellID <- rownames(meta_df)

meta_merge <- merge(meta_df, all_pheno, 
                    by.x = "gsm_id", by.y = "geo_accession", 
                    all.x = T)

rownames(meta_merge) <- meta_merge$cellID
meta_add <- meta_merge[,c("treatment:ch1", 
                          "patientID", 
                          "characteristics_ch1.1", 
                          "characteristics_ch1.2", 
                          "characteristics_ch1.3", 
                          "characteristics_ch1.4", 
                          "characteristics_ch1.5")]
meta_add$characteristics_ch1.5
colnames(meta_add) <- c("treatment:ch1", "patientID", "age", "gender", "location", "diagnosis", "treatment")
meta_add

obj <- AddMetaData(obj, metadata = meta_add)

################################################################################
# Panobinostat TCS integration
################################################################################

tcs_df <- read.csv(file = "TCS_out/phase2_compound_TCS_v0.1_20230914_112256.csv", row.names = 1)

# regress cell cycle and scale
DefaultAssay(obj) <- "RNA"
obj <- ScaleData(obj, features = rownames(obj))

total.transpose <- t(obj@assays$RNA@scale.data) # which assay to integrate with?

# counter <- 1
# Final_Matrix <- data.frame()
# for (i in 1:length(rownames(LINCS.ResponseSigs))){
# cmpdToOverlay <- rownames(LINCS.ResponseSigs)[i]
cmpdToOverlay <- "panobinostat"
# progress6$inc(1/length(rownames(LINCS.ResponseSigs)), detail = paste("Calculating Correlations Against ", cmpdToOverlay))
tcs_df <- as.data.frame(t(tcs_df))
head(rownames(tcs_df))
tcs_df$compound <- rownames(tcs_df)
cmpd <- subset(tcs_df, tcs_df$compound == cmpdToOverlay)
cmpd
cmpd$Genes <- NULL
cmpdGenes <- colnames(cmpd)
cmpd_overlap <- colnames(total.transpose)[which(colnames(total.transpose) %in% cmpdGenes)]
total.transpose.cmpd <- total.transpose[,cmpd_overlap] # FileA

cmpd_ordered <- as.numeric(as.vector(t(cmpd)))
names(cmpd_ordered) <- colnames(cmpd)
cmpd_ordered2 <- cmpd_ordered[cmpd_overlap] # Character list in brackets orders to fit that character list...
cmpd_ordered2[cmpd_overlap]
SC <- cor(cmpd_ordered2, t(total.transpose.cmpd), method = "spearman", use = "complete.obs")
sc_df <- as.data.frame(t(SC))
colnames(sc_df) <- cmpdToOverlay
obj <- AddMetaData(obj, metadata = sc_df)
obj$panobinostat_inverseCor <- obj$panobinostat*-1

library(viridis)
FeaturePlot(obj, features = c("panobinostat_inverseCor"), order = T, cols = viridis(256, option = "turbo"))

VlnPlot(obj, features = "panobinostat", group.by = "treatment")
# Final_Matrix <- rbind(Final_Matrix, SC)
# rownames(Final_Matrix) <- rownames(LINCS.ResponseSigs)[1:counter]
# counter <- counter + 1
# }

DimPlot(obj, group.by = c("patientID", "treatment"))

################################################################################
# save annotated RDS obj
################################################################################

saveRDS(obj, file = "zhao_integrated_TCSadded.RDS")
