# Can we use lme4 to identify drugs with sensitization changes?
# Then use the results to ID overlap between PDX and patient data.

library(lme4)
library(edgeR)
library(EnhancedVolcano)
library(Seurat)
library(dplyr)
library(ggpubr)
library(Hmisc)
library(readxl)
library(cowplot)
library(stringr)
library(viridis)
library(ComplexHeatmap)

# Patient Data
patient <- readRDS(file = "GBM_CombinationScoring_ISOSCELES_RV/isosceles_suterJohnsonMerge_fourStates_v5.RDS") # isosceles_suterJohnsonMerge_fourStates.RDS

Idents(patient) <- patient$Neftel_State
tumorCells <- WhichCells(patient, idents = c("AC", "MES", "NPC", "OPC"))

patient <- subset(patient, cells = tumorCells)

patientid <- data.frame(patient$patientID)
# rm(patient)

tcs_2020 <- read.csv(file = "GBM_CombinationScoring_ISOSCELES_RV/phase2_compound_TCS_v0.1_20230914_112256.csv", row.names = 1)
head(rownames(tcs_2020))
head(colnames(tcs_2020))

# read in combination summary data
wman <- read_excel(path = "GBM_CombinationScoring_ISOSCELES_RV/Synergy Scores.xlsx", sheet = 2)
# drugs to search for
inScreen <- gsub(x = tolower(unique(c(wman$`Drug 1`, wman$`Drug 2`))),
                 pattern = "-", replacement = ".", fixed = T)
screenIntersect <- intersect(tolower(colnames(tcs_2020)), inScreen)
screenIntersect

# what comparisons do we have then?

wman_sub <- wman[which(tolower(wman$`Drug 1`) %in% screenIntersect & tolower(wman$`Drug 2`) %in% screenIntersect),]
wman_sub

toScore <- unique(c(wman_sub$`Drug 1`, wman_sub$`Drug 2`))
toScore
# inscreendf <- as.data.frame(inScreen)
# colnamedf <- as.data.frame(colnames(tcs_2020))

colnames(tcs_2020) <- tolower(colnames(tcs_2020))

tcs_2020$erlotinib
class(toScore)
class(tcs_2020)
sub_tcs_2020 <- tcs_2020[,tolower(toScore)]
sub_tcs_2020

################################################################################
# use this subset of TCSs and score cell-drug connectivity
################################################################################

a <- as.data.frame(t(sub_tcs_2020)) # LINCS.ResponseSigs
a$compound <- rownames(a)

total.transpose <- t(patient@assays$RNA@scale.data) # which assay to integrate with?

# loop through compounds
counter <- 1
Final_Matrix <- data.frame()
for (i in 1:length(rownames(a))){
  cmpdToOverlay <- rownames(a)[i]
  # progress6$inc(1/length(rownames(LINCS.ResponseSigs)), detail = paste("Calculating Correlations Against ", cmpdToOverlay))
  cmpd <- subset(a, a$compound == cmpdToOverlay)
  cmpd$Genes <- NULL
  cmpdGenes <- colnames(cmpd)
  cmpd_overlap <- colnames(total.transpose)[which(colnames(total.transpose) %in% cmpdGenes)]
  total.transpose.cmpd <- total.transpose[,cmpd_overlap] # FileA
  cmpd_ordered <- as.numeric(as.vector(t(cmpd)))
  names(cmpd_ordered) <- colnames(cmpd)
  cmpd_ordered2 <- cmpd_ordered[cmpd_overlap] # Character list in brackets orders to fit that character list...
  # head(colnames(total.transpose.cmpd) == names(cmpd_ordered2))
  SC <- cor(cmpd_ordered2, t(total.transpose.cmpd), method = "spearman", use = "complete.obs")
  Final_Matrix <- rbind(Final_Matrix, SC)
  rownames(Final_Matrix) <- rownames(a)[1:counter]
  counter <- counter + 1
}

Final_Matrix


################################################################################
#
################################################################################

colnames(Final_Matrix)
rownames(patientid)

patient_drugdat <- Final_Matrix
patient_drugdat <- patient_drugdat + 1
d0 <- DGEList(patient_drugdat)
d0 <- calcNormFactors(d0)

snames <- colnames(patient_drugdat)
patientid <- subset(patientid, rownames(patientid) %in% colnames(patient_drugdat))


################################################################################
### loop through included compounds...
# need to generate sensitivity dataframe...

rownames(patient_drugdat) <- tolower(rownames(patient_drugdat))
screenIntersect
forSensitivity <- as.data.frame(t(patient_drugdat))
forSensitivity <- forSensitivity - 1
colnames(forSensitivity)

sens_df <- data.frame(row.names = rownames(forSensitivity))
# rownames(sens_df) <- rownames(forSensitivity)
for (drug in colnames(forSensitivity)){
  new_df <- as.data.frame(eval(parse(text = paste0("forSensitivity$", drug))))
  colnames(new_df) <- c(drug)
  rownames(new_df) <- rownames(forSensitivity)
  cutOff <- mean(eval(parse(text = paste0("new_df$", drug))))
  for (i in 1:length(rownames(new_df))){
    if (eval(parse(text = paste0("new_df$", drug)))[i] > cutOff){ # should the cut-off be set dynamically? i.e. geometric mean of all cells connectivities?
      new_df$sensitivity[i] <- "resistant"
    } else {
      new_df$sensitivity[i] <- "sensitive"
      # eval(parse(text = paste0("new_df$", drug, "_sensitivity")))[i] <- "sensitive"
    }
  }
  print(colnames(new_df))
  colnames(new_df) <- c(drug, paste0(drug, "_sensitivity"))
  print(dim(new_df))
  sens_df <- cbind(sens_df, new_df)
}
sens_df

dim(sens_df)
dim(patient)

################################################################################
#
################################################################################

# Visualization of connectivity

patient <- AddMetaData(patient, metadata = sens_df)

library(dittoSeq)

sensitivityBarplotList <- list()
for (col in grep(pattern = "*_sensitivity", x = colnames(sens_df), value = T)){
  sensitivityBarplotList[[col]] <- dittoSeq::dittoBarPlot(patient, var = col, group.by = "patientID")
}

pdf(file = "sensitivity_barPlots.pdf", width = 17.5)
plot_grid(plotlist = sensitivityBarplotList, ncol = 5)
dev.off()

################################################################################
#
################################################################################

# since np only how do I plot hierarchy here?

# erlotinib_featureUmap <- FeaturePlot(patient, features = c("erlotinib"), cols = viridis(200, option = "C", alpha = 0.2), order = T, reduction = "umap_int_merge") + theme_void()
# pazopanib_featureUmap <- FeaturePlot(patient, features = c("pazopanib"), cols = viridis(200, option = "C", alpha = 0.2), order = T, reduction = "umap_int_merge") + theme_void()
# lapatanib_featureUmap <- FeaturePlot(patient, features = c("lapatinib"), cols = viridis(200, option = "C", alpha = 0.2), order = T, reduction = "umap_int_merge") + theme_void()
# gemcitabine_featureUmap <- FeaturePlot(patient, features = c("gemcitabine"), cols = viridis(200, option = "C", alpha = 0.2), order = T, reduction = "umap_int_merge") + theme_void()
# thapsigargin_featureUmap <- FeaturePlot(patient, features = c("thapsigargin"), cols = viridis(200, option = "C", alpha = 0.2), order = T, reduction = "umap_int_merge") + theme_void()
# tipifarnib_featureUmap <- FeaturePlot(patient, features = c("tipifarnib"), cols = viridis(200, option = "C", alpha = 0.2), order = T, reduction = "umap_int_merge") + theme_void()
# obatoclax_featureUmap <- FeaturePlot(patient, features = c("obatoclax"), cols = viridis(200, option = "C", alpha = 0.2), order = T, reduction = "umap_int_merge") + theme_void()
# 
# pdf(file = "featureUMAPs.pdf")
# plot_grid(plotlist = list(erlotinib_featureUmap, pazopanib_featureUmap, lapatanib_featureUmap,
#                           gemcitabine_featureUmap, thapsigargin_featureUmap,
#                           tipifarnib_featureUmap, obatoclax_featureUmap),
#           ncol = 2,
#           labels = "AUTO")
# dev.off()

sens_df
################################################################################

## Loop through limma pipeline using sens_df sensitivity columns...

colsSens <- grep(x = colnames(sens_df), pattern = "*_sensitivity", value = T)
patient_mm_list <- list()
patient_y_list <- list()
patient_fit_list <- list()
for (comparison in colsSens){
  print(comparison)
  sensitivity <- as.data.frame(eval(parse(text = paste0("sens_df$", comparison))))
  rownames(sensitivity) <- rownames(sens_df)
  colnames(sensitivity) <- c(comparison)
  patient_group <- interaction(eval(parse(text = paste0("sensitivity$", comparison))),
                               patientid$patient.patientID)
  patient_mm_list[[comparison]] <- model.matrix(~0 + patient_group)
  pdf(file = paste0("mean_variance_voom_plots", comparison, ".pdf"))
  patient_y_list[[comparison]] <- voom(patient_drugdat, patient_mm_list[[comparison]], plot = T)
  dev.off()
  patient_fit_list[[comparison]] <- lmFit(patient_y_list[[comparison]], patient_mm_list[[comparison]])
}

# saveRDS(patient_mm_list, file = "patient_mm_list.RDS")
# saveRDS(patient_y_list, file = "patient_y_list.RDS")
# saveRDS(patient_fit_list, file = "patient_fit_list.RDS")
################################################################################

patient_mm_list <- readRDS(file = "GBM_CombinationScoring_ISOSCELES_RV/patient_mm_list.RDS")
patient_y_list <- readRDS(file = "GBM_CombinationScoring_ISOSCELES_RV/patient_y_list.RDS")
patient_fit_list <- readRDS(file = "GBM_CombinationScoring_ISOSCELES_RV/patient_fit_list.RDS")


patientContrasts <- c(
  "patient_groupresistant.GBM21-patient_groupsensitive.GBM21",
  "patient_groupresistant.GBM41-patient_groupsensitive.GBM41",
  "patient_groupresistant.GBM47-patient_groupsensitive.GBM47",
  "patient_groupresistant.GBM49-patient_groupsensitive.GBM49",
  "patient_groupresistant.GBM51-patient_groupsensitive.GBM51",
  "patient_groupresistant.GBM53-patient_groupsensitive.GBM53",
  "patient_groupresistant.SM006-patient_groupsensitive.SM006",
  "patient_groupresistant.SM011-patient_groupsensitive.SM011",
  "patient_groupresistant.SM012-patient_groupsensitive.SM012",
  "patient_groupresistant.SM017-patient_groupsensitive.SM017",
  "patient_groupresistant.SM018-patient_groupsensitive.SM018"
)

names(patient_fit_list)
colnames(patient_fit_list[["erlotinib_sensitivity"]])

compound_patient_resdf <- data.frame()
for (comparison in colsSens){
  patient_resdf <- data.frame()
  for(i in 1:length(patientContrasts)){
    print(i)
    patient_contr_loop <- makeContrasts(contrasts = patientContrasts[i], levels = colnames(coef(patient_fit_list[[comparison]])))
    patient_tmp_loop <- contrasts.fit(patient_fit_list[[comparison]], patient_contr_loop)
    patient_tmp_loop <- eBayes(patient_tmp_loop)
    patient_top.table_loop <- topTable(patient_tmp_loop, sort.by = "P", n = Inf)
    patient_top.table_loop$contrast <- patientContrasts[i]
    patient_top.table_loop$compound <- rownames(patient_top.table_loop)
    patient_resdf <- rbind(patient_resdf, patient_top.table_loop)
  }
  patient_resdf$reference_drug <- comparison
  compound_patient_resdf <- rbind(compound_patient_resdf, patient_resdf)
}

compound_patient_resdf

############################
############################
# write.csv(compound_patient_resdf, file = "GBM_combinationscreen_patient_DrugDiscordancelimma_results.csv", row.names = F)

# Now, loop through results and aggregate over individual patients for each reference drug...
unique(compound_patient_resdf$reference_drug)

limma_agg_list <- list()
for (comparison in colsSens){
  print(comparison)
  reference_sub <- subset(compound_patient_resdf, compound_patient_resdf$reference_drug == comparison)
  print(dim(reference_sub))
  cols <- colnames(reference_sub[which(colnames(reference_sub) %nin% c("compound", "contrast", "reference_drug"))])
  agg <- reference_sub %>%
    group_by(compound) %>%
    summarise(across(all_of(cols), mean), .groups = 'drop') %>%
    distinct() %>%
    as.data.frame()
  limma_agg_list[[comparison]] <- agg
}
limma_agg_list
# generic volcano of each ref drug...

volc_list_01 <- list()
volc_list_02 <- list()
for (comparison in colsSens){
  res <- limma_agg_list[[comparison]]
  print(res$compound[1:6])
  # pdf(file = paste0(comparison, "_limmaDiscordance_volc.pdf"))
  volc_list_01[[comparison]] <- EnhancedVolcano(toptable = res,
                  lab = res$compound,
                  x = "logFC", y = "adj.P.Val",
                  pCutoff = 0.05,
                  xlim = c(min(res[["logFC"]], na.rm = TRUE) - .05, max(res[["logFC"]], na.rm = TRUE) + .05),
                  FCcutoff = 0.01,
                  title = paste0(comparison, " Induced Sensitization & Desensitization"),
                  subtitle = "ISOSCELES + Limma")

  # what about a volcano with only these 10 compounds in it...?
  res_sub <- subset(res, res$compound %in% screenIntersect)
  volc_list_02[[comparison]] <- EnhancedVolcano(toptable = res_sub,
                   lab = res_sub$compound,
                   x = "logFC", y = "adj.P.Val",
                   pCutoff = 0.05,
                   xlim = c(min(res_sub[["logFC"]], na.rm = TRUE) - .05, max(res_sub[["logFC"]], na.rm = TRUE) + .05),
                   FCcutoff = 0.01,
                   title = paste0(comparison, " Induced Sensitization & Desensitization"),
                   subtitle = "ISOSCELES + Limma")
  # dev.off()

}
volc_list_01[[10]]
volc_list_02[[10]]

################################################################################
# Merge Limma Agg Lists w/ Resistant Cell Means
################################################################################


agg_rcm_merge_list <- list()
for (reference in names(limma_agg_list)){
  print(reference)
  res <- limma_agg_list[[reference]]
  # merge res with mean resistant population... ################################
  # pull sensitivity groupings
  sens <- as.data.frame(eval(parse(text = paste0("sens_df$", reference))))
  rownames(sens) <- rownames(sens_df)
  colnames(sens) <- c(reference)
  # print(head(rownames(sens)))
  # print(eval(parse(text = paste0("sens$", reference))))
  resistant_drugdat <- patient_drugdat[,rownames(subset(sens,
                       eval(parse(text = paste0("sens$", reference))) == "resistant"))]
  # print(head(resistant_drugdat))
  mean_resistant_connectivity <- as.data.frame(rowMeans(resistant_drugdat))
  colnames(mean_resistant_connectivity) <- c("Resistant_Cell_Connectivity")
  rcm <- merge(res, mean_resistant_connectivity, by.x = "compound", by.y = "row.names")
  rcm <- subset(rcm, rcm$logFC < 0)
  # print(rcm$Resistant_Cell_Connectivity[1:5])
  rcm$Resistant_Cell_Connectivity <- rcm$Resistant_Cell_Connectivity - 1
  disc_compounds <- rcm[which(rcm$Resistant_Cell_Connectivity < 0),]$compound
  conc_compounds <- rcm[which(rcm$Resistant_Cell_Connectivity > 0),]$compound
  rcm <- rcm %>% mutate(discordant = if_else(rcm$compound %in% disc_compounds, "discordant", "concordant"))
  rcm <- rcm %>% mutate(combinationScore = logFC * Resistant_Cell_Connectivity)
  # for (i in 1:length(rownames(rcm))){
  #   if (rcm$compound[i] %in% disc_compounds){
  #     rcm$discordant[i] <- "discordant"
  #   } else {
  #     rcm$discordant[i] <- "concordant"
  #   }
  # }
  agg_rcm_merge_list[[reference]] <- rcm
}

rcm_volc_plotlist <- list()
for (reference in names(agg_rcm_merge_list)){
  print(reference)
  res <- agg_rcm_merge_list[[reference]]
  # print(colnames(res))
  rcm_volc_plotlist[[reference]] <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = Resistant_Cell_Connectivity)) +
    geom_point() +
    geom_text(aes(label = compound),vjust = 0, hjust = 0) +
    viridis::scale_color_viridis(option = "D") +
    # scale_color_gradient(low = "blue", high = "red") + # Customize the color scale
    labs(x = "Connectivity Shift (logFC)", y = "-log10(P-Adj)") +
    theme_minimal()
}

pdf(file = "GBM_combinations_volcanoes.pdf")
plot_grid(plotlist = rcm_volc_plotlist, labels = gsub(pattern = "*_sensitivity", x = names(rcm_volc_plotlist), replacement = ""), ncol = 2)
dev.off()

rcmy_volc_plotlist <- list()
for (reference in names(agg_rcm_merge_list)){
  print(reference)
  res <- agg_rcm_merge_list[[reference]]
  # print(colnames(res))
  rcmy_volc_plotlist[[reference]] <- ggplot(res, aes(x = logFC, y = Resistant_Cell_Connectivity, color = -log10(adj.P.Val))) +
    geom_point() +
    geom_text(aes(label = compound),vjust = 0, hjust = 0) +
    viridis::scale_color_viridis(option = "D") +
    # scale_color_gradient(low = "blue", high = "red") + # Customize the color scale
    labs(x = "Connectivity Shift (logFC)", y = "MRC") +
    theme_minimal()
}

pdf(file = "GBM_combinations_volcanoes_02.pdf", width = 14, height = 10)
plot_grid(plotlist = rcmy_volc_plotlist, labels = gsub(pattern = "*_sensitivity", x = names(rcm_volc_plotlist), replacement = ""), ncol = 2)
dev.off()

combScore_volc_plotlist <- list()
for (reference in names(agg_rcm_merge_list)){
  print(reference)
  res <- agg_rcm_merge_list[[reference]]
  # print(colnames(res))
  combScore_volc_plotlist[[reference]] <- ggplot(res, aes(x = logFC, y = combinationScore, color = -log10(adj.P.Val))) +
    geom_point() +
    geom_text(aes(label = compound),vjust = 0, hjust = 0) +
    viridis::scale_color_viridis(option = "D") +
    # scale_color_gradient(low = "blue", high = "red") + # Customize the color scale
    labs(x = "Connectivity Shift (logFC)", y = "Combination Score") +
    theme_minimal()
}

pdf(file = "GBM_combinations_volcanoes_03.pdf", width = 14, height = 10)
plot_grid(plotlist = combScore_volc_plotlist, labels = gsub(pattern = "*_sensitivity", x = names(rcm_volc_plotlist), replacement = ""), ncol = 2)
dev.off()

# lapatinib & pazopanib
combScore_volc_plotlist <- list()
combScore_barplotlist <- list()
for (reference in names(agg_rcm_merge_list)[2:3]){
  print(reference)
  reference2 <- (gsub("*_sensitivity", x = reference, replacement = ""))
  res <- agg_rcm_merge_list[[reference]]
  res <- res[which(res$compound %in% c("gemcitabine", "obatoclax", "thapsigargin", "tipifarnib")),]
  # print(colnames(res))
  combScore_volc_plotlist[[reference]] <- ggplot(res, aes(x = logFC, y = Resistant_Cell_Connectivity, color = combinationScore)) +
    geom_point(aes(size = combinationScore)) +
    geom_text(aes(label = compound),vjust = 0, hjust = 0.2, nudge_y = 0.002, nudge_x = 0.0008) +
    viridis::scale_color_viridis(option = "D", end = 0.9) +
    # scale_color_gradient(low = "blue", high = "red") + # Customize the color scale
    labs(x = "Connectivity Shift (logFC)", y = "Mean Resistant Cell Connectivity", size = 6) +
    theme_minimal() + ggtitle(label = paste0("Reference Compound: ", reference2))

  res$compound <- reorder(res$compound, res$combinationScore)
  combScore_barplotlist[[reference]] <- ggplot(res, aes(x = compound, y = combinationScore, fill = combinationScore)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_c(option = "magma", end = 0.9) +
    labs(x = "Compound", y = "Combination Score", size = 6) +
    theme_minimal() + ggtitle(label = paste0("Reference Compound: ", reference2))

}

pdf(file = "GBM_combinations_volcanoes_select.pdf", width = 7, height = 14)
plot_grid(plotlist = combScore_volc_plotlist,
          labels = gsub(pattern = "*_sensitivity", x = names(combScore_volc_plotlist), replacement = ""),
          ncol = 1)
dev.off()

pdf(file = "GBM_combinations_barplots_select.pdf", width = 7, height = 14)
plot_grid(plotlist = combScore_barplotlist,
          labels = gsub(pattern = "*_sensitivity", x = names(combScore_barplotlist), replacement = ""),
          ncol = 1)
dev.off()


################################################################################
#
################################################################################

# Create a dataframe of reference drugs (colnames) with each partner drug's Combination Score as values

combinationScoreList <- lapply(agg_rcm_merge_list, function(x) data.frame(x$combinationScore, row.names = x$compound))
combinationScore_df <- bind_rows(combinationScoreList, .id = "Source")
combinationScore_df$compound <- gsub("^(.*?)\\.\\.\\..*", "\\1", x = rownames(combinationScore_df))
combinationScore_df <- as.data.frame(tidyr::pivot_wider(combinationScore_df, names_from = Source, values_from = x.combinationScore))
rownames(combinationScore_df) <- combinationScore_df$compound
combinationScore_df


wman_sub
# test case erlotinib
# for (reference in unique(c(wman$`Drug 1`, wman$`Drug 2`))){
combinationPredix_scatterlist <- list()
for (reference in names(agg_rcm_merge_list)){
  reference2 <- (gsub("*_sensitivity", x = reference, replacement = ""))
  print(reference2)
  wman_sub2 <- subset(wman_sub, tolower(wman_sub$`Drug 1`) == reference2 | tolower(wman_sub$`Drug 2`) == reference2)
  wman_sub2$`Drug 1` <- tolower(wman_sub2$`Drug 1`)
  wman_sub2$`Drug 2` <- tolower(wman_sub2$`Drug 2`)
  sub_cs <- as.data.frame(eval(parse(text = paste0("combinationScore_df$",
                                                   reference))),
                          row.names = rownames(combinationScore_df))
  colnames(sub_cs) <- tolower(reference2)
  print(sub_cs)


  forCorr <- merge(wman_sub2, sub_cs, by.x = 'Drug 2', by.y = "row.names")
  forCorr_8 <- na.omit(subset(forCorr, forCorr$Supp == 8))
  forCorr_9 <- na.omit(subset(forCorr, forCorr$Supp == 9))
  # print(forCorr)

  #   # 8
  eight <- cor(x = forCorr_8$`Average LOEWL Scores`, y = eval(parse(text = paste0("forCorr_8$", reference2))))
  print(eight)
  # # 9
  nine <- cor(x = forCorr_9$`Average BLISS Score`, y = eval(parse(text = paste0("forCorr_9$", reference2))))
  print(nine)

  # pdf(file = paste0(reference2, "_predictionScatter_test.pdf")
  combinationPredix_scatterlist[[reference]] <- ggscatter(data = forCorr_9, x = "Average BLISS Score", y = reference2, add = "reg.line", conf.int = T, cor.coef = T, cor.method = "spearman")
  # dev.off()
}

pdf(file = "predictionScatter_test_02.pdf")
plot(combinationPredix_scatterlist[[2]])
dev.off()
#

################################################################################
#
################################################################################

# Calculate correlations across individual cell lines...
library(dplyr)
bliss_scores <- as.data.frame(t(read.csv(file = "GBM_CombinationScoring_ISOSCELES_RV/supp_9_bliss.csv", row.names = 1)))
bliss_scores$combination <- tolower(rownames(bliss_scores))
bliss_scores <- bliss_scores %>%
  mutate(Drug1 = sapply(strsplit(combination, "...", fixed=T), function(x) x[1]))
bliss_scores <- bliss_scores %>%
  mutate(Drug2 = sapply(strsplit(combination, "...", fixed=T), function(x) x[2]))
bliss_scores$Drug2 <- sub("^\\.", "", bliss_scores$Drug2)
colnames(bliss_scores) <- trimws(colnames(bliss_scores), "right")
colnames(bliss_scores) <- gsub(x = colnames(bliss_scores), pattern = "-", replacement = "", fixed = T)

################################################################################
# lapatinib
################################################################################

multi_line_corPlot_list <- list()
multi_line_corPlot_list2 <- list()
overallCorr_list <- list()
names(agg_rcm_merge_list)
for (reference in names(agg_rcm_merge_list[2])){ # lapatinib

  # subset results for single reference compound
  ##############################################################################

  reference2 <- (gsub("*_sensitivity", x = reference, replacement = ""))
  print(reference2)
  bliss_sub <- bliss_scores[bliss_scores$Drug1 == reference2 | bliss_scores$Drug2 == reference2, ]
  bliss_sub$combination <- NULL
  sub_cs <- as.data.frame(eval(parse(text = paste0("combinationScore_df$",
                                                   reference))),
                          row.names = rownames(combinationScore_df))

  colnames(sub_cs) <- tolower(reference2)

  # insert ifelse here based on ref in drug1 or drug2
  # if the reference compound is in drug2, swap drug1 and drug2

  bliss_sub$new_Drug1 <- ifelse(bliss_sub$Drug2 == reference2, bliss_sub$Drug2, bliss_sub$Drug1)
  bliss_sub$new_Drug2 <- ifelse(bliss_sub$Drug2 == reference2, bliss_sub$Drug1, bliss_sub$Drug2)
  bliss_sub$Drug1 <- bliss_sub$new_Drug1
  bliss_sub$Drug2 <- bliss_sub$new_Drug2
  bliss_sub <- bliss_sub[,!(colnames(bliss_sub) %in% c("new_Drug1", "new_Drug2"))]

  forCorr <- na.omit(merge(bliss_sub, sub_cs, by.x = 'Drug2', by.y = "row.names"))
  forCorr

  if (length(unique(forCorr$Drug2)) != 1){

    print("True and plotting...")

    mlt_df <- reshape2::melt(forCorr)
    mlt_df <- merge(mlt_df, forCorr[c("Drug2", reference2)], by.x = "Drug2", by.y = "Drug2")
    mlt_df <- mlt_df[which(mlt_df$variable != reference2),]
    ##############################################################################
    # Dplyr spearman correlations
    ##############################################################################

    # dplyr spearman corr
    correlations <- mlt_df %>%
      group_by(variable) %>%
      summarise(correlation = cor(eval(parse(text = reference2)),
                                  value,
                                  method = "spearman"),
                p_value = cor.test(eval(parse(text = reference2)),
                                   value,
                                   method = "spearman")$p.value)

    overallCorr <- cor.test(eval(parse(text = paste0("mlt_df$",reference2))),
                            mlt_df$value,
                            method = "spearman")
    overallCorr_list[[reference]] <- overallCorr
    ##############################################################################
    # plot 1
    ##############################################################################

    # Create a scatter plot with lines connecting points within the same group
    multi_line_corPlot_list[[reference]] <- ggplot(mlt_df,
                                                   aes(x = eval(parse(text = reference2)),
                                                       y = value,
                                                       color = variable,
                                                       group = variable)) +
      # geom_point() +
      geom_smooth(method = "glm", se = FALSE) +
      labs(x = "Combination Score", y = "BLISS Raw Sum Synergy") +
      theme_minimal()

    ##############################################################################
    # Plot 2
    ##############################################################################

    multi_line_corPlot_list2[[reference]] <- ggplot(mlt_df,
                                                    aes(x = eval(parse(text = reference2)),
                                                        y = value,
                                                        color = variable,
                                                        group = variable
                                                    )) +
      # geom_point() +
      geom_smooth(method = "lm", se = FALSE, aes(alpha = 0.2), size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, aes(alpha = 0.2, group = 1), color = "red") +
      scale_color_viridis_d(alpha = 0.3) +
      annotate("text", x = min(eval(parse(text = paste0("mlt_df$", reference2)))), y = max(mlt_df$value) - 1,
               label = paste("Spearman rho =", round(overallCorr$estimate, 2),
                             "p-value =", format(overallCorr$p.value, digits = 2)),
               hjust = -0.1, vjust = 1, size = 5) +
      labs(x = "ISOSCELES Combination Score", y = "BLISS Raw Sum Synergy (Houweling et al, 2023 Screen)") +
      ggtitle(label = paste0("Ref. Compound: ", reference2)) +
      theme_minimal() # + NoLegend()

    ##############################################################################
  }
}

pdf(file = "lapatinib_ref.pdf")
plot_grid(plotlist = multi_line_corPlot_list2)
dev.off()

################################################################################
# pazopanib
################################################################################

multi_line_corPlot_list <- list()
multi_line_corPlot_list2 <- list()
overallCorr_list <- list()
names(agg_rcm_merge_list)
for (reference in names(agg_rcm_merge_list[3])){ # pazopanib

  # subset results for single reference compound
  ##############################################################################

  reference2 <- (gsub("*_sensitivity", x = reference, replacement = ""))
  print(reference2)
  bliss_sub <- bliss_scores[bliss_scores$Drug1 == reference2 | bliss_scores$Drug2 == reference2, ]
  bliss_sub$combination <- NULL
  sub_cs <- as.data.frame(eval(parse(text = paste0("combinationScore_df$",
                                                   reference))),
                          row.names = rownames(combinationScore_df))

  colnames(sub_cs) <- tolower(reference2)

  # insert ifelse here based on ref in drug1 or drug2
  # if the reference compound is in drug2, swap drug1 and drug2

  bliss_sub$new_Drug1 <- ifelse(bliss_sub$Drug2 == reference2, bliss_sub$Drug2, bliss_sub$Drug1)
  bliss_sub$new_Drug2 <- ifelse(bliss_sub$Drug2 == reference2, bliss_sub$Drug1, bliss_sub$Drug2)
  bliss_sub$Drug1 <- bliss_sub$new_Drug1
  bliss_sub$Drug2 <- bliss_sub$new_Drug2
  bliss_sub <- bliss_sub[,!(colnames(bliss_sub) %in% c("new_Drug1", "new_Drug2"))]

  forCorr <- na.omit(merge(bliss_sub, sub_cs, by.x = 'Drug2', by.y = "row.names"))
  forCorr
  if (length(unique(forCorr$Drug2)) != 1){

    print("True and plotting...")

    mlt_df <- reshape2::melt(forCorr)
    mlt_df <- merge(mlt_df, forCorr[c("Drug2", reference2)], by.x = "Drug2", by.y = "Drug2")
    mlt_df <- mlt_df[which(mlt_df$variable != reference2),]
    ##############################################################################
    # Dplyr spearman correlations
    ##############################################################################

    # dplyr spearman corr
    correlations <- mlt_df %>%
      group_by(variable) %>%
      summarise(correlation = cor(eval(parse(text = reference2)),
                                  value,
                                  method = "spearman"),
                p_value = cor.test(eval(parse(text = reference2)),
                                   value,
                                   method = "spearman")$p.value)

    overallCorr <- cor.test(eval(parse(text = paste0("mlt_df$",reference2))),
                            mlt_df$value,
                            method = "spearman")
    overallCorr_list[[reference]] <- overallCorr
    ##############################################################################
    # plot 1
    ##############################################################################

    # Create a scatter plot with lines connecting points within the same group
    multi_line_corPlot_list[[reference]] <- ggplot(mlt_df,
                                                   aes(x = eval(parse(text = reference2)),
                                                       y = value,
                                                       color = variable,
                                                       group = variable)) +
      # geom_point() +
      geom_smooth(method = "glm", se = FALSE) +
      labs(x = "Combination Score", y = "BLISS Raw Sum Synergy") +
      theme_minimal()

    ##############################################################################
    # Plot 2
    ##############################################################################

    multi_line_corPlot_list2[[reference]] <- ggplot(mlt_df,
                                                    aes(x = eval(parse(text = reference2)),
                                                        y = value,
                                                        color = variable,
                                                        group = variable
                                                    )) +
      # geom_point() +
      geom_smooth(method = "lm", se = FALSE, aes(alpha = 0.2), size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, aes(alpha = 0.2, group = 1), color = "red") +
      scale_color_viridis_d(alpha = 0.3) +
      annotate("text", x = min(eval(parse(text = paste0("mlt_df$", reference2)))), y = max(mlt_df$value) - 1,
               label = paste("Spearman rho =", round(overallCorr$estimate, 2),
                             "p-value =", format(overallCorr$p.value, digits = 2)),
               hjust = -0.1, vjust = 1, size = 5) +
      labs(x = "ISOSCELES Combination Score", y = "BLISS Raw Sum Synergy (Houweling et al, 2023 Screen)") +
      ggtitle(label = paste0("Ref. Compound: ", reference2)) +
      theme_minimal() # + NoLegend()

    ##############################################################################
  }
}

pdf(file = "pazopanib_ref.pdf")
plot_grid(plotlist = multi_line_corPlot_list2)
dev.off()


################################################################################
################################################################################
################################################################################

################################################################################
# complex heatmap for combination figure
################################################################################


# res <- agg_rcm_merge_list[[2]]
# names(res)

heatmap_mergelist <- list()
for (reference in names(agg_rcm_merge_list)){ # pazopanib

  # subset results for single reference compound
  ##############################################################################

  reference2 <- (gsub("*_sensitivity", x = reference, replacement = ""))
  print(reference2)
  bliss_sub <- bliss_scores[bliss_scores$Drug1 == reference2 | bliss_scores$Drug2 == reference2, ]
  bliss_sub$combination <- NULL
  sub_cs <- as.data.frame(eval(parse(text = paste0("combinationScore_df$",
                                                   reference))),
                          row.names = rownames(combinationScore_df))

  colnames(sub_cs) <- tolower(reference2)

  # insert ifelse here based on ref in drug1 or drug2
  # if the reference compound is in drug2, swap drug1 and drug2

  bliss_sub$new_Drug1 <- ifelse(bliss_sub$Drug2 == reference2, bliss_sub$Drug2, bliss_sub$Drug1)
  bliss_sub$new_Drug2 <- ifelse(bliss_sub$Drug2 == reference2, bliss_sub$Drug1, bliss_sub$Drug2)
  bliss_sub$Drug1 <- bliss_sub$new_Drug1
  bliss_sub$Drug2 <- bliss_sub$new_Drug2
  bliss_sub <- bliss_sub[,!(colnames(bliss_sub) %in% c("new_Drug1", "new_Drug2"))]

  forCorr <- na.omit(merge(bliss_sub, sub_cs, by.x = 'Drug2', by.y = "row.names"))
  heatmap_mergelist[[reference2]] <- forCorr
}

total_hm_df <- bind_rows(heatmap_mergelist)
# coerce combination score values to a single column...
total_hm_df <- total_hm_df %>%
  rowwise() %>%
  mutate(CombinationScore = coalesce(erlotinib,
                                     lapatinib,
                                     pazopanib,
                                     sunitinib,
                                     thapsigargin, #
                                     tipifarnib,
                                     gemcitabine,
                                     obatoclax, #
                                     vinorelbine, #
                                     midostaurin)) #

total_hm_df <- total_hm_df[which(total_hm_df$Drug1 %nin% c("thapsigargin", "obatoclax", "vinorelbine", "midostaurin")),]

total_hm_df <- total_hm_df %>% mutate(colAnnotation = CombinationScore) %>% as.data.frame()
total_hm_df$colAnnotation
rownames(total_hm_df) <- paste0(total_hm_df$Drug1, "_", total_hm_df$Drug2)
total_hm_df

mat_cols <- grep(pattern = "^G", x = colnames(total_hm_df), value = T)
mat_cols

total_hm_df
total_hm_df <- total_hm_df[-which(rownames(total_hm_df) %in% c("gemcitabine_erlotinib", "gemcitabine_pazopanib", "tipifarnib_sunitinib", "tipifarnib_pazopanib", "tipifarnib_lapatinib")),]
total_hm_df

combScore_min <- min(total_hm_df$colAnnotation)
combScore_max <- max(total_hm_df$colAnnotation)

total_hm_df <- total_hm_df %>% mutate(colAnnotation2 = (colAnnotation - combScore_min) / (combScore_max - combScore_min))

unique(total_hm_df)
colMeans(t(total_hm_df[,mat_cols]))



column_ha <- HeatmapAnnotation(ScaledCombinationScore = anno_barplot(total_hm_df$colAnnotation2, gp = gpar(fill = 1:10)),
                               # meanBLISS = anno_simple(colMeans(t(total_hm_df[,mat_cols]))),
                               # sumBLISS = anno_barplot(colSums(t(total_hm_df[,mat_cols]))),
                               height = unit(2, "cm"))
column_ha

overall_corr <- cor.test(total_hm_df$colAnnotation2, colMeans(t(total_hm_df[,mat_cols])), method = "spearman")

combinationVsBLISS_allLines <- data.frame(CombinationScore=total_hm_df$colAnnotation2,
                                          meanBLISS=colMeans(t(total_hm_df[,mat_cols])),
                                          combination = rownames(total_hm_df),
                                          referenceDrug = total_hm_df$Drug1)

# ggscatter(combinationVsBLISS_allLines, x = "CombinationScore", y = "meanBLISS", cor.coef = T, cor.method = "spearman", add = c("reg.line"), conf.int = T, label = "combination", color = "referenceDrug")


write.csv(combinationVsBLISS_allLines, file = "CombinationIndex_vs_InVitroBLISSS.csv")

overall_corr_plot <- ggplot(data = combinationVsBLISS_allLines, aes(x = CombinationScore, y = meanBLISS, color = referenceDrug)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = c("#440154FF",
                                "#31688EFF",
                                "#35B779FF",
                                "#FDE725FF")) +
  labs(title = "Scatter Plot with Colored Points and Regression Line") +
  geom_text(aes(x = min(CombinationScore),
                y = max(meanBLISS),
                label = paste("Spearman's rho = ", round(overall_corr$estimate, 2), "\n", "p-value = ", format(overall_corr$p.value, scientific = TRUE)), hjust = 0, vjust = 1)) +
  geom_text(aes(label = combination), vjust = -0.5, hjust = 1)

pdf(file = "overall_corr_plot.pdf")
overall_corr_plot
dev.off()

corr_res_df <- data.frame()
for (col in mat_cols){
  print(col)
  print(cor.test(total_hm_df$colAnnotation2, eval(parse(text = paste0("total_hm_df$", col))), method = "spearman"))
  cor_res <- cor.test(total_hm_df$colAnnotation2, eval(parse(text = paste0("total_hm_df$", col))), method = "spearman")
  cor_res_df1 <- data.frame(rho = cor_res$estimate, p.val = cor_res$p.value)
  print(cor_res_df1)
  cor_res_df1$cellLine <- col
  print(cor_res_df1)
  corr_res_df <- rbind(corr_res_df, cor_res_df1)
}

corr_res_df

pvalue = corr_res_df$p.val
is_sig = pvalue < 0.05
pch = rep("*", length(pvalue))
pch[!is_sig] = NA

# color mapping
# pvalue_col_fun = viridis(length(pvalue), option = "B")
# col_fun = circlize::colorRamp2(seq(0,1,length = 256), viridis(256, option = "D"))
col_fun <- circlize::colorRamp2(seq(0,1, length.out = 256), viridis(256, option = "B"))
col_fun2 <- circlize::colorRamp2(seq(0,1, length.out = 256), viridis(256, option = "C"))
row_ha <- HeatmapAnnotation(which = "row", SpearmanRho = anno_simple(corr_res_df$rho, col = col_fun2),
                        pvalue = anno_simple(-log10(pvalue), col = col_fun, pch = pch))

row_ha

ref_anno <- total_hm_df$Drug1
ref_color_palette <- viridis_pal(option = "D")(length(unique(ref_anno)))
# ref_annotation_colors <- ref_color_palette[match(ref_anno, unique(ref_anno))]


bottom_ha <- HeatmapAnnotation(which = "column",
                               reference = anno_simple(total_hm_df$Drug1,
                                                       col = c('erlotinib' = "#440154FF",
                                                               'lapatinib' = "#31688EFF",
                                                               'pazopanib' = "#35B779FF",
                                                               'sunitinib' = "#FDE725FF")))

heatmap <- Heatmap(
  matrix = t(as.matrix(total_hm_df[,mat_cols])),
  name = "BLISS raw sum synergy",
  # name = "Value1",  # Specify the column to display in the heatmap
  column_title = "Heatmap Title",
  row_title = "Row Annotation",
  # row_names_gp = gpar(fontsize = 12, col = "blue"),  # Customize row annotation
  top_annotation = column_ha,
  right_annotation = row_ha,
  bottom_annotation = bottom_ha,
  width = unit(8, "cm"),
  height = unit(10, "cm"),
  # column_km = 2,
  # ColumnAnnotation = anno_barplot(hm_df$colAnnotation,
  #                                 col = "red",
  #                                 border = "black")
  # )
  col = viridis(100, option = "D"),
  cluster_columns = F
  # col = colorRamp2(c(0, 10), c("white", "red")),  # Customize color scale
  # width = unit(3, "cm"),  # Specify the width of the heatmap
  # cluster_rows = FALSE,  # Disable row clustering
  # show_row_names = TRUE,  # Show row names
  # show_column_names = TRUE  # Show column names
)

# now we generate two legends, one for the p-value
# see how we define the legend for pvalue
lgd_pvalue = Legend(title = "p-value", col_fun = col_fun, at = c(0, 0.25, 0.5, 0.75, 1),
                    labels = rev(c(round(quantile(pvalue), digits = 2))))

lgd_rho = Legend(title = "Spearman Rho", col_fun = col_fun2, at = c(0, 0.25, 0.5, 0.75),
                 labels = c(round(quantile(corr_res_df$rho), digits = 2)))

# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.05")

pdf(file = "heatmaptest.pdf", height = 10, width = 20)
draw(heatmap, annotation_legend_list = list(lgd_pvalue, lgd_sig, lgd_rho))
dev.off()

write.csv(corr_res_df, file = "IndividualCellLine_scFOCAL_CombinationIndexVsBLISS_Spearman.csv")

write.csv(total_hm_df[,mat_cols], file = "CombinationInVitro_HeatmapSourceData.csv")
################################################################################


