# Can we use lme4 to identify drugs with sensitization changes? 
# Then use the results to ID overlap between PDX and patient data. 

library(lme4)
library(edgeR)
library(EnhancedVolcano)
library(Seurat)
library(dplyr)
library(ggpubr)

# Need to fix cell IDs in the corrmat... gsub periods for -'s

# Patient Data
patient <- readRDS(file = "isosceles_suterJohnsonMerge_fourStates.RDS")
Idents(patient) <- patient$Neftel_State
tumorCells <- WhichCells(patient, idents = c("AC", "MES", "NPC", "OPC"))

# newrownames <- gsub(x = rownames(patient_corrmat), pattern = "\\.", replacement = "-")
# rownames(patient_corrmat) <- newrownames

# patient_drugdat <- t(patient_corrmat + 1)
# patient_drugdat <- patient_drugdat[,colnames(patient_drugdat) %in% tumorCells]
# saveRDS(patient_drugdat, file = "limma/patient_drugdat_prefilt.rds")
patient_drugdat <- readRDS(file = "limma/patient_drugdat_prefilt.rds")
################################################################################
# Filter patient drugdat here for compounds with atleast n number genes

tcs_size <- read.csv(file = "L1000_2017_TCS_size.csv", row.names = 1)
colnames(tcs_size)

# need to match strings in tcs size file and in drugdat...
dim(tcs_size)
newRownames <- gsub(x = rownames(tcs_size), pattern = "-", replacement = ".", fixed = T)
rownames(tcs_size) <- newRownames
length(intersect(rownames(tcs_size), rownames(patient_drugdat)))

tcs_threshold = 0
keep <- rownames(subset(tcs_size, tcs_size$nonzero > tcs_threshold))
dim(tcs_size)
length(keep)
keep

################################################################################

patient_drugdat <- patient_drugdat[rownames(patient_drugdat) %in% keep,]

d0 <- DGEList(patient_drugdat)
d0 <- calcNormFactors(d0)

snames <- colnames(patient_drugdat)
patientid <- data.frame(patient$patientID)
patientid <- subset(patientid, rownames(patientid) %in% colnames(patient_drugdat))
sensitivity <- data.frame(patient$sensitivity)
sensitivity <- subset(sensitivity, rownames(sensitivity) %in% tumorCells)
patient_group <- interaction(sensitivity$patient.sensitivity, patientid$patient.patientID)

patient_mm <- model.matrix(~0 + patient_group)
patient_mm

patient_y <- voom(patient_drugdat, patient_mm, plot = T)

# save group, y, patient_drugdat, mm as rds files...
saveRDS(patient_drugdat, file = "limma/patient_drugdat.rds")
saveRDS(patient_group, file = "limma/patient_group.rds")
saveRDS(patient_mm, file = "limma/patient_mm.rds")
saveRDS(patient_y, file = "limma/patienty.rds")

################################################################################
# PDX data

dmso_mat2 <- readRDS(file = "pdx_patient_comparison/dmso_mat2.RDS")
alisertib_mat2 <- readRDS(file = "pdx_patient_comparison/alisertib_mat2.RDS")

pdx_drugdat <- rbind(dmso_mat2, alisertib_mat2)
pdx_drugdat <- t(pdx_drugdat + 1)
pdx_drugdat <- pdx_drugdat[rownames(pdx_drugdat) %in% keep,]

pdx_d0 <- DGEList(pdx_drugdat)
pdx_d0 <- calcNormFactors(pdx_d0)

dmso <- readRDS(file = "pdx_patient_comparison/dmso_isosceled.RDS")
alisertib <- readRDS(file = "pdx_patient_comparison/alisertib_isosceled.RDS")
pdxArmDf <- rbind(data.frame(Arm = dmso$Arm), data.frame(Arm = alisertib$Arm))
# coerce order of drugdat to match grouping 
pdx_drugdat <- pdx_drugdat[,rownames(pdxArmDf)]
pdx_group <- interaction(pdxArmDf$Arm)

pdx_mm <- model.matrix(~0 + pdx_group)
pdx_mm

pdx_y <- voom(pdx_drugdat, pdx_mm, plot = T)

# save group, y, patient_drugdat, mm as rds files...
saveRDS(pdx_drugdat, file = "limma/pdx_drugdat.rds")
saveRDS(pdx_group, file = "limma/pdx_group.rds")
saveRDS(pdx_mm, file = "limma/pdx_mm.rds")
saveRDS(pdx_y, file = "limma/pdxy.rds")

################################################################################

# to filter low-expressed genes normally, let's proceed without this...
# probably replace with filter for signatures above a certain threshold

# cutoff <- 1
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# d <- d0[-drop,] 
# dim(d) # number of genes left

################################################################################



################################################################################
# Load precomputed models
################################################################################
patient_drugdat <- readRDS(file = "limma/patient_drugdat.rds")
patient_group <- readRDS(file = "limma/patient_group.rds")
patient_mm <- readRDS(file = "limma/patient_mm.rds")
patient_y <- readRDS(file = "limma/patienty.rds")

#PDX
pdx_drugdat <- readRDS(file = "limma/pdx_drugdat.rds")
pdx_group <- readRDS(file = "limma/pdx_group.rds")
pdx_mm <- readRDS(file = "limma/pdx_mm.rds")
pdx_y <- readRDS(file = "limma/pdxy.rds")

################################################################################

patient_fit <- lmFit(patient_y, patient_mm)
head(coef(patient_fit))

# comparisons

colnames(coef(patient_fit))

# PDX
pdx_fit <- lmFit(pdx_y, pdx_mm)
head(coef(pdx_fit))

# comparisons

colnames(coef(pdx_fit))

saveRDS(patient_fit, file = "limma/patient_fit.rds")
saveRDS(pdx_fit, file = "limma/pdx_fit.rds")

# patientContrastList <- list(
#   c("groupresistant.GBM21", "groupsensitive.GBM21"),
#   c("groupresistant.GBM41", "groupsensitive.GBM41"),
#   c("groupresistant.GBM47", "groupsensitive.GBM47"),
#   c("groupresistant.GBM49", "groupsensitive.GBM49"),
#   c("groupresistant.GBM51", "groupsensitive.GBM51"),
#   c("groupresistant.GBM53", "groupsensitive.GBM53"),
#   c("groupresistant.SM006", "groupsensitive.SM006"),
#   c("groupresistant.SM011", "groupsensitive.SM011"),
#   c("groupresistant.SM012", "groupsensitive.SM012"),
#   c("groupresistant.SM017", "groupsensitive.SM017"),
#   c("groupresistant.SM018", "groupsensitive.SM018")
# )

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
) # how do we automate this from a selection of variables in shiny?

pdxContrasts <- c(
  "pdx_groupAlisertib-pdx_groupDMSO"
)

# Patient
patient_resdf <- data.frame()
for(i in 1:length(patientContrasts)){
  print(i)
  patient_contr_loop <- makeContrasts(contrasts = patientContrasts[i], levels = colnames(coef(patient_fit)))
  patient_tmp_loop <- contrasts.fit(patient_fit, patient_contr_loop)
  patient_tmp_loop <- eBayes(patient_tmp_loop)
  patient_top.table_loop <- topTable(patient_tmp_loop, sort.by = "P", n = Inf)
  patient_top.table_loop$contrast <- patientContrasts[i]
  patient_top.table_loop$compound <- rownames(patient_top.table_loop)
  patient_resdf <- rbind(patient_resdf, patient_top.table_loop)
}

# PDX
pdx_resdf <- data.frame()
for(i in 1:length(pdxContrasts)){
  print(i)
  pdx_contr_loop <- makeContrasts(contrasts = pdxContrasts[i], levels = colnames(coef(pdx_fit)))
  pdx_tmp_loop <- contrasts.fit(pdx_fit, pdx_contr_loop)
  pdx_tmp_loop <- eBayes(pdx_tmp_loop)
  pdx_top.table_loop <- topTable(pdx_tmp_loop, sort.by = "P", n = Inf)
  pdx_top.table_loop$contrast <- pdxContrasts[i]
  pdx_top.table_loop$compound <- rownames(pdx_top.table_loop)
  pdx_resdf <- rbind(pdx_resdf, pdx_top.table_loop)
}

write.csv(patient_resdf, file = "patient_DrugDiscordancelimma_results.csv", row.names = F)
write.csv(pdx_resdf, file = "pdx_DrugDiscordancelimma_results.csv", row.names = F)

# goal: identify compounds consistently made sensitive to across all patients...
################################################################################
# Patient
################################################################################
# What are the top 5 desensitized? "alisertib"   "GSK.1070916" "NVP.TAE684"  "palbociclib" "AZD.8330"    "GSK.2126458"
patient_desensdf <- patient_resdf %>% group_by(contrast) %>% top_n(wt = logFC, n = 10) 
dim(patient_desensdf)
length(patient_desensdf$compound[duplicated(patient_desensdf$compound)])
unique(patient_desensdf$compound[duplicated(patient_desensdf$compound)])

# What are the top 5 sensitized? 
patient_sensdf <- patient_resdf %>% group_by(contrast) %>% top_n(wt = -logFC, n = 10)
dim(patient_sensdf)
length(patient_sensdf$compound[duplicated(patient_sensdf$compound)])
unique(patient_sensdf$compound[duplicated(patient_sensdf$compound)])

patient_desensdf <- patient_resdf %>% group_by(contrast) %>% dplyr::filter(logFC > 0.025)
patient_sensdf <- patient_resdf %>% group_by(contrast) %>% dplyr::filter(logFC < -0.025)
unique(patient_desensdf$compound)
unique(patient_sensdf$compound)
# For each comparison, what are the overlapping compounds across x number patients
# below a threshold p value?

################################################################################
# PDX
################################################################################
# What are the top 5 desensitized? "alisertib"   "GSK.1070916" "NVP.TAE684"  "palbociclib" "AZD.8330"    "GSK.2126458"
pdx_desensdf <- pdx_resdf %>% group_by(contrast) %>% top_n(wt = logFC, n = 10) 
dim(pdx_desensdf)
length(pdx_desensdf$compound)
unique(pdx_desensdf$compound)

# What are the top 5 sensitized? 
pdx_sensdf <- pdx_resdf %>% group_by(contrast) %>% top_n(wt = -logFC, n = 10)
dim(pdx_sensdf)
length(pdx_sensdf$compound)
unique(pdx_sensdf$compound)

hist(pdx_resdf$logFC)
pdx_desensdf <- pdx_resdf %>% group_by(contrast) %>% dplyr::filter(logFC > 0.005)
pdx_sensdf <- pdx_resdf %>% group_by(contrast) %>% dplyr::filter(logFC < -0.005)
unique(pdx_desensdf$compound)
unique(pdx_sensdf$compound)

intersect(unique(pdx_sensdf$compound), patient_sensdf$compound)

EnhancedVolcano(pdx_resdf, lab = rownames(pdx_resdf), x = "logFC", y = "adj.P.Val", 
                                xlim = c(min(pdx_resdf[["logFC"]], na.rm = TRUE) - .01,
                                         max(pdx_resdf[["logFC"]], na.rm = TRUE) + .01), FCcutoff = 0.005, pCutoff = 5*10^-50)

# For each comparison, what are the overlapping compounds across x number patients
# below a threshold p value?

rownames(patient_resdf)

patient_resdf3 <- patient_resdf[patient_resdf$adj.P.Val < 0.05,]
pdx_resdf3 <- pdx_resdf[pdx_resdf$adj.P.Val < 0.05,]

EnhancedVolcano(toptable = pdx_resdf3, 
                lab = rownames(pdx_resdf3), 
                x = "logFC", y = "adj.P.Val", FCcutoff = 0.01,
                xlim = c(min(pdx_resdf3[["logFC"]], na.rm = TRUE) - .01, max(pdx_resdf3[["logFC"]], na.rm = TRUE) + .01))

EnhancedVolcano(toptable = patient_resdf3, 
                lab = rownames(patient_resdf3), 
                x = "logFC", y = "adj.P.Val", FCcutoff = 0.01,
                xlim = c(min(patient_resdf3[["logFC"]], na.rm = TRUE) - .01, max(patient_resdf3[["logFC"]], na.rm = TRUE) + .01))




patient_resdf2 <- patient_resdf[c("logFC", "compound", "P.Value", "adj.P.Val")]
collapsed_patient_res <- patient_resdf2 %>% group_by(compound) %>% mutate(across(.fns = mean)) %>% distinct()
collapsed_patient_res <- as.data.frame(collapsed_patient_res)
rownames(collapsed_patient_res) <- collapsed_patient_res$compound
pdx_resdf <- pdx_resdf[rownames(collapsed_patient_res),]

EnhancedVolcano(toptable = collapsed_patient_res, 
                lab = rownames(collapsed_patient_res), 
                x = "logFC", y = "adj.P.Val", FCcutoff = 0.01,
                xlim = c(min(patient_resdf3[["logFC"]], na.rm = TRUE) - .01, max(patient_resdf3[["logFC"]], na.rm = TRUE) + .01))


patient_resdf4 <- patient_resdf3[c("logFC", "compound")]
collapsed_patient_res4 <- patient_resdf4 %>% group_by(compound) %>% mutate(across(.fns = mean)) %>% distinct()
collapsed_patient_res4 <- as.data.frame(collapsed_patient_res4)
rownames(collapsed_patient_res4) <- collapsed_patient_res4$compound
pdx_resdf4 <- pdx_resdf3[rownames(collapsed_patient_res),]

comp_df <- data.frame(row.names = rownames(pdx_resdf), pdx_logFC = pdx_resdf$logFC, patient_logFC = collapsed_patient_res$logFC)
comp_df2 <- data.frame(row.names = rownames(pdx_resdf4), pdx_logFC = pdx_resdf4$logFC, patient_logFC = collapsed_patient_res4$logFC)


pdf(file = "limma_discordanceShiftScatter.pdf")
ggscatter(data = comp_df, x = "pdx_logFC", y = "patient_logFC", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson")
dev.off()

neg_comp_df <- comp_df2[comp_df2$pdx_logFC < 0 & comp_df2$patient_logFC < 0,]
neg_comp_df$compound <- rownames(neg_comp_df)
ggscatter(data = neg_comp_df, x = "pdx_logFC", y = "patient_logFC", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson", label = "compound")

pos_comp_df <- comp_df2[comp_df2$pdx_logFC > 0 & comp_df2$patient_logFC > 0,]
pos_comp_df$compound <- rownames(pos_comp_df)
ggscatter(data = pos_comp_df, x = "pdx_logFC", y = "patient_logFC", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson", label = "compound")


pdf(file = "limma_discordanceShiftScatter2.pdf")
ggscatter(data = comp_df2, x = "pdx_logFC", y = "patient_logFC", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson")
dev.off()

################################################################################

# MOA subplotting

rep <- read.delim(file = "repurposing_drugs_20200324.txt", skip = 9, sep = "\t")
rep$new_iname <- gsub(pattern = "-", replacement = ".", x = rep$pert_iname, fixed = T)


pi3k <- subset(rep, rep$moa == "PI3K inhibitor")

pi3k_subdf <- subset(comp_df, rownames(comp_df) %in% pi3k$new_iname)
colnames(pi3k_subdf)
ggscatter(data = pi3k_subdf, x = "pdx_logFC", y = "patient_logFC", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson")

mtor <- subset(rep, rep$moa == "mTOR inhibitor")
mtor_subdf <- subset(comp_df, rownames(comp_df) %in% mtor$new_iname)
ggscatter(data = mtor_subdf, x = "pdx_logFC", y = "patient_logFC", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson")


aurora <- subset(rep, rep$moa == "Aurora kinase inhibitor")
aurora_subdf <- subset(comp_df, rownames(comp_df) %in% aurora$new_iname)
colnames(aurora_subdf)
ggscatter(data = aurora_subdf, x = "pdx_logFC", y = "patient_logFC", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson")

BET <- subset(rep, rep$moa == "bromodomain inhibitor")
BET_subdf <- subset(comp_df, rownames(comp_df) %in% BET$new_iname)
colnames(BET_subdf)
ggscatter(data = BET_subdf, x = "pdx_logFC", y = "patient_logFC", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson")

MEK <- subset(rep, rep$moa == "MEK inhibitor")
MEK_subdf <- subset(comp_df, rownames(comp_df) %in% MEK$new_iname)
colnames(MEK_subdf)
ggscatter(data = MEK_subdf, x = "pdx_logFC", y = "patient_logFC", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson")


combo <- subset(rep, rep$moa %in% c("Aurora kinase inhbitor", "bromodomain inhibitor"))
combo
combo_subdf <- subset(comp_df, rownames(comp_df) %in% combo$new_iname)
colnames(combo_subdf)
ggscatter(data = combo_subdf, x = "pdx_logFC", y = "patient_logFC", add = "reg.line", conf.int = T, cor.coef = T, cor.method = "pearson")

################################################################################
#
################################################################################

# Which compound classes show the highest correlations between ISO and PDX?

scatterList <- list()
corList <- list()
corDF <- data.frame(row.names = unique(rep$moa))
for (i in 1:length(unique(rep$moa))){
  # print(unique(rep$moa)[i])
  rep2 <- subset(rep, rep$moa == unique(rep$moa)[i])
  subdf <- subset(comp_df, rownames(comp_df) %in% rep2$new_iname)
  ggsc <- ggscatter(data = subdf, 
                    x = "pdx_logFC", 
                    y = "patient_logFC", 
                    add = "reg.line", 
                    conf.int = T, 
                    cor.coef = T, 
                    cor.method = "spearman", 
                    title = unique(rep$moa)[i], color = c("darkblue"), size = 4 
                    # label = rownames(subdf)
                    ) + theme(text = element_text(size = 32))
  scatterList[[i]] <- ggsc
  corRes <- cor(x = subdf$patient_logFC, y = subdf$pdx_logFC, method = "spearman")
  # corRes2 <- cor.test(subdf$patient_logFC, y = subdf$pdx_logFC, method = "spearman", use = "complete.obs")
  # rownames(corDF)[i] <- rep$moa[i]
  corDF$moa[i] <- unique(rep$moa)[i]
  corDF$res[i] <- corRes
}

topMOA <- subset(corDF, corDF$res > 0.5)
names(scatterList) <- unique(rep$moa)
rownames(topMOA)
topScatterList <- scatterList[rownames(topMOA)]
library(cowplot)

filteredClasses <- rownames(topMOA)[c(3, 7, 14, 22, 38, 41, 43, 74)]
filteredClasses
filteredScatterList <- scatterList[filteredClasses]

pdf(file = "moa_scatters_spearman.pdf", height = 100, width = 7)
plot_grid(plotlist = topScatterList, ncol = 2)
dev.off()

pdf(file = "moa_scatters_spearman_filtered.pdf", width = 56, height = 28)
plot_grid(plotlist = filteredScatterList, ncol = 4, labels = "AUTO", label_size = 32)
dev.off()






scatterList <- list()
corList <- list()
corDF <- data.frame(row.names = unique(rep$moa))
for (i in 1:length(unique(rep$moa))){
  print(unique(rep$moa)[i])
  rep2 <- subset(rep, rep$moa == unique(rep$moa)[i])
  subdf <- subset(comp_df, rownames(comp_df) %in% rep2$new_iname)
  ggsc <- ggscatter(data = subdf, 
                    x = "pdx_logFC", 
                    y = "patient_logFC", 
                    add = "reg.line", 
                    conf.int = T, 
                    cor.coef = T, 
                    cor.method = "pearson", 
                    title = unique(rep$moa)[i] 
                    # label = rownames(subdf)
                    )
  scatterList[[i]] <- ggsc
  corRes <- cor(x = subdf$patient_logFC, y = subdf$pdx_logFC, method = "pearson")
  # rownames(corDF)[i] <- rep$moa[i]
  corDF$moa[i] <- unique(rep$moa)[i]
  corDF$res[i] <- corRes
}

topMOA <- subset(corDF, corDF$res > 0.5)
names(scatterList) <- unique(rep$moa)
topScatterList <- scatterList[rownames(topMOA)]
library(cowplot)

pdf(file = "moa_scatters_pearson.pdf", height = 100, width = 7)
plot_grid(plotlist = topScatterList, ncol = 2)
dev.off()

# res1 <- resdf[resdf$contrast == patientContrasts[4],]
# res1
# 
# 
# EnhancedVolcano(toptable = resdf, lab = rownames(resdf), x = "logFC", y = "adj.P.Val", 
#                 xlim = c(min(resdf[["logFC"]], na.rm = TRUE) - .01, 
#                          max(resdf[["logFC"]], na.rm = TRUE) + .01), FCcutoff = 0.01, pCutoff = 0.05)
# 
# EnhancedVolcano(toptable = res1, lab = rownames(res1), x = "logFC", y = "adj.P.Val", 
#                 xlim = c(min(res1[["logFC"]], na.rm = TRUE) - .01, 
#                          max(res1[["logFC"]], na.rm = TRUE) + .01), FCcutoff = 0.01, pCutoff = 5*10^-80)


# # contr <- makeContrasts(groupresistant.SM018 - groupsensitive.SM018, levels = colnames(coef(fit)))
# # contr
# 
# contr <- makeContrasts(contrasts = patientContrasts, levels = colnames(coef(fit)))
# contr
# 
# # estimate contrast for each drug (gene)
# tmp <- contrasts.fit(fit, contr)
# 
# tmp <- eBayes(tmp)
# 
# top.table <- topTable(tmp, sort.by = "P", n = Inf)
# head(top.table)
