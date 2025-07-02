library(Seurat)
library(ISOSCELES)
library(limma)
library(EnhancedVolcano)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(Hmisc)
library(dittoSeq)

ctresp <- read.csv(file = "/data/ISOupdates/lfc_t_sigres.csv")

corr_mat <- read.csv(file = "~/Downloads/20250225_RDS_upload_L1000_consensus_corrmat.csv", row.names = 1)
head(rownames(corr_mat))

rephub <- read.delim(file = "/data/GBM_CombinationScoring_ISOSCELES_RV/repurposing_drugs_20200324.txt",
                     row.names = 1, skip = 9, header = T)
rephub$pert_iname <- rownames(rephub)
unique(rephub$clinical_phase)
clinical_rephub <- subset(rephub, rephub$clinical_phase %nin% c("Preclinical", "Withdrawn"))# & rephub$disease_area == "oncology")
dim(clinical_rephub)
clinical_rephub$new_iname <- gsub(x = clinical_rephub$pert_iname,
                                  pattern = "-", replacement = ".",
                                  fixed = T)


sensitivity <- read.csv(file = "~/Downloads/CT179_sensitivityByTCS_2025-03-02.csv", row.names = 1)
resistantCells <- rownames(sensitivity)[which(sensitivity$sensitivity == "resistant")]
sensitiveCells <- rownames(sensitivity)[which(sensitivity$sensitivity == "sensitive")]
length(c(resistantCells, sensitiveCells))

ct179_connectivities <- read.csv(file = "~/Downloads/CT179_connectivities_2025-03-02.csv", row.names = 1)

RDSseurat <- readRDS(file = "/data/GBM_CombinationScoring_ISOSCELES_RV/isosceles_suterJohnsonMerge_fourStates_v5.RDS")
RDSseurat <- AddMetaData(RDSseurat, metadata = ct179_connectivities)

    RDSseurat <- SetIdent(RDSseurat, value = "Neftel_State")
    corr_mat_2 <- as.data.frame(t(corr_mat))
    tumorCells <- WhichCells(RDSseurat, idents = c("NPC", "OPC", "AC", "MES"))
    tumorCells <- gsub(x = tumorCells, pattern = "-", replacement = ".", fixed = T)
    compoundSpearmans <- RDSseurat@meta.data["CT179"]
    dim(compoundSpearmans)
    compoundSpearmans <- subset(compoundSpearmans, rownames(compoundSpearmans) %in% tumorCells)
    dim(compoundSpearmans)
    corr_mat_3 <- corr_mat_2[which(rownames(corr_mat_2) %in% tumorCells),which(colnames(corr_mat_2) %in% clinical_rephub$new_iname)]
    z_matrix <- 0.5 * log((1 + corr_mat_3) / (1 - corr_mat_3)) # apply Fisher z-transformation

    RDSseurat <- AddMetaData(RDSseurat, metadata = sensitivity)
    np <- readRDS(file = "/data/ISOupdates/isoscles_suterJohnsonMerge_fourStates_np.RDS")
    np <- AddMetaData(np, metadata = sensitivity)

    # set up groupings and construct design matrix...
    snames <- colnames(corr_mat_3)

    subjectID <- data.frame(subjectID = RDSseurat@meta.data[,"patientID"],
                            row.names = rownames(RDSseurat@meta.data))
    rownames(subjectID) <- gsub(x = rownames(subjectID),
                                pattern = "-", replacement = ".",
                                fixed = TRUE)
    subjectID <- subjectID[which(rownames(subjectID) %in% rownames(corr_mat_3)),]

    rownames(sensitivity) <- gsub(x = rownames(sensitivity),
                                  pattern = "-", replacement = ".",
                                  fixed = TRUE)
    length(sensitivity)
    sensitivity <- sensitivity[which(rownames(sensitivity) %in% rownames(corr_mat_3)),]
    # subject_group <- interaction(sensitivity$sensitivity, subjectID$subjectID)
    sensitivity <- factor(sensitivity)
    subjectID <- factor(subjectID)
    design <- model.matrix(~ subjectID + sensitivity)
    z_matrix <- t(z_matrix)
    fit <- lmFit(z_matrix, design)
    fit <- eBayes(fit)
    results <- topTable(fit, coef = 2, adjust.method = "fdr", number = Inf)
    print(head(results))
    print(head(colnames(results)))
    print(head(colnames(z_matrix)))
    max(-log10(results[["adj.P.Val"]]))

    Idents(RDSseurat) <- RDSseurat$sensitivity
    sens_bp <- dittoBarPlot(RDSseurat, var = "sensitivity", group.by = "patientID",
                 cells.use = WhichCells(RDSseurat,
                                        idents = c("resistant","sensitive")),
                 color.panel = c("darkred", "darkblue"), main = "CT-179 Predicted Sensitivity")

    # dittoDimPlot(RDSseurat, var = "CT179", order = "randomize")
    hier_sen <- dittoDimPlot(np, var = "sensitivity", order = "randomize",
                 reduction = "hierarchy",
                 color.panel = c("darkred", "darkblue"), size = 0.5,
                 main = "CT-179 Predicted Sensitivity", legend.show = F)

    RDSseurat$Neftel_State <- factor(RDSseurat$Neftel_State,
                                         levels = c("AC", "MES", "NPC", "OPC", "Non-Neoplastic"))
    connectbox <- dittoPlot(RDSseurat, var = "CT179", group.by = "Neftel_State",
              jitter.size = 0.5, plots = c("vlnplot", "boxplot"), vlnplot.width = 1.2, boxplot.width = 0.8,
              ylab = c("CT-179 Connectivity"),
              color.panel = c("forestgreen", "red","blue", "purple", "gray"),
              theme = theme_minimal())

    ############################################################################

    volc <- EnhancedVolcano(toptable = results, # title = paste0(input$Custom_TCS_nameInput, " induced alterations"),
                            # subtitle = "fisher's Z transformation of raw connectivity values + limma",
                            lab = rownames(results), selectLab = c("ARRY.334543", "indibulin", "docetaxel"), drawConnectors = T,
                            x = "logFC",   xlab = bquote("Connectivity" ~Log[2] ~ "fold change"),
                            y = "adj.P.Val", FCcutoff = 0.01,
                            xlim = c(min(results[["logFC"]], na.rm = TRUE) - 0.1, max(results[["logFC"]], na.rm = TRUE) + 0.1),
                            # ylim = c(0, max(-log10(results[["adj.P.Val"]]), na.rm = TRUE))
                            )
    volc_facet <- volc + ggforce::facet_zoom(xlim = c(min(results[["logFC"]], na.rm = TRUE) - 0.05,0))
    volc_facet <- volc_facet + theme(panel.spacing.y = unit(0.1, "cm"))

    pdf(file = "/data/ISOupdates/full_CT179_connectivityLimmaVolc.pdf", height = 10.5)
    volc_facet
    dev.off()
    
    dim(results) # 1382
    dim(results[which(results$logFC < 0 & results$adj.P.Val < 0.05),]) # 721



################################################################################

################################################################################

        resistMat <- z_matrix[,which(colnames(z_matrix) %in% resistantCells)]
        # resistMat(resistMat)
        mrc <- data.frame(mrc = rowMeans(resistMat), row.names = rownames(resistMat))
        mrc$Group <- rownames(mrc)

        mrc

mrc <- mrc %>%
  arrange(mrc) %>%
  mutate(Group = factor(Group, levels = Group),
         Group_num = row_number())  # creates a numeric index 1, 2, ..., n

dim(mrc)
dim(mrc[which(mrc$mrc < 0),])

sensitized <- results[which(results$logFC < 0 & results$adj.P.Val < 0.05),]
discordant <- mrc[which(mrc$mrc < 0),]

head(sensitized)
head(discordant)

length(intersect(rownames(sensitized), rownames(discordant)))
# Create the bar plot using the numeric index for the x-axis.
g <- ggbarplot(mrc, x = "Group_num", y = "mrc", fill = "mrc", color = "mrc",
               ggtheme = theme_minimal_grid(),
               title = "Mean connectivities to CT-179 resistant cells") +
  scale_fill_viridis_c(option = "turbo", direction = -1) +
  scale_color_viridis_c(option = "turbo", direction = -1) +
  # Map the numeric index back to the original group labels on the x-axis.
  scale_x_continuous(breaks = mrc$Group_num, labels = levels(mrc$Group)) +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x = element_blank()) +  # Remove x-axis labels globally
  # Zoom in on the first 30 groups (using the numeric index)
  ggforce::facet_zoom(x = Group_num <= 824, ylim = c(min(mrc$mrc),0), zoom.size = 3)

# Optionally, add annotations for the zoomed-in area.
anndat <- subset(mrc, Group %in% c("ARRY.334543", "indibulin", "docetaxel"))
g2 <- g + geom_text_repel(data = anndat, force = 5, max.time = 10, nudge_x = 20, nudge_y = 0.005,
                    aes(x = Group_num, y = mrc, label = Group),
                    color = "black", size = 12, max.overlaps = 30)

pdf(file = "/data/ISOupdates/mrctest.pdf")
g2
dev.off()

################################################################################
#                                                                              #
################################################################################

mer <- merge(results, mrc, by.x = "row.names", by.y = "Group")
mer <- subset(mer, mer$logFC < 0 & mer$mrc < 0)
mer$combination <- mer$logFC * mer$mrc
mer$Group <- rownames(mer)
mer <- mer %>%
  arrange(-combination) %>%
  mutate(Group = factor(Group, levels = Group),
         Group_num2 = row_number())  # creates a numeric index 1, 2, ..., n

gs <- ggscatter(mer, x = "mrc", y = "logFC", size = 4,star.plot = T,
                color = "combination", rug = T,  # font.label = c(14, "bold", "black"),
                xlab = "Mean CT-179 Resistant Cell Connectivity",
                ylab = "CT-179 Resistant vs Sensitive Connectivity logFC") +
  scale_color_viridis_c(option = "turbo", direction = 1, begin = 0.5) +
  scale_x_continuous(position = "top")#  + coord_fixed()
  anndat <- subset(mer, Group_num2 <= 30)
  anndat2 <- subset(mer, Row.names %in% c("ARRY.334543", "indibulin", "docetaxel"))
  gs2 <- gs + geom_text_repel(data = anndat2, force = 10, force_pull = 0, max.time = 10, nudge_x = 0.005, nudge_y = 0.005,
                              fontface = "bold",
                              arrow = arrow(angle = 30, length = unit(0.01, "inches"), ends = "last", type = "open"),
                              aes(x = mrc, y = logFC, label = Row.names),
                              color = "black", size = 12, max.overlaps = 30)
pdf(file = "/data/ISOupdates/CT-179_connectivity_scatter.pdf")
gs2
dev.off()

gr1 <- plot_grid(plotlist = list(sens_bp, hier_sen, connectbox), ncol = 3, labels = "AUTO")
# gr1 <- plot_grid(plotlist = list(sens_bp, volc_facet, g2), ncol = 3, labels = "AUTO")
gr2 <- plot_grid(plotlist = list(volc_facet, g2), ncol = 2, labels = c("C", "D"), rel_widths = c(1,2))
gr3 <- plot_grid(plotlist = list(gs2), ncol = 3, labels = c("E", "F", "G"))
gr_full <- plot_grid(plotlist = list(gr1, gr2, gr3), nrow = 3, rel_heights = c(1,1,1), labels = c("", "D"))

pdf("/data/ISOupdates/ct179_scoring_figures.pdf", height = 20, width = 20)
gr_full
dev.off()


pdf(file = "/data/ISOupdates/EGFR_OLIG2_vln.pdf", width = 14)
dittoPlot(object = RDSseurat, var = c("EGFR", "OLIG2"),
          group.by = "Neftel_State",
          slot = "scale.data",
          plots = c("vlnplot", "boxplot"), boxplot.width = 0.8)
dev.off()
