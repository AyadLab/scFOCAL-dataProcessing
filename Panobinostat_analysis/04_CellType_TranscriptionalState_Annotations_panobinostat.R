library(ggplot2)
library(cowplot)
library(GEOquery)
library(Seurat)
library(Hmisc)
library(dplyr)
library(scrabble)
library(singscore)
library(viridis)
library(GSEABase)

set.seed(1)

obj <- readRDS(file = "zhao_integrated_TCSadded.RDS")

################################################################################
# Separate out discrete cell types
################################################################################

DefaultAssay(obj) <- "integrated"
obj <- RunPCA(obj, verbose = FALSE)
ElbowPlot(obj)
# obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.25)
obj <- RunUMAP(obj, dims = 1:15, n.neighbors = 50, min.dist = 0.5, reduction.name = "umap_int")
treatmentUMAP <-DimPlot(obj, reduction = "umap_int", 
                        group.by = c("treatment"), 
                        cols = viridis(2, option = "turbo", begin = 0.1, end = 0.6, alpha = 0.5), 
                        shuffle = T, pt.size = 0.25) + coord_fixed() + theme_void()
pdf(file = "panobinostat_out/panobinostat_treatment_UMAP.pdf")
treatmentUMAP
dev.off()


# DimPlot(obj, reduction = "umap_int", 
#         group.by = c("integrated_snn_res.0.25"), 
#         # cols = viridis(18, option = "turbo", begin = 0.1, end = 0.6, alpha = 0.5), 
#         shuffle = T, pt.size = 0.25, label = T) + coord_fixed() + theme_void()
Idents(obj) <- obj$integrated_snn_res.0.25
# cl11_markers <- FindMarkers(obj, ident.1 = "11", only.pos = T)
# cll11_markers %>% top_n(10, avg_log2FC)
# DotPlot(obj, features = c("EGFR", "MBP", "PLP1", "CD74", "CD3E", "PTPRC", "PDGFRA", "PDGFRB"), col.min = 0, group.by = "integrated_snn_res.0.2")
# FeaturePlot(obj, features = c("EGFR", "PDGFRB", "CD3E", "PTPRC", "CD74", "MBP", "PLP1"), min.cutoff = 0)
# FeaturePlot(obj, features = c("PLP1", "CD14", "TRAC"))
# VlnPlot(obj, group.by = "seurat_clusters", features = c("PDGFRA", "PDGFRB", "CD3E", "PTPRC", "CD74", "MBP", "PLP1"))

Idents(obj) <- obj$integrated_snn_res.0.25
obj <- RenameIdents(obj, '0' = 'Neoplastic',
                    '1' = 'Oligodendrocytes',
                    '2' = 'Myeloid',
                    '3' = 'Neoplastic',
                    '4' = 'Neoplastic',
                    '5' = 'Neoplastic',
                    '6' = 'Myeloid',
                    '7' = 'Neoplastic',
                    '8' = 'T-Cells',
                    '9' = 'Neoplastic',
                    '10' = 'Myeloid',
                    '11' = 'Pericytes', 
                    '12' = "Myeloid")

obj$cellType <- Idents(obj)
# DotPlot(obj, features = c("EGFR", "MBP", "PLP1", "CD74", "CD3E", "PTPRC", "CLDN5"),
#         col.min = 0, group.by = "cellType")
cellTypeUMAP <- DimPlot(obj, reduction = "umap_int", 
                        group.by = c("cellType"), 
                        # cols = viridis(18, option = "turbo", begin = 0.1, end = 0.6, alpha = 0.5), 
                        shuffle = T, pt.size = 0.25) + coord_fixed() + theme_void()

pdf(file = "panobinostat_out/panobinostat_cellType_umap.pdf")
cellTypeUMAP
dev.off()


# SNNmarkers <- FindAllMarkers(obj, test.use = "MAST", )

################################################################################
#
################################################################################

# Annotate tumor and TME cell types
Idents(obj) <- obj$cellType 
NeoplasticCells <- WhichCells(obj, idents = c("Neoplastic"))
MyeloidCells <- WhichCells(obj, idents = c("Myeloid"))
OligodendrocyteCells <- WhichCells(obj, idents = c("Oligodendrocytes"))

SingscoreExprs <- GetAssayData(obj[["RNA"]], slot = "scale.data")
np <- SingscoreExprs[,NeoplasticCells]

### from patient data preprocessing...

# load in cell type signatures file
suva.sigs <- read.csv(file = 'SuvaMetaSignatures.csv', skip = 4)

# set up ranked data based on scRNA signatures

# Idents(obj.int) <- obj.int$ClassificationNPvsNon
# obj.np <- subset(obj.int, idents = c("Neoplastic"))
# DefaultAssay(obj.np) <- "RNA"
# obj.seurat.data <- GetAssayData(obj.np, assay = "RNA", slot = "scale.data")
np <- as.data.frame(np)
rank.data <- singscore::rankGenes(np)

# create singscore for each cell state
for (i in 1:length(colnames(suva.sigs))){
  # capture unique gene signatures (removes duplicate cells)
  uni <- unique(suva.sigs)[i]
  
  # capture name of cell state
  state <- colnames(suva.sigs)[i]
  
  # removes empty cells
  uni <- subset(uni, uni[,1] != "")
  uni <- as.vector(uni[,1])
  
  # scoring with singscore
  set <- GeneSet()
  set@geneIds <- as.character(uni)
  scored <- singscore::simpleScore(rankData = rank.data, upSet = set)
  
  # renaming "$Sig" to appropriate cell state
  scored$Sig <- as.character(state)
  
  assign(paste0("scored.", state), scored)
}

scored.AC
scored.MES1
# append singscores to seurat object
singscores.list <- c("scored.AC", "scored.G1.S", "scored.G2.M", "scored.MES1", "scored.MES2", "scored.NPC1", "scored.NPC2", "scored.OPC")

for (n in singscores.list){
  state <- unlist(strsplit(as.character(n), split='.', fixed = TRUE))[2]
  scored <- eval(parse(text = n))
  toAdd <- as.data.frame(scored$TotalScore)
  rownames(toAdd) <- rownames(scored)
  obj <- AddMetaData(obj, metadata = toAdd, col.name = paste0(state, "singscore"))
}

FeaturePlot(obj, features = c("ACsingscore",
                              "MES1singscore",
                              "MES2singscore",
                              "NPC1singscore",
                              "NPC2singscore",
                              "OPCsingscore",
                              "G1singscore",
                              "G2singscore"), ncol = 4, reduction = "umap_int", order = T, cols = viridis(256, option = "D"))

# pull enrichment scores for each state out of the seurat object to a new data frame
singscoredf <- data.frame(row.names = rownames(obj@meta.data),
                          NPC1 = obj@meta.data$NPC1singscore,
                          NPC2 = obj@meta.data$NPC2singscore,
                          OPC = obj@meta.data$OPCsingscore,
                          AC = obj@meta.data$ACsingscore,
                          MES1 = obj@meta.data$MES1singscore,
                          MES2 = obj@meta.data$MES2singscore
)

singscoredf <- singscoredf[NeoplasticCells,]
### then loop through and just label each cell based on what signature is highest - here, you’ll want to add in the g1/s signatures you have also

singscoredf$assignment <- "Ambiguous"
for (i in 1:length(rownames(singscoredf))){
  if (singscoredf$NPC1[i] > singscoredf$NPC2[i]){
    if (singscoredf$NPC1[i] > singscoredf$OPC[i]){
      if (singscoredf$NPC1[i] > singscoredf$AC[i]){
        if (singscoredf$NPC1[i] > singscoredf$MES1[i]){
          if (singscoredf$NPC1[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "NPC1"
          }
        }
      }
    }
  }
  if (singscoredf$NPC2[i] > singscoredf$NPC1[i]){
    if (singscoredf$NPC2[i] > singscoredf$OPC[i]){
      if (singscoredf$NPC2[i] > singscoredf$AC[i]){
        if (singscoredf$NPC2[i] > singscoredf$MES1[i]){
          if (singscoredf$NPC2[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "NPC2"
          }
        }
      }
    }
  }
  if (singscoredf$OPC[i] > singscoredf$NPC1[i]){
    if (singscoredf$OPC[i] > singscoredf$NPC2[i]){
      if (singscoredf$OPC[i] > singscoredf$AC[i]){
        if (singscoredf$OPC[i] > singscoredf$MES1[i]){
          if (singscoredf$OPC[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "OPC"
          }
        }
      }
    }
  }
  if (singscoredf$AC[i] > singscoredf$NPC1[i]){
    if (singscoredf$AC[i] > singscoredf$NPC2[i]){
      if (singscoredf$AC[i] > singscoredf$OPC[i]){
        if (singscoredf$AC[i] > singscoredf$MES1[i]){
          if (singscoredf$AC[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "AC"
          }
        }
      }
    }
  }
  if (singscoredf$MES1[i] > singscoredf$NPC1[i]){
    if (singscoredf$MES1[i] > singscoredf$NPC2[i]){
      if (singscoredf$MES1[i] > singscoredf$OPC[i]){
        if (singscoredf$MES1[i] > singscoredf$AC[i]){
          if (singscoredf$MES1[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "MES1"
          }
        }
      }
    }
  }
  if (singscoredf$MES2[i] > singscoredf$NPC1[i]){
    if (singscoredf$MES2[i] > singscoredf$NPC2[i]){
      if (singscoredf$MES2[i] > singscoredf$OPC[i]){
        if (singscoredf$MES2[i] > singscoredf$AC[i]){
          if (singscoredf$MES2[i] > singscoredf$MES1[i]){
            singscoredf$assignment[i] <- "MES2"
          }
        }
      }
    }
  }
}

singscoredf
singscoredf$cell <- rownames(singscoredf)

# relabel substates to parent (i.e. Mes1 and Mes2 = MES) (4 total states)

for (i in 1:length(rownames(singscoredf))){
  print(singscoredf$assignment[i])
  if (singscoredf$assignment[i] %in% c("NPC1", "NPC2")){
    singscoredf$assignment2[i] <- "NPC"
  } else if (singscoredf$assignment[i] %in% c("MES1", "MES2")){
    singscoredf$assignment2[i] <- "MES"
  } else {
    singscoredf$assignment2[i] <- as.character(singscoredf$assignment[i])
  }
}

singscoredf

singscorebin <- singscoredf[,c("cell", "assignment", "assignment2")]
colnames(singscorebin) <- c("cell", "Neftel_State_6", "Neftel_State_4")
singscorebin$cell <- NULL
singscorebin
obj <- AddMetaData(obj, metadata = singscorebin)

################################################################################
# Aggregate to 4 signatures...

singscoredf$MES <- rowMeans(singscoredf[ , c("MES1", "MES2")], na.rm = T)
singscoredf$NPC <- rowMeans(singscoredf[ , c("NPC1", "NPC2")], na.rm = T)

inputdf <- singscoredf[c("AC", "MES", "NPC", "OPC")]
inputdf

obj <- AddMetaData(obj, metadata = inputdf)

hierarchy <- scrabble::hierarchy(m = inputdf, quadrants = NULL, log.scale = T)
pdf(file = "panobinostat_prelimHierarchy.pdf")
scrabble::plot_hierarchy(hierarchy)
dev.off()



################################################################################


saveRDS(obj, file = "zhao_integrated_TCSadded_cellTyped.RDS")
obj <- readRDS(file = "zhao_integrated_TCSadded_cellTyped.RDS")
stateUMAP <- DimPlot(obj, reduction = "umap_int", 
                     group.by = c("Neftel_State_4"), 
                     # cols = viridis(18, option = "turbo", begin = 0.1, end = 0.6, alpha = 0.5), 
                     shuffle = T, pt.size = 0.25) + coord_fixed() + theme_void()

# stateUMAP

# VlnPlot(obj, group.by = "treatment", features = "panobinostat", pt.size = 0)
np <- subset(obj, cells = NeoplasticCells)

# add hierarchy data to np object - 
hier <- as.matrix(hierarchy)
colnames(hier) <- c("hierarchy_1", "hierarchy_2")
np[["hierarchy"]] <- CreateDimReducObject(embeddings = hier, key = "hierarchy_", assay = "RNA")

np_hier <- DimPlot(np, reduction = "hierarchy", group.by = "Neftel_State_4") + 
  xlim(c(-0.6,0.6)) + ylim(c(-0.6,0.6)) +
  xlab(label = "Relative Meta-Module Score [log2(| SC1 - SC2 | + 1)]") + 
  ylab(label = "Relative Meta-Module Score [log2(| SC1 - SC2 | + 1)]") + 
  theme_minimal()

pdf(file = "panobinostat_hier_dimplot_01.pdf")
np_hier
dev.off()

np_treatment_hier <- DimPlot(np, reduction = "hierarchy", group.by = "treatment") + 
  xlim(c(-0.6,0.6)) + ylim(c(-0.6,0.6)) +
  xlab(label = "Relative Meta-Module Score [log2(| SC1 - SC2 | + 1)]") + 
  ylab(label = "Relative Meta-Module Score [log2(| SC1 - SC2 | + 1)]") + 
  theme_minimal()

pdf(file = "panobinostat_treatment_hier_dimplot_01.pdf")
np_treatment_hier
dev.off()

np_pano_vln <- VlnPlot(np, group.by = "treatment", features = "panobinostat", pt.size = 0, cols = viridis(2, option = "turbo", begin = 0.1, end = 0.6, alpha = 0.5)) + theme_minimal()
my <- subset(obj, cells = MyeloidCells)
my_pano_vln <- VlnPlot(my, group.by = "treatment", features = "panobinostat", pt.size = 0, cols = viridis(2, option = "turbo", begin = 0.1, end = 0.6, alpha = 0.5)) + theme_minimal()

# wilcox test...

unique(np$treatment)
pano_np_meta <- np@meta.data[which(np@meta.data$treatment == c("treatment: 0.2 uM panobinostat")),]
dmso_np_meta <- np@meta.data[which(np@meta.data$treatment == c("treatment: vehicle (DMSO)")),]

np_wil_res <- wilcox.test(dmso_np_meta$panobinostat, pano_np_meta$panobinostat, paired = F)

unique(my$treatment)
pano_my_meta <- my@meta.data[which(my@meta.data$treatment == c("treatment: 0.2 uM panobinostat")),]
dmso_my_meta <- my@meta.data[which(my@meta.data$treatment == c("treatment: vehicle (DMSO)")),]

my_wil_res <- wilcox.test(dmso_my_meta$panobinostat, pano_my_meta$panobinostat, paired = F)
my_wil_res


pdf(file = "panobinostat_out/panobinostat_tcs_neoplastic_myeloid.pdf", width = 14)
plot_grid(plotlist = list(np_pano_vln, my_pano_vln))
dev.off()

# dittoSeq::dittoBarPlot(obj, var = "singscorebin", group.by = "treatment", scale = "percent", cells.use = NeoplasticCells)

pano_tcs_UMAP <- FeaturePlot(obj, features = "panobinostat_inverseCor", cols = viridis(option = "turbo", begin = 0.1, end = 1, n = 256, alpha = 0.7), order = T, min.cutoff = 0, pt.size = 0.25) + coord_fixed() + theme_void()

pdf(file = "panobinostat_out/panobinostat_tcs_umap.pdf")
pano_tcs_UMAP
dev.off()
################################################################################
################################################################################
################################################################################

# Panobinostat TCS and tumor cell population changes...

np_meta <- np@meta.data

# Overall KDE Panobinostat TCS

pdf(file = "np_panobinostatTreatment_MPI_normKDE_02.pdf", height = 1.75)
ggplot(np_meta, aes(x = panobinostat, fill = treatment)) +
  stat_density(geom = "area", position = "fill", alpha = 0.5) +
  labs(
    title = "Area-Normalized KDE Plot",
    x = "Panobinostat TCS Connectivity",
    y = "Density"
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal()
dev.off()


pdf(file = "np_panobinostatTreatment_neftelStates_normKDE_01.pdf")
ggplot(np_meta, aes(x = ACsingscore, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Normalized KDE Plot",
    x = "AC enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
ggplot(np_meta, aes(x = MES1singscore, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Normalized KDE Plot",
    x = "MES1 enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
ggplot(np_meta, aes(x = MES2singscore, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Normalized KDE Plot",
    x = "MES2 enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
ggplot(np_meta, aes(x = NPC1singscore, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Normalized KDE Plot",
    x = "NPC1 enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
ggplot(np_meta, aes(x = NPC2singscore, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Normalized KDE Plot",
    x = "NPC2 enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
ggplot(np_meta, aes(x = OPCsingscore, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Normalized KDE Plot",
    x = "OPC enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
dev.off()

################################################################################
# Reformat using cowplot...

np_kde_plotlist <- list()
for (state in unique(np_meta$Neftel_State_4)){
  print(state)
  np_kde_plotlist[[state]] <- ggplot(np_meta, aes_string(x = state, color = "treatment", fill = "treatment")) +
    geom_density(alpha = 0.5) +
    labs(
      title = "Normalized KDE Plot",
      x = paste(state, "enrichment", sep = " "),
      y = "Density"
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal()
}

pdf(file = "np_kde_plotgrid.pdf", height = 1.75*4)
plot_grid(plotlist = np_kde_plotlist, ncol = 1)
dev.off()

## 4 states
pdf(file = "np_panobinostatTreatment_neftelStates_normKDE_02_fourStates.pdf")
ggplot(np_meta, aes(x = AC, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Normalized KDE Plot",
    x = "AC enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
ggplot(np_meta, aes(x = MES, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Normalized KDE Plot",
    x = "MES enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
ggplot(np_meta, aes(x = NPC, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Normalized KDE Plot",
    x = "NPC enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
ggplot(np_meta, aes(x = OPC, color = treatment, fill = treatment)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Normalized KDE Plot",
    x = "OPC enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
dev.off()

# normalized KDE...

# pdf(file = "macRes_panobinostatTreatment_TCS_normKDE_02.pdf", height = 1.75)
# ggplot(macres_tcs_merge, aes(x = panobinostat, fill = Feature)) +
#   stat_density(geom = "area", position = "fill", alpha = 0.5) +
#   labs(
#     title = "Area-Normalized KDE Plot",
#     x = "Panobinostat TCS Connectivity",
#     y = "Density"
#   ) +
#   scale_fill_brewer(palette = "Set1") +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme_minimal()
# dev.off()

pdf(file = "np_panobinostatTreatment_neftelStates_normKDE_02_fourStates.pdf")
ggplot(np_meta, aes(x = AC, color = treatment, fill = treatment)) +
  stat_density(geom = "area", position = "fill", alpha = 0.5) +
  labs(
    title = "Area-Normalized KDE Plot",
    x = "AC enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
ggplot(np_meta, aes(x = MES, color = treatment, fill = treatment)) +
  stat_density(geom = "area", position = "fill", alpha = 0.5) +
  labs(
    title = "Area-Normalized KDE Plot",
    x = "MES enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
ggplot(np_meta, aes(x = NPC, color = treatment, fill = treatment)) +
  stat_density(geom = "area", position = "fill", alpha = 0.5) +
  labs(
    title = "Area-Normalized KDE Plot",
    x = "NPC enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
ggplot(np_meta, aes(x = OPC, color = treatment, fill = treatment)) +
  stat_density(geom = "area", position = "fill", alpha = 0.5) +
  labs(
    title = "Area-Normalized KDE Plot",
    x = "OPC enrichment",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
#
dev.off()

# Reformat using cowplot...

np_areaNorm_kde_plotlist <- list()
for (state in unique(np_meta$Neftel_State_4)){
  print(state)
  np_areaNorm_kde_plotlist[[state]] <- ggplot(np_meta, aes_string(x = state, color = "treatment", fill = "treatment")) +
    stat_density(geom = "area", position = "fill", alpha = 0.5) +
    labs(
      title = "Area-Normalized KDE Plot",
      x = paste(state, "module enrichment", sep = " "),
      y = "Density"
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal()
}

pdf(file = "np_areaNorm_kde_plotgrid.pdf", height = 1.75*4)
plot_grid(plotlist = np_areaNorm_kde_plotlist, ncol = 1)
dev.off()

# Hierarchy Plotting of NP cells...
# grab code from other script

################################################################################
# Neoplastic predictions in treatment naive (DMSO) samples
################################################################################




Idents(np) <- np$treatment
Idents(np)
dmso_np <- subset(np, idents = c("treatment: vehicle (DMSO)"))
# need to rescale... re-score for state and connectivity...

dmso_np <- CellCycleScoring(dmso_np, g2m.features = cc.genes.updated.2019$g2m.genes, s.features = cc.genes.updated.2019$s.genes)

DefaultAssay(dmso_np) <- "RNA"
dmso_np <- ScaleData(dmso_np, features = rownames(dmso_np), vars.to.regress = c("G2M.Score", "S.Score"))
saveRDS(dmso_np, file = "panobinostat_dmso_np_ccreg.RDS")


# re-singscore...

dmso_np_SingscoreExprs <- GetAssayData(dmso_np[["RNA"]], slot = "scale.data")

### from patient data preprocessing...

# load in cell type signatures file
suva.sigs <- read.csv(file = 'SuvaMetaSignatures.csv', skip = 4)

# set up ranked data based on scRNA signatures

# Idents(dmso_np.int) <- dmso_np.int$ClassificationNPvsNon
# dmso_np.np <- subset(dmso_np.int, idents = c("Neoplastic"))
# DefaultAssay(dmso_np.np) <- "RNA"
# dmso_np.seurat.data <- GetAssayData(dmso_np.np, assay = "RNA", slot = "scale.data")
dmso_np_SingscoreExprs <- as.data.frame(dmso_np_SingscoreExprs)
rank.data <- singscore::rankGenes(dmso_np_SingscoreExprs)

# create singscore for each cell state
for (i in 1:length(colnames(suva.sigs))){
  # capture unique gene signatures (removes duplicate cells)
  uni <- unique(suva.sigs)[i]
  
  # capture name of cell state
  state <- colnames(suva.sigs)[i]
  
  # removes empty cells
  uni <- subset(uni, uni[,1] != "")
  uni <- as.vector(uni[,1])
  
  # scoring with singscore
  set <- GeneSet()
  set@geneIds <- as.character(uni)
  scored <- singscore::simpleScore(rankData = rank.data, upSet = set)
  
  # renaming "$Sig" to appropriate cell state
  scored$Sig <- as.character(state)
  
  assign(paste0("scored.", state), scored)
}

# append singscores to seurat dmso_npect
singscores.list <- c("scored.AC", "scored.G1.S", "scored.G2.M", "scored.MES1", "scored.MES2", "scored.NPC1", "scored.NPC2", "scored.OPC")

for (n in singscores.list){
  state <- unlist(strsplit(as.character(n), split='.', fixed = TRUE))[2]
  scored <- eval(parse(text = n))
  toAdd <- as.data.frame(scored$TotalScore)
  rownames(toAdd) <- rownames(scored)
  dmso_np <- AddMetaData(dmso_np, metadata = toAdd, col.name = paste0(state, "singscore_dmso"))
}

# FeaturePlot(dmso_np, features = c("ACsingscore",
#                               "MES1singscore",
#                               "MES2singscore",
#                               "NPC1singscore",
#                               "NPC2singscore",
#                               "OPCsingscore",
#                               "G1singscore",
#                               "G2singscore"), ncol = 4, reduction = "umap_int", order = T, cols = viridis(256, option = "D"))

# pull enrichment scores for each state out of the seurat dmso_npect to a new data frame
singscoredf <- data.frame(row.names = rownames(dmso_np@meta.data),
                          NPC1 = dmso_np@meta.data$NPC1singscore_dmso,
                          NPC2 = dmso_np@meta.data$NPC2singscore_dmso,
                          OPC = dmso_np@meta.data$OPCsingscore_dmso,
                          AC = dmso_np@meta.data$ACsingscore_dmso,
                          MES1 = dmso_np@meta.data$MES1singscore_dmso,
                          MES2 = dmso_np@meta.data$MES2singscore_dmso
)

# singscoredf <- singscoredf[NeoplasticCells,]
### then loop through and just label each cell based on what signature is highest - here, you’ll want to add in the g1/s signatures you have also

singscoredf$assignment <- "Ambiguous"
for (i in 1:length(rownames(singscoredf))){
  if (singscoredf$NPC1[i] > singscoredf$NPC2[i]){
    if (singscoredf$NPC1[i] > singscoredf$OPC[i]){
      if (singscoredf$NPC1[i] > singscoredf$AC[i]){
        if (singscoredf$NPC1[i] > singscoredf$MES1[i]){
          if (singscoredf$NPC1[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "NPC1"
          }
        }
      }
    }
  }
  if (singscoredf$NPC2[i] > singscoredf$NPC1[i]){
    if (singscoredf$NPC2[i] > singscoredf$OPC[i]){
      if (singscoredf$NPC2[i] > singscoredf$AC[i]){
        if (singscoredf$NPC2[i] > singscoredf$MES1[i]){
          if (singscoredf$NPC2[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "NPC2"
          }
        }
      }
    }
  }
  if (singscoredf$OPC[i] > singscoredf$NPC1[i]){
    if (singscoredf$OPC[i] > singscoredf$NPC2[i]){
      if (singscoredf$OPC[i] > singscoredf$AC[i]){
        if (singscoredf$OPC[i] > singscoredf$MES1[i]){
          if (singscoredf$OPC[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "OPC"
          }
        }
      }
    }
  }
  if (singscoredf$AC[i] > singscoredf$NPC1[i]){
    if (singscoredf$AC[i] > singscoredf$NPC2[i]){
      if (singscoredf$AC[i] > singscoredf$OPC[i]){
        if (singscoredf$AC[i] > singscoredf$MES1[i]){
          if (singscoredf$AC[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "AC"
          }
        }
      }
    }
  }
  if (singscoredf$MES1[i] > singscoredf$NPC1[i]){
    if (singscoredf$MES1[i] > singscoredf$NPC2[i]){
      if (singscoredf$MES1[i] > singscoredf$OPC[i]){
        if (singscoredf$MES1[i] > singscoredf$AC[i]){
          if (singscoredf$MES1[i] > singscoredf$MES2[i]){
            singscoredf$assignment[i] <- "MES1"
          }
        }
      }
    }
  }
  if (singscoredf$MES2[i] > singscoredf$NPC1[i]){
    if (singscoredf$MES2[i] > singscoredf$NPC2[i]){
      if (singscoredf$MES2[i] > singscoredf$OPC[i]){
        if (singscoredf$MES2[i] > singscoredf$AC[i]){
          if (singscoredf$MES2[i] > singscoredf$MES1[i]){
            singscoredf$assignment[i] <- "MES2"
          }
        }
      }
    }
  }
}

singscoredf
singscoredf$cell <- rownames(singscoredf)

# relabel substates to parent (i.e. Mes1 and Mes2 = MES) (4 total states)

for (i in 1:length(rownames(singscoredf))){
  print(singscoredf$assignment[i])
  if (singscoredf$assignment[i] %in% c("NPC1", "NPC2")){
    singscoredf$assignment2[i] <- "NPC"
  } else if (singscoredf$assignment[i] %in% c("MES1", "MES2")){
    singscoredf$assignment2[i] <- "MES"
  } else {
    singscoredf$assignment2[i] <- as.character(singscoredf$assignment[i])
  }
}

singscoredf

singscorebin <- singscoredf[,c("cell", "assignment", "assignment2")]
colnames(singscorebin) <- c("cell", "Neftel_State_6", "Neftel_State_4")
singscorebin$cell <- NULL
singscorebin
dmso_np <- AddMetaData(dmso_np, metadata = singscorebin)

################################################################################
# Aggregate to 4 signatures...

singscoredf$MES <- rowMeans(singscoredf[ , c("MES1", "MES2")], na.rm = T)
singscoredf$NPC <- rowMeans(singscoredf[ , c("NPC1", "NPC2")], na.rm = T)
singscoredf

inputdf <- singscoredf[c("AC", "MES", "NPC", "OPC")]
inputdf

dmso_np <- AddMetaData(dmso_np, metadata = inputdf)
dmso_np <- AddMetaData(dmso_np, metadata = singscoredf)
################################################################################
# Recalculate panobinostat TCS
################################################################################

tcs_df <- read.csv(file = "TCS_out/phase2_compound_TCS_v0.1_20230914_112256.csv", row.names = 1)

# # regress cell cycle and scale
# DefaultAssay(obj) <- "RNA"
# obj <- ScaleData(obj, features = rownames(obj))

total.transpose <- t(dmso_np@assays$RNA@scale.data) # which assay to integrate with?

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
dmso_np <- AddMetaData(dmso_np, metadata = sc_df)
dmso_np$panobinostat_inverseCor <- dmso_np$panobinostat*-1


################################################################################
#
################################################################################

dmso_np_meta <- dmso_np@meta.data
dmso_np_meta <- dmso_np_meta %>%
  mutate(panobinostat_predicted_sensitivity = ifelse(panobinostat > mean(dmso_np_meta$panobinostat), "Resistant", "Sensitive"))

dmso_np_areaNorm_kde_plotlist <- list()
for (state in unique(dmso_np_meta$Neftel_State_6)){
  print(state)
  dmso_np_areaNorm_kde_plotlist[[state]] <- ggplot(dmso_np_meta, 
    aes_string(x = state, 
               color = "panobinostat_predicted_sensitivity", 
               fill = "panobinostat_predicted_sensitivity")) +
    stat_density(geom = "area", position = "fill", alpha = 0.5) +
    labs(
      title = "Area-Normalized KDE Plot",
      x = paste(state, "module enrichment", sep = " "),
      y = "Density"
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal()
}

pdf(file = "dmso_np_predictions_areaNorm_kde_plotgrid.pdf", height = 1.75*6)
plot_grid(plotlist = dmso_np_areaNorm_kde_plotlist, ncol = 1)
dev.off()

pdf(file = "dmso_panobinostat_TCS_by_state.pdf")
VlnPlot(dmso_np, group.by = "patientID", split.by = "Neftel_State_4", features = c("panobinostat"))
dev.off()

################################################################################
#
################################################################################

# olig <- SingscoreExprs[,OligodendrocyteCells]


# rescale in context of myeloid cells only

Idents(obj) <- obj$cellType
myeloid <- subset(obj, idents = c("Myeloid"))
myeloid <- ScaleData(myeloid, features = rownames(myeloid))

my <- GetAssayData(myeloid[["RNA"]], slot = "scale.data")
my <- as.data.frame(my)


toConvert <- as.data.frame(rownames(my))
colnames(toConvert) <- c("gene")
rawFileForENSEMBLids <- read.delim(file = "set_up_counts/cts_files/GSM4483741_PW029-701.cts.txt", row.names = 1)
conversionTable <- as.data.frame(rawFileForENSEMBLids$gene)
colnames(conversionTable) <- c("gene")
rownames(conversionTable) <- rownames(rawFileForENSEMBLids)
conversionTable$ensembl <- rownames(conversionTable)
conversionTable_ordered <- conversionTable[match(toConvert$gene, conversionTable$gene),]
my_merge <- merge(my, conversionTable, by.x = "row.names", by.y = "gene")
rownames(my_merge) <- my_merge$ensembl
my_merge$Row.names <- NULL
my_merge$geneid <- rownames(my_merge)

# prepare for macSpectrum...
# Column to be moved to the first position
column_to_move <- my_merge$geneid

# Create a new dataframe with the selected column as the first column
macSpec_inputMtx <- data.frame(geneid = column_to_move)

# Add the remaining columns from the original dataframe to the new dataframe
macSpec_inputMtx <- cbind(macSpec_inputMtx, my_merge[, -which(names(my_merge) == "geneid")])
macSpec_inputMtx$ensembl <- NULL

class(macSpec_inputMtx)
tail(colnames(macSpec_inputMtx))

library(macSpectrum)
# toConvert <- as.data.frame(rownames(macSpec_inputMtx))
# for (i in 1:length(rownames(toConvert))){
#   toConvert$newRownames[i] <- unlist(strsplit(toConvert$`rownames(macSpec_inputMtx)`, split = ".", fixed = T))[1]
# }
# rownames(macSpec_inputMtx) <- toConvert$newRownames

row_names <- rownames(macSpec_inputMtx)
row_names_clean <- sub("\\..*$", "", row_names)
rownames(macSpec_inputMtx) <- row_names_clean
macSpec_inputMtx$geneid <- rownames(macSpec_inputMtx)
macres <- macspec(mac_mtx = macSpec_inputMtx, feature = myeloid$treatment, select_hu_mo = "hum")

# is panobinostat changing macrophage polarization by macSpectrum?
colnames(macres)
q <- quantile(macres$MPI, 0.95)
pdf(file = "macRes_panobinostatTreatment_MPI_violin_01.pdf")
ggplot(macres, aes(x = Feature, y = MPI)) +
  geom_violin() +
  geom_boxplot(width = 1, outlier.shape = NA, coef = 0) + 
  geom_point(data = subset(macres, MPI < q),
             aes(x = Feature, y = MPI), color = "red", size = 3, 
             position = position_jitter(width = 0.2)) + 
  labs(
    title = "Violin Plot Example",
    x = "Treatment",
    y = "macSpectrum MPI"
  ) +
  theme_minimal()
dev.off()

# normalized KDE plot...

pdf(file = "macRes_panobinostatTreatment_MPI_normKDE_01.pdf")
ggplot(macres, aes(x = MPI, color = Feature, fill = Feature)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Normalized KDE Plot",
    x = "Value",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()
dev.off()



## what about predix in DMSO sample using panobinostat TCS?

myeloid_meta <- myeloid@meta.data
myeloid_meta
macres

macres_tcs_merge <- merge(macres, myeloid_meta, by.x = "Samples", by.y = "row.names")
macres_tcs_merge

pdf(file = "macRes_panobinostatTreatment_TCS_normKDE_02.pdf", height = 1.75)
ggplot(macres_tcs_merge, aes(x = panobinostat, fill = Feature)) +
  stat_density(geom = "area", position = "fill", alpha = 0.5) +
  labs(
    title = "Area-Normalized KDE Plot",
    x = "Panobinostat TCS Connectivity",
    y = "Density"
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal()
dev.off()

################################################################################
# What does Panobinostat TCS predict in Panobinostat-Naive slice culture?
################################################################################

Idents(myeloid) <- myeloid$treatment
Idents(myeloid)
dmso_myeloid <- subset(myeloid, idents = c("treatment: vehicle (DMSO)"))

################################################################################
# Need to re-run MacSpectrum on treatment naive myeloid population...
################################################################################

dmso_macSpec_inputMtx <- macSpec_inputMtx[c("geneid", colnames(dmso_myeloid))]
dmso_macSpec_inputMtx
dmso_macres <- macspec(mac_mtx = dmso_macSpec_inputMtx, feature = dmso_myeloid$treatment, select_hu_mo = "hum")
dmso_macres


# How mixed of a population is myeloid?
# add in MPI metadata to dmso_myeloid...
# toAdd <- subset(macres_tcs_merge, macres_tcs_merge$Samples %in% colnames(dmso_myeloid))
toAdd <- subset(dmso_macres, dmso_macres$Samples %in% colnames(dmso_myeloid))
intersect(toAdd$Samples, colnames(dmso_myeloid))

addMetaDF <- data.frame(row.names = toAdd$Samples, MPIscore = toAdd$MPI)
addMetaDF
dmso_myeloid <- AddMetaData(dmso_myeloid, metadata = addMetaDF)
dmso_myeloid$MPIscore

DefaultAssay(dmso_myeloid) <- "integrated"
dmso_myeloid <- FindVariableFeatures(dmso_myeloid)
dmso_myeloid <- ScaleData(dmso_myeloid, features = rownames(dmso_myeloid))
dmso_myeloid <- RunPCA(dmso_myeloid)
pdf(file = "dmso_myeloid_elbow.pdf")
ElbowPlot(dmso_myeloid)
dev.off()
dmso_myeloid <- FindNeighbors(dmso_myeloid, dims = 1:15)
dmso_myeloid <- FindClusters(dmso_myeloid, resolution = 0.2)
dmso_myeloid <- RunUMAP(dmso_myeloid, dims = 1:15, reduction.name = "dmso_umap")

# Visualize MPI score...
pdf(file = "DMSO_MPI_Score_DimPlot.pdf")
FeaturePlot(dmso_myeloid, features = "MPIscore", cols = viridis(256, option = "turbo"), reduction = "dmso_umap")
FeaturePlot(dmso_myeloid, features = "panobinostat", cols = viridis(256, option = "turbo"), reduction = "dmso_umap")
FeatureScatter(dmso_myeloid, feature1 = "MPIscore", feature2 = "panobinostat")
VlnPlot(dmso_myeloid, group.by = "seurat_clusters", features = c("MPIscore", "panobinostat"), pt.size = 0)
dev.off()


pdf(file = "dmso_myeloid_umap_01.pdf")
DimPlot(dmso_myeloid, group.by = c("seurat_clusters", "patientID"), reduction = "dmso_umap")
FeaturePlot(dmso_myeloid, features = c("CD74", "panobinostat", "CD36", "FCGR3A", "CD14", "CD8A", "NKG2A"), reduction = "dmso_umap", min.cutoff = 0, order = T)
dev.off()

myeloid_markers <- FindAllMarkers(dmso_myeloid, only.pos = T)
top_myeloid_markers <- myeloid_markers %>% 
  group_by(cluster) %>% 
  top_n(5, wt = avg_log2FC)

pdf(file = "myeloid_snn_markers.pdf")
DoHeatmap(dmso_myeloid, features = top_myeloid_markers$gene, assay = "RNA")
dev.off()

# FeatureScatter(dmso_myeloid, feature1 = "MPI", feature2 = "panobinostat")

dmso_macres_tcs_merge <- subset(macres_tcs_merge, macres_tcs_merge$treatment == "treatment: vehicle (DMSO)")

# Create the scatterplot
p <- ggplot(data = dmso_macres_tcs_merge, aes(x = MPI, y = panobinostat)) +
  geom_point() +                  # Add scatterplot points
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Add linear regression line
  labs(
    title = "Scatterplot with Correlation Statistics",
    x = "MPI",
    y = "panobinostat"
  )

# Calculate and display the correlation coefficient
cor_coef <- cor(dmso_macres_tcs_merge$MPI, dmso_macres_tcs_merge$panobinostat, method = "spearman")
cor_coef
p <- p + 
  annotate("text", x = max(dmso_macres_tcs_merge$MPI), y = max(dmso_macres_tcs_merge$panobinostat), 
           label = paste("Correlation =", c(round(cor_coef, 2))))

# Display the plot
pdf(file = "dmso_MPIvsPanobinostatTCSscatter.pdf")
print(p)
dev.off()

################################################################################
# Hexbin...
################################################################################

# library(hexbin)
# 
# # Create hexbin data
# hb <- hexbin(dmso_macres_tcs_merge$MPI, dmso_macres_tcs_merge$panobinostat)
# 
# # Create the hexbin plot with variable opacity
# dmso_macres_tcs_merge <- data.frame(x = hb@xcm, y = hb@ycm, hex_alpha = hb@count)
# 
# p <- ggplot(data = dmso_macres_tcs_merge, aes(x = x, y = y)) +
#   geom_hex(aes(fill = hex_alpha, alpha = hex_alpha)) +  # Add hexbin plot with variable opacity
#   labs(
#     title = "Hexbin Plot with Correlation Statistics",
#     x = "MPI",
#     y = "panobinostat"
#   ) +
#   scale_fill_gradient(low = "white", high = "blue") +  # Specify color palette
#   scale_alpha_continuous(range = c(0.1, 1))  # Adjust the range of alpha values
# 
# # Calculate and display the correlation coefficient
# # dmso_macres_tcs_merge$
# cor_coef <- cor(x = dmso_macres_tcs_merge$x, y = dmso_macres_tcs_merge$y)
# p <- p + 
#   annotate("text", x = max(dmso_macres_tcs_merge$x), y = max(dmso_macres_tcs_merge$y), 
#            label = paste("Correlation =", round(cor_coef, 2)))
#            
# # Display the plot
# pdf(file = "dmso_MPIvsPanobinostatTCSscatter_hex.pdmso_macres_tcs_merge")
# print(p)
# dev.off()


################################################################################

# rescale...
# pano TCS - does this need recalculated? 
# MPI vs TCS scatterplot?




macSpectrum::m.hum_mou_map
## The function is currently defined as
# function (mac_mtx, feature, select_hu_mo = "mou")
# {
#   rownames(M1_mean) <- M1_mean$GeneID
#   rownames(M2_mean) <- M2_mean$GeneID
#   rownames(M0_mean) <- M0_mean$GeneID
#   inFile <- mac_mtx
#   if (is.null(inFile)) {
#     return(NULL)
#   }
#   if (select_hu_mo == "hum") {
#     inFile_gene_id <- 1:nrow(inFile)
#     for (i in 1:nrow(inFile)) {
#       inFile_gene_id[i] <- as.character((hum_mou_map[,
#                                                      2])[inFile[i, 1] == hum_mou_map[, 1]])[1]
#     }
#     inFile_gene_id[is.na(inFile_gene_id)] <- "no_match"
#     inFile[, 1] <- inFile_gene_id
#   }
#   inFile <- inFile[!duplicated(inFile[, 1]), ]
#   rownames(inFile) <- inFile[, 1]
#   inFile_feature <- feature
#   if (is.null(inFile_feature)) {
#     return(NULL)
#   }
#   inFile <- inFile[, 2:ncol(inFile)]
#   inFile <- inFile - rowMeans(inFile)
#   MPI_genes <- intersect(M1_mean$GeneID, rownames(inFile))
#   M1_mean <- M1_mean[MPI_genes, ]
#   M2_mean <- M2_mean[MPI_genes, ]
#   inFile_bak <- inFile
#   inFile <- inFile[MPI_genes, ]
#   AMDI_genes <- intersect(M0_mean$GeneID, rownames(inFile_bak))
#   M0_mean <- M0_mean[AMDI_genes, ]
#   inFile_m0 <- inFile_bak[AMDI_genes, ]
#   inFile_sigma <- 1:ncol(inFile)
#   total_gene_number <- nrow(inFile)
#   for (i in 1:ncol(inFile)) {
#     options(digits = 9)
#     inFile_sigma[i] <- (sum(inFile[, i]^2)/total_gene_number)^0.5
#   }
#   inFile_sigma_m0 <- 1:ncol(inFile_m0)
#   total_gene_number <- nrow(inFile_m0)
#   for (i in 1:ncol(inFile_m0)) {
#     options(digits = 9)
#     inFile_sigma_m0[i] <- (sum(inFile_m0[, i]^2)/total_gene_number)^0.5
#   }
#   total_gene_number <- nrow(M0_mean)
#   M0_sigma <- (sum(M0_mean$value^2)/total_gene_number)^0.5
#   total_gene_number <- nrow(M1_mean)
#   M1_sigma <- (sum(M1_mean$value^2)/total_gene_number)^0.5
#   total_gene_number <- nrow(M2_mean)
#   M2_sigma <- (sum(M2_mean$value^2)/total_gene_number)^0.5
#   total_gene_number <- nrow(inFile_m0)
#   inFile_Pearson_per_cell_m0 <- 1:ncol(inFile_m0)
#   for (j in 1:ncol(inFile_m0)) {
#     inFile_Pearson_per_cell_m0[j] <- sum((inFile_m0[, j]/inFile_sigma_m0[j]) *
#                                            (M0_mean$value/M0_sigma))/total_gene_number
#   }
#   total_gene_number <- nrow(M2_mean)
#   inFile_Pearson_per_cell_m1 <- 1:ncol(inFile)
#   for (j in 1:ncol(inFile)) {
#     inFile_Pearson_per_cell_m1[j] <- sum((inFile[, j]/inFile_sigma[j]) *
#                                            (M1_mean$value/M1_sigma))/total_gene_number
#   }
#   total_gene_number <- nrow(M2_mean)
#   inFile_Pearson_per_cell_m2 <- 1:ncol(inFile)
#   for (j in 1:ncol(inFile)) {
#     inFile_Pearson_per_cell_m2[j] <- sum((inFile[, j]/inFile_sigma[j]) *
#                                            (M2_mean$value/M2_sigma))/total_gene_number
#   }
#   a <- 0.991414467
#   b <- 1
#   c <- -0.0185412856
#   x0 <- inFile_Pearson_per_cell_m1
#   y0 <- inFile_Pearson_per_cell_m2
#   d_sqr <- (a * x0 + b * y0 + c)^2/(a^2 + b^2)
#   x_start <- -1
#   y_start <- (-a) * x_start + (-c)
#   x_end <- 1
#   y_end <- (-a) * x_end + (-c)
#   l <- ((x0 - x_start)^2 + (y0 - y_start)^2 - d_sqr)^0.5
#   l_max <- ((x_end - x_start)^2 + (y_end - y_start)^2 - d_sqr)^0.5
#   MPI <- (l - 0)/(l_max - 0) * 100 - 50
#   AMDI <- -inFile_Pearson_per_cell_m0 * 50
#   mac_output <- data.frame(colnames(inFile), inFile_Pearson_per_cell_m1,
#                            inFile_Pearson_per_cell_m2, inFile_Pearson_per_cell_m0,
#                            inFile_feature[, 1], l, MPI, AMDI, row.names = colnames(inFile),
#                            stringsAsFactors = F)
#   colnames(mac_output) <- c("Samples", "rm1", "rm2", "r_m0",
#                             "Feature", "l", "MPI", "AMDI")
#   return(mac_output)
# }

# read in macrophage signatures file...

library(readxl)

ghosh <- as.data.frame(read_xlsx(path = "ghosh_et_al_mmc2_supp.xlsx", sheet = 1))
c13 <- data.frame(c13 = ghosh$`Gene Signature: C13`)
colnames(c13) <- "genes"
c13 <- c13[which(is.na(c13$genes) == FALSE),]
c13 <- data.frame(c13 = c13[-1])

#myeloid rank-data
my.rank.data <- singscore::rankGenes(my)

for (i in 1:length(colnames(c13))){
  # capture unique gene signatures (removes duplicate cells)
  uni <- unique(c13)[i]
  
  # capture name of cell state
  state <- colnames(c13)[i]
  
  # removes empty cells
  uni <- subset(uni, uni[,1] != "")
  uni <- as.vector(uni[,1])
  
  # scoring with singscore
  set <- GeneSet()
  set@geneIds <- as.character(uni)
  scored <- singscore::simpleScore(rankData = my.rank.data, upSet = set)
  
  # renaming "$Sig" to appropriate cell state
  scored$Sig <- as.character(state)
  
  assign(paste0("scored.", state), scored)
}

dim(scored)
length(myeloid$treatment)
pdf(file = "c13_dispersion_panobinostat.pdf")
singscore::plotDispersion(scored, annot = myeloid$treatment)
hist(scored$TotalScore)
dev.off()



################################################################################
# Neftel Scoring
################################################################################



