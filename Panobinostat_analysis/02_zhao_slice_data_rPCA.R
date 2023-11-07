##----
# load libraries

library(ggplot2)
library(cowplot)
library(GEOquery)
library(Seurat)
library(Hmisc)
library(dplyr)

################################################################################

# saveRDS(obj.list, file = "zhao_individual_object_list.RDS")
obj.list <- readRDS(file = "zhao_individual_object_list.RDS")
##----
# reciprocal pca integration to ID discrete cell types

# integrate post-qc filtering

# obj.list <- lapply(X = obj.list, FUN = function(x) {
#   x <- NormalizeData(x, verbose = FALSE)
#   x <- FindVariableFeatures(x, verbose = FALSE)
# })

features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = obj.list, reference = c(1, 2), reduction = "rpca",
                                  dims = 1:50)

# length(unlist(lapply(obj.list, FUN = "rownames")))
to_integrate <- unique(unlist(lapply(obj.list, FUN = "rownames")))

obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)#, features.to.integrate = to_integrate)
rm(obj.list)
rm(anchors)
gc()

obj.integrated <- ScaleData(obj.integrated, verbose = TRUE)
obj.integrated <- RunPCA(obj.integrated, verbose = FALSE)
obj.integrated <- FindNeighbors(obj.integrated, dims = 1:50)
obj.integrated <- FindClusters(obj.integrated)

obj.integrated <- Seurat::RunUMAP(obj.integrated, dims = 1:50, reduction.name = "umap_int")

# FeaturePlot(obj.integrated, features = c("AURKA", "CD74", "CDK4", "APOE"), min.cutoff = 0, slot = "scale.data", cols = c("grey", "red"), order = T, reduction = "umap_int")

saveRDS(obj.integrated, file = "zhao_integrated_varFeatures.RDS")

# DimPlot(obj.integrated, group.by = c("gsm_id", "seurat_clusters"), reduction = "umap_int")

dim(obj.integrated)
length(unique(obj.integrated$cts_file))
