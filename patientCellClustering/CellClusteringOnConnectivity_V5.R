################################################################################
#
# Full Analysis Pipeline: Per-Patient Clustering & Cross-Patient Metaclustering
#
################################################################################

# --- Load All Necessary Libraries ---
# Ensure all these packages are installed before running the script.
# install.packages(c("Seurat", "ComplexHeatmap", "circlize", "dplyr", "purrr", "factoextra", "dittoSeq"))

library(Seurat)
library(ComplexHeatmap)
library(colorRamp2)
library(scales)
library(circlize)
library(tibble)
library(dplyr)
library(corrplot)
library(factoextra)
library(cluster)
library(dittoSeq)
library(purrr)


################################################################################
# PART 1: PER-PATIENT ANALYSIS FUNCTION
################################################################################

#' Analyze Drug Correlation and Generate Plots for Each Patient
#'
#' This function iterates through each unique patient in a Seurat object,
#' subsets the drug correlation data, performs clustering, and generates
#' plots and cluster data files.
#'
#' @param seurat_obj A Seurat object containing the single-cell data.
#' @param full_corr_matrix A data frame or matrix of the full drug correlation data.
#' @param fda_onc_corr_matrix A data frame or matrix of the FDA-approved oncology drug correlation data.
#' @param output_dir A string specifying the directory path to save the output.
#' @param k_for_heatmap An integer specifying the number of clusters (k) to use for heatmap annotations. Default is 7.
#' @param save_cluster_data A logical value indicating whether to save the metadata with appended
#' cluster information for each patient as a CSV and RDS file. Default is TRUE.
#'
#' @return This function does not return a value but saves plots and data files
#' to the specified output directory.
#'
analyze_patient_drug_correlations <- function(seurat_obj,
                                              full_corr_matrix,
                                              fda_onc_corr_matrix,
                                              output_dir = "patient_analysis_output",
                                              k_for_heatmap = 7,
                                              save_cluster_data = TRUE) {
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get the complete metadata from the Seurat object
  meta <- seurat_obj@meta.data
  
  # Get a list of unique patient IDs and exclude 'GBM41'
  # patient_ids <- setdiff(unique(meta$patientID), "GBM41")
  patient_ids <- "GBM41"
  # Check if patientID column exists
  if (is.null(unique(meta$patientID))) {
    stop("Error: 'patientID' column not found in the Seurat object metadata.")
  }
  
  # Loop through each patient ID
  for (patient_id in patient_ids) {
    
    message(paste("Processing patient:", patient_id))
    
    # Subset metadata for the current patient and filter out non-neoplastic cells
    patient_meta <- subset(meta, patientID == patient_id)
    patient_meta <- subset(patient_meta, subset = Neftel_State != "Non-Neoplastic")
    
    if (nrow(patient_meta) < 10) {
      warning(paste("Skipping patient", patient_id, "due to insufficient neoplastic cells (< 10)."))
      next
    }
    
    patient_cells <- rownames(patient_meta)
    
    # Subset the correlation matrices
    patient_subset_fda_onc <- fda_onc_corr_matrix[rownames(fda_onc_corr_matrix) %in% patient_cells,
                                                  colnames(fda_onc_corr_matrix) %in% patient_cells]
    
    # Perform clustering
    clust_fda_onc <- patient_subset_fda_onc %>%
      dist(method = "euclidean") %>%
      hclust(method = "ward.D2")
    
    clusters_fda_onc <- as.data.frame(cutree(clust_fda_onc, k = seq(2, 30, 1)))
    colnames(clusters_fda_onc) <- paste0("FDA_Onc_k_", colnames(clusters_fda_onc))
    
    # Generate Elbow Plot
    calc <- patient_subset_fda_onc
    diag(calc) <- 0
    elbow_plot_path <- file.path(output_dir, paste0(patient_id, "_fda_onc_elbow.pdf"))
    
    pdf(file = elbow_plot_path)
    print(
      fviz_nbclust(calc, kmeans, method = "wss", k.max = 24) +
        theme_minimal() +
        ggtitle(paste(patient_id, "Elbow Plot"))
    )
    dev.off()
    
    # Combine and Save Data
    patient_meta <- cbind(patient_meta, clusters_fda_onc[rownames(patient_meta), , drop = FALSE])
    
    if (save_cluster_data) {
      csv_path <- file.path(output_dir, paste0(patient_id, "_metadata_with_clusters.csv"))
      rds_path <- file.path(output_dir, paste0(patient_id, "_metadata_with_clusters.rds"))
      
      write.csv(patient_meta, file = csv_path, row.names = TRUE)
      saveRDS(patient_meta, file = rds_path)
      
      message(paste("Saved clustering data for", patient_id))
    }
    
    # Prepare Annotations for Heatmap
    annotation_criteria <- select(patient_meta, c('Neftel_State', 'Phase', 'results.CytoTRACE'))
    row_clusters <- cutree(clust_fda_onc, k = k_for_heatmap)
    row_annotation <- data.frame(cluster = as.character(row_clusters))
    
    # --- Generate and Save Heatmaps ---
    
    # Define color palette and breaks
    custpal <- colorRampPalette(c("steelblue3", "white", "firebrick3"))(100)
    max_abs_val <- max(abs(patient_subset_fda_onc), na.rm = TRUE)
    breaks <- seq(-max_abs_val, max_abs_val, length.out = length(custpal) + 1)
    
    # Dynamically create colors for the k clusters
    cluster_colors <- scales::hue_pal()(k_for_heatmap)
    names(cluster_colors) <- 1:k_for_heatmap
    
    annotation_colors_list <- list(
      cluster = cluster_colors,
      Neftel_State = c(AC = "green", MES = "red", NPC = "blue", OPC = "purple", "Non-Neoplastic" = "grey"),
      Phase = c(G1 = "aquamarine", G2M = "maroon", S = "violet"),
      results.CytoTRACE = c("white", "firebrick")
    )
    
    # Generate the heatmap object (non-rasterized)
    chm <- pheatmap::pheatmap(
      patient_subset_fda_onc, use_raster = FALSE, color = custpal, breaks = breaks,
      annotation_col = annotation_criteria, annotation_row = row_annotation,
      clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
      annotation_colors = annotation_colors_list, clustering_method = "ward.D2",
      show_rownames = FALSE, show_colnames = FALSE, silent = TRUE,
      main = paste(patient_id, "FDA Oncology Drug Correlations")
    )
    
    # Generate the rasterized heatmap object
    chm_rast <- pheatmap::pheatmap(
      patient_subset_fda_onc, use_raster = TRUE, color = custpal, breaks = breaks,
      annotation_col = annotation_criteria, annotation_row = row_annotation,
      clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
      annotation_colors = annotation_colors_list, clustering_method = "ward.D2",
      show_rownames = FALSE, show_colnames = FALSE, silent = TRUE,
      main = paste(patient_id, "FDA Oncology Drug Correlations (Rasterized)")
    )
    
    # Save the non-rasterized heatmaps (PDF and PNG)
    heatmap_path_pdf <- file.path(output_dir, paste0(patient_id, "_FDAOnc_chm.pdf"))
    heatmap_path_png <- file.path(output_dir, paste0(patient_id, "_FDAOnc_chm.png"))
    pdf(file = heatmap_path_pdf, width = 12, height = 10); print(chm); dev.off()
    png(file = heatmap_path_png, width = 12, height = 10, units = "in", res = 300); print(chm); dev.off()
    
    # Save the rasterized heatmaps (PDF and PNG)
    heatmap_rast_path_pdf <- file.path(output_dir, paste0(patient_id, "_FDAOnc_chm_rast.pdf"))
    heatmap_rast_path_png <- file.path(output_dir, paste0(patient_id, "_FDAOnc_chm_rast.png"))
    pdf(file = heatmap_rast_path_pdf, width = 12, height = 10); print(chm_rast); dev.off()
    png(file = heatmap_rast_path_png, width = 12, height = 10, units = "in", res = 300); print(chm_rast); dev.off()
    
    message(paste("Finished processing for patient:", patient_id))
  }
}


################################################################################
# SCRIPT EXECUTION WORKFLOW
################################################################################

# --- 1. Load Primary Data ---
# This section should be edited with the correct paths to your files.
message("Loading required data objects...")
obj <- readRDS(file = "isosceles_suterJohnsonMerge_fourStates.RDS")
comp_cluster_corres_full2017 <- readRDS(file = "patient_2017_drugDat_CorrMatrix.RDS")
comp_cluster_corres_full2017_FDAOnc <- readRDS(file = "patient_2017_drugDat_CorrMatrix_FDAOnc.RDS")


# --- 2. Run Per-Patient Analysis ---
# This will generate the individual patient output files needed for the next step.
output_directory <- "GBM_CombinationScoring_ISOSCELES_RV/Patient_Analysis_Output_2024"
# message("Starting analysis for all patients...")
analyze_patient_drug_correlations(
  seurat_obj = obj,
  full_corr_matrix = comp_cluster_corres_full2017,
  fda_onc_corr_matrix = comp_cluster_corres_full2017_FDAOnc,
  output_dir = output_directory,
  k_for_heatmap = 7,
  save_cluster_data = TRUE
)
message("Per-patient analysis complete.")


################################################################################
# PART 2: META-CLUSTERING ANALYSIS
################################################################################

# --- 3. Aggregate All Patient Cluster Data ---
message("Aggregating results for metaclustering...")
rds_files <- list.files(path = output_directory, pattern = "_metadata_with_clusters.rds$", full.names = TRUE)

# Check if files were found
if (length(rds_files) == 0) {
  stop(paste("No RDS files found in directory:", output_directory, ". Please run Part 1 first."))
}

all_patients_meta <- map_dfr(rds_files, readRDS)


# --- 4. Create a Profile for Each Cluster ---
# We summarize each cluster by its biological properties.
# Here, we'll use the k=7 clustering results as an example.
message("Creating summary profiles for each cluster...")
cluster_profiles <- all_patients_meta %>%
  group_by(patientID, FDA_Onc_k_7) %>%
  summarise(
    # Calculate proportions for Neftel States
    prop_AC = mean(Neftel_State == "AC", na.rm = TRUE),
    prop_MES = mean(Neftel_State == "MES", na.rm = TRUE),
    prop_NPC = mean(Neftel_State == "NPC", na.rm = TRUE),
    prop_OPC = mean(Neftel_State == "OPC", na.rm = TRUE),
    
    # Calculate proportions for Cell Cycle Phases
    prop_G1 = mean(Phase == "G1", na.rm = TRUE),
    prop_S = mean(Phase == "S", na.rm = TRUE),
    prop_G2M = mean(Phase == "G2M", na.rm = TRUE),
    
    # Calculate mean CytoTRACE score
    mean_cytoTRACE = mean(results.CytoTRACE, na.rm = TRUE),
    
    # Count number of cells in the cluster
    n_cells = n(),
    .groups = 'drop' # Ungroup after summarising
  ) %>%
  # Create a unique ID for each patient-specific cluster
  mutate(cluster_id = paste(patientID, FDA_Onc_k_7, sep = "_k"))

message("Cluster profiles created:")
print(head(cluster_profiles))


# --- 5. Cluster the Cluster Profiles (Metaclustering) ---
message("Performing metaclustering on cluster profiles...")
# Prepare the data for clustering (use the profiles, not the IDs)
profiles_for_clustering <- cluster_profiles %>%
  select(prop_AC, prop_MES, prop_NPC, prop_OPC, prop_G1, prop_S, prop_G2M, mean_cytoTRACE) %>%
  as.matrix()

# Assign the unique cluster IDs as row names
rownames(profiles_for_clustering) <- cluster_profiles$cluster_id

# Scale the data so that features with large values don't dominate
profiles_for_clustering_scaled <- scale(profiles_for_clustering)

# Perform hierarchical clustering
metacluster_hclust <- hclust(dist(profiles_for_clustering_scaled), method = "ward.D2")


# --- 6. Visualize and Interpret the Metaclusters ---
message("Generating metacluster heatmap...")
metacluster_heatmap_path <- file.path(output_directory, "Metacluster_Heatmap.pdf")

# Determine a good number of metaclusters (e.g., k=5)
num_metaclusters <- 5 

metacluster_colors <- colorRampPalette(c("steelblue", "white", "firebrick"))(100)

# Create symmetric breaks for the color scale centered at zero
max_abs_val <- max(abs(profiles_for_clustering_scaled), na.rm = TRUE)
heatmap_breaks <- seq(-max_abs_val, max_abs_val, length.out = 101) # length.out = number of colors + 1

pdf(metacluster_heatmap_path, width = 14, height = 3)
pheatmap(
  t(profiles_for_clustering_scaled), # Transpose the matrix
  color = metacluster_colors,
  breaks = heatmap_breaks, # Add the breaks for a zero-centered scale
  cluster_rows = FALSE, # Features are now rows, keep them in original order
  cluster_cols = metacluster_hclust, # Clusters are now columns
  show_rownames = TRUE,
  show_colnames = TRUE, # Hide individual cluster names to avoid clutter
  fontsize_row = 10,
  cutree_cols = num_metaclusters, # Cut the column tree to show metaclusters
  main = "Metaclusters of Patient Drug-Response Profiles"
)
dev.off()

message(paste("Metaclustering complete. Heatmap saved to:", metacluster_heatmap_path))

write.csv(profiles_for_clustering, file = "metacluster_heatmap_sourceData.csv")
write.csv(t(profiles_for_clustering_scaled), file = "metacluster_heatmap_scaled.csv")
