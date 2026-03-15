#!/usr/bin/env Rscript
# Convert DLBCL CosMx RDS objects to h5ad for Python analysis
# Outputs:
#   results/cosmx_prt.h5ad  — protein panel (67 markers, 687K cells, has final_cell_type)
#   results/cosmx_rna.h5ad  — RNA panel (1000 genes, 464K cells, has cluster labels)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(Matrix)
})

DATA_DIR <- "projects/cosmx_1k/lymph_dlbcl/data"
OUT_DIR  <- "projects/cosmx_1k/lymph_dlbcl/results"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

convert_seurat_to_h5ad <- function(rds_path, out_stem, celltype_col = NULL,
                                   umap_key = NULL) {
  cat("Loading:", rds_path, "\n")
  obj <- readRDS(rds_path)
  cat("  Dims:", dim(obj), "\n")

  # --- rename metadata columns to clean names ---
  meta <- obj@meta.data

  # Find cell type column
  if (!is.null(celltype_col) && celltype_col %in% colnames(meta)) {
    meta$celltype <- meta[[celltype_col]]
    cat("  Cell types:", paste(sort(unique(meta$celltype)), collapse=", "), "\n")
  }

  # Spatial coordinates
  spatial_cols <- c("x_slide_mm", "y_slide_mm", "x_FOV_px", "y_FOV_px")
  for (col in spatial_cols) {
    if (col %in% colnames(meta)) meta[[col]] <- meta[[col]]
  }

  # Sample
  if ("Run_Tissue_name" %in% colnames(meta)) {
    meta$sample <- meta$Run_Tissue_name
  }
  if ("fov" %in% colnames(meta)) meta$fov <- meta$fov

  # Replace UUID column names with clean names
  nn_cols <- grep("^nn_|^RNA_nbclust_|^RNA_spatialclust_|^spatialclust_|^c77|^c88", colnames(meta))
  if (length(nn_cols) > 0) {
    meta <- meta[, -nn_cols, drop = FALSE]
  }

  obj@meta.data <- meta

  # Find best UMAP reduction
  if (!is.null(umap_key) && umap_key %in% names(obj@reductions)) {
    # rename to "umap" for standard handling
    obj@reductions[["umap"]] <- obj@reductions[[umap_key]]
  } else {
    umap_reds <- grep("umap", names(obj@reductions), ignore.case = TRUE, value = TRUE)
    if (length(umap_reds) > 0) {
      obj@reductions[["umap"]] <- obj@reductions[[umap_reds[1]]]
      cat("  Using UMAP reduction:", umap_reds[1], "\n")
    }
  }

  # Save as .h5seurat then convert to .h5ad
  h5seurat_path <- file.path(OUT_DIR, paste0(out_stem, ".h5seurat"))
  h5ad_path     <- file.path(OUT_DIR, paste0(out_stem, ".h5ad"))

  cat("  Saving h5seurat:", h5seurat_path, "\n")
  SaveH5Seurat(obj, filename = h5seurat_path, overwrite = TRUE)

  cat("  Converting to h5ad:", h5ad_path, "\n")
  Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)

  # Clean up intermediate file
  file.remove(h5seurat_path)

  cat("  Done:", h5ad_path, "\n")
  invisible(h5ad_path)
}

# --- Convert PRT object (protein panel, has cell type labels) ---
prt_celltype_col <- "c77ade37f.4ab3.44de.bb38.e7593c07ef6c_1_final_cell_type"
convert_seurat_to_h5ad(
  file.path(DATA_DIR, "DLBCL_cosmx_run1_CTMA_PRT.RDS"),
  out_stem = "cosmx_prt",
  celltype_col = prt_celltype_col,
)

# --- Convert RNA object (1000-plex, has cluster labels) ---
rna_cluster_col <- "RNA_nbclust_60580b7c.7ccc.4bc3.afd2.a28eaa025e28_1_clusters"
convert_seurat_to_h5ad(
  file.path(DATA_DIR, "DLBCL_cosmx_run1_CTMA.RDS"),
  out_stem = "cosmx_rna",
  celltype_col = rna_cluster_col,
)

cat("Conversion complete.\n")
