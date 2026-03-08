#!/usr/bin/env Rscript
# Convert Seurat v5 RDS files to AnnData h5ad for sc_tools pipeline
#
# Usage:
#   Rscript --vanilla convert_rds_to_h5ad.R <rds_path> <output_h5ad> <sample_num> <panel> [csv_path]
#
# Example:
#   Rscript --vanilla convert_rds_to_h5ad.R \
#     /athena/project-saha/data_IBD/CosMx_1k_16/1_smi_so_qcnorm.RDS \
#     data/cosmx_1k_01/adata.p0.h5ad \
#     1 cosmx_1k \
#     /athena/project-saha/data_IBD/NC_compare_11122025_meta.csv
#
# Handles:
#   - Seurat v5 layer structure (@layers$counts, not @counts)
#   - Assay names: Nanostring (CosMx) or Xenium
#   - Spatial coordinates from GetTissueCoordinates()
#   - Cell type labels (ct_minor, ct_major, ct_minor_new where available)
#   - Disease metadata joined from CSV by sample number + panel
#   - PCA/UMAP reductions preserved in obsm
#
# Output: h5ad with X = raw counts (sparse), obsm['spatial'], obs metadata

.libPaths(c("~/R/libs_R441", .libPaths()))

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})

# We write intermediate files (MTX + CSVs) and call Python to create h5ad
# This avoids issues with R anndata package var_df conversion

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: convert_rds_to_h5ad.R <rds_path> <output_h5ad> <sample_num> <panel> [csv_path]")
}

rds_path   <- args[1]
output_h5ad <- args[2]
sample_num <- as.integer(args[3])
panel      <- args[4]
csv_path   <- if (length(args) >= 5) args[5] else NULL

cat(sprintf("Converting: %s -> %s (sample=%d, panel=%s)\n", rds_path, output_h5ad, sample_num, panel))

# --- Load Seurat object ---
obj <- readRDS(rds_path)
cat(sprintf("Loaded: %d cells x %d genes, class=%s\n", ncol(obj), nrow(obj), class(obj)[1]))

# --- Identify primary assay ---
assay_names <- names(obj@assays)
primary_assay <- NULL
for (candidate in c("Nanostring", "Xenium", "RNA")) {
  if (candidate %in% assay_names) {
    primary_assay <- candidate
    break
  }
}
if (is.null(primary_assay)) {
  stop(sprintf("No known assay found. Available: %s", paste(assay_names, collapse = ", ")))
}
cat(sprintf("Primary assay: %s\n", primary_assay))

# --- Extract raw counts (Seurat v5: @layers$counts) ---
assay_obj <- obj@assays[[primary_assay]]
if ("counts" %in% names(assay_obj@layers)) {
  # Seurat v5 format
  X <- assay_obj@layers$counts
  cat(sprintf("Raw counts from @layers$counts: %dx%d\n", nrow(X), ncol(X)))
} else if (!is.null(tryCatch(assay_obj@counts, error = function(e) NULL))) {
  # Seurat v4 fallback
  X <- assay_obj@counts
  cat(sprintf("Raw counts from @counts: %dx%d\n", nrow(X), ncol(X)))
} else {
  stop("Cannot find raw counts in assay")
}

# Get gene and cell names from the assay object (Seurat v5 layers don't have dimnames)
gene_names <- as.character(rownames(assay_obj))
cell_names <- as.character(colnames(assay_obj))
# Fallback: try the Seurat object itself
if (length(gene_names) == 0) gene_names <- as.character(rownames(obj))
if (length(cell_names) == 0) cell_names <- as.character(colnames(obj))
cat(sprintf("Gene names: %d, Cell names: %d\n", length(gene_names), length(cell_names)))
X_t <- t(X)
rownames(X_t) <- cell_names
colnames(X_t) <- gene_names
cat(sprintf("Transposed: %dx%d (cells x genes)\n", nrow(X_t), ncol(X_t)))

# --- Build obs (cell metadata) ---
obs_df <- obj@meta.data

# Add sample identifiers
obs_df$sample <- sprintf("%s_%02d", panel, sample_num)
obs_df$library_id <- obs_df$sample
obs_df$raw_data_dir <- rds_path
obs_df$panel <- panel

# Rename cell type columns for consistency
if ("ct_minor" %in% colnames(obs_df))     obs_df$celltype <- obs_df$ct_minor
if ("ct_major" %in% colnames(obs_df))     obs_df$celltype_broad <- obs_df$ct_major
if ("ct_minor_new" %in% colnames(obs_df)) obs_df$celltype_new <- obs_df$ct_minor_new

# --- Join disease metadata from CSV ---
if (!is.null(csv_path) && file.exists(csv_path)) {
  cat("Joining metadata from CSV...\n")
  meta_csv <- read.csv(csv_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  # Clean column names (remove BOM if present)
  colnames(meta_csv) <- trimws(gsub("^\uFEFF", "", colnames(meta_csv)))

  # Map panel name to CSV panel column values
  panel_map <- list(
    cosmx_1k       = "CosMx 1K",
    cosmx_6k       = "CosMx 6k",
    xenium_mt      = "Xenium 377 MT",
    xenium_5k      = "Xenium 5K",
    xenium_noseg   = "Xenium 377 MT (no segmentation marker)",
    xenium_withseg = "Xenium 377 MT (with segmentation marker)",
    xenium_colon   = "Xenium colon panel"
  )

  csv_panel <- panel_map[[panel]]
  if (!is.null(csv_panel)) {
    row <- meta_csv[meta_csv$panel == csv_panel & meta_csv$sample_id == sample_num, ]
    if (nrow(row) == 1) {
      obs_df$patient_id    <- row$block
      obs_df$disease       <- row$disease
      obs_df$disease_state <- row$disease_state
      obs_df$tissue_type   <- row$tissue_type
      cat(sprintf("  Matched: patient=%s, disease=%s, state=%s, tissue=%s\n",
                  row$block, row$disease, row$disease_state, row$tissue_type))
    } else {
      cat(sprintf("  WARNING: %d rows matched panel=%s sample=%d in CSV\n",
                  nrow(row), csv_panel, sample_num))
    }
  } else {
    cat(sprintf("  WARNING: panel '%s' not in panel_map\n", panel))
  }
} else if (!is.null(csv_path)) {
  cat(sprintf("  WARNING: CSV not found: %s\n", csv_path))
}

# --- Spatial coordinates ---
obsm <- list()
coords <- tryCatch({
  tc <- GetTissueCoordinates(obj)
  as.matrix(tc[, c("x", "y")])
}, error = function(e) {
  # Fallback: try meta.data columns
  if (all(c("x_centroid", "y_centroid") %in% colnames(obs_df))) {
    as.matrix(obs_df[, c("x_centroid", "y_centroid")])
  } else if (all(c("CenterX_global_px", "CenterY_global_px") %in% colnames(obs_df))) {
    as.matrix(obs_df[, c("CenterX_global_px", "CenterY_global_px")])
  } else {
    cat("  WARNING: no spatial coordinates found\n")
    NULL
  }
})

if (!is.null(coords)) {
  obsm[["spatial"]] <- coords
  cat(sprintf("Spatial coords: %dx%d, x=[%.1f,%.1f] y=[%.1f,%.1f]\n",
              nrow(coords), ncol(coords),
              min(coords[, 1]), max(coords[, 1]),
              min(coords[, 2]), max(coords[, 2])))
}

# --- Reductions (PCA, UMAP) ---
for (red_name in names(obj@reductions)) {
  emb <- Embeddings(obj@reductions[[red_name]])
  obsm_key <- sprintf("X_%s", red_name)
  obsm[[obsm_key]] <- emb
  cat(sprintf("Reduction %s -> obsm['%s']: %dx%d\n", red_name, obsm_key, nrow(emb), ncol(emb)))
}

# --- Write intermediate files and call Python to build h5ad ---
# Clean obs_df: convert factors to character, remove problematic columns
for (col in colnames(obs_df)) {
  if (is.factor(obs_df[[col]])) {
    obs_df[[col]] <- as.character(obs_df[[col]])
  }
}
keep_cols <- sapply(obs_df, function(x) is.atomic(x) && !is.list(x))
obs_df <- obs_df[, keep_cols, drop = FALSE]

# Create output directory and temp dir
dir.create(dirname(output_h5ad), recursive = TRUE, showWarnings = FALSE)
tmp_dir <- tempdir()

# Write sparse matrix (cells x genes) in MatrixMarket format
mtx_path <- file.path(tmp_dir, "matrix.mtx")
writeMM(X_t, mtx_path)

# Write gene names and cell barcodes as single-column CSVs
gene_df <- data.frame(gene = gene_names, stringsAsFactors = FALSE)
write.csv(gene_df, file.path(tmp_dir, "genes.csv"), row.names = FALSE)
barcode_df <- data.frame(barcode = cell_names, stringsAsFactors = FALSE)
write.csv(barcode_df, file.path(tmp_dir, "barcodes.csv"), row.names = FALSE)
cat(sprintf("Wrote %d genes, %d barcodes\n", length(gene_names), length(cell_names)))

# Write obs metadata
obs_path <- file.path(tmp_dir, "obs.csv")
write.csv(obs_df, obs_path, row.names = TRUE)

# Write spatial coordinates
if (!is.null(coords)) {
  coord_path <- file.path(tmp_dir, "spatial.csv")
  # coords may have different row names than cell_names; use coords as-is
  coord_df <- data.frame(x = coords[, 1], y = coords[, 2])
  rownames(coord_df) <- rownames(coords)
  write.csv(coord_df, coord_path)
}

# Write embeddings
for (emb_name in names(obsm)) {
  emb_path <- file.path(tmp_dir, sprintf("%s.csv", emb_name))
  write.csv(obsm[[emb_name]], emb_path, row.names = TRUE)
}

cat(sprintf("Intermediate files written to %s\n", tmp_dir))

# Call Python to assemble h5ad
py_script <- sprintf('
import scipy.io, pandas as pd, numpy as np, anndata as ad
from pathlib import Path
from scipy.sparse import csr_matrix

tmp = "%s"
out = "%s"

# R writeMM writes cells x genes (already transposed in R)
X = scipy.io.mmread(f"{tmp}/matrix.mtx")
X = csr_matrix(X)

import csv
with open(f"{tmp}/genes.csv") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    genes = [row[0] for row in reader]
with open(f"{tmp}/barcodes.csv") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    barcodes = [row[0] for row in reader]

print(f"X shape: {X.shape}, genes: {len(genes)}, barcodes: {len(barcodes)}")

# Ensure X is cells x genes
if X.shape[0] == len(genes) and X.shape[1] == len(barcodes):
    X = X.T  # was genes x cells, transpose to cells x genes
    print(f"Transposed X to: {X.shape}")

obs = pd.read_csv(f"{tmp}/obs.csv", index_col=0)
obs.index = obs.index.astype(str)

var = pd.DataFrame(index=genes)

print(f"obs: {obs.shape}, var: {var.shape}, X: {X.shape}")

adata = ad.AnnData(X=X, obs=obs, var=var)

# Spatial coords
sp = Path(f"{tmp}/spatial.csv")
if sp.exists():
    coords = pd.read_csv(sp, index_col=0)
    adata.obsm["spatial"] = coords[["x","y"]].values.astype(np.float32)

# Embeddings
for f in Path(tmp).glob("X_*.csv"):
    key = f.stem
    emb = pd.read_csv(f, index_col=0).values.astype(np.float32)
    adata.obsm[key] = emb

adata.write_h5ad(out)
print(f"Written: {out} ({adata.n_obs} cells x {adata.n_vars} genes, {Path(out).stat().st_size/1024**2:.1f} MB)")
', tmp_dir, output_h5ad)

# Find Python - prefer conda sc_tools env, fall back to system
py_bin <- Sys.which("python3")
if (nchar(py_bin) == 0) py_bin <- Sys.which("python")
if (nchar(py_bin) == 0) py_bin <- "/usr/bin/python3"

py_file <- file.path(tmp_dir, "assemble_h5ad.py")
writeLines(py_script, py_file)

exit_code <- system2(py_bin, py_file)
if (exit_code != 0) {
  # Try with conda env
  conda_py <- "/home/fs01/juk4007/.conda/envs/sc_tools/bin/python"
  if (file.exists(conda_py)) {
    cat("Retrying with conda sc_tools env...\n")
    exit_code <- system2(conda_py, py_file)
  }
}

if (exit_code != 0) {
  stop("Python h5ad assembly failed")
}

cat("DONE\n")
