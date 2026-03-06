#!/usr/bin/env Rscript
# convert_rds_to_h5ad.R — Convert Seurat RDS objects to h5ad on cayuga
#
# Usage: Rscript scripts/convert_rds_to_h5ad.R
#
# Reads RDS from backup path, converts via SeuratDisk, writes to results/seurat_converted/

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

BACKUP_DIR <- "/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2"
OUTPUT_DIR <- "results/seurat_converted"

# Map of RDS files to convert (relative to BACKUP_DIR -> relative output dir)
files_to_convert <- list(
  # Essential files
  list(rds = "tcell_2_preprocessing/t2_SO_seurat.rds",
       out = "tcell_2_preprocessing/t2_SO_seurat.h5ad"),
  list(rds = "stroma_2_preprocessing/S1_seurat_SO.rds",
       out = "stroma_2_preprocessing/S1_seurat_SO.h5ad"),
  list(rds = "tcell_merged/SO2_seurat_bcell.rds",
       out = "tcell_merged/SO2_seurat_bcell.h5ad"),
  list(rds = "stroma_merged/SO2_seurat_bcell.rds",
       out = "stroma_merged/SO2_seurat_bcell.h5ad"),
  list(rds = "stroma_spatial/SO_k30_community_cluster.rds",
       out = "stroma_spatial/SO_k30_community_cluster.h5ad"),
  list(rds = "s_tme_seurat.rds",
       out = "s_tme_seurat.h5ad"),

  # B-cell subsets (supp_fig2)
  list(rds = "tcell_2_preprocessing/t2_bcell2_seurat.rds",
       out = "tcell_2_preprocessing/t2_bcell2_seurat.h5ad"),
  list(rds = "stroma_2_preprocessing/S1_seurat_bcell.rds",
       out = "stroma_2_preprocessing/S1_seurat_bcell.h5ad"),

  # T-cell/myeloid subsets (supp_fig3)
  list(rds = "tcell_2_preprocessing/1.29.23.tcell_seurat.rds",
       out = "tcell_2_preprocessing/1.29.23.tcell_seurat.h5ad"),
  list(rds = "tcell_2_preprocessing/t2_tcell2_seurat.rds",
       out = "tcell_2_preprocessing/t2_tcell2_seurat.h5ad"),
  list(rds = "tcell_2_preprocessing/t2_myeloid2_seurat.rds",
       out = "tcell_2_preprocessing/t2_myeloid2_seurat.h5ad"),
  list(rds = "stroma_1_preprocessing/S1_seurat_myeloid2.rds",
       out = "stroma_1_preprocessing/S1_seurat_myeloid2.h5ad"),

  # Stroma subset (supp_fig4 vessel)
  list(rds = "stroma_2_preprocessing/S1_seurat_stroma.rds",
       out = "stroma_2_preprocessing/S1_seurat_stroma.h5ad")
)

cat("=== RDS to h5ad Conversion ===\n")
cat(sprintf("Backup dir: %s\n", BACKUP_DIR))
cat(sprintf("Output dir: %s\n", OUTPUT_DIR))
cat(sprintf("Files to convert: %d\n\n", length(files_to_convert)))

# Create output directories
for (f in files_to_convert) {
  dir.create(file.path(OUTPUT_DIR, dirname(f$out)), recursive = TRUE, showWarnings = FALSE)
}

success <- 0
failed <- 0

for (f in files_to_convert) {
  rds_path <- file.path(BACKUP_DIR, f$rds)
  h5ad_path <- file.path(OUTPUT_DIR, f$out)
  h5seurat_path <- sub("\\.h5ad$", ".h5seurat", h5ad_path)

  # Skip if already converted
  if (file.exists(h5ad_path)) {
    cat(sprintf("SKIP (exists): %s\n", f$out))
    success <- success + 1
    next
  }

  if (!file.exists(rds_path)) {
    cat(sprintf("MISSING: %s\n", rds_path))
    failed <- failed + 1
    next
  }

  cat(sprintf("Converting: %s -> %s\n", f$rds, f$out))
  tryCatch({
    obj <- readRDS(rds_path)

    # Ensure default assay
    if (inherits(obj, "Seurat")) {
      # Save as h5Seurat first, then convert to h5ad
      SaveH5Seurat(obj, filename = h5seurat_path, overwrite = TRUE)
      Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)

      # Clean up intermediate h5seurat
      if (file.exists(h5seurat_path)) file.remove(h5seurat_path)

      cat(sprintf("  OK: %s (%.1f MB)\n", f$out,
                  file.info(h5ad_path)$size / 1e6))
      success <- success + 1
    } else {
      cat(sprintf("  WARN: Not a Seurat object: %s (class: %s)\n",
                  f$rds, class(obj)[1]))
      failed <- failed + 1
    }
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n  %s\n", f$rds, conditionMessage(e)))
    failed <<- failed + 1
  })

  # Free memory
  gc(verbose = FALSE)
}

cat(sprintf("\n=== Done: %d success, %d failed ===\n", success, failed))
