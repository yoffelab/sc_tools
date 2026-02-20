#!/usr/bin/env Rscript
# Convert a Seurat RDS file to AnnData h5ad using SeuratDisk.
# Usage: Rscript seurat_rds_to_h5ad.R <input.rds> <output.h5ad> [work_dir]
# Exits with 0 on success, non-zero on error.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  message("Usage: Rscript seurat_rds_to_h5ad.R <input.rds> <output.h5ad> [work_dir]")
  quit(save = "no", status = 1)
}

input_rds <- args[1]
output_h5ad <- args[2]
work_dir <- if (length(args) >= 3) args[3] else NULL

if (nzchar(work_dir)) {
  if (!dir.exists(work_dir)) {
    message("Error: work_dir does not exist: ", work_dir)
    quit(save = "no", status = 1)
  }
  setwd(work_dir)
}

if (!file.exists(input_rds)) {
  message("Error: input RDS not found: ", input_rds)
  quit(save = "no", status = 1)
}

# Resolve paths (after setwd if work_dir was set)
input_rds <- normalizePath(input_rds, mustWork = TRUE)
output_h5ad <- normalizePath(output_h5ad, mustWork = FALSE)
out_dir <- dirname(output_h5ad)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Temp h5Seurat in same dir as output so Convert() writes .h5ad alongside
tmp_base <- file.path(out_dir, paste0("_tmp_", basename(tools::file_path_sans_ext(input_rds))))
tmp_h5seurat <- paste0(tmp_base, ".h5Seurat")
tmp_h5ad_auto <- paste0(tmp_base, ".h5ad")

tryCatch(
  {
    library(Seurat, quietly = TRUE)
    library(SeuratDisk, quietly = TRUE)

    seurat_obj <- readRDS(input_rds)
    SaveH5Seurat(seurat_obj, filename = tmp_h5seurat, overwrite = TRUE)
    Convert(tmp_h5seurat, dest = "h5ad", overwrite = TRUE)

    # SeuratDisk writes .h5ad with same base name as .h5Seurat
    if (!file.exists(tmp_h5ad_auto)) {
      message("Error: Convert did not produce expected file: ", tmp_h5ad_auto)
      quit(save = "no", status = 1)
    }
    if (normalizePath(tmp_h5ad_auto, mustWork = TRUE) != normalizePath(output_h5ad, mustWork = FALSE)) {
      file.rename(tmp_h5ad_auto, output_h5ad)
    }
  },
  error = function(e) {
    message("Error: ", conditionMessage(e))
    quit(save = "no", status = 1)
  },
  finally = {
    if (file.exists(tmp_h5seurat)) {
      unlink(tmp_h5seurat)
    }
  }
)

quit(save = "no", status = 0)
