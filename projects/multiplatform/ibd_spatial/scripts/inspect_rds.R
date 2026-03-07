#!/usr/bin/env Rscript
# Step 0 metadata inspection script
# Run on cayuga: Rscript --vanilla scripts/inspect_rds.R
# Output: ~/elementolab/projects/ibd_spatial/logs/inspect_<date>.out
#
# Confirms: raw counts, cell type columns, spatial coords, noseg cell status

panels <- list(
  cosmx_1k       = "/athena/project-saha/data_IBD/CosMx_1k_16/1_smi_so_qcnorm.RDS",
  cosmx_6k       = "/athena/project-saha/data_IBD/CosMx_6k_4/1_smi_so_qcnorm.RDS",
  xenium_5k      = "/athena/project-saha/data_IBD/Xenium_5K_4/1_xen_so_qcnorm.RDS",
  xenium_mt      = "/athena/project-saha/data_IBD/Xenium_MT_16/1_xen_so_qcnorm.RDS",
  xenium_noseg   = "/athena/project-saha/data_IBD/Xenium_MT_noseg_4/1_xen_377_no_seg_so_qcnorm.RDS",
  xenium_withseg = "/athena/project-saha/data_IBD/Xenium_MT_withseg_4/1_xen_377_with_seg_so_qcnorm.RDS",
  xenium_colon   = "/athena/project-saha/data_IBD/Xenium_colon_4/1_xen_colon_with_seg_so_qcnorm.RDS"
)

load_obj <- function(path) {
  # Try readRDS first; fall back to load() for RData/save() format
  tryCatch(readRDS(path), error = function(e) {
    env <- new.env()
    tryCatch({
      load(path, envir = env)
      vars <- ls(env)
      if (length(vars) == 0) stop("load() found no objects")
      get(vars[1], envir = env)
    }, error = function(e2) {
      stop(sprintf("readRDS failed (%s) AND load() failed (%s)", conditionMessage(e), conditionMessage(e2)))
    })
  })
}

inspect_panel <- function(name, path) {
  cat(sprintf("\n========== %s ==========\n", name))
  cat(sprintf("Path: %s\n", path))
  cat(sprintf("Exists: %s\n", file.exists(path)))
  if (!file.exists(path)) return(invisible(NULL))

  sobj <- tryCatch(load_obj(path), error = function(e) {
    cat("LOAD ERROR:", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(sobj)) return(invisible(NULL))

  cat("Class:", class(sobj), "\n")
  cat("n_cells:", ncol(sobj), " n_genes:", nrow(sobj), "\n")
  cat("Assays:", paste(names(sobj@assays), collapse = ", "), "\n")

  # S0.2: Raw counts check
  rna_assay <- sobj@assays$RNA
  if (!is.null(rna_assay)) {
    counts_mat <- rna_assay@counts
    data_mat   <- rna_assay@data
    n_rows <- min(100, nrow(sobj))
    n_cols <- min(10, ncol(sobj))
    cnt_max <- tryCatch(max(counts_mat[seq_len(n_rows), seq_len(n_cols)]), error = function(e) NA)
    dat_max <- tryCatch(max(data_mat[seq_len(n_rows), seq_len(n_cols)]), error = function(e) NA)
    cnt_int  <- tryCatch(all(counts_mat[seq_len(n_rows), seq_len(n_cols)] == floor(counts_mat[seq_len(n_rows), seq_len(n_cols)])), error = function(e) NA)
    cat(sprintf("S0.2 @counts max=%.3f (integers=%s) | @data max=%.3f\n", cnt_max, cnt_int, dat_max))
    if (!is.na(cnt_max) && cnt_max > 10) {
      cat("     => RAW COUNTS AVAILABLE (scVI usable)\n")
    } else {
      cat("     => WARNING: may be normalized only (scVI NOT usable)\n")
    }
  }

  # Metadata columns
  meta_cols <- colnames(sobj@meta.data)
  cat("meta.data cols:", paste(meta_cols, collapse = ", "), "\n")

  # Diagnosis columns (should match what CSV provides)
  diag_pat <- "diag|disease|condition|IBD|UC|CD|status|phenotype|tissue|inflam|sample|patient|block|slide|fov"
  diag_cols <- grep(diag_pat, meta_cols, value = TRUE, ignore.case = TRUE)
  if (length(diag_cols) > 0) {
    cat("DIAGNOSIS columns found:\n")
    for (col in diag_cols) {
      vals <- sort(unique(as.character(sobj@meta.data[[col]])))
      cat(sprintf("  [%s] (%d unique): %s\n", col, length(vals),
                  paste(head(vals, 10), collapse = " | ")))
    }
  } else {
    cat("DIAGNOSIS columns: NONE found -- check meta_cols above\n")
  }

  # Cell type columns
  ct_pat <- "cell.?type|celltype|ct_|cluster|annot|label|ident"
  ct_cols <- grep(ct_pat, meta_cols, value = TRUE, ignore.case = TRUE)
  if (length(ct_cols) > 0) {
    cat("CELLTYPE columns:\n")
    for (col in ct_cols[seq_len(min(4, length(ct_cols)))]) {
      vals <- sort(unique(as.character(sobj@meta.data[[col]])))
      cat(sprintf("  [%s] (%d types): %s\n", col, length(vals),
                  paste(head(vals, 10), collapse = " | ")))
    }
  } else {
    cat("CELLTYPE columns: NONE found\n")
  }

  # S0.4: Spatial coordinates
  cat("S0.4 Spatial coordinates:\n")
  has_xy <- all(c("x", "y") %in% meta_cols)
  cat(sprintf("  x/y in meta.data: %s\n", has_xy))
  if (has_xy) {
    x_range <- range(sobj@meta.data$x, na.rm = TRUE)
    y_range <- range(sobj@meta.data$y, na.rm = TRUE)
    cat(sprintf("  x range: [%.1f, %.1f] y range: [%.1f, %.1f]\n",
                x_range[1], x_range[2], y_range[1], y_range[2]))
  }
  coord_cols <- grep("centroid|coord|cell_x|cell_y|CenterX|CenterY|x_loc|y_loc",
                     meta_cols, value = TRUE, ignore.case = TRUE)
  if (length(coord_cols) > 0) cat("  Other coord cols:", paste(coord_cols, collapse = ", "), "\n")
  cat("  @images slots:", paste(names(sobj@images), collapse = ", "), "\n")
  cat("  @reductions:", paste(names(sobj@reductions), collapse = ", "), "\n")
  tryCatch({
    coords <- SeuratObject::GetTissueCoordinates(sobj)
    cat(sprintf("  GetTissueCoordinates: OK dims=%s cols=%s\n",
                paste(dim(coords), collapse = "x"),
                paste(colnames(coords), collapse = ", ")))
    cat(sprintf("  coord ranges: col1=[%.1f,%.1f] col2=[%.1f,%.1f]\n",
                min(coords[, 1]), max(coords[, 1]),
                min(coords[, 2]), max(coords[, 2])))
  }, error = function(e) {
    cat("  GetTissueCoordinates: ERROR -", conditionMessage(e), "\n")
  })

  # S0.3: noseg cell status (for noseg panel)
  if (grepl("noseg", name)) {
    cat("S0.3 (noseg) n_cells:", ncol(sobj), "(>0 means cells present -- keep for integration)\n")
  }

  rm(sobj)
  gc(verbose = FALSE)
  invisible(NULL)
}

for (nm in names(panels)) {
  inspect_panel(nm, panels[[nm]])
}

cat("\n========== DONE ==========\n")
