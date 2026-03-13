# Spatial Omics QC — Agent Skill

**Format:** K-Dense-AI claude-scientific-skills
**Version:** 1.0
**Scope:** sc_tools pipeline, phase `qc_filter`

---

## 1. Overview

Spatial omics QC differs fundamentally from single-cell RNA-seq (scRNA-seq) QC because observations are spatially registered. Each observation (a spot, a bin, a cell, or an ROI) carries a physical coordinate that must be preserved through filtering. Aggressive filtering that is acceptable in scRNA-seq — where cells are randomly distributed — can destroy spatial structure in tissue sections, removing entire tissue compartments rather than just low-quality cells.

The core goals of spatial QC are:

- Remove technically failed observations (empty spots, debris, zero-count cells, segmentation artifacts).
- Preserve biologically valid low-count observations that may represent real tissue biology (e.g., sparse stromal areas, acellular zones adjacent to tumors).
- Flag or remove samples that fail platform-level quality thresholds before concatenation.
- Produce a `results/adata.filtered.h5ad` that satisfies the metadata contract for the `qc_filter` checkpoint (see Architecture.md §2.2).

Unlike scRNA-seq, spatial QC also involves spatially variable gene (SVG) detection, cross-sample batch-effect visualization anchored in tissue space, and platform-specific normalization choices. The choice of normalization transform must be locked at QC time and carried forward consistently: transcriptomic platforms use `log1p` (or raw for VAE-based models), while IMC uses `arcsinh(x/5)`.

---

## 2. When to Apply This Skill

Apply this skill during the `qc_filter` pipeline phase. This phase runs after `ingest_load` (per-sample AnnData objects exist at `data/{sample_id}/adata.ingested.h5ad`) and before `metadata_attach` or `preprocess`.

Supported modalities and their platforms:

| Modality slug | Platform | Cell / observation unit |
|---------------|----------|------------------------|
| `visium` | 10x Visium | ~55 µm spot (~10-50 cells) |
| `visium_hd` | 10x Visium HD | 8 µm bin (1-few cells) |
| `visium_hd_cell` | 10x Visium HD + SpaceRanger 4 | Single cell (segmented) |
| `xenium` | 10x Xenium | Single cell (segmented) |
| `imc` | Imaging Mass Cytometry (ElementoLab/steinbock) | Single cell (segmented) |
| `cosmx_1k` / `cosmx_6k` / `cosmx_full_library` | NanoString CosMx | Single cell (segmented) |

Do NOT apply deconvolution QC steps to `visium_hd_cell`, Xenium, CosMx, or IMC — those modalities produce single-cell resolution data directly and deconvolution is skipped.

---

## 3. Key QC Metrics by Modality

These are the primary metrics computed by `sc_tools.qc.calculate_qc_metrics()`. Thresholds are defaults; project-specific values override them when documented in `metadata/sample_metadata.csv` or the project Mission.md.

### 3.1 Visium (and Visium HD bin-level)

| Metric | Key | Default threshold | Action if violated |
|--------|-----|-------------------|--------------------|
| Total UMI counts per spot | `total_counts` | > 500 (Visium); > 100 (Visium HD) | Remove spot |
| Number of genes detected | `n_genes_by_counts` | > 200 | Remove spot |
| Mitochondrial fraction | `pct_counts_mt` | < 20% | Remove spot if > 20% |
| Number of spots per sample | `n_spots` | > 1000 (Visium); project-specific for HD | Flag sample |

For Visium HD, the 8 µm bin size means many bins will have lower counts than a 55 µm Visium spot; lower `total_counts` thresholds are expected and appropriate.

Mitochondrial genes are identified by the pattern `^MT-` (human) or `^mt-` (mouse). Filter patterns for ribosomal (`^RP[SL]`) and hemoglobin (`^HB[^P]`) genes from the HVG set but do not remove them from `adata.X` by default.

### 3.2 Xenium

| Metric | Key | Default threshold | Action if violated |
|--------|-----|-------------------|--------------------|
| Transcript count per cell | `transcript_count` | > 5 | Remove cell |
| Cell area (µm²) | `cell_area` | 20–2000 µm² | Remove if outside range |
| Nucleus-to-cell area ratio | `nucleus_ratio` | > 0.1 | Flag as potential segmentation failure |
| Mitochondrial fraction | `pct_counts_mt` | < 15% | Remove cell if > 15% |

Xenium panels are targeted (hundreds to low thousands of genes), so absolute transcript counts are lower than scRNA-seq. The 15% MT threshold (rather than 20%) reflects the panel design: MT transcripts are usually not panel targets, so a high MT fraction signals cytoplasmic contamination or poor segmentation more than metabolic state.

Cell segmentation quality is determined upstream by Xenium Ranger before QC runs. If segmentation is poor (large median `cell_area`, many cells at boundary thresholds), flag the sample for human review rather than silently dropping cells.

### 3.3 IMC (Imaging Mass Cytometry)

| Metric | Key | Default threshold | Action if violated |
|--------|-----|-------------------|--------------------|
| Number of proteins detected | `n_proteins_detected` | > 5 (of panel) | Remove cell |
| Mean protein intensity (arcsinh) | `mean_intensity` | > 0.1 (arcsinh scale) | Flag |
| Cell area (µm²) | `cell_area` | 10–500 µm² | Remove if outside range |
| Signal-to-noise ratio | `signal_to_noise` | > 1.5 | Flag ROI if median SNR below threshold |

IMC uses `arcsinh(x/5)` normalization — not `log1p`. This is non-negotiable. Never apply `normalize_total` or `log1p` to IMC protein intensity data. The cofactor of 5 is standard for CyTOF / IMC mass spectrometry data. Intensities below 0 (instrument noise) are valid after arcsinh and should not be clipped.

IMC segmentation quality is the most critical upstream factor. Poor ilastik probability maps or steinbock segmentation produce cells with unrealistic shapes. Visual inspection of `adata.uns['spatial'][sample_id]['images']['mask']` is recommended before proceeding.

### 3.4 CosMx (1k, 6k, full_library)

| Metric | Key | Default threshold | Action if violated |
|--------|-----|-------------------|--------------------|
| Transcript count per cell | `transcript_count` | > 10 | Remove cell |
| Cell area (µm²) | `cell_area` | 20–5000 µm² | Remove if outside range |
| Number of genes detected | `n_genes` | > 5 | Remove cell |

CosMx thresholds vary by panel tier. The 1k panel detects fewer transcripts per cell than the 6k or full_library panels; set lower absolute count thresholds for 1k data. The NanoString detection algorithm changed between 2023 and 2024 software versions; verify which version produced the data before comparing across datasets.

For `cosmx_full_library` (~18k genes), `adata.X` may be stored as a backed sparse array on disk. Keep `adata.X` sparse throughout QC; do not densify.

---

## 4. Spatially Variable Gene (SVG) Detection

SVG detection identifies genes whose expression patterns are spatially structured in tissue — i.e., genes that are not randomly distributed but are enriched in particular spatial domains (tumor cores, stromal margins, immune infiltrates, etc.).

When to run SVG detection:
- After initial QC filtering, on the per-sample `adata.ingested.h5ad` objects (before concatenation), OR
- After concatenation on `results/adata.filtered.h5ad` with `sample` as a grouping key.
- Always run before the `preprocess` phase to inform HVG selection.

sc_tools function: `sc_tools.qc.spatially_variable_genes(adata, method="moran", n_jobs=-1)`

This wraps `squidpy.gr.spatial_autocorr(mode="moran")` with spatial graph construction via `squidpy.gr.spatial_neighbors()`. The output adds `adata.var['spatially_variable']` (bool) and `adata.var['moranI']` (float). Genes with high Moran I scores are spatially structured.

For the `preprocess` phase HVG selection, intersect HVGs with SVGs to prioritize genes that are both variable and spatially informative. This is modality-specific:
- Visium / Visium HD: always intersect HVG with SVG.
- Xenium / CosMx: SVG detection is optional; targeted panels already restrict the feature space.
- IMC: SVG detection applies to proteins, not genes. With 30-50 proteins, all markers are typically included; no HVG/SVG intersection step.

Do not run SVG detection on data that has been normalized or log-transformed — run on raw counts for transcriptomic modalities.

---

## 5. Multi-Sample QC Comparison and Batch Effect Detection

After per-sample QC, concatenate all samples and perform cross-sample comparison before proceeding.

Key checks to run on the concatenated `adata.filtered.h5ad`:

**Distribution comparison:** Plot `total_counts`, `n_genes_by_counts`, and `pct_counts_mt` as violin or box plots grouped by `obs['sample']`. Samples with dramatically different distributions (median counts 2x or more above/below cohort median) are potential batch outliers.

**Spatial %MT visualization:** Plot `pct_counts_mt` in spatial coordinates (`obsm['spatial']`) for each sample. High MT% concentrated at tissue edges (not uniformly distributed) indicates mechanical damage or poor tissue quality rather than high metabolic state — do not remove these spots blindly; flag for human review.

**Cell/spot count per sample:** Flag samples with fewer than 1000 spots/cells as potentially failed. Do not automatically exclude; let `sample_qc.classify_samples()` produce a `pass`/`flag`/`fail` label for human review.

**Cross-sample gene detection rate:** For targeted panels (Xenium, CosMx), check that the same genes are detected at comparable rates across samples. A gene detected in <10% of cells in most samples but >50% in one sample suggests a platform artifact.

When to flag vs. remove samples:
- `flag`: sample metrics are outside normal range but within 2x of cohort median; keep but note.
- `fail`: sample metrics are extreme (>5x deviation from median), or the sample has <500 cells/spots; remove from downstream analysis.
- Document every removed sample in `Journal.md` with the reason.

---

## 6. QC Report Outputs

The pipeline produces four date-versioned HTML reports saved to `figures/QC/`. All filenames use the `YYYYMMDD` format for auditability. Reports are generated by `sc_tools.qc.report.generate_qc_report()`.

| Report | Filename | Phase trigger | Primary content |
|--------|----------|---------------|-----------------|
| 1 | `pre_filter_qc_YYYYMMDD.html` | `qc_filter` entry | Per-sample spot/cell counts, gene/protein detection histograms, %MT distributions; raw data before any filtering |
| 2 | `post_filter_qc_YYYYMMDD.html` | `qc_filter` exit / `metadata_attach` | Pre-vs-post comparison, HVG/SVG summary, per-sample pass/flag/fail labels |
| 3 | `post_integration_qc_YYYYMMDD.html` | `preprocess` exit | Batch score (primary metric), UMAP grid, cluster distributions, integration benchmark; bio metrics informational only at this stage |
| 4 | `post_celltyping_qc_YYYYMMDD.html` | `celltype_manual` exit | Full bio and batch scores with validated celltypes; re-evaluation of all candidate integration methods |

Key principle: Report 3 uses batch score as the primary integration selection criterion. Bio metrics (ARI, NMI, ASW celltype) are reported but must not drive method selection at this stage because celltype labels are preliminary leiden clusters, not validated annotations. Bio metrics become the primary evaluation criterion only after celltyping is complete (Report 4). Selecting an integration method based on preliminary bio metrics introduces circular reasoning.

---

## 7. sc_tools Functions Reference

All functions below are in the `sc_tools.qc` subpackage unless otherwise noted.

**`sc_tools.qc.calculate_qc_metrics(adata, modality, mt_pattern="^MT-")`**
Computes per-observation QC metrics appropriate for the modality. Adds columns to `adata.obs`: `total_counts`, `n_genes_by_counts`, `pct_counts_mt` (transcriptomic) or `n_proteins_detected`, `mean_intensity`, `cell_area` (IMC / CosMx). Wraps `scanpy.pp.calculate_qc_metrics()` with modality-aware defaults.

**`sc_tools.qc.filter_cells(adata, modality, thresholds=None)`**
Applies per-observation filtering. `thresholds` dict overrides defaults per metric. Returns filtered AnnData. Does not modify in place by default (returns a copy). Always preserves `obsm['spatial']` alignment.

**`sc_tools.qc.filter_genes(adata, min_cells=10, exclude_patterns=None)`**
Removes genes detected in fewer than `min_cells` observations. `exclude_patterns` is a list of regex patterns (e.g. `["^MT-", "^RP[SL]", "^HB[^P]"]`) to exclude from HVG candidates but NOT from the feature matrix (genes are flagged in `adata.var`, not dropped).

**`sc_tools.qc.sample_qc.classify_samples(adata, groupby="sample", thresholds=None)`**
Computes per-sample aggregate QC metrics and assigns a `qc_status` label (`pass`, `flag`, `fail`) to `adata.obs`. Uses `adata.obs['sample']` as the grouping key. Outputs a summary table to `outputs/sample_qc_summary.csv`.

**`sc_tools.qc.spatially_variable_genes(adata, method="moran", n_jobs=-1)`**
Detects SVGs using spatial autocorrelation. Requires `obsm['spatial']` to be set. Adds `adata.var['spatially_variable']` and `adata.var['moranI']`. Must be run on raw counts.

**`sc_tools.qc.report.generate_qc_report(adata, phase, output_path, metadata=None)`**
Generates the HTML QC report for the given phase (`pre_filter`, `post_filter`, `post_integration`, `post_celltyping`). Output path should follow the naming convention `figures/QC/{phase}_qc_YYYYMMDD.html`.

**`sc_tools.qc.highly_variable_genes(adata, n_top_genes=3000, batch_key="sample", flavor="seurat_v3")`**
Wraps `scanpy.pp.highly_variable_genes()` with sc_tools defaults. Use `flavor="seurat_v3"` for raw counts input (scVI path); use `flavor="seurat"` on log-normalized data. For IMC, this function is not applicable — all proteins are retained.

---

## 8. HPC / SLURM Notes

Per-sample QC runs as a SLURM job array — one task per sample — before the concatenation step. Concatenation runs as a single downstream job with `--dependency=afterok` on the array job.

Resource guidelines per modality:

| Modality | CPUs | RAM | Walltime |
|----------|------|-----|----------|
| Visium | 4-8 | 32-64 GB | 1h |
| Visium HD | 8-16 | 64-128 GB | 2h |
| Xenium | 4-8 | 32-64 GB | 1h |
| IMC | 4 | 16-32 GB | 30m |
| CosMx 1k/6k | 8 | 64 GB | 2h |
| CosMx full_library | 16 | 128 GB | 4h |

MCP tools for generating platform-specific SLURM scripts:

- `mcp__sc-tools__generate_sbatch_spaceranger` — generates Space Ranger SLURM scripts for Visium / Visium HD `ingest_raw` jobs
- `mcp__sc-tools__generate_sbatch_xenium` — generates Xenium Ranger SLURM scripts for Xenium `ingest_raw` jobs
- `mcp__sc-tools__generate_sbatch_imc` — generates IMC pipeline SLURM scripts for IMC `ingest_raw` jobs

These MCP tools produce ready-to-submit sbatch files with correct `--array`, `--cpus-per-task`, `--mem`, and `--gres` parameters for the registered HPC servers (brb, cayuga). After generating, review the `--gres` line: on cayuga, GPU type must be explicit (e.g. `--gres=gpu:a40:1`); generic `--gres=gpu:1` fails on cayuga.

Always write SLURM outputs to scratch (`/athena/elementolab/scratch/juk4007/`), never to login node home directories.

---

## 9. Common Pitfalls

### IMC: Never use log1p

IMC protein intensities must be transformed with `arcsinh(x/5)`, not `log1p` or `normalize_total`. Applying `log1p` to mass spectrometry intensities produces biologically meaningless results and breaks downstream CytoVI integration. The cofactor 5 is standard for IMC / CyTOF; do not change it unless the project explicitly documents a different cofactor with justification.

### Xenium: Upstream segmentation quality drives downstream QC

If Xenium Ranger cell segmentation is poor, no amount of per-cell QC filtering will rescue the data. Check segmentation metrics (median cell area, fraction of cells near area bounds) before filtering. A bimodal `cell_area` distribution usually indicates mixed segmentation quality (some cells segmented correctly, others merged or split). Flag the sample for visual inspection rather than filtering blindly.

### Visium: Do not remove low-count spots near tissue boundaries

The tissue edge and regions near adipose, necrotic, or acellular areas often have lower counts. These are biologically real. Apply a minimal absolute threshold (e.g. `total_counts > 100`) to exclude truly empty spots (background), but do not use a high threshold that removes real tissue observations. Avoid aggressive automated filtering.

### Visium HD: Expect many zero-count bins

At 8 µm resolution, a large fraction of bins will overlap extracellular matrix or empty space between cells and will have zero or near-zero counts. This is expected and normal. Filter at a very low threshold (e.g. `total_counts > 10`) rather than the Visium default of 500.

### CosMx: Detection algorithm version matters

The NanoString/10x CosMx detection algorithm changed between 2023 and 2024 software versions. Datasets processed with different versions are not directly comparable and may require different QC thresholds. Always check the `software_version` field in the CosMx metadata before setting thresholds.

### Spatial coordinates must survive filtering

Every call to `filter_cells()` or `adata[mask]` must preserve `obsm['spatial']` alignment. AnnData handles this automatically when subsetting with boolean masks on `obs_names`, but explicit index manipulation (e.g. modifying `adata.obs_names` directly) can break alignment. Always subset via `adata[mask, :]` syntax, never by direct matrix slicing.

### Do not normalize before QC

QC filtering must run on raw counts (`adata.X` = integer counts). Normalization (`normalize_total`, `log1p`, `arcsinh`) is applied in the `preprocess` phase, after `qc_filter` completes. Running QC on normalized data produces incorrect MT fraction estimates and breaks the metadata contract for the `qc_filter` checkpoint.

### Multi-sample batch effects at QC time are expected

Cross-sample differences in `total_counts` distributions observed at QC time do not necessarily indicate technical failure. Some batch effect is expected; integration methods in `preprocess` will address it. At QC time, only flag samples that are extreme outliers (>5x deviation from cohort median) or that clearly failed platform processing.
