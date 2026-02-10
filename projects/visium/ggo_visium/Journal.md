# Research Journal & Decision Log: GGO Visium Project

This journal documents the technical evolution, parameters, and rationale behind **project-specific** analytical decisions for the Lung Tumor Evolution and TLS transcriptomics study. Repository-level decisions (architecture, script sanity check, scalable layout) are in the root `Journal.md`.

**Project:** `projects/visium/ggo_visium`  
**Reference:** Root `Mission.md` (toolkit); this project's `Mission.md` (study aims).

---

## Log Entries

### [2025-12-30] - Initialization of Modern Frameworks
- **Action:** Transitioned from standard log-normalization to probabilistic modeling via `scvi-tools`.
- **Rationale:** Traditional normalization was masking biological variation in low-count spots, particularly in GGO (Ground Glass Opacity) regions.
- **Decision:** Use `scvi.model.SCVI` for all future integration steps.
- **Result:** Improved separation of tumor vs. stromal clusters in the integrated UMAP.

### [2025-12-28] - Statistical Visualization Standards
- **Action:** Enforced mandatory significance bars and numerical p-values for all manuscript figures.
- **Rationale:** Maximum transparency and scientific rigor for computational oncology publications.
- **Decision:** Default to **1-vs-rest** comparisons for tumor progression stages (Normal, GGO, Part-Solid, Solid) to highlight stage-wise programs. (Dodging and comparison labeling later resolved.)

### [2025-01-02] - Phase III Cell Type Deconvolution Implementation
- **Action:** Created `celltype_deconvolution_phase3.py` for Tangram-based deconvolution; proportions integrated into spatial AnnData.
- **Decision:** Tangram cluster mode, 500 epochs per batch; process by `library_id`; store proportions in `obsm['cell_type_proportions']` and `obs['prop_*']`.
- **Output:** `results/adata.deconvolution.h5ad`.

### [2025-01-02] - Phase III Memory Optimization and Segmentation Fault Fix
- **Issue:** Segmentation fault from loading full spatial dataset at once.
- **Solution:** Process one `library_id` at a time; load spatial per library with backed mode; explicit cleanup and minimal copies.
- **Result:** Script processes libraries individually with proper cleanup.

### [2025-01-05] - Phase III Modularization: Multi-Method Deconvolution Support
- **Action:** Refactored deconvolution into modular functions with fallback (DestVI → Cell2location → Tangram).
- **Result:** Pipeline robust to method failures; method usage tracked.

### [2025-01-05] - Memory Profiling and Optimization for Deconvolution
- **Action:** Memory profiling (psutil, tracemalloc), logging to `output/deconvolution/logs/`, adaptive gene/parameter reduction, aggressive cleanup.
- **Result:** Detailed memory diagnostics; optimizations to reduce peak usage.

### [2025-01-05] - Efficiency Optimization: Remove QC Filtering from Single-Cell Data
- **Action:** Removed QC filtering step from single-cell reference loading in deconvolution script.
- **Result:** Faster loading, fewer copies; QC labels still excluded from cell type list.

### [2025-01-05] - Tumor Differences Plotting: Significance Bar Overlap and Comparison Labeling
- **Action:** Dodged significance bars; all comparisons in text box format "Group A - Rest: p*-value (p=X.XXe-XX)".
- **Result:** No bar overlap; all p-values displayed.

### [2025-01-05] - Tumor Differences Plotting: Final Refinements
- **Action:** One line per comparison; whisker protection; dynamic bottom margin.
- **Result:** Maximally readable boxplots.

### [2025-01-05] - Code Cleanup: Gene Signature Scoring Refactoring
- **Action:** Extracted scoring from `signature_cleaned.py` into `score_gene_signatures.py`; archived old file to `scripts/old_code/`.
- **Result:** Clean scoring script; scores all signatures including macrophage states for Phase IV.

### [2025-01-05] - Tumor Differences: Added Liron Myeloid and T Cell Signatures
- **Action:** Added 11 Liron signatures (6 T cell, 5 myeloid) to `tumor_differences.py`.
- **Result:** 23 total immune programs tested across tumor types.

### [2025-01-06] - Macrophage Localization: Colocalization with Proliferative Program
- **Action:** Created `macrophage_localization.py`; scatterplots (proliferative vs macrophage scores); stratification by Normal/Non-Solid/Solid; Pearson per type; annotations.
- **Output:** `figures/manuscript/macrophage_localization/`, `figures/stats/`.

### [2025-01-06] - Macrophage Localization: Regression Lines and Filtering
- **Action:** Filter proliferative and macrophage >= 0; reference lines x=0, y=0; regression lines with confidence intervals; layering fixes.
- **Result:** Clear trend lines and excluded points in gray.

### [2025-01-06] - Macrophage Localization: Confidence Interval Fix
- **Action:** Dual filtering (proliferative and macrophage >= 0); robust regression validation; confidence intervals for all tumor types.
- **Result:** All types show regression with CI when valid.

### [2025-01-05] - Macrophage Localization: CI Formula and Visibility
- **Fix:** Corrected CI formula (confidence vs prediction interval); increased alpha for visibility.
- **Result:** CI bands correct and visible for Normal, Non-Solid, Solid.

### [2025-01-05] - Macrophage Localization: Removed Confidence Intervals
- **Action:** Removed CI bands; regression lines only for clarity.
- **Result:** Cleaner plots.

### [2025-01-05] - Gene Signature Scoring: Transition to Seurat-based Method
- **Action:** Switched from simple mean to Scanpy `score_genes` (Seurat AddModuleScore) in `score_gene_signatures.py`.
- **Rationale:** Control for technical variation and expression magnitude; standard practice.
- **Result:** More robust signature scores.

### [2025-01-20] - Process Colocalization Analysis: Modular Implementation
- **Action:** Created `process_colocalization.py` with Pearson correlation, Moran's I (manual), optional neighborhood enrichment, volcano plots (Normal vs Non-Solid vs Solid).
- **Output:** Correlation heatmap, Moran's I heatmap, volcano plots, CSVs in `figures/process_colocalization/` and `results/process_colocalization/`.

### [2025-01-20] - Process Colocalization: Moran's I Fix and Volcano Enhancements
- **Fix:** Manual Moran's I for signature scores in `obs`; volcano annotations for significant processes.
- **Result:** Moran's I correct; volcano plots with labeled significant processes.

### [2025-01-20] - Process Colocalization: Faceted Volcano and Signature Filtering
- **Action:** Single faceted figure (3 subplots); VOLCANO_SIGNATURES_INCLUDE/EXCLUDE; correlation clustermaps.
- **Result:** Flexible signature selection; clustermaps for process groupings.

### [2025-01-20] - Signature Heatmap Visualization
- **Action:** Created `signature_heatmap.py`; rows = spots (patient/solidity sorted), columns = signatures (clustered); 3-level hierarchy; 4 figures (2 heatmaps, 2 clustermaps).
- **Output:** `figures/manuscript/signature_heatmaps/`.

### [2025-01-20] - Macrophage Localization: Correlation Heatmaps and Clustermaps
- **Action:** Correlation matrix (macrophage signatures × tumor types); heatmap and clustermap with significance stars.
- **Output:** `correlation_heatmap.pdf`, `correlation_clustermap.pdf` in `figures/manuscript/macrophage_localization/`.

---

## Future Rationale & Pending Decisions

### TLS Transcriptional Heterogeneity
- **Pending:** Cell2location vs DestVI for final TLS deconvolution. Tangram batches as baseline; Visium HD may require SpatialData for memory.

### Gene Signature Refinement
- **Sanity Check:** Verify signatures do not overlap heavily with housekeeping genes; consider HLCA, MSigDB, TCGA LUAD.
