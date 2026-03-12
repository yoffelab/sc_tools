# Mission: Lung Tumor Evolution and TLS Transcriptomics (GGO Visium)

**Project:** `projects/visium/ggo_visium`  
**Current Status:** Testing setup — unit tests first, then sc_tools tests, then function implementation  
**Author:** Junbum Kim  
**Last Updated:** 2025-02-09

This file holds **project-specific** goals. Repository-level pipeline and phase definitions are in the `docs/Mission.md` and `README.md` (Pipeline Workflow).

---

## 1. Objective

Identify transcriptional and spatial differences across the lung tumor progression spectrum (Normal, Non-Solid, and Solid tumors) with a specific focus on the heterogeneity of Tertiary Lymphoid Structures (TLS).

---

## 2. Phase Alignment (sc_tools Phasing Scheme)

The project follows the repo phasing scheme (see `docs/Mission.md` and `README.md`). `preprocess` = preprocessing and clustering only; `scoring` = gene scoring, automated cell typing, deconvolution (separate branch from `preprocess`, parallel to `demographics`); `celltype_manual` = manual cell typing (skippable if `scoring` adequate).

| Phase | GGO Visium Status | Key Tasks |
|-------|-------------------|-----------|
| **`qc_filter`** | Done | Ingestion (cloupe→AnnData), QC raw |
| **`metadata_attach`** | Done | Clinical metadata joined (solidity, patient) |
| **`preprocess`** | Done | scVI, clustering (no cell typing in `preprocess`) |
| **`demographics`** | Partial | Demographics / Figure 1 cohort description |
| **`scoring`** | Done | Gene scoring, automated cell typing, optional deconvolution |
| **`celltype_manual`** | Done | Manual cell typing (phenotyped); was used for refinement |
| **`biology`** | Active | tumor_differences, macrophage localization, neutrophil–cytotoxic T-cell colocalization, process_colocalization, TLS, signature heatmaps, manuscript figures |
| **`meta_analysis`** | Pending | ROI/patient aggregation, meta analysis |

---

## 3. Completed Tasks

- [x] **`qc_filter`:** Data ingestion; AnnData with `sample`, spatial coords, H&E images.
- [x] **`metadata_attach`:** Grouping by pathology (Normal, Non-Solid, Solid) in `adata.obs`.
- [x] **`preprocess`:** scVI integration, Leiden clustering (preprocessing only).
- [x] **`scoring`:** Gene scoring (Seurat-based), automated cell typing, optional deconvolution; phenotyped AnnData and `adata.normalized.scored.p35.h5ad`.
- [x] **`celltype_manual`:** Manual cell typing refinement; phenotyped AnnData.
- [x] **`biology` (partial):** Differential program analysis, macrophage localization, neutrophil–cytotoxic T-cell colocalization (SLC16A3+ neutrophil vs cytotoxic T-cell), process colocalization, signature heatmaps (versioned), TLS B-cell/T-cell, ligand-receptor, manuscript spatial plots.

---

## 4. Implementation Roadmap (Current Priority)

**Order of work:**
1. **Unit tests for ggo_visium** — Implement first. Validate Makefile and scripts (`qc_filter`, `metadata_attach`, `preprocess`, `scoring`, `celltype_manual`, `biology`) run correctly with fixtures. Establish baseline before refactoring.
2. **Unit tests for sc_tools** — Implement second. Ensure guaranteed behavior of sc_tools functions.
3. **Implement functions** — Refactor scripts to use sc_tools; new code must compile and pass both test layers.

Pipeline phasing is defined in `docs/Mission.md` (Phasing scheme) and `README.md` (diagram).

---

## 5. Active Tasks (Roadmap)

### Unit Tests (Next)
- [ ] **ggo_visium project tests:** Create `projects/visium/ggo_visium/tests/`. Integration/smoke tests for Makefile targets and scripts. Use fixtures (minimal h5ad, metadata). `qc_filter` through `preprocess` not fully tested yet; `celltype_manual` and `biology` should work (or with few fixes) with legacy code.
- [ ] **sc_tools package tests:** After ggo_visium tests; create `sc_tools/tests/`. Unit tests for pl, tl, qc with synthetic AnnData.

### `demographics`
- [ ] **Figure 1:** Cohort stats (piechart, violin, bar, heatmap) for manuscript.

### `biology`: Downstream Biology
- [x] **Deconvolution (Tangram):** Tangram run successful (29,952 spots x 31 cell types) via `sc_tools.tl.deconvolution()`. Output: `results/adata.deconvolution.tangram.h5ad`.
- [x] **Deconvolution (Cell2location):** Completed on CPU using reference_profiles shortcut (29,952 spots x 31 cell types). Output: `results/adata.deconvolution.cell2location.h5ad`. Spatial plots in `figures/deconvolution/{method}/` (PDF + per-library PNGs at 300 DPI).
- [ ] **Gene Signature Refinement:** Update `projects/visium/ggo_visium/metadata/gene_signatures.json`; validate against HLCA, MSigDB, TCGA LUAD.
- [ ] **Spatially Variable Genes (SVG):** Moran's I / `squidpy.gr.spatial_autocorr`; spatial plots.
- [ ] **Cell Type Colocalization:** Pearson, Moran's I, neighborhood enrichment on deconvolution proportions.
- [ ] **Spatial Transition Areas:** Transcriptional changes (Normal↔Non-Solid↔Solid; tumor core vs TLS).
- [ ] **Macrophage State Comparison:** 1-vs-all testing across tumor types.
- [x] **Neutrophil–cytotoxic T-cell colocalization:** SLC16A3+ neutrophil vs Liron cytotoxic T-cell score; scatter and Pearson per tumor type (script: neutrophil_cytotoxic_tcell_localization.py).
- [ ] **Differential Cell Type Proportions:** 1-vs-all on deconvolution proportions.
- [ ] **TLS Niche Extraction:** Subset lymphoid-rich neighborhoods; spatial TLS distribution.

### `meta_analysis` (Optional)
- [ ] Aggregate to ROI and patient level.
- [ ] Downstream analysis on aggregated data.

---

## 6. Blockers and Sanity Checks

- **Sanity Check:** Do not over-filter low-count Non-Solid (GGO) spots.
- **Statistical Correction:** Benjamini-Hochberg (FDR) for all multiple comparisons.
- **Memory Management:** Deconvolution batch per `library_id` to avoid segmentation faults.

---

## 7. Technical Decision Log (Reference this project's Journal.md)

- **scVI Latent Dimensions:** n=30.
- **Comparison Mode:** Default 1-vs-all for tumor types.
- **Gene Signature Scoring:** Seurat-based (Scanpy `score_genes`).
- **Deconvolution:** Fallback DestVI → Cell2location → Tangram; batch per library_id.
- **Process Colocalization:** Pearson, Moran's I, faceted volcano; signature filtering.
- **Signature Heatmaps:** 3-level hierarchy; solidity order Normal → Non-Solid → Solid.
