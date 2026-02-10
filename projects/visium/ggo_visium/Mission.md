# Mission: Lung Tumor Evolution and TLS Transcriptomics (GGO Visium)

**Project:** `projects/visium/ggo_visium`  
**Current Status:** Testing setup — unit tests first, then sc_tools tests, then function implementation  
**Author:** Junbum Kim  
**Last Updated:** 2025-02-09

This file holds **project-specific** goals. Repository-level pipeline and phase definitions are in the root `Mission.md` and `WORKFLOW.md`.

---

## 1. Objective

Identify transcriptional and spatial differences across the lung tumor progression spectrum (Normal, Non-Solid, and Solid tumors) with a specific focus on the heterogeneity of Tertiary Lymphoid Structures (TLS).

---

## 2. Phase Alignment (New Schema)

The project uses the 7-phase workflow. Entry point for GGO Visium: Phases 1–3 largely complete; focus on Phases 4–5 and optional 6–7.

| New Phase | GGO Visium Status | Key Tasks |
|-----------|-------------------|-----------|
| **1** | Done | Ingestion (cloupe→AnnData), QC raw |
| **2** | Done | Clinical metadata joined (solidity, patient) |
| **3** | Done | scVI, clustering, automated typing |
| **3.5** | Partial | Demographics / Figure 1 cohort description |
| **4** | Done | Manual cell typing (phenotyped) |
| **5** | Active | Gene scoring, deconvolution, tumor_differences, process_colocalization, TLS, macrophage |
| **6–7** | Pending | ROI/patient aggregation, meta analysis |

---

## 3. Completed Tasks

- [x] **Phase 1:** Data ingestion; AnnData with `sample`, spatial coords, H&E images.
- [x] **Phase 2:** Grouping by pathology (Normal, Non-Solid, Solid) in `adata.obs`.
- [x] **Phase 3:** scVI integration, Leiden clustering, automated cell typing.
- [x] **Phase 4:** Manual cell typing; phenotyped AnnData.
- [x] **Phase 5 (partial):** Gene scoring (Seurat-based), differential program analysis, macrophage localization, process colocalization, signature heatmaps, TLS B-cell/T-cell, ligand-receptor.

---

## 4. Implementation Roadmap (Current Priority)

**Order of work:**
1. **Unit tests for ggo_visium** — Implement first. Validate Makefile and scripts (Phases 1–5) run correctly with fixtures. Establish baseline before refactoring.
2. **Unit tests for sc_tools** — Implement second. Ensure guaranteed behavior of sc_tools functions.
3. **Implement functions** — Refactor scripts to use sc_tools; new code must compile and pass both test layers.

---

## 5. Active Tasks (Roadmap)

### Unit Tests (Next)
- [ ] **ggo_visium project tests:** Create `projects/visium/ggo_visium/tests/`. Integration/smoke tests for Makefile targets and scripts. Use fixtures (minimal h5ad, metadata). Phases 1–3 not fully tested yet; Phases 4–5 should work (or with few fixes) with legacy code.
- [ ] **sc_tools package tests:** After ggo_visium tests; create `sc_tools/tests/`. Unit tests for pl, tl, qc with synthetic AnnData.

### Phase 3.5: Demographics
- [ ] **Figure 1:** Cohort stats (piechart, violin, bar, heatmap) for manuscript.

### Phase 5: Downstream Biology
- [ ] **Deconvolution:** Cell2location, DestVI with batch processing; identify memory limits.
- [ ] **Gene Signature Refinement:** Update `projects/visium/ggo_visium/metadata/gene_signatures.json`; validate against HLCA, MSigDB, TCGA LUAD.
- [ ] **Spatially Variable Genes (SVG):** Moran's I / `squidpy.gr.spatial_autocorr`; spatial plots.
- [ ] **Cell Type Colocalization:** Pearson, Moran's I, neighborhood enrichment on deconvolution proportions.
- [ ] **Spatial Transition Areas:** Transcriptional changes (Normal↔Non-Solid↔Solid; tumor core vs TLS).
- [ ] **Macrophage State Comparison:** 1-vs-all testing across tumor types.
- [ ] **Differential Cell Type Proportions:** 1-vs-all on deconvolution proportions.
- [ ] **TLS Niche Extraction:** Subset lymphoid-rich neighborhoods; spatial TLS distribution.

### Phase 6–7: Meta Analysis (Optional)
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
