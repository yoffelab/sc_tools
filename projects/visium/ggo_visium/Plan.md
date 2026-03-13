# Plan: GGO Visium — Lung Tumor Evolution and TLS Transcriptomics

Identify transcriptional and spatial differences across the lung tumor progression spectrum (Normal, Non-Solid, Solid) with focus on Tertiary Lymphoid Structure (TLS) heterogeneity. Visium spatial transcriptomics, 29,952 spots, 31 cell types.

---

## Phase Status

- [x] `qc_filter` — Ingestion (cloupe to AnnData), QC raw
- [x] `metadata_attach` — Clinical metadata joined (solidity, patient)
- [x] `preprocess` — scVI integration, Leiden clustering
- [ ] `demographics` — Partial; Figure 1 cohort description pending
- [x] `scoring` — Gene scoring, automated cell typing, deconvolution
- [x] `celltype_manual` — Manual cell typing refinement
- [ ] `biology` — Active; deconvolution done, downstream analyses in progress
- [ ] `meta_analysis` — Pending

---

## Tasks

### Testing (priority order)

- [ ] ggo_visium project integration tests (`tests/`): fixtures, smoke tests for Makefile targets and scripts across all phases
- [ ] sc_tools package unit tests (`sc_tools/tests/`): pl, tl, qc with synthetic AnnData
- [ ] Refactor scripts to use sc_tools; all code must pass both test layers

### Demographics

- [ ] Figure 1: cohort stats (piechart, violin, bar, heatmap) for manuscript

### Biology (active)

- [x] Deconvolution — Tangram (29,952 x 31 cell types): `results/adata.deconvolution.tangram.h5ad`
- [x] Deconvolution — Cell2location (CPU, reference_profiles shortcut): `results/adata.deconvolution.cell2location.h5ad`
- [x] Neutrophil-cytotoxic T-cell colocalization (SLC16A3+ neutrophil vs cytotoxic T-cell score)
- [ ] Gene signature refinement: validate against HLCA, MSigDB, TCGA LUAD
- [ ] Spatially variable genes (SVG): Morans I / squidpy spatial_autocorr
- [ ] Cell type colocalization: Pearson, Morans I, neighborhood enrichment on deconvolution proportions
- [ ] Spatial transition areas: transcriptional changes (Normal to Non-Solid to Solid; tumor core vs TLS)
- [ ] Macrophage state comparison: 1-vs-all testing across tumor types
- [ ] Differential cell type proportions: 1-vs-all on deconvolution proportions
- [ ] TLS niche extraction: subset lymphoid-rich neighborhoods, spatial TLS distribution

### Meta-Analysis (optional)

- [ ] Aggregate to ROI and patient level
- [ ] Downstream analysis on aggregated data

---

## Blocked / Sanity Checks

- Do not over-filter low-count Non-Solid (GGO) spots
- Benjamini-Hochberg (FDR) for all multiple comparisons
- Deconvolution batch per `library_id` to avoid OOM / segfaults

---

## Technical Decisions

- scVI latent dimensions: n=30
- Comparison mode: 1-vs-all for tumor types (default)
- Gene signature scoring: Seurat-based (Scanpy `score_genes`)
- Deconvolution fallback chain: DestVI > Cell2location > Tangram; batch per library_id
- Process colocalization: Pearson, Morans I, faceted volcano; signature filtering
- Signature heatmaps: 3-level hierarchy; solidity order Normal > Non-Solid > Solid
- Statistical output: significance text box in top-right corner; dodged bars for multi-group

---

## Key Files

| Script | Output |
|--------|--------|
| `scripts/tumor_differences.py` | `figures/manuscript/tumor_differences/` |
| `scripts/macrophage_localization.py` | `figures/manuscript/macrophage_localization/` |
| `scripts/neutrophil_cytotoxic_tcell_localization.py` | `figures/manuscript/neutrophil_cytotoxic_tcell_localization/` |
| `scripts/process_colocalization.py` | `figures/process_colocalization/`, `results/process_colocalization/` |
| `scripts/signature_heatmap_versioned.py` | `figures/manuscript/signature_heatmaps_versioned.done` |
| `scripts/tls_analysis.py` | `figures/manuscript/tls_analysis/` |

**Snakemake targets:** `phase35b`, `localization_only`, `signature_heatmaps_versioned`
