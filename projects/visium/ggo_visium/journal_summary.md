# Journal Summary: GGO Visium Project

Condensed summary of `Journal.md` for this project. Full entries are in `Journal.md`.

## Project scope

- Lung tumor evolution and TLS transcriptomics (Visium). GGO (Ground Glass Opacity) and tumor progression stages: Normal, GGO, Part-Solid, Solid. Uses scVI for integration; significance bars and numerical p-values on all manuscript figures.

## Recent phase (2025-01–02)

- **Deconvolution (2026-02-27):** Generic `sc_tools.tl.deconvolution()` with backend registry. Both Tangram (29,952 spots x 31 cell types) and Cell2location completed successfully (CPU mode with reference_profiles shortcut). Outputs: `results/adata.deconvolution.{tangram,cell2location}.h5ad`. Snakefile Phase 3.5b updated with configurable `deconv_method` and figure rules. Spatial plots: per-library PNGs at 300 DPI in `figures/deconvolution/{method}/`; improved spot sizing and contrast (98th percentile vmax).
- **Gene signatures:** Seurat-style scoring in `score_gene_signatures.py`; Liron myeloid/T cell signatures in tumor_differences; signature heatmaps and versioned figures.
- **Spatial/process:** Macrophage localization (proliferative vs macrophage scores, regression, correlation heatmaps); neutrophil–cytotoxic T-cell colocalization (SLC16A3+ neutrophil vs Liron cytotoxic T-cell, scatter and Pearson per tumor type); process colocalization (Pearson, Moran's I, volcano by Normal/Non-Solid/Solid). Outputs under `figures/manuscript/` and `results/process_colocalization/`.

## Key conventions

- Comparisons: 1-vs-rest for tumor stages; significance text box format; dodged bars.
- Checkpoint names: semantic slug names adopted (adata.annotated.h5ad, adata.scored.h5ad, adata.celltyped.h5ad). Backwards compat with old p-code names via `_old` fallback vars in all scripts.
