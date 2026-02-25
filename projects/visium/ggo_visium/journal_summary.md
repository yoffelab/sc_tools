# Journal Summary: GGO Visium Project

Condensed summary of `Journal.md` for this project. Full entries are in `Journal.md`.

## Project scope

- Lung tumor evolution and TLS transcriptomics (Visium). GGO (Ground Glass Opacity) and tumor progression stages: Normal, GGO, Part-Solid, Solid. Uses scVI for integration; significance bars and numerical p-values on all manuscript figures.

## Recent phase (2025-01–02)

- **Deconvolution:** Tangram-based (cluster mode, 500 epochs per library_id); fallback DestVI → Cell2location → Tangram. Proportions in `obsm['cell_type_proportions']` and `obs['prop_*']`. Process per library with backed mode to avoid segfaults; memory profiling and cleanup in place.
- **Gene signatures:** Seurat-style scoring in `score_gene_signatures.py`; Liron myeloid/T cell signatures in tumor_differences; signature heatmaps and versioned figures.
- **Spatial/process:** Macrophage localization (proliferative vs macrophage scores, regression, correlation heatmaps); neutrophil–cytotoxic T-cell colocalization (SLC16A3+ neutrophil vs Liron cytotoxic T-cell, scatter and Pearson per tumor type); process colocalization (Pearson, Moran's I, volcano by Normal/Non-Solid/Solid). Outputs under `figures/manuscript/` and `results/process_colocalization/`.

## Key conventions

- Comparisons: 1-vs-rest for tumor stages; significance text box format; dodged bars.
- Checkpoint names aligned with root Architecture (e.g. adata.normalized.scored.p35.h5ad after checkpoint migration).
