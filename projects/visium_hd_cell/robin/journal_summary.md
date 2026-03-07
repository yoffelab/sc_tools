# Journal Summary: robin (visium_hd_cell)

Condensed summary of `Journal.md`. Full entries in Journal.md.

## Project scope

Single-cell Visium HD analysis (Robin) using SpaceRanger 4 cell segmentation. Parallel to bin-level project at `projects/visium_hd/robin/`. 15 samples, 1,574,694 cells x 18,132 genes.

## Recent phase (2026-03-06)

- **Phase 0b ingest COMPLETE:** All 15 samples loaded. Fixed `load_visium_hd_cell_sample` to handle SR4 directory layout (`segmented_outputs/`, `filtered_feature_cell_matrix.h5`) and extract cell centroids from `cell_segmentations.geojson` with `cellid_XXXXXXXXX-1` index format. Pat5_Samp3 excluded (QC fail). Per-sample checkpoints at `data/{sample_id}/adata.h5ad`. Concatenated at `results/adata.raw.h5ad` (4.1 GB).

## Key conventions

- Checkpoint names: `adata.raw.h5ad`, `adata.normalized.h5ad`, `adata.celltyped.h5ad`
- Uses xenium preprocessing recipe (single-cell resolution)
- Shares SR4 runs with `projects/visium_hd/robin/` (same raw data, different output modality)
