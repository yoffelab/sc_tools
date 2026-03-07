# Research Journal: robin (visium_hd_cell)

**Project:** `projects/visium_hd_cell/robin`

---

## Log Entries

### 2026-03-06 — Project created and Phase 0b ingest complete

- **Action:** Created project via `create_project.sh robin visium_hd_cell`. Wrote `scripts/ingest_phase0b.py` and `scripts/ingest_phase0b.sbatch` for cell segmentation ingest.
- **Loader fixes (`sc_tools/ingest/loaders.py`):**
  1. Added search for `segmented_outputs/` (SR4 default) in addition to `cell_segmentation/` (legacy).
  2. Added search for `filtered_feature_cell_matrix.h5` (SR4) in addition to `filtered_feature_bc_matrix.h5` (legacy).
  3. Added `_extract_centroids_from_geojson()` to compute cell polygon centroids from `cell_segmentations.geojson`.
  4. Fixed index mismatch: SR4 h5 barcodes use `cellid_XXXXXXXXX-1` format while geojson has integer `cell_id`. Auto-detects format from `obs_names`.
  5. Also fixed `load_visium_hd_sample()` to search `outs/binned_outputs/{bin_size}` (SR4 path).
- **Results:** 15/15 samples loaded (Pat5_Samp3 excluded, QC fail). 1,574,694 cells x 18,132 genes. Per-sample `data/{sample_id}/adata.h5ad` + `results/adata.raw.h5ad` (4.1 GB). SLURM job 13971144 on brb, 2 min runtime.
- **Binned ingest (parallel):** Also completed `projects/visium_hd/robin/` Phase 0b: 1,818,487 spots x 18,132 genes, `results/adata.raw.h5ad` (5.1 GB). Job 13971137.
