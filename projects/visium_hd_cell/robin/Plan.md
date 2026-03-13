# Plan: robin (visium_hd_cell)

Single-cell resolution spatial transcriptomics of Visium HD data using SpaceRanger 4 cell segmentation. 15 samples, 1,574,694 cells x 18,132 genes. Parallel to bin-level analysis in `projects/visium_hd/robin/`.

---

## Status

- [x] `ingest_load` — 15/15 samples loaded, Pat5_Samp3 excluded (QC fail)
- [ ] `qc_filter` — script ready, awaiting SLURM submission
- [ ] `metadata_attach` — script ready, clinical metadata from sibling project
- [ ] `preprocess` — script ready, xenium recipe + Harmony integration
- [ ] `scoring` — script ready, project sigs + Hallmark
- [ ] `celltype_manual` — pending, decide after scoring evaluation
- [ ] `biology` — pending
- [ ] `meta_analysis` — pending

## Tasks

- [ ] Submit `scripts/run_pipeline.sbatch` on brb (`qc_filter` through `scoring`)
- [ ] Validate checkpoints after pipeline completes
- [ ] Generate QC reports (pre-filter, post-filter, post-integration)
- [ ] Evaluate clusters and decide manual celltyping vs automated
- [ ] If manual celltyping: create `metadata/celltype_map.json`, run iteratively
- [ ] Downstream biology analysis

## Blockers

- Pipeline submission pending (scripts written, not yet run)
- Same SR4 runs as `projects/visium_hd/robin/` but uses `segmented_outputs/` instead of `binned_outputs/square_008um/`

## Technical Decisions

- **Preprocessing recipe:** Xenium-like (single-cell resolution, not bin-based)
- **Cell coordinates:** Centroids from mean of polygon exterior in `cell_segmentations.geojson`; not all geojson cells appear in filtered matrix
- **Index format:** SR4 `cellid_{cell_id:09d}-1`; loader auto-detects and converts from geojson integer `cell_id`
- **No deconvolution:** Cell segmentation provides single-cell resolution directly
