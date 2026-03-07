# Mission: robin (visium_hd_cell)

**Project:** `projects/visium_hd_cell/robin`
**Current Status:** Phase 0b COMPLETE. Ready for QC/filtering.
**Last Updated:** 2026-03-06

Project-specific goals. Repository-level pipeline and toolkit goals are in the root `Mission.md`.

---

## 1. Objective

Single-cell resolution spatial transcriptomics analysis of Visium HD data (Robin project) using SpaceRanger 4 cell segmentation output. Parallel to the bin-level (008um) analysis in `projects/visium_hd/robin/`. Cell segmentation provides true single-cell resolution from nuclei detection on H&E images.

---

## 2. Phase Alignment

| Slug | Status | Key Tasks |
|------|--------|-----------|
| **ingest_load** | Done | 15/15 samples loaded (1,574,694 cells x 18,132 genes). Cell centroids from geojson. |
| **qc_filter** | Pending | Per-sample QC, filtering, concatenation QC report |
| **metadata_attach** | Pending | Metadata attachment |
| **preprocess** | Pending | Preprocessing, normalization, clustering |
| **demographics** | Pending | Demographics (optional branch) |
| **scoring** | Pending | Gene scoring, cell typing, deconvolution |
| **celltype_manual** | Pending | Manual cell typing (optional, iterative) |
| **biology** | Pending | Downstream biology |
| **meta_analysis** | Pending | Meta analysis (optional) |

---

## 3. Completed Tasks

- [x] **Phase 0b ingest:** 15 samples loaded from SR4 `segmented_outputs/`. Cell centroids extracted from `cell_segmentations.geojson`. Per-sample `data/{sample_id}/adata.h5ad` + concatenated `results/adata.raw.h5ad` (4.1 GB). Pat5_Samp3 excluded (QC fail). Job 13971144 on brb.

---

## 4. Active Tasks

- [ ] QC filtering: run per-sample QC, generate pre-filter report

---

## 5. Blockers and Sanity Checks

- Data shares same SR4 runs as `projects/visium_hd/robin/` but uses `segmented_outputs/` instead of `binned_outputs/square_008um/`.
- Uses xenium-like preprocessing recipe (single-cell resolution, not bin-based).
- Obs names are `cellid_XXXXXXXXX-1` format (SR4 convention).

---

## 6. Technical Decisions

- **Cell coordinate extraction:** Centroids computed as mean of polygon exterior coordinates from `cell_segmentations.geojson`. Not all geojson cells are in the filtered matrix (54,450 in geojson vs 53,105 in filtered h5 for PT01-1_NAT).
- **Index format:** SR4 uses `cellid_{cell_id:09d}-1`; geojson has integer `cell_id`. Loader auto-detects and converts.
