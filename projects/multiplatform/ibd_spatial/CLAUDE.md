# ibd_spatial — Claude Code Configuration

## Sync Before Work

1. @Mission.md — current todo list and phase status
2. @Journal.md — recent decisions and conventions

For repo-wide rules (container, conventions, testing): see repo root CLAUDE.md.

---

## Project Context

**Goal:** Patient-matched cross-platform spatial integration validation. The Saha lab measured
the SAME tissue blocks on multiple spatial platforms (CosMx 1k, CosMx 6k, Xenium MT, Xenium 5K,
Xenium colon). Evaluate whether cross-platform batch correction recovers patient identity and
IBD-relevant biology (CD vs UC, Inflamed vs Non-inflamed).

**Platform:** Multiplatform (CosMx 1k + CosMx 6k + Xenium MT + Xenium 5K + Xenium colon)

**Collaboration:** Saha lab (Cornell/Weill Cornell)

**Runtime:** cayuga (`/athena/elementolab/scratch/juk4007/sc_tools/projects/multiplatform/ibd_spatial`)

**Data source:** `/athena/project-saha/data_IBD/` (read-only; 51 samples across 7 panels)

**Phase status** (see Mission.md for full detail):
- `ingest_load` (`S0.5`): Complete — all 52 RDS files converted to h5ad via SLURM array job
- **`qc_filter` / `preprocess` (Active):** Cross-platform integration benchmarking (Milestone 1)

---

## Running This Project

```bash
# From repo root:
./scripts/run_container.sh projects/multiplatform/ibd_spatial python projects/multiplatform/ibd_spatial/scripts/<script>.py
snakemake -d projects/multiplatform/ibd_spatial -s projects/multiplatform/ibd_spatial/Snakefile <target>

# HPC (cayuga):
sbatch projects/multiplatform/ibd_spatial/scripts/run_m1.sh
```

---

## Key Files and Checkpoints

| Path | Description |
|------|-------------|
| `results/adata.filtered.h5ad` | Concatenated, QC-filtered AnnData across all 7 panels |
| `results/adata.annotated.h5ad` | With clinical metadata (disease, disease_state, tissue, patient) |
| `results/adata.normalized.h5ad` | Integrated across panels; selected integration method |
| `metadata/NC_compare_11122025_meta.csv` | 51-sample metadata CSV (patient block IDs, disease labels) |
| `metadata/phase0/` | Per-panel batch manifests (7 TSV files, one per panel) |
| `data/` | 34 Xenium Ranger and CosMx flat-file output dirs (`SAHA_{CP|CR|XR}_*`) |
| `figures/QC/` | Pre/post-filter QC reports; integration benchmark |

---

## Project-Specific Conventions

- **Matched design:** Same patient block measured on multiple platforms — this is the ground truth for integration validation
- **Assay names:** `Nanostring` (CosMx), `Xenium` (Xenium) — NOT `RNA`; use assay-aware loading
- **Raw counts:** Seurat v5 `@layers$counts` (integer sparse); scVI usable on raw counts
- **Spatial coords:** CosMx = `CenterX/Y_global_px`; Xenium = `x_centroid`/`y_centroid`
- **Cell types:** CosMx and most Xenium panels have `ct_minor`/`ct_major`; Xenium noseg/withseg/colon have NO annotations (k-means only)
- **Primary integration metric:** Patient recovery (batch = patient ID) is the key metric since matched design provides ground truth
- **R packages:** Install under `~/R/libs_R441` (R/4.4.1 module on cayuga); SeuratObject, Matrix, sp, Rcpp
