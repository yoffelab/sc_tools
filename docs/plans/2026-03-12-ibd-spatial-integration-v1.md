---
status: complete
created: 2026-03-07
project: ibd_spatial
summary: "IBD Spatial: cross-platform integration benchmark (CosMx + Xenium), resolVI+scANVI completion, final report"
---

# IBD Spatial Cross-Platform Integration

## Context

Patient-matched cross-platform spatial transcriptomics integration (Saha lab collab).
Same tissue blocks measured on CosMx (1k, 6k) and Xenium (MT, 5K, colon, noseg, withseg).
51 samples across 7 panels. Goal: recover patient identity and IBD biology after
cross-platform batch correction.

Project dir: `projects/multiplatform/ibd_spatial/`
HPC: cayuga (`/home/fs01/juk4007/elementolab/projects/ibd_spatial/`)

## Completed Steps

### Infrastructure (2026-03-07)
- [x] S0 metadata verification (diagnosis, raw counts, spatial coords, noseg status)
- [x] R package installation on cayuga (SeuratObject, Matrix, sp, Rcpp under R/4.4.1)
- [x] RDS-to-h5ad conversion pipeline (`convert_rds_to_h5ad.R`)
- [x] Batch conversion of all 52 RDS files (SLURM 2700745)
- [x] Batch manifests from NC_compare CSV (7 TSV files)

### Milestone 0: Technical replicate baseline (2026-03-08)
- [x] Xenium noseg vs withseg (8 samples, 4 patients, 377 genes)
- [x] All methods batch_score > 0.99 (near-zero batch effect)
- [x] resolVI completed (batch_score=0.994, entropy=0.962) — SLURM 2702309

### Milestone 1: Cross-platform matched patients (2026-03-08)
- [x] CosMx 1k + Xenium MT (32 samples, 16 patients, 119 shared genes)
- [x] PCA, Harmony, scVI, BBKNN benchmarked; Harmony best (0.971)
- [x] IBD biology: 21 cell types show CD vs UC differences (Fisher p<0.05)
- [x] Bio evaluation + scib reports generated (2026-03-09)
- [x] resolVI completed (batch_score=0.916) — SLURM 2702309
- [x] scANVI: trained but prediction failed in job 2702309; fix submitted in job 2702344 — completed

### Milestone 2: High-plex cross-platform (2026-03-09)
- [x] CosMx 6k + Xenium 5K (8 samples, 4 patients, 2552 shared genes)
- [x] scVI best (batch_score=0.992); overtakes Harmony with more genes
- [x] M2 followup (SLURM 2702091): scANVI, Scanorama, IBD markers, entropy, criteria
- [x] All 6 base methods benchmarked; scVI recommended, scANVI best bio conservation
- [x] 8/8 IBD canonical markers present; platform entropy 0.632 (PASS)
- [x] resolVI: was mid-training when job 2702309 cancelled; resubmitted in 2702344 — M2 resolVI completed (batch_score=0.907)
- [x] Final report with resolVI included (job 2702436, 74 MB HTML)

## Completed Steps (2026-03-12)

### Step 1: Check job 2702344 (resolVI+scANVI resubmit) — DONE
- M2 resolVI: 100 epochs on A40 (25 min), batch_score=0.907, entropy=0.539
- M0/M1 resolVI: skipped (already cached)
- scANVI M1/M2: skipped (already cached)

### Step 2: Update report script with resolVI — DONE
- Added `load_resolvi_benchmark()` + `merge_resolvi_into_benchmark()` helpers
- resolVI UMAPs added to all 3 milestones
- resolVI added to scib embeddings dict for all milestones

### Step 3: Regenerate project report — DONE
- SLURM 2702436: M0 (4 methods), M1 (5 methods), M2 (6 methods in scib)
- Report: `figures/QC/ibd_spatial_integration_report.html` (74 MB)

### Step 4: Update Mission.md and Journal.md — DONE
- resolVI results recorded per milestone
- M2 "Regenerate project report" marked complete
- Plan status set to complete

### Step 5 (future): Milestone 3 — full cross-platform
- All 44 samples, all platforms, ~100-300 shared genes
- n_latent=8, n_hidden=32 for scVI (small gene set)
- Full IBD biology: CD vs UC spatial patterns

## Key Technical Notes

- **resolVI** uses Pyro backend (not PyTorch Lightning like scVI/scANVI)
  - Requires `obsm['X_spatial']` (not `obsm['spatial']`)
  - `train()` does NOT accept `train_size`, `early_stopping`, `early_stopping_patience`
  - Cannot chain with scANVI via `from_scvi_model()` — must run fresh scVI->scANVI separately
- **Scanorama** `integrate_scanpy()` returns None (modifies in-place)
- **scANVI** categorical fillna: use `.astype(str).replace("nan", "Unknown")`
- **gimVI**: dropped (removed from scvi-tools 1.4.2)
- **LLOKI**: dropped (Python 3.8 incompatible)

## Final Results Summary

| Milestone | Best Method | Batch Score | Bio (ct_broad ASW) | resolVI | Notes |
|-----------|-------------|-------------|---------------------|---------|-------|
| M0 | scVI | 0.997 | N/A (no labels) | 0.994 | Near-zero batch effect |
| M1 | Harmony | 0.971 | -0.008 | 0.916 | 119 genes; CD+UC biology preserved |
| M2 | scVI | 0.992 | 0.009 | 0.907 | 2552 genes; scANVI best bio (0.074) |

resolVI results (complete):
- M0: batch_score=0.994, entropy=0.962 (2nd best, excellent)
- M1: batch_score=0.916, entropy=0.356 (5th, below Harmony/scVI)
- M2: batch_score=0.907, entropy=0.539, ct_broad_ASW=0.025 (4th, passes entropy threshold)
