# Mission: IBD Spatial Integration Project (ibd_spatial)

**Platform:** Multiplatform (CosMx + Xenium)
**Collaboration:** Saha lab (Cornell/Weill Cornell)
**Last Updated:** 2026-03-07

---

## Scientific Context

**Key design insight (confirmed from NC_compare_11122025_meta.csv):**
This is a **patient-matched cross-platform evaluation**. The Saha lab measured the SAME
tissue blocks on multiple spatial platforms. This is NOT a random cross-platform collection
-- the matched design provides biological ground truth for integration validation:
cells from the same patient block should cluster together after integration.

**Scientific question:** Can we recover patient identity and IBD-relevant biology
(CD vs UC, Inflamed vs Non-inflamed) after cross-platform batch correction? How much
does gene-panel divergence degrade integration quality?

**S0 status (from metadata CSV inspection 2026-03-07):**
- [x] Diagnosis labels CONFIRMED: `disease` (CD/UC/Healthy), `disease_state` (Inflamed/Non-inflamed)
- [x] Patient-matched design CONFIRMED: same block IDs across panels
- [x] Actual filenames CONFIRMED (corrected from plan; RDS in `data_IBD/` subdir)
- [x] Tissue types CONFIRMED: Ileum + Rectum (16-patient group); Rectum only (4-patient group)
- [x] Raw counts CONFIRMED: Seurat v5 `@layers$counts` (integer sparse); `@layers$data` (log-norm); `@layers$scale.data`. scVI usable.
- [x] Cell type columns CONFIRMED: CosMx has `ct_minor`/`ct_major`/`ct_minor_new`; Xenium 5K/MT have `ct_minor`/`ct_major`; Xenium noseg/withseg/colon have NO cell type annotations (only k-means clusters).
- [x] Spatial coordinates CONFIRMED: `GetTissueCoordinates()` works on all panels. CosMx: `CenterX/Y_global_px`. Xenium: `x_centroid`/`y_centroid`.
- [x] noseg cell status CONFIRMED: 28,858 cells present (Xenium default segmentation; usable for integration).

**Note on data path:** All RDS files are under `/athena/project-saha/data_IBD/`
(NOT `/athena/project-saha/` directly).

---

## Data Inventory (confirmed from CSV)

### Panel structure (51 samples total across 7 panels)

| Panel | N | File pattern | Tissue | Disease | Notes |
|-------|---|-------------|--------|---------|-------|
| CosMx 1k | 16 | `data_IBD/CosMx_1k_16/{1-16}_smi_so_qcnorm.RDS` | Ileum (8) + Rectum (8) | CD + UC | Matched to Xenium MT 16 |
| Xenium MT 16 | 16 | `data_IBD/Xenium_MT_16/{1-16}_xen_so_qcnorm.RDS` | Ileum (8) + Rectum (8) | CD + UC | Matched to CosMx 1k |
| CosMx 6k | 4 | `data_IBD/CosMx_6k_4/{1-4}_smi_so_qcnorm.RDS` | Rectum | 1 Healthy + 3 UC | Matched to Xenium 5K |
| Xenium 5K | 4 | `data_IBD/Xenium_5K_4/{1-4}_xen_so_qcnorm.RDS` | Rectum | 1 Healthy + 3 UC | Matched to CosMx 6k |
| Xenium MT withseg | 4 | `data_IBD/Xenium_MT_withseg_4/{1-4}_xen_377_with_seg_so_qcnorm.RDS` | Rectum | 1 Healthy + 3 UC | Same 4 patients as 6k/5K |
| Xenium MT noseg | 4 | `data_IBD/Xenium_MT_noseg_4/{1-4}_xen_377_no_seg_so_qcnorm.RDS` | Rectum | 1 Healthy + 3 UC | Same 4 patients; no seg stain |
| Xenium colon | 4 | `data_IBD/Xenium_colon_4/{1-4}_xen_colon_with_seg_so_qcnorm.RDS` | Rectum | 1 Healthy + 3 UC | Same 4 patients; colon panel |

**Disease composition (16-patient group):**
- CD: I0262, I0275, I0278, I0294, I0303 (5 patients)
- UC: I0276, I0277, I0284, I0286, I0310, I0321 (6 patients)
- (5 remaining unaccounted -- check CSV vs RDS patient IDs)
- Inflamed: subset of CD + UC; Non-inflamed: subset of both

**Disease composition (4-patient group -- high-plex panels):**
- Healthy: I0400 (1 patient, Rectum)
- UC: I0355, I0380, I0387 (3 patients, Rectum)
- **NO CD in this group** -- limits IBD biology for Milestones 2-3

**Additional resources:**
- `data_IBD/SAHA_IBD_RNA.h5ad` (~11 GB) -- scRNA-seq reference; `obs['ct_major_new']`, `obs['ct_minor_new']`
- `data_IBD/NC_compare_11122025_meta.csv` -- sample metadata (51 rows); confirmed readable
- S3 paths in `aws_folder` column of CSV for each sample
- `data/` directory: 34 Xenium Ranger and CosMx flat-file output dirs (`SAHA_{CP|CR|XR}_*`)

**Matched patient IDs (from block column in CSV):**
```
Patient I0294: CosMx_1k #1 (Ileum,CD,Inflamed) + Xenium_MT #1 (same)
Patient I0278: CosMx_1k #2 (Ileum,CD,Non-inf)  + Xenium_MT #2
Patient I0262: CosMx_1k #3 (Ileum,CD,Non-inf)  + Xenium_MT #3
Patient I0321: CosMx_1k #4 (Ileum,UC,Non-inf)  + Xenium_MT #4
...
Patient I0400: CosMx_6k #1 + Xenium_5K #1 + Xenium_MT_withseg #1 + Xenium_MT_noseg #1 + Xenium_colon #1
Patient I0355: CosMx_6k #2 + Xenium_5K #2 + ...  (same for #3=I0380, #4=I0387)
```

---

## Step 0: Pre-work (blocking verification) — STATUS

### S0.1 Diagnosis metadata -- CONFIRMED (2026-03-07)
- `disease`: CD, UC, Healthy (confirmed from NC_compare_11122025_meta.csv)
- `disease_state`: Inflamed, Non-inflamed, Healthy
- Tissue: Ileum, Rectum (16-patient group); Rectum (4-patient group)
- CSV fully readable; 51 rows; patient block IDs confirmed matched across panels

### S0 BLOCKER: RESOLVED (2026-03-07)

**Previous issue:** 6 of 7 panels had zero-filled RDS files from failed S3-to-Lustre sync.
**Resolution (2026-03-07):** Re-sync completed by `jip2007`. All 7 panels now have valid
gzip RDS files (magic bytes `1f8b` confirmed via `xxd`). All milestones UNBLOCKED.

| Panel | Status | Files | Sizes |
|-------|--------|-------|-------|
| CosMx_1k | VALID | 16 | 2.8M-178M |
| CosMx_6k | VALID | 4 | ~430M each |
| Xenium_MT_16 | VALID | 16 | varies |
| Xenium_5K | VALID | 4 | varies |
| Xenium_MT_withseg | VALID | 4 | varies |
| Xenium_MT_noseg | VALID | 4 | varies |
| Xenium_colon | VALID | 4 | varies |

**Note on SeuratObject:** System R on cayuga login node (`/usr/bin/Rscript`) does NOT
have `SeuratObject` package. Need to install before full RDS inspection.

### S0.2 Raw counts -- CONFIRMED (2026-03-07)
- All panels: Seurat v5 `@layers$counts` = integer sparse (raw); `@layers$data` = log-norm
- Assay names: `Nanostring` (CosMx), `Xenium` (Xenium) -- NOT `RNA`
- scVI usable on raw counts

### S0.3 Xenium noseg -- CONFIRMED (2026-03-07)
- 28,858 cells present with spatial coordinates (Xenium default segmentation)
- No cell type annotations (only k-means clusters)

### S0.4 Spatial coordinates -- CONFIRMED (2026-03-07)
- `GetTissueCoordinates()` works on all 7 panels

### S0.5 Batch conversion -- COMPLETE (2026-03-07)
- All 52 RDS files converted to h5ad (SLURM job 2700745, ~24s each)
- R packages reinstalled under R/4.4.1 module in `~/R/libs_R441` (SeuratObject, Matrix, sp, Rcpp)
- Validated: raw counts (int), spatial coords, disease metadata joined from CSV

---

## Todo List

### Infrastructure

- [x] ~~BLOCKER: Contact Saha lab~~ — RESOLVED: all 7 panels re-synced and valid
- [x] Install SeuratObject in R on cayuga (~/R/libs; 2026-03-07)
- [x] Create HPC working directory structure
- [x] Write `inspect_rds.R` (scripts/inspect_rds.R in this repo)
- [x] Run full S0 inspection on ALL 7 panels (2026-03-07) — Seurat v5 layers, raw counts confirmed, spatial OK
- [x] Write and test `convert_rds_to_h5ad.R` — tested on Xenium MT sample 1 (5138 cells x 377 genes, 2.9MB h5ad)
- [x] Batch conversion SLURM array job 2700745 — all 52 tasks completed (2026-03-07)
- [x] Build batch manifests from NC_compare CSV (all 7 TSVs written)
- [x] Reinstall R packages under R/4.4.1 module (`~/R/libs_R441`; SeuratObject, Matrix, sp, Rcpp)

### Milestone 0: Technical replicate baseline (upper bound) -- UNBLOCKED

**Scope:** Xenium MT noseg vs Xenium MT withseg (same 4 patients, Rectum, identical ~377-gene panel)
- **Status: UNBLOCKED** -- all RDS files valid (re-synced 2026-03-07)
- [x] Convert 4 noseg + 4 withseg RDS files (done in batch job 2700745)
- [x] Run integration benchmark (SLURM 2700824): PCA, Harmony, scVI, BBKNN — all batch_score > 0.99
- [x] Confirm: near-zero batch effect (batch ASW ~ 0, batch_score 0.992-0.997)
- [x] Disease signal preserved: Healthy vs UC separates on UMAP
- Note: Scanorama failed silently; no celltype annotations available for these panels

### Milestone 0.5 (UNBLOCKED): CosMx_6k single-platform QC (4 samples)

**Scope:** CosMx_6k alone (4 samples: 1 Healthy + 3 UC, all Rectum)
- **Purpose:** Validate conversion pipeline; establish single-platform baseline
- **Status: UNBLOCKED** (CosMx_6k files are valid gzip RDS)
- [x] Install SeuratObject in R on cayuga
- [x] Run `inspect_rds.R` — all checks passed
- [x] Write + test `convert_rds_to_h5ad.R` on CosMx_6k sample 1
- [x] Convert all 4 CosMx_6k samples (done in batch job 2700745)
- [ ] QC report for CosMx_6k
- This validates the full pipeline before the matched cross-platform milestones

### Milestone 1: Cross-platform, matched patients, same plex (KEY MILESTONE) -- COMPLETE

**Scope:** CosMx 1k (16 samples) + Xenium MT 16 (same 16 patients), 119 shared genes
- **Status: COMPLETE** (SLURM 2700857, 2h7m on GPU A40)
- [x] Convert all 16 CosMx_1k + 16 Xenium_MT RDS files (done in batch job 2700745)
- [x] Metadata from CSV joined during conversion: `disease`, `disease_state`, `tissue_type`, `patient_id`
- [x] Run integration benchmark (PCA, Harmony, scVI, BBKNN; Scanorama failed)
- [x] Best method: Harmony (batch_score=0.971); all methods > 0.93
- [x] IBD biology: 21 cell types show significant CD vs UC differences (Fisher p < 0.05)
- [x] Cell type proportions differ: DC (CD-enriched), Macrophage/Inflammatory fibroblast (UC-enriched)
- Note: Only 119 shared genes (CosMx 950 ∩ Xenium 377); celltype ASW negative (insufficient genes for fine separation)
- Note: First attempt (job 2700834) timed out on CPU; fixed with GPU + reduced epochs

### Milestone 2: Cross-platform, high-plex, matched patients -- COMPLETE

**Scope:** CosMx 6k (4 samples) + Xenium 5K (same 4 patients), 2,552 shared genes
- **Status: COMPLETE** (SLURM 2701598, 54 min on GPU A40)
- [x] Convert CosMx_6k + Xenium_5K RDS files (done in batch job 2700745)
- [x] HVG selection: 2,000 / 2,552 genes (batch-aware)
- [x] Run integration benchmark (PCA, Harmony, scVI, BBKNN; Scanorama failed)
- [x] Best method: scVI (batch_score=0.992) -- dominates with higher gene count + model capacity
- [x] Bio eval: per-platform ASW, kNN transfer, disease composition
- Note: scVI overtakes Harmony as best method (vs M1 where Harmony won)
- Note: Fine celltype ASW still negative (-0.070) even with 2,552 genes
- Note: PCA batch ASW higher (0.195) than M1 (0.065) -- more genes reveal more platform effects
- [x] IBD marker check: 8/8 canonical markers present in 2,552-gene intersection — PASS
- [x] Platform entropy: scVI median=0.632 (>0.5 threshold) — PASS
- [x] Formal success criteria: 4/4 applicable criteria PASS for scVI
- [x] scANVI evaluation: batch_score=0.881, ct_broad_ASW=0.074 (best bio conservation), prediction accuracy 94.9%
- [x] Scanorama fix: batch_score=0.837, ct_broad_ASW=-0.089 (poor; worst bio conservation)
- ~~gimVI~~: dropped (removed from scvi-tools 1.4.2)
- ~~LLOKI~~: dropped (Python 3.8 incompatible with sc_tools env)
- [ ] Add Xenium colon 4 (same 4 patients, colon panel ~400 genes) as extension
- [ ] Regenerate project report with all 6 methods

### Milestone 3: Full cross-platform integration

**Scope:** All 44 samples (exclude noseg from primary analysis), all platforms
- **Patient-level batch key:** use `block` (patient ID from CSV) as grouping variable
  to separate technical (platform) from biological (disease) variation
- **Gene intersection:** ~100-300 genes across all panels
- **Imputation strategy:** gimVI with SAHA_IBD_RNA.h5ad as reference
- [x] All samples converted (52/52, batch job 2700745)
- [ ] Run with n_latent=8, n_hidden=32 for scVI (small gene set)
- [ ] Full IBD biology: CD vs UC spatial patterns (ileum only for CD)

---

## HPC Setup (`~/elementolab/projects/ibd_spatial/`)

```bash
# Working dir (do once when connectivity allows)
mkdir -p ~/elementolab/projects/ibd_spatial/{data/{raw/{cosmx_1k,cosmx_rds,xenium,xenium_rds},h5ad,reference},results/tmp/integration_test,figures,metadata/phase0,outputs,logs}

# Soft links -- CosMx flat-file dirs (SAHA_CP_* format)
for d in /athena/project-saha/data/SAHA_CP_*; do ln -sf "$d" ~/elementolab/projects/ibd_spatial/data/raw/cosmx_1k/; done

# Soft links -- Xenium Ranger output dirs (SAHA_XR_* format)
for d in /athena/project-saha/data/SAHA_XR_*; do ln -sf "$d" ~/elementolab/projects/ibd_spatial/data/raw/xenium/; done

# RDS directories (in data_IBD subdir -- CORRECTED from original plan)
ln -sf /athena/project-saha/data_IBD/CosMx_1k_16         ~/elementolab/projects/ibd_spatial/data/raw/cosmx_rds/cosmx_1k
ln -sf /athena/project-saha/data_IBD/CosMx_6k_4          ~/elementolab/projects/ibd_spatial/data/raw/cosmx_rds/cosmx_6k
ln -sf /athena/project-saha/data_IBD/Xenium_5K_4         ~/elementolab/projects/ibd_spatial/data/raw/xenium_rds/xenium_5k
ln -sf /athena/project-saha/data_IBD/Xenium_MT_16        ~/elementolab/projects/ibd_spatial/data/raw/xenium_rds/xenium_mt
ln -sf /athena/project-saha/data_IBD/Xenium_colon_4      ~/elementolab/projects/ibd_spatial/data/raw/xenium_rds/xenium_colon
ln -sf /athena/project-saha/data_IBD/Xenium_MT_withseg_4 ~/elementolab/projects/ibd_spatial/data/raw/xenium_rds/xenium_mt_withseg
ln -sf /athena/project-saha/data_IBD/Xenium_MT_noseg_4   ~/elementolab/projects/ibd_spatial/data/raw/xenium_rds/xenium_mt_noseg
ln -sf /athena/project-saha/data_IBD/SAHA_IBD_RNA.h5ad   ~/elementolab/projects/ibd_spatial/data/reference/SAHA_IBD_RNA.h5ad
```

---

## Batch Manifests (from NC_compare CSV)

Manifests are populated from `NC_compare_11122025_meta.csv`. Column mapping:
- `sample_id` from CSV `sample_id` + panel suffix
- `patient_id` from CSV `block` column
- `disease` from CSV `disease` (CD/UC/Healthy)
- `disease_state` from CSV `disease_state` (Inflamed/Non-inflamed/Healthy)
- `tissue` from CSV `tissue_type` (Ileum/Rectum)
- `panel` from CSV `panel` column
- `rds_path` constructed from panel dir + sample number

See `scripts/build_batch_manifests.py` for automated population.

---

## Integration Strategy Notes (revised)

### Batch key hierarchy (corrected from original plan)
- `platform` = primary batch (CosMx vs Xenium) -- the main technical effect
- `patient_id` (block) = **biological anchor** (NOT a batch to correct)
- `tissue_type` = stratification variable (analyze Ileum/Rectum separately first)
- `disease` / `disease_state` = biological target variables
- `institution` = confounded with panel (not independently testable)

### Method selection for small gene sets (< 500 genes)
```python
# For M0 and M1 (~377 genes): use ALL shared genes (no HVG selection)
# Reduce scVI dimensions
run_scvi(adata, n_latent=10, n_hidden=64)  # reduced from defaults 30/128
# For M2 (~1500 genes): HVG+SVG intersection
# For M3 (~100-300 genes): n_latent=8, n_hidden=32
```

### Patient-matched validation metric (new -- not in original plan)
After integration, compute:
- Per-patient silhouette score across platforms (should be HIGH -- patients should NOT separate)
- Per-platform silhouette score (should be LOW -- platforms should mix)
- This is more informative than generic ASW_batch given the matched design

---

## Success Criteria (revised)

| Metric | M0 (threshold) | M1 (threshold) | M2-M3 (threshold) |
|--------|----------------|----------------|-------------------|
| ASW_platform (lower=better integration) | < 0.2 | < 0.35 | < 0.5 |
| Per-patient cross-platform mixing | > 0.8 | > 0.6 | > 0.4 |
| Bio conservation (ASW_celltype) | vs unintegrated | vs unintegrated | vs unintegrated |
| IBD markers in shared set | N/A | >= 5/8 canonical | >= 3/8 |
| CD vs UC proportions differ | N/A | Fisher p < 0.05 | Limited (UC only) |

---

## Key Decisions and Rationale

| Decision | Rationale |
|----------|-----------|
| Use patient-matched metric, not generic ASW_batch | Matched design gives ground truth; patient ID is the anchor |
| Analyze Ileum and Rectum separately first | Tissue type confounds platform comparisons; CD is Ileum-enriched |
| Start M1 before M0 if noseg status unclear | M1 has more samples (16+16) and CD biology |
| Do not attempt CD biology in M2-M3 | 4-patient group has NO CD; UC only |
| Keep noseg samples for technical analysis only | Segmentation stain may affect transcript detection |
| Use `block` column (not sample_id) as patient key | CSV `block` = tissue block ID = patient identifier |

---

## Files

| File | Status |
|------|--------|
| `Mission.md` (this file) | Done |
| `Journal.md` | Done |
| `journal_summary.md` | Done |
| `config.yaml` | Done |
| `scripts/inspect_rds.R` | Done (on cayuga) |
| `scripts/convert_rds_to_h5ad.R` | Done |
| `scripts/run_convert_array.sh` | Done |
| `scripts/build_batch_manifests.py` | TODO |
| `scripts/run_integration_benchmark.py` | TODO |
| `scripts/run_batch_factor_analysis.py` | TODO |
| `metadata/phase0/cosmx_1k_samples.tsv` | Done (from CSV) |
| `metadata/phase0/cosmx_6k_samples.tsv` | Done |
| `metadata/phase0/xenium_5k_samples.tsv` | Done |
| `metadata/phase0/xenium_mt_samples.tsv` | Done |
