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
- [ ] Raw counts in @assays$RNA@counts: **PENDING** (SLURM job submitted; connectivity issues)
- [ ] Cell type columns in RDS meta.data: **PENDING**
- [ ] Spatial coordinate extraction method: **PENDING**
- [ ] noseg cell-level segmentation status: **PENDING**

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

### S0 BLOCKER: Most RDS files are zero-filled (corrupted transfer) -- CRITICAL

**Finding (2026-03-07):** All panels EXCEPT CosMx_6k have RDS files that are
entirely zero-filled (null bytes from start to end). These files are unreadable.

| Panel | File status | Evidence |
|-------|------------|---------|
| **CosMx_6k** | VALID -- gzip RDS (magic: `1f8b`) | readRDS succeeds; Seurat object |
| CosMx_1k | CORRUPT -- all null bytes | xxd confirms; readRDS fails |
| Xenium_MT_16 | CORRUPT -- all null bytes | xxd confirms; readRDS fails |
| Xenium_5K | CORRUPT -- all null bytes | xxd confirms |
| Xenium_MT_withseg | CORRUPT -- all null bytes | xxd confirms |
| Xenium_MT_noseg | CORRUPT -- all null bytes | xxd confirms |
| Xenium_colon | CORRUPT -- all null bytes | xxd confirms |

**Root cause:** Bulk S3-to-Lustre sync ran on 2026-03-06 18:49 (`jip2007` user).
All files show identical modification timestamps, suggesting a single sync job. Only
CosMx_6k (430 MB) appears to have transferred fully. The others (36-180 MB each)
were pre-allocated on Lustre but the data was never written (or the sync failed).

**Action required:**
1. **Contact Saha lab (jip2007)** to re-run the S3 sync for the 6 failed panels, OR
2. **Get AWS credentials** to download from S3 directly:
   `s3://saha-dulai-collaboration/NC_compare_11122025/{CosMx_1k_16,Xenium_MT_16,...}/`
3. **In the meantime:** proceed with CosMx_6k only (4 samples, UC + Healthy)

**Note on SeuratObject:** System R on cayuga login node (`/usr/bin/Rscript`) does NOT
have `SeuratObject` package. Need to either install it or find a conda env with Seurat.
When inspecting CosMx_6k: class=Seurat confirmed, but further inspection halted on
missing SeuratObject for method dispatch. Must install before full inspection.

### S0.2 Raw counts (PARTIAL -- CosMx_6k only)
- CosMx_6k: `readRDS` returns Seurat object; full inspection pending (SeuratObject needed)
- Expected: `@assays$RNA@counts` = integer counts; `@assays$RNA@data` = log-normalized
- All other panels: blocked until re-sync

### S0.3 Xenium noseg (BLOCKED -- zero-filled file)

### S0.4 Spatial coordinates (BLOCKED -- zero-filled files)

---

## Todo List

### Infrastructure

- [ ] **BLOCKER: Contact Saha lab** to re-sync the 6 failed panels from S3, OR get AWS credentials
- [ ] Install SeuratObject in R on cayuga (for CosMx_6k inspection and conversion)
- [x] Create HPC working directory structure
- [x] Write `inspect_rds.R` (scripts/inspect_rds.R in this repo)
- [ ] Run full S0 inspection once files are available and SeuratObject is installed
- [ ] Write and test `convert_rds_to_h5ad.R` on CosMx_6k first (only usable panel)
- [x] Build batch manifests from NC_compare CSV (all 7 TSVs written)
- [ ] Batch SLURM job array for all-sample conversion (blocked for most panels)

### Milestone 0: Technical replicate baseline (upper bound) -- BLOCKED

**Scope:** Xenium MT noseg vs Xenium MT withseg (same 4 patients, Rectum, identical ~377-gene panel)
- **Status: BLOCKED** -- both withseg and noseg RDS files are zero-filled
- [ ] (BLOCKED) Unblock when files re-synced from S3
- [ ] Convert 4 noseg + 4 withseg RDS files
- [ ] Run `run_integration_benchmark(modality="xenium", batch_key="panel_variant")`
- [ ] Confirm: ASW_batch >> 0.7 (nearly perfect integration expected)
- [ ] Check celltype concordance between the two runs

### Milestone 0.5 (UNBLOCKED): CosMx_6k single-platform QC (4 samples)

**Scope:** CosMx_6k alone (4 samples: 1 Healthy + 3 UC, all Rectum)
- **Purpose:** Validate conversion pipeline; establish single-platform baseline
- **Status: UNBLOCKED** (CosMx_6k files are valid gzip RDS)
- [ ] Install SeuratObject in R on cayuga
- [ ] Run `inspect_rds.R` on CosMx_6k to confirm raw counts, celltype cols, coords
- [ ] Write + test `convert_rds_to_h5ad.R` on CosMx_6k sample 1
- [ ] Convert all 4 CosMx_6k samples via SLURM job array
- [ ] QC report for CosMx_6k
- This validates the full pipeline before the matched cross-platform milestones

### Milestone 1: Cross-platform, matched patients, same plex (KEY MILESTONE) -- BLOCKED

**Scope:** CosMx 1k (16 samples) + Xenium MT 16 (same 16 patients), ~377 shared genes
- **Status: BLOCKED** -- CosMx_1k and Xenium_MT_16 files are zero-filled
- **This is the scientifically most important milestone** -- CD + UC, Ileum + Rectum, matched
- [ ] (BLOCKED) Unblock when files re-synced from S3
- [ ] Convert all 16 CosMx_1k + 16 Xenium_MT RDS files
- [ ] Add metadata from CSV: `disease`, `disease_state`, `tissue_type`, `block` (patient ID)
- [ ] Run `run_integration_benchmark(batch_key="platform", celltype_key="cell_type")`
  - Methods: Harmony, scVI (if raw counts OK), ComBat, BBKNN, Scanorama
  - Add `patient_id` as additional batch key in secondary analysis
- [ ] Evaluate: patient clustering (patients should mix within celltypes but separate between)
- [ ] IBD biology: do CD vs UC cell proportions differ? (Fisher p < 0.05)
- [ ] Extend to Ileum + Rectum after Ileum validates

### Milestone 2: Cross-platform, high-plex, matched patients

**Scope:** CosMx 6k (4 samples) + Xenium 5K (same 4 patients), ~1,500-2,000 shared genes
- **Same 4 patients as M0** -- can use M0 integration as scaffold
- **Disease limitation:** UC only + 1 Healthy (no CD) -- IBD biology limited
- [ ] Convert CosMx_6k + Xenium_5K RDS files
- [ ] HVG + SVG gene selection for high-plex panels
- [ ] gimVI evaluation: impute to reference `SAHA_IBD_RNA.h5ad`
- [ ] Add Xenium colon 4 (same 4 patients, colon panel ~400 genes) as extension
- [ ] LLOKI evaluation (best-effort; Python 3.8 conda env)

### Milestone 3: Full cross-platform integration

**Scope:** All 44 samples (exclude noseg from primary analysis), all platforms
- **Patient-level batch key:** use `block` (patient ID from CSV) as grouping variable
  to separate technical (platform) from biological (disease) variation
- **Gene intersection:** ~100-300 genes across all panels
- **Imputation strategy:** gimVI with SAHA_IBD_RNA.h5ad as reference
- [ ] All samples converted
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
| `Journal.md` | TODO |
| `journal_summary.md` | TODO |
| `config.yaml` | TODO |
| `scripts/inspect_rds.R` | TODO (write to cayuga directly) |
| `scripts/convert_rds_to_h5ad.R` | TODO |
| `scripts/build_batch_manifests.py` | TODO |
| `scripts/run_integration_benchmark.py` | TODO |
| `scripts/run_batch_factor_analysis.py` | TODO |
| `metadata/phase0/cosmx_1k_samples.tsv` | TODO (populate from CSV) |
| `metadata/phase0/cosmx_6k_samples.tsv` | TODO |
| `metadata/phase0/xenium_5k_samples.tsv` | TODO |
| `metadata/phase0/xenium_mt_samples.tsv` | TODO |
