# Plan: IBD Spatial Integration Project (ibd_spatial)

## Context

**Platform:** Multiplatform (CosMx + Xenium)
**Collaboration:** Saha lab (Cornell/Weill Cornell)
**Scientific question:** Can we recover patient identity and IBD-relevant biology (CD vs UC, Inflamed vs Non-inflamed) after cross-platform batch correction? How much does gene-panel divergence degrade integration quality?
**Key design insight:** Patient-matched cross-platform evaluation -- same tissue blocks measured on multiple spatial platforms. Matched design provides biological ground truth for integration validation.

---

## Data Inventory

### Panel structure (54 samples, 8 panels)

52 RDS files across 7 panels (51 unique biological; noseg is technical replicate of withseg) + 2 CosMx Ileum IBD Seurat objects converted separately.

| Panel | N | File pattern | Tissue | Disease | Notes |
|-------|---|-------------|--------|---------|-------|
| CosMx 1k | 16 | `data_IBD/CosMx_1k_16/{1-16}_smi_so_qcnorm.RDS` | Ileum (8) + Rectum (8) | CD + UC | Matched to Xenium MT 16 |
| Xenium MT 16 | 16 | `data_IBD/Xenium_MT_16/{1-16}_xen_so_qcnorm.RDS` | Ileum (8) + Rectum (8) | CD + UC | Matched to CosMx 1k |
| CosMx 6k | 4 | `data_IBD/CosMx_6k_4/{1-4}_smi_so_qcnorm.RDS` | Rectum | 1 Healthy + 3 UC | Matched to Xenium 5K |
| Xenium 5K | 4 | `data_IBD/Xenium_5K_4/{1-4}_xen_so_qcnorm.RDS` | Rectum | 1 Healthy + 3 UC | Matched to CosMx 6k |
| Xenium MT withseg | 4 | `data_IBD/Xenium_MT_withseg_4/{1-4}_xen_377_with_seg_so_qcnorm.RDS` | Rectum | 1 Healthy + 3 UC | Same 4 patients as 6k/5K |
| Xenium MT noseg | 4 | `data_IBD/Xenium_MT_noseg_4/{1-4}_xen_377_no_seg_so_qcnorm.RDS` | Rectum | 1 Healthy + 3 UC | No seg stain |
| Xenium colon | 4 | `data_IBD/Xenium_colon_4/{1-4}_xen_colon_with_seg_so_qcnorm.RDS` | Rectum | 1 Healthy + 3 UC | Colon panel |
| CosMx 6k ILE IBD | 2 | `data/cosmx_ile_ibd_{01,02}/` | Ileum | IBD | NEW: 387K + 54K cells, 6175 genes, from Saha lab raw Seurat objects |

**Matched patient IDs (from block column in CSV):**
- 16-patient group: CosMx 1k + Xenium MT 16 (CD: I0262, I0275, I0278, I0294, I0303; UC: I0276, I0277, I0284, I0286, I0310, I0321; 5 remaining TBD)
- 4-patient group: CosMx 6k + Xenium 5K + Xenium MT withseg/noseg + Xenium colon (Healthy: I0400; UC: I0355, I0380, I0387)
- **NO CD in 4-patient group** -- limits IBD biology for M2-M3

**Additional resources:**
- `data_IBD/SAHA_IBD_RNA.h5ad` (~11 GB) -- scRNA-seq reference (`obs['ct_major_new']`, `obs['ct_minor_new']`)
- `data_IBD/NC_compare_11122025_meta.csv` -- sample metadata (51 rows)
- `data/` directory: 34 Xenium Ranger and CosMx flat-file output dirs (`SAHA_{CP|CR|XR}_*`)

---

## Phase Status

- [x] **S0: Pre-work verification** -- COMPLETE (all 7 panels valid, RDS re-synced, 52/52 converted to h5ad)
- [x] **M0: Technical replicate baseline** -- COMPLETE (noseg vs withseg, batch_score 0.992-0.997)
- [ ] **M0.5: CosMx 6k single-platform QC** -- QC report pending
- [x] **M1: Cross-platform matched, same plex** -- COMPLETE (Harmony best, batch_score=0.971)
- [x] **M2: Cross-platform high-plex matched** -- COMPLETE (scVI best, batch_score=0.992)
- [ ] **M3: Full cross-platform integration** -- NOT STARTED

---

## Todo

### Infrastructure (COMPLETE)

- [x] Install SeuratObject in R on cayuga (`~/R/libs_R441`)
- [x] Create HPC working directory structure
- [x] Write and test `inspect_rds.R` and `convert_rds_to_h5ad.R`
- [x] Batch conversion SLURM array job 2700745 -- all 52 tasks completed
- [x] Build batch manifests from NC_compare CSV (all 7 TSVs written)

### M0: Technical replicate baseline (COMPLETE)

- [x] Convert 4 noseg + 4 withseg RDS files
- [x] Run integration benchmark (PCA, Harmony, scVI, BBKNN) -- all batch_score > 0.99
- [x] Confirm near-zero batch effect (batch ASW ~ 0)
- [x] Disease signal preserved: Healthy vs UC separates on UMAP
- Note: Scanorama failed silently; no celltype annotations for these panels

### M0.5: CosMx 6k single-platform QC

- [x] Convert all 4 CosMx_6k samples
- [ ] Generate QC report for CosMx_6k (validates pipeline before cross-platform milestones)

### M1: Cross-platform matched, same plex (COMPLETE)

- [x] Convert 16 CosMx_1k + 16 Xenium_MT RDS files
- [x] Metadata joined during conversion: `disease`, `disease_state`, `tissue_type`, `patient_id`
- [x] Integration benchmark (PCA, Harmony, scVI, BBKNN; Scanorama failed)
- [x] Best method: Harmony (batch_score=0.971); all methods > 0.93
- [x] IBD biology: 21 cell types show significant CD vs UC differences (Fisher p < 0.05)
- [x] Cell type proportions differ: DC (CD-enriched), Macrophage/Inflammatory fibroblast (UC-enriched)
- Note: 119 shared genes (CosMx 950 intersect Xenium 377); celltype ASW negative (insufficient genes for fine separation)

### M2: Cross-platform high-plex matched (COMPLETE)

- [x] HVG selection: 2,000 / 2,552 shared genes (batch-aware)
- [x] Integration benchmark (PCA, Harmony, scVI, BBKNN; Scanorama fixed)
- [x] Best method: scVI (batch_score=0.992) -- dominates with higher gene count
- [x] IBD marker check: 8/8 canonical markers present -- PASS
- [x] Platform entropy: scVI median=0.632 (>0.5 threshold) -- PASS
- [x] Formal success criteria: 4/4 applicable criteria PASS for scVI
- [x] scANVI evaluation: batch_score=0.881, ct_broad_ASW=0.074 (best bio conservation), prediction accuracy 94.9%
- [x] Scanorama fix: batch_score=0.837, ct_broad_ASW=-0.089 (worst bio conservation)
- [x] resolVI completed for all milestones (SLURM 2702344): M0=0.994, M1=0.916, M2=0.907
- [x] Project report regenerated with all 7 methods including resolVI (SLURM 2702436)
- [x] **scANVI hyperparameter sweep** (SLURM 2703115, 3h12m A100): 12 configs, 11 succeeded; best = A6_3layer_genebatch (ct_broad_asw=0.189, batch=0.903, entropy=0.576) -- 2.5x over vanilla
- [x] **resolVI-SS hyperparameter sweep** (SLURM 2703116, 5h3m A100): 10 configs, 6 succeeded; best = B8_per_sample (ct_broad_asw=0.104, batch=0.836, entropy=0.296) -- 2.4x over vanilla
- [x] Added `sc_tools/pp/integration_configs.py` with `get_scanvi_config()` and `get_resolvi_ss_config()` helpers
- [x] **Updated best method: scANVI A6_3layer_genebatch** -- Pareto-dominant on bio, batch, entropy (replaces scVI as top recommendation)
- [x] **Random hyperparameter search** (SLURM 2703212, 40-config array, A100): 29/40 completed; **R023 achieves ct_broad_asw=0.396** (2.1x over A6). Key: 4 layers + cls_ratio=192 + dropout=0.20
- [x] Cancelled retry jobs (2703394, 2703420) -- not needed; 29 configs sufficient for search
- [x] Update integration_configs.py with R023 findings (4-layer architecture, cls_ratio=192, dropout=0.20)
- [x] **M2 final benchmark job (SLURM 2703646, cayuga A100) -- COMPLETE**
  - 10 methods benchmarked; R023 best (ct_broad_asw=0.378, batch_score=0.909)
  - Output files: `m2_final_benchmark.csv`, `figures/m2/m2_umap_grid.png`, `figures/QC/m2_final_benchmark_report.html`
  - Scripts: `scripts/run_m2_final_benchmark.py` + `.sh`
- [ ] Add Xenium colon 4 (same 4 patients, colon panel ~400 genes) as extension
- [ ] Ablation experiments: test R023 config on M1 (119-gene) to confirm generalization
- [ ] Final M2 integration with best config on full dataset (not just benchmark subsample)

### M3: Full cross-platform integration

- [x] All 52 samples converted
- [ ] Run with n_latent=8, n_hidden=32 for scVI (small gene set ~100-300 shared genes)
- [ ] Full IBD biology: CD vs UC spatial patterns (ileum only for CD)
- [ ] Patient-level batch key: use `block` (patient ID) as grouping variable

---

## Blocked

- **M0.5 QC report:** Low priority; pipeline validated through M1/M2 completion
- **M3 gene intersection:** Expect ~100-300 genes across all panels; may need imputation (gimVI dropped -- removed from scvi-tools 1.4.2; LLOKI dropped -- Python 3.8 incompatible)

---

## Technical Decisions

### Integration strategy

| Decision | Rationale |
|----------|-----------|
| Use patient-matched metric, not generic ASW_batch | Matched design gives ground truth; patient ID is the anchor |
| Analyze Ileum and Rectum separately first | Tissue type confounds platform comparisons; CD is Ileum-enriched |
| Do not attempt CD biology in M2-M3 | 4-patient group has NO CD; UC only |
| Keep noseg samples for technical analysis only | Segmentation stain may affect transcript detection |
| Use `block` column as patient key | CSV `block` = tissue block ID = patient identifier |

### Batch key hierarchy

- `platform` = primary batch (CosMx vs Xenium) -- main technical effect
- `patient_id` (block) = biological anchor (NOT a batch to correct)
- `tissue_type` = stratification variable
- `disease` / `disease_state` = biological target variables

### Method selection for small gene sets

- M0, M1 (~377 genes): use ALL shared genes, no HVG; scVI n_latent=10, n_hidden=64
- M2 (~2,552 genes): HVG selection (2,000 genes, batch-aware); **scANVI R023** recommended (n_latent=53, n_hidden=512, n_layers=4, NB likelihood, gene-batch dispersion, classification_ratio=192, dropout=0.20, 125+30 epochs; ct_broad_asw=0.396)
- M3 (~100-300 genes): n_latent=8, n_hidden=32

### Hyperparameter sweep findings (2026-03-13)

- Gene-batch dispersion is the single most impactful parameter (per-platform noise modeling frees latent space for biology)
- Fine-grained celltype labels hurt cross-platform -- use celltype_broad for semi-supervised training
- **4 layers is critical** -- all top random search configs use n_layers=4; 2-layer configs cluster below 0.18
- **Higher classification_ratio (120-190) helps** -- pushes model to respect celltype structure more aggressively
- **Higher dropout (0.15-0.20) improves generalization** -- prevents platform-specific memorization
- Hidden dim matters less than depth -- R026 (h=128) beats R007 (h=512) at same layer count
- Longer training (200+80 ep) causes batch overcorrection -- 125+30 is sufficient
- resolVI-SS benefits from batch_key=sample (not platform), but still underperforms scANVI

### Success criteria

| Metric | M0 | M1 | M2-M3 |
|--------|----|----|-------|
| ASW_platform (lower=better) | < 0.2 | < 0.35 | < 0.5 |
| Per-patient cross-platform mixing | > 0.8 | > 0.6 | > 0.4 |
| Bio conservation (ASW_celltype) | vs unintegrated | vs unintegrated | vs unintegrated |
| IBD markers in shared set | N/A | >= 5/8 canonical | >= 3/8 |
| CD vs UC proportions differ | N/A | Fisher p < 0.05 | Limited (UC only) |

---

## HPC Setup (`~/elementolab/projects/ibd_spatial/`)

```bash
# Working dir
mkdir -p ~/elementolab/projects/ibd_spatial/{data/{raw/{cosmx_1k,cosmx_rds,xenium,xenium_rds},h5ad,reference},results/tmp/integration_test,figures,metadata/phase0,outputs,logs}

# Soft links -- CosMx flat-file dirs
for d in /athena/project-saha/data/SAHA_CP_*; do ln -sf "$d" ~/elementolab/projects/ibd_spatial/data/raw/cosmx_1k/; done

# Soft links -- Xenium Ranger output dirs
for d in /athena/project-saha/data/SAHA_XR_*; do ln -sf "$d" ~/elementolab/projects/ibd_spatial/data/raw/xenium/; done

# RDS directories (in data_IBD subdir)
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

## Conventions

- **Matched design:** Same patient block measured on multiple platforms -- ground truth for integration validation
- **Assay names:** `Nanostring` (CosMx), `Xenium` (Xenium) -- NOT `RNA`; use assay-aware loading
- **Raw counts:** Seurat v5 `@layers$counts` (integer sparse); scVI usable on raw counts
- **Spatial coords:** CosMx = `CenterX/Y_global_px`; Xenium = `x_centroid`/`y_centroid`
- **Cell types:** CosMx and most Xenium panels have `ct_minor`/`ct_major`; Xenium noseg/withseg/colon have NO annotations (k-means only)
- **Primary integration metric:** Patient recovery (batch = patient ID) is the key metric since matched design provides ground truth
- **R packages:** Install under `~/R/libs_R441` (R/4.4.1 module on cayuga); SeuratObject, Matrix, sp, Rcpp
