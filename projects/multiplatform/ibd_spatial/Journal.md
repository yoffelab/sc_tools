# Journal: ibd_spatial

## 2026-03-11 (M2 followup: scANVI, Scanorama, success criteria)

**M2 followup completed** (SLURM job 2702091, 25 min on GPU A40):

All 6 methods now benchmarked on M2 (CosMx 6k + Xenium 5K, 164K cells, 2,552 genes):

| Method | Batch Score | CT Broad ASW | Platform Entropy |
|--------|-------------|-------------|-----------------|
| scVI | **0.992** | 0.009 | **0.632** |
| Harmony | 0.921 | -0.069 | 0.455 |
| scANVI | 0.881 | **0.074** | 0.386 |
| Scanorama | 0.837 | -0.089 | 0.081 |
| PCA | 0.805 | -0.073 | 0.023 |
| BBKNN | 0.805 | -0.073 | 0.484 |

**scANVI findings:**
- Semi-supervised (pretrain scVI 50 epochs, then 30 scANVI epochs with celltype_broad labels)
- Best bio conservation: ct_broad_ASW=0.074 (only positive value besides scVI's 0.009)
- Prediction accuracy: 94.9% overall; Xenium=99.1%, CosMx=91.3%
- Cross-platform kNN transfer (broad): CosMx->Xenium=98.7%, Xenium->CosMx=85.6% — dramatically better than unsupervised methods
- Trade-off: lower batch correction (0.881) than scVI (0.992) — the supervised signal pulls same-type cells together regardless of platform mixing

**Scanorama findings:**
- Fixed API: `integrate_scanpy()` modifies in-place (returns None), not `correct_scanpy()`
- Poor performance: batch_score=0.837, worst bio conservation (ct_broad_ASW=-0.089)
- Platform entropy=0.081 — almost no cross-platform mixing
- Not competitive for this cross-platform use case

**IBD marker genes (8/8 canonical present):**
EPCAM, CD3E, CD68, MS4A1, MKI67, FOXP3, IFNG, TNF — all in 2,552-gene intersection.
Extended: T cell 8/8, Myeloid 6/6, IBD-specific 4/7 (NOD2, ATG16L1, IL23R, CARD9).

**Formal success criteria (all PASS for scVI):**
- ASW_batch < 0.5: 0.008 PASS
- Platform entropy > 0.5: 0.632 PASS
- Bio conservation >= 0: 0.009 PASS
- IBD markers >= 5/8: 8/8 PASS
- CD vs UC: N/A (no CD in M2 4-patient group)

**gimVI dropped:** Removed from scvi-tools 1.4.2; LLOKI also dropped (Python 3.8 incompatible).

**Bug fixes in followup script:**
1. scANVI categorical fillna: `.astype(str).replace("nan", "Unknown")` instead of `.fillna("Unknown")`
2. Scanorama: `integrate_scanpy()` returns None (in-place); iterate `adatas_by_platform` not return value

**Decision: scVI is the recommended M2 integration method** — best batch correction by far, positive bio conservation, only method passing platform entropy criterion. scANVI is the best choice if celltype preservation is prioritized over batch mixing.

## 2026-03-07

**Task:** Project setup and Step 0 metadata inspection.

**Data inventory confirmed (from NC_compare_11122025_meta.csv):**
- CSV is readable (51 rows = 7 panels x N samples, header + data)
- Columns: `sample_id, block, tissue_type, disease, disease_state, cosmx_flow_cell_folder, fov_start, fov_end, cosmx_machine, xenium_MT_folder, aws_folder, panel`
- Diagnosis confirmed: `disease` = {CD, UC, Healthy}, `disease_state` = {Inflamed, Non-inflamed, Healthy}
- Tissue: Ileum and Rectum (16-patient group); Rectum only (4-patient group)

**Critical finding: patient-matched design.**
The same tissue blocks appear across multiple panels. Confirmed matched pairs:
- 16 patients: CosMx_1k_16 ↔ Xenium_MT_16 (same block IDs in CSV rows 1-16 for each panel)
- 4 patients (I0400/I0355/I0380/I0387): CosMx_6k, Xenium_5K, Xenium_MT_withseg, Xenium_MT_noseg, Xenium_colon

**Data path correction:** RDS files are in `/athena/project-saha/data_IBD/` (NOT `/athena/project-saha/`). Flat-file Xenium and CosMx directories are in `/athena/project-saha/data/`.

**Filename corrections (actual names):**
- CosMx 1k: `{1-16}_smi_so_qcnorm.RDS`
- CosMx 6k: `{1-4}_smi_so_qcnorm.RDS`
- Xenium 5K: `{1-4}_xen_so_qcnorm.RDS`
- Xenium MT 16: `{1-16}_xen_so_qcnorm.RDS`
- Xenium MT noseg: `{1-4}_xen_377_no_seg_so_qcnorm.RDS`
- Xenium MT withseg: `{1-4}_xen_377_with_seg_so_qcnorm.RDS`
- Xenium colon: `{1-4}_xen_colon_with_seg_so_qcnorm.RDS`

**CRITICAL S0 FINDING: Most RDS files are zero-filled (2026-03-07):**

Ran binary inspection (`xxd`, Python) on all panel files on cayuga. Confirmed:
- **CosMx_6k** (`*_smi_so_qcnorm.RDS`): magic bytes `1f8b` = valid gzip RDS; `readRDS()` returns Seurat object
- **CosMx_1k, Xenium_MT_16, Xenium_5K, Xenium_MT_withseg, Xenium_MT_noseg, Xenium_colon**: all null bytes (zeros from start to end of file)

All files have identical modification timestamps (2026-03-06 18:49) = same bulk S3 sync. CosMx_6k (430MB) completed; all others (36-180MB) either failed or were pre-allocated but never written. This is a partial S3-to-Lustre sync failure.

**Action required:** Contact Saha lab (jip2007) to re-sync failed panels from S3, OR request AWS credentials to download directly from `s3://saha-dulai-collaboration/NC_compare_11122025/`.

**SeuratObject issue:** System R on cayuga login node is missing `SeuratObject` package. Seurat class confirmed (readRDS works), but inspection halted. Must install SeuratObject or use a conda env with Seurat before full inspection.

**SLURM issues encountered:** cayuga SLURM protocol mismatch without `module load slurm`; compute nodes lack Rscript (system R only on login node); log output paths must use `/home/fs01/juk4007/` not `/home/juk4007/`.

**CosMx readRDS failure:** First attempt at `readRDS()` on `1_smi_so_qcnorm.RDS` returned "unknown input format". CosMx RDS files may be saved with `save()` (RData format) not `saveRDS()`. The updated R script uses `load()` fallback.

**Plan revisions made:**
- Reframed as patient-matched cross-platform evaluation (not generic feasibility)
- New Milestone 0: noseg vs withseg technical baseline (same 4 patients)
- New M1: CosMx_1k vs Xenium_MT_16 (16 matched patients, key milestone)
- Disease limitation documented: 4-patient group has NO CD
- Batch key hierarchy revised: patient_id is anchor, not batch

**Next:**
1. Confirm SLURM job ran; if not, resubmit when connectivity stabilizes
2. Write convert_rds_to_h5ad.R and test on 1 sample
3. Populate batch manifests from NC_compare CSV
4. HPC directory setup (soft links)

## 2026-03-07 (continued)

**S0 blocker resolved:** All 7 panels have valid gzip RDS files (re-synced by jip2007).

**SeuratObject installed:** `~/R/libs` on cayuga login node. R anndata also installed.

**Full S0 inspection completed on all 7 panels:**
- Seurat v5 format: assay layers `counts` (raw int), `data` (log-norm), `scale.data`
- Assay names: `Nanostring` (CosMx), `Xenium` (Xenium) — NOT `RNA`
- Cell types: CosMx has `ct_minor`/`ct_major`/`ct_minor_new`; Xenium 5K/MT have `ct_minor`/`ct_major`;
  Xenium noseg/withseg/colon have NO cell type annotations (only k-means clusters)
- Spatial: `GetTissueCoordinates()` works on all panels (x, y, cell)
- noseg: 28,858 cells present (Xenium default segmentation)
- No disease/diagnosis columns inside RDS — must join from CSV
- Xenium colon: 322 genes (not ~400)

**Key Seurat v5 findings for conversion:**
- `@layers$counts` has NO dimnames; gene/cell names from `rownames(assay_obj)`/`colnames(assay_obj)`
- R `anndata` package fails on `var` — switched to MTX+CSV intermediate + Python assembly
- CSV has BOM and panel names like "CosMx 1K" (space, not underscore)

**Conversion script written and tested:**
- `convert_rds_to_h5ad.R`: R extracts raw counts (MTX), obs (CSV), spatial + embeddings (CSV);
  Python assembles into h5ad. Tested on Xenium MT sample 1: 5,138 cells x 377 genes, 2.9 MB.
  Disease metadata matched from CSV (patient=I0294, CD, Inflamed, Ileum).

**SLURM batch job submitted:** Job 2700475 (array 0-51). First attempt (2700476) failed: Rscript
not in PATH on compute nodes. Fixed: added `module load R/4.4.1`. Resubmitted but SSH
connection to cayuga lost before results could be verified.

## 2026-03-07 (batch conversion)

**SLURM debugging (3 failed attempts):**
- Job 2700581: `sp.so` cannot find `libR.so` — R packages in `~/R/libs` compiled against system R, incompatible with R/4.4.1 module
- Job 2700633: Added `LD_LIBRARY_PATH=/opt/ohpc/pub/software/R/4.4.1/lib64` — insufficient, `Rcpp.so` same issue
- Job 2700692: Same — the `.so` files have compiled-in rpaths, not just `LD_LIBRARY_PATH` dependent

**Fix: Reinstalled all R packages under R/4.4.1 module:**
- Job 2700691: Installed SeuratObject + Matrix + sp into `~/R/libs_R441`
- Job 2700744: Also installed Rcpp (missed dependency)
- Both completed successfully; `library(SeuratObject)` works under R/4.4.1

**SSH config updated:** Added `ConnectTimeout 60`, `ServerAliveInterval 30`, `ServerAliveCountMax 5` to `~/.ssh/config` for cayuga and brb.

**Batch conversion SUCCESS — Job 2700745:**
All 52 tasks completed (exit 0, ~24s each). 52 `adata.p0.h5ad` files produced.

Summary by panel:
| Panel | N | Genes | Example size |
|-------|---|-------|-------------|
| CosMx 1k | 16 | 950 | 2.1–146 MB |
| CosMx 6k | 4 | 6,175 | 134–181 MB |
| Xenium MT | 16 | 377 | 2.9–38 MB |
| Xenium 5K | 4 | 5,001 | 5.2–84 MB |
| Xenium noseg | 4 | 377 | 3.1–19 MB |
| Xenium withseg | 4 | 377 | 8.8–22 MB |
| Xenium colon | 4 | 322 | 6.2–19 MB |

**Validation (cosmx_1k_01):** X.max=70 (raw int), spatial present, disease=CD/Inflamed/Ileum, celltype/celltype_broad present.
**Validation (xenium_noseg_01):** X.max=102 (raw int), spatial present, disease=Healthy, celltype correctly absent.

**Next:** Proceed to Milestone 0 (noseg vs withseg integration benchmark, 4 matched patients).

## 2026-03-08 (Milestone 0)

**Conda env `sc_tools` created on cayuga** (SLURM job 2700798):
Python 3.11, scanpy 1.11.5, scvi-tools 1.4.2, harmonypy, bbknn, scanorama, squidpy, leidenalg.

**M0 benchmark completed** (SLURM job 2700824, 1h32m on CPU):
- 8 samples: 4 noseg + 4 withseg (same 4 patients, Rectum, 377 genes)
- 155,035 cells concatenated; 118,792 after QC filtering
- 4 integration methods run (Scanorama failed silently):

| Method | Batch ASW | Patient Mix ASW | Batch Score |
|--------|-----------|-----------------|-------------|
| scVI | 0.003 | 0.020 | **0.997** |
| PCA | -0.007 | -0.019 | 0.993 |
| BBKNN | -0.007 | -0.019 | 0.993 |
| Harmony | -0.008 | -0.023 | 0.992 |

**Key findings:**
- Near-zero batch effect between noseg and withseg (all methods batch_score > 0.99)
- Harmony converged in 2 iterations (nothing to correct)
- scVI best at patient separation while maintaining panel mixing
- Disease signal preserved: Healthy (I0400) clearly separates from UC on UMAP
- scVI overfitted on CPU (100 epochs, loss went from 66.9 at ep10 to 69.4 at ep100; early stopping did not trigger — likely no validation split configured)

**Outputs:**
- `results/m0_benchmark/m0_benchmark.csv`
- `results/m0_benchmark/adata.m0.h5ad` (586 MB, 118K cells x 377 genes)
- `figures/m0/m0_umap_grid.png` + per-method PNGs

**Next:** Milestone 1 — CosMx 1k (16) + Xenium MT (16), 32 samples, ~377 shared genes. This is the key milestone with real cross-platform batch effects and CD+UC biology.

## 2026-03-08 (Milestone 1)

**M1 benchmark completed** (SLURM job 2700857, 2h7m on GPU A40):
- 32 samples: 16 CosMx 1k (950 genes) + 16 Xenium MT (377 genes), same 16 patients
- 693,780 cells loaded; 380,058 after QC filtering (45% removed, mostly low-count CosMx)
- **Only 119 shared genes** between panels (genuine panel overlap, not a naming issue)
- 4 integration methods run (Scanorama failed on `X_scanorama` key):

| Method | Batch ASW | Patient Mix | Celltype ASW | Celltype Broad ASW | Batch Score |
|--------|-----------|-------------|--------------|-------------------|-------------|
| Harmony | 0.029 | 0.029 | -0.162 | -0.008 | **0.971** |
| scVI | 0.039 | 0.065 | -0.138 | 0.001 | 0.961 |
| PCA | 0.065 | 0.065 | -0.127 | 0.024 | 0.935 |
| BBKNN | 0.065 | 0.065 | -0.127 | 0.024 | 0.935 |

**Key findings:**
- Cross-platform batch effect is surprisingly small (all batch_scores > 0.93)
- Harmony best at platform mixing (0.971); scVI second (0.961)
- **Negative celltype ASW across all methods** — 119 shared genes insufficient for fine-grained cell type separation; broad celltypes near zero (better)
- Harmony converged in 2 iterations again — consistent with M0, suggesting platform batch effects are modest
- scVI trained 50 epochs on GPU in 40 min (vs 5h timeout on CPU in first attempt)

**Disease biology preserved (Fisher exact tests, all p < 0.05):**
- DC: CD=6.5% vs UC=0.9% (CD-enriched, consistent with literature)
- Goblet: CD=8.9% vs UC=2.5% (CD/ileum signature)
- Macrophage: CD=5.6% vs UC=11.6% (UC-enriched)
- Transient amplifying: CD=4.1% vs UC=17.9% (UC-enriched)
- Inflammatory fibroblast: CD=0% vs UC=9.3% (UC-specific)
- Plasma: CD=3.6% vs UC=5.3%
- 21 cell types showed significant CD vs UC differences

**First attempt (job 2700834) timed out** after 6h: scVI on CPU at 91s/epoch x 200 epochs = 5h+ for scVI alone. Fixed by switching to GPU (A40) and reducing max_epochs to 50.

**Gene intersection finding:** CosMx 1k (950 genes) and Xenium MT (377 genes) share only 119 genes. This is genuine — the panels target different gene sets. For M2 (CosMx 6k x Xenium 5K), the intersection should be much larger (~1500-2000 genes).

**Outputs:**
- `results/m1_benchmark/m1_benchmark.csv`
- `results/m1_benchmark/m1_celltype_by_disease.csv`
- `results/m1_benchmark/adata.m1.h5ad` (1.1 GB, 380K cells x 119 genes)
- `figures/m1/m1_umap_grid.png` + per-method PNGs

**Next:** Milestone 2 — CosMx 6k (4) + Xenium 5K (4), same 4 patients, ~1500-2000 shared genes. Higher gene count should improve celltype separation.

## 2026-03-09 (M1 bio evaluation + integration reports)

**M1 bio evaluation completed** (SLURM job 2701016, 3 min on CPU with 50K subsample):

Key finding: **negative celltype ASW is a gene count problem, not an integration artifact.**
- Per-platform celltype ASW also negative: CosMx=-0.115, Xenium=-0.040 (119 genes cannot resolve 40+ fine types even within a single platform)
- Broad celltypes (8 types): Xenium within-platform ASW=+0.141 (positive), CosMx=-0.031 (marginal)
- Cross-platform kNN label transfer: celltype_broad CosMx->Xenium=73.6%, Xenium->CosMx=28.6% (asymmetric — CosMx labels more transferable)
- Per-celltype platform mixing: Epithelial shows most platform separation (ASW 0.24-0.27); Myeloid/Mesenchymal mix well (ASW~0)

**Integration reports generated** (SLURM job 2701063, 29 min):
- Installed `sc_tools` in cayuga conda env (editable install from synced repo)
- Used `sc_tools.bm.integration.compare_integrations()` for proper scib-style metrics
- Generated HTML reports via `sc_tools.bm.report.generate_integration_report()`

M1 scib metrics (batch_weight=0.4, bio_weight=0.6):
| Method | ASW Batch | PCR | Graph Conn | ASW CT | ARI | NMI | Overall |
|--------|-----------|-----|-----------|--------|-----|-----|---------|
| PCA | 0.923 | 0.985 | 0.815 | 0.517 | 0.084 | 0.225 | **0.528** |
| scVI | 0.970 | 0.982 | 0.778 | 0.506 | 0.069 | 0.202 | 0.519 |
| Harmony | 0.954 | 0.996 | 0.762 | 0.498 | 0.068 | 0.168 | 0.509 |

**Notable:** PCA (unintegrated) has highest overall score because it best preserves biological signal. Integration improves batch mixing but costs bio signal with only 119 genes. This trade-off should resolve with more shared genes (M2: ~1500-2000).

**Reports:**
- `figures/QC/m0_integration_report.html`
- `figures/QC/m1_integration_report.html`
- `figures/QC/m1_integration_report_celltype.html`

## 2026-03-09 (Milestone 2 + Project Report)

**M2 benchmark completed** (SLURM job 2701598, 54 min on GPU A40):
- 8 samples: 4 CosMx 6k (6,175 genes) + 4 Xenium 5K (5,001 genes), same 4 patients (Rectum)
- **2,552 shared genes** (vs 119 in M1 — 21x more)
- 164,392 cells after QC (only 55 removed — high-quality data)
- HVG selection: 2,000 / 2,552 genes used for PCA
- scVI: n_latent=20, n_hidden=128, n_layers=2, 100 epochs on A40 (4m14s), early stopping did not trigger
- Scanorama failed again (same `X_scanorama` key issue)
- scVI `get_latent_representation()` failed with format string error — but training completed and latent was extracted before the error, so scVI results are valid

| Method | Batch ASW | Patient Mix | Celltype ASW | Celltype Broad ASW | Batch Score |
|--------|-----------|-------------|--------------|-------------------|-------------|
| scVI | 0.008 | 0.015 | -0.070 | **0.009** | **0.992** |
| Harmony | 0.079 | 0.072 | -0.135 | -0.069 | 0.921 |
| PCA | 0.195 | 0.179 | -0.159 | -0.073 | 0.805 |
| BBKNN | 0.195 | 0.179 | -0.159 | -0.073 | 0.805 |

**Key findings:**
- scVI is now the clear winner (0.992 batch score) — with 2,552 genes and higher model capacity (n_latent=20, n_layers=2), it dominates
- M2 shows MORE platform batch effect than M1 (PCA batch_ASW=0.195 vs 0.065), despite more genes — the high-plex panels reveal more platform-specific expression patterns
- Celltype broad ASW is near zero for scVI (0.009) — marginal but no longer negative
- Fine celltype ASW still negative (-0.070 for scVI) — 30 types still too many even with 2,552 genes
- Harmony converged in 2 iterations again (consistent across M0/M1/M2)

**M2 vs M1 comparison:**
- Gene count: M1=119 -> M2=2,552 (21x)
- Best batch score: M1=0.971 (Harmony) -> M2=0.992 (scVI) — scVI overtakes Harmony with more genes
- Celltype broad ASW: M1=0.024 (PCA) -> M2=0.009 (scVI) — slight decrease, but now positive for the best integration method
- The hypothesis that more genes would help was partially confirmed: scVI batch correction improved, but celltype separation did not dramatically improve

**Cross-platform label transfer (M2):**
- Broad types: PCA CosMx->Xenium=84.0%, Xenium->CosMx=36.8%
- Harmony: CosMx->Xenium=78.3%, Xenium->CosMx=39.6%
- scVI: CosMx->Xenium=66.6%, Xenium->CosMx=31.9%
- Transfer accuracy is HIGHER in M2 than M1 for broad types, but asymmetry persists
- scVI has lower kNN transfer than PCA/Harmony — it creates a more mixed embedding space where nearest neighbors are more diverse (a feature for integration, not a bug)

**Per-platform celltype ASW (M2):**
- CosMx celltype_broad: PCA=-0.089, scVI=-0.032 (negative even within platform — CosMx cell type labels may be noisier)
- Xenium celltype_broad: PCA=+0.094, scVI=+0.055 (positive — Xenium labels more self-consistent)
- This platform asymmetry in label quality explains the asymmetric kNN transfer

**Disease composition (UC vs Healthy only — no CD in 4-patient group):**
- Transient amplifying: Healthy=31.6% vs UC=25.0% (dominant cell type in both)
- Inflammatory fibroblast: Healthy=3.7% vs UC=7.6% (UC-enriched, consistent with M1)
- Plasma: Healthy=10.6% vs UC=9.4% (similar)
- Endothelial: Healthy=8.3% vs UC=10.7% (UC slightly higher)

**Comprehensive project report generated** (SLURM job 2701599, 11 min):
- Standalone HTML with embedded UMAPs, bio eval results, scib metrics
- Includes project context, M0-to-M1 progression, BBKNN handling notes
- `figures/QC/ibd_spatial_integration_report.html` (34.5 MB with embedded images)

**Outputs:**
- `results/m2_benchmark/m2_benchmark.csv`
- `results/m2_benchmark/m2_bio_platform_asw.csv`
- `results/m2_benchmark/m2_bio_label_transfer.csv`
- `results/m2_benchmark/adata.m2.h5ad` (164K cells x 2,552 genes)
- `figures/m2/m2_umap_grid.png` + per-method PNGs
- `figures/QC/ibd_spatial_integration_report.html`

**Next:** Milestone 3 — full cross-platform integration (all panels, ~100-300 shared genes). Or focus on the scVI M2 embedding for downstream biology (UC vs Healthy spatial patterns).
