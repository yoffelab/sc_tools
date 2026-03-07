# Journal: ibd_spatial

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
