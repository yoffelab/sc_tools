# Journal Summary: ibd_spatial

## Project scope

Patient-matched cross-platform spatial transcriptomics integration study.
Saha lab data: same tissue blocks measured by CosMx (1k and 6k) and Xenium (5K, 377-MT, colon panel).
Scientific question: can cross-platform batch correction recover patient identity and IBD biology?

## Key facts (confirmed 2026-03-07)

- **Data path:** RDS files in `/athena/project-saha/data_IBD/`; flat files in `/athena/project-saha/data/`
- **Diagnosis confirmed:** `disease` (CD/UC/Healthy), `disease_state` (Inflamed/Non-inflamed/Healthy)
- **Matched design:** 16 patients in both CosMx_1k AND Xenium_MT_16; same 4 patients in 5 high-plex panels
- **Disease caveat:** 4-patient high-plex group = UC only + 1 Healthy (no CD)
- **Tissue:** 16-patient group: Ileum + Rectum; 4-patient group: Rectum only
- **CosMx RDS format:** May use `save()` not `saveRDS()` -- use `load()` fallback in R

## BLOCKER: Most RDS files zero-filled (2026-03-07)

- **CosMx_6k only**: valid gzip RDS (magic `1f8b`); readRDS works; SeuratObject class confirmed
- **All other panels** (CosMx_1k, Xenium_MT_16, Xenium_5K, withseg, noseg, colon): entirely zero bytes -- partial S3-to-Lustre sync failure on 2026-03-06 18:49
- **Action**: Contact Saha lab (jip2007) to re-sync, OR get AWS creds for `s3://saha-dulai-collaboration/NC_compare_11122025/`
- **SeuratObject** not installed in system R on cayuga login node -- install before full inspection

## HPC notes (cayuga)

- Home dir: `/home/fs01/juk4007/` (NOT `/home/juk4007/`)
- SLURM requires `module load slurm` before sbatch/squeue
- System R at `/usr/bin/Rscript` on login node only (not on compute nodes)
- Use `--wrap` with full path to Rscript in SLURM jobs

## Milestone order

- M0: Xenium noseg vs withseg (same 4 patients, technical replicate)
- M1: CosMx_1k vs Xenium_MT_16 (16 matched patients, main cross-platform test)
- M2: CosMx_6k vs Xenium_5K (4 matched patients, high-plex)
- M3: Full integration (all 44 samples, imputation via gimVI)
