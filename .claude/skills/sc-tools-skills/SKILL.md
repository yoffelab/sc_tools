---
name: sc-tools-skills
description: Apply sc_tools coding and analysis standards. Use when writing or reviewing single-cell/spatial analysis code, pipelines, statistical plots, CI, or sandbox runs.
---

# sc_tools Skills and Standards

Read and follow the full standards in @skills.md.

## When to use

- Writing/reviewing analysis code (Scanpy, scvi-tools, Squidpy, SpatialData)
- Defining or running workflows (Snakemake, Apptainer)
- Statistical figures with significance bars
- QC, normalization, clustering, spatial analysis

## Defaults

- **Container:** Apptainer (Linux/HPC, primary); Docker (macOS/Windows, fallback). Auto-detect via `scripts/run_container.sh`.
- **Workflow:** Snakemake (all environments — sandbox, local, production, CI)
- **Run script:** `scripts/run_container.sh <project_path> [command]`
- **Snakemake:** `--use-apptainer` flag (or `--use-singularity`)

## Quick reference

| Topic | Rule |
|-------|------|
| Data containers | AnnData (single), SpatialData (multi-modal) |
| Large data | Dask-backed arrays + Zarr |
| Modeling | scvi-tools VAEs for integration (raw counts, not log-normed) |
| FDR | Always Benjamini–Hochberg |
| Significance | `*` < 0.05, `**` < 0.01, `***` < 0.001, `****` < 0.0001 |
| 2 groups | Pairwise comparison |
| >2 groups | 1-vs-rest (default); all-pairwise as option |
| Colors | `adata.uns[f'{name}_colors']`; do not overwrite if length matches |
| Linting | Ruff; never commit failing lint |
| Tests | Unit + integration; fail-proof with empty/sub/full fixtures |

$ARGUMENTS
