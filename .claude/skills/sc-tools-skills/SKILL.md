---
name: sc-tools-skills
tier: 1
description: "[Tier 1 — Core] Apply sc_tools coding and analysis standards. Use when writing or reviewing single-cell/spatial analysis code, pipelines, statistical plots, CI, or sandbox runs."
allowed-tools: Read Glob Grep
---

# sc_tools Skills and Standards

Read and follow the full standards in @docs/skills.md.

## When to use

- Writing/reviewing analysis code (Scanpy, scvi-tools, Squidpy, SpatialData)
- Defining or running workflows (Snakemake, Apptainer)
- Statistical figures with significance bars
- QC, normalization, clustering, spatial analysis
- HPC cluster work (brb, cayuga)

## Defaults

- **Container:** Apptainer (Linux/HPC, primary); Docker (macOS/Windows, fallback). Auto-detect via `scripts/run_container.sh`.
- **Workflow:** Snakemake (all environments — sandbox, local, production, CI)
- **Run script:** `scripts/run_container.sh <project_path> [command]`
- **Snakemake:** `--use-apptainer` flag (or `--use-singularity`)
- **GPU clustering:** rapids-singlecell preferred (falls back to scanpy when GPU unavailable)

## Quick reference

| Topic | Rule |
|-------|------|
| Data containers | AnnData (single), SpatialData (multi-modal) |
| Large data | Dask-backed arrays + Zarr |
| Modeling | scvi-tools VAEs for integration (raw counts, not log-normed) |
| FDR | Always Benjamini-Hochberg |
| Significance | `*` < 0.05, `**` < 0.01, `***` < 0.001, `****` < 0.0001 |
| 2 groups | Pairwise comparison |
| >2 groups | 1-vs-rest (default); all-pairwise as option |
| Colors | `adata.uns[f'{name}_colors']`; do not overwrite if length matches |
| Linting | Ruff; never commit failing lint |
| Tests | Unit + integration; fail-proof with empty/sub/full fixtures |
| Checkpoints | Standard names only — see docs/Architecture.md §2.1 |
| Signature scores | `obsm['signature_score']` / `obsm['signature_score_z']`; not in obs |
