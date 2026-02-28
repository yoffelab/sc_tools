---
name: sc-tools-skills
description: Apply sc_tools coding and analysis standards from repository skills.md. Use when writing or running single-cell/spatial analysis code, pipelines, statistical plots, CI, or sandbox runs in sc_tools and subprojects.
---

# sc_tools Skills and Standards

When working in this repository on single-cell or spatial omics analysis, code generation, pipelines, or sandbox execution, **read and follow the standards in the repository root file `skills.md`**.

## When to Use This Skill

- Writing or reviewing analysis code (Scanpy, scvi-tools, Squidpy, etc.)
- Defining or running workflows (Snakemake, containers)
- Producing statistical figures (significance bars, FDR, group comparisons)
- Setting up or running sandbox or local pipeline runs
- Implementing QC, normalization, clustering, or spatial analysis

## Defaults for Sandbox and Local Runs

- **Container:** Use **Apptainer** for all sandbox, local, and remote HPC runs .
- **Workflow:** Use **Snakemake** as the workflow engine for all environments (sandbox, local, production, CI). Each project has its own Snakefile.

## Quick Reference

- **Data:** AnnData (single-modality), SpatialData (multi-modal); Dask/Zarr for large data.
- **Modeling:** Prefer scvi-tools (VAEs) over log-normalization for integration.
- **Stats:** Benjamini–Hochberg (FDR); significance bars and asterisks on plots (* < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001); 2 groups = pairwise, >2 groups = 1-vs-rest by default.
- **Colors:** Use `adata.uns[f'{name}_colors']` for categorical obs; do not overwrite if already set and length matches.
- **Containers:** Default to Apptainer.

## Full Reference

For complete sections (data ingestion, QC, normalization, DE, reproducibility, testing, CI), see **`skills.md`** at the repository root.
