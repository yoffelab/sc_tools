# sc_tools — Scientific analysis tool for spatial and single-cell multiomics

**sc_tools** is a **scientific analysis tool** for **spatial and single-cell multiomics**: a Python package and pipeline for spatial transcriptomics (e.g. Visium) and single-cell RNA-seq. It provides reusable utilities (plotting, statistics, colocalization, I/O) in a scanpy-style API and is used here for **lung tumor progression and Tertiary Lymphoid Structure (TLS) transcriptomics**—identifying transcriptional and spatial differences across **Normal**, **Non-Solid (GGO)**, and **Solid** lung tumors, with a focus on TLS heterogeneity and tumor–immune niches.

The project is **buildable as an installable package** via `pyproject.toml` for use as a dependency or as a base for custom spatial/single-cell workflows.

---

## Installation

From the repository root:

```bash
# Editable install (recommended for development)
pip install -e .

# Full pipeline (deconvolution: scvi-tools, tangram)
pip install -e ".[deconvolution]"

# With optional GPU support
pip install -e ".[deconvolution,gpu]"
```

Conda users can still use `conda env create -f environment.yml` and then `pip install -e .` inside that environment. See `ENVIRONMENT_SETUP.md` for details.

---

## What this project does

- **Data:** Converts raw Visium and scRNA-seq into standardized AnnData; integrates with **scVI** for batch-corrected latent space; scores gene signatures (Seurat-based via Scanpy) from `metadata/gene_signatures.json`.
- **Spatial analysis:** Compares tumor and immune programs across tumor types (1-vs-rest); analyzes **process colocalization** (Pearson, Moran’s I, neighborhood enrichment); visualizes signature heatmaps and volcano plots; studies **macrophage–proliferation** colocalization and TLS B-cell/T-cell states.
- **Deconvolution:** Estimates cell-type proportions in spots via **Tangram** (with optional Cell2location/DestVI fallback), in batches per `library_id` for memory efficiency.
- **Outputs:** Processed AnnData objects (`results/`), manuscript-ready figures (`figures/manuscript/`), and CSV statistics. All comparisons use **FDR (Benjamini–Hochberg)** and follow the statistical and plotting rules in `skills.md`.

---

## What’s important

- **Phased pipeline:** Ingestion → QC & scVI integration → Gene scoring → Deconvolution → Spatial & TLS analysis → Visualization. See `Architecture.md` for data flow and script roles.
- **Reusable library:** The **`sc_tools`** Python package (scanpy-style API) provides generic plotting and tools so analysis scripts stay thin and reproducible:
  - **`st.pl`** — spatial plots, heatmaps/clustermaps, statistical annotations, volcano plots, versioned figure saving (PDF+PNG, dpi=300).
  - **`st.tl`** — Mann–Whitney, FDR, colocalization (Pearson, Moran’s I, neighborhood enrichment), deconvolution helpers, versioned `write_h5ad`.
  - **`st.data`** — caching, I/O. **`st.memory`** — profiling, GPU detection.
- **Standards:** Statistical rigor (significance bars, asterisks, adjusted p-values), no over-filtering of low-count GGO spots, and batch processing for deconvolution are mandatory. See `skills.md` and `Mission.md`.

---

## What’s still being worked on

- **Deconvolution:** Pinpointing minimum data size and memory limits for Cell2location/DestVI; refining gene signatures and validating against HLCA/MSigDB/TCGA LUAD.
- **Spatial:** Spatially variable genes (SVG), cell-type colocalization from deconvolution, spatial transition areas (Normal ↔ Non-Solid ↔ Solid, tumor–TLS), and macrophage state comparison across tumor types.
- **TLS:** TLS niche extraction and lymphoid-rich neighborhood analysis.

---

## Repository layout

| Path | Purpose |
|------|--------|
| `sc_tools/` | Reusable Python package (pl, tl, data, memory, utils) |
| `sc_analysis/` | Project-specific analysis modules (placeholder) |
| `scripts/` | Analysis scripts (preprocessing, deconvolution, spatial, figures) |
| `metadata/` | Gene signatures (JSON), sample metadata |
| `results/` | Processed AnnData (e.g. `adata.img.genescores.h5ad`), CSVs |
| `figures/` | Output plots (manuscript, spatial, TLS, process_colocalization, etc.) |

Key docs: **`Mission.md`** (objectives and task list), **`Architecture.md`** (phases and data flow), **`Journal.md`** (decisions and rationale), **`skills.md`** (statistical and coding standards).

---

## Setup (environments)

- **Conda:** `conda env create -f environment.yml`, then `conda activate ggo_visium` and `pip install -e .` to install the package.
- **pip only:** `python -m venv venv`, activate, then `pip install -e ".[deconvolution]"` (or `pip install -r requirements.txt` for a non-package install).

See `ENVIRONMENT_SETUP.md` for details. Core stack: **scanpy**, **squidpy**, **anndata**; optional **scvi-tools**, **tangram**, **cell2location** for deconvolution.

---

## Running the pipeline

The **Makefile** defines targets in execution order (ingestion → preprocessing → scoring → deconvolution → spatial/figures). Run `make` or `make help` and use the annotated targets for reproducibility.

---

## License and contact

Internal research project. For questions or collaboration, contact the repository maintainers or the Yoffe Lab.
