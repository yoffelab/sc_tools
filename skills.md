# Skills for Single-Cell and Spatial Omics Analysis

Author: Junbum Kim

This document defines actionable, reproducible skills and best-practice guidelines for agentic scientific programming in computational oncology. It serves as a standard for AI-assisted analysis, code reviews, and repository maintenance, bridging legacy Scanpy/Squidpy workflows with modern probabilistic modeling and high-resolution spatial frameworks.

---

## 1. Data Ingestion, Integrity, and Multi-Modal Containers

### Core Skills
- Load data into standardized containers: **AnnData** for single-modality, **MuData** for multi-modal (multi-omics) in-memory, and **SpatialData** for multi-modal datasets (images, masks, points, and expression).
- **Scalability:** Utilize **Dask-backed** arrays and **Zarr** storage to handle "out-of-memory" issues common with Visium HD and Xenium. For GPU-accelerated preprocessing and PCA, consider **rapids-singlecell** (RAPIDS) where compatible.
- Maintain coordinate systems (physical microns vs. pixel space) for high-resolution platforms.
- Ensure consistent cell/spot identifiers across all modalities and verify presence of expression matrices, metadata, and spatial coordinates.

### Required Checks
- `adata.X` is non-empty, numeric, and preferably raw counts for probabilistic models.
- `adata.obs` and `adata.var` are fully populated with relevant metadata.
- Spatial datasets include valid coordinate columns in `adata.obsm['spatial']`.

---

## 2. Quality Control and Filtering

### Single-Cell RNA-seq
- Compute per-cell QC metrics (counts, genes, mitochondrial fraction).
- Evaluate doublets using robust methods such as `scDblFinder` or `Solo` (within `scvi-tools`).
- Filter low-quality cells using explicit, documented thresholds. Avoid fixed cutoffs that might bias against rare, high-metabolic states.

### Spatial Transcriptomics
- Assess tissue coverage and spot/cell quality.
- Avoid aggressive filtering that removes biologically informative low-count areas in sparse tissues.
- Record QC metrics before and after filtering to ensure transparency.

---

## 3. Normalization and Transformation

### Traditional Guidelines
- Apply library-size normalization followed by `log1p` transformation for exploratory analysis.
- Scale data only when required for specific non-generative downstream models.

### Modern Generative Practices
- **scVI-tools:** Use raw counts for training Variational Autoencoders (VAEs) to learn a latent space that accounts for technical noise and library size without manual scaling.
- Ensure normalization logic does not introduce spatial artifacts or over-correct biological variance.

---

## 4. Feature Selection

### Core Skills
- Identify highly variable genes (HVGs) prior to dimensionality reduction.
- Tune feature counts based on dataset size and biological complexity.
- Utilize spatially aware methods (e.g., `squidpy.gr.spatial_autocorr`) for discovering spatially variable genes.

---

## 5. Dimensionality Reduction and Latent Space

### Core Steps
- **scVI:** Prioritize deep generative modeling for dimensionality reduction and batch integration.
- **PCA:** Use for initial visualization or as input for traditional graph construction if VAEs are not applicable.
- Choose components based on variance explained (elbow plots) or reconstruction error.

### Validation
- Inspect latent space mixing using metrics like ASW (Average Silhouette Width) or iLISI.
- Confirm graph connectivity is biologically reasonable and not driven by technical batches.

---

## 6. Clustering and Annotation

### Core Skills
- Use graph-based methods (Leiden or Louvain).
- **Spatial Clustering:** Incorporate spatial proximity into clustering using **Banksy** or **GraphST** to identify distinct tissue domains.
- **Automated Validation:** Cross-reference identified marker genes against established cell-type databases (e.g., CellTypist) to assign biological meaning.

---

## 7. Integration of Single-Cell and Spatial Data

### Strategies
- **Deconvolution:** Estimate cell-type mixtures in spatial spots using **cell2location**, **tangram**, or **DestVI**.
- **Optimal Transport:** Align multi-omic or temporal datasets using **Moscot**.
- **Mapping:** Project single-cell identities into spatial contexts while preserving original embeddings for comparison.

---

## 8. Spatial Statistics and Neighborhood Analysis

### Core Skills
- Build spatial neighbor graphs using physical distances, Delaunay triangulation, or grid topology.
- Quantify spatial enrichment/depletion using permutation-based methods.
- **Signaling:** Infer spatially constrained ligand-receptor interactions via `squidpy.gr.ligprec`.

---

## 9. Differential and Spatially Variable Features

### Differential Expression (DE)
- Use non-parametric tests or model-based DE (e.g., `model.differential_expression()` in `scvi-tools`).
- **Strict Correction:** Always apply Benjamini-Hochberg (FDR) adjustment for multiple hypothesis testing.
- Report effect sizes (log-fold change) alongside adjusted p-values.

### Robust Statistical Visualization
For all statistical comparisons (boxplots, violin plots, strip plots), the following must be enforced:
- **Visual Annotation:** Draw significance bars and asterisks strictly following:
  - $* < 0.05$
  - $** < 0.01$
  - $*** < 0.001$
  - $**** < 0.0001$
- **Numerical Reporting:** Write numerical adjusted p-values in a clear, non-interfering area of the figure (e.g., top-right corner or bottom margin). If there is more than one comparison, list Group A - Group B: p*-val.
- **Comparison Modes:** If there is more than one comparison, stack dodging bars
  - **Two groups:** Perform pairwise comparison.
  - **More than two groups:** Default to **1-vs-rest** comparisons.
  - **Option:** Allow a toggle for **all-pairwise** comparisons when required.

---

## 10. Imaging-Based Spatial Modalities (IMC, CosMx, Xenium)

### Core Skills
- Validate segmentation masks/quality prior to expression analysis.
- Handle "binless" or subcellular data by linking molecular coordinates to cell identifiers or neighborhood grids.
- Perform biological sanity checks using known tissue structures (e.g., Keratin in epithelium).

---

## 11. Reproducibility and Workflow Standards

### Sandbox and local runs (defaults)
- **Always use Snakemake as the workflow engine for all runs** (sandbox, local, and production). Each project has its own Snakefile. Run rules inside the project container image so that execution is reproducible and consistent.
- Prefer Snakemake for one-off and development runs; use a project Snakefile unless the user or project explicitly overrides (e.g. bare Python or Makefile).

### Required Practices
- Use version control (Git) for all analysis code and notebook checkpoints.
- **Workflow Manager:** **Snakemake** is the workflow engine for all environments (dev, sandbox, production, CI). Each project has a Snakefile implementing the phase-dependent pipeline (Phases 1-7 documented in Mission.md, Architecture.md, README). Makefile may coexist for convenience targets. Run Snakemake dry-run in CI to validate the workflow.
- **Containerization:** Use **Apptainer/Singularity** as the primary container runtime (Linux/HPC). Use **Docker** as the fallback for macOS and Windows where Apptainer is not natively available. The pipeline **auto-configures** via `scripts/run_container.sh`: detect platform and select the appropriate runtime. Define one container image per pipeline; publish to a registry (e.g. Docker Hub, GHCR) so HPC can pull via Apptainer. Ensure environment parity across local and HPC runs.
- **Packaging and distribution:** Make package compatible with pip, uv, and poetry; support **PyPI deployment** so the package can be installed via `pip install <package_name>`. Use `pyproject.toml` with build-backend and versioning; publish via GitHub Actions on release (trusted publishing or secrets), never store PyPI credentials in the repo.
- **Linting:** Apply a consistent linter (e.g. **Ruff** or flake8) and optionally a formatter (e.g. Ruff format or Black); config in repo (`pyproject.toml` or `ruff.toml`). Do not commit code that fails the configured lint check; fix or explicitly ignore with justification.
- **API documentation:** Generate and maintain API docs with **Sphinx** (autodoc from docstrings); public modules and functions should have docstrings that appear in the built docs. Document build command (e.g. `make docs`) and where built docs live (`docs/` source, `_build/html`, or deployed URL).
- Capture all software versions and seed parameters to ensure exact reproducibility.
- Produce and edit a Makefile for python scripts in the order in which all scripts should be processed for unit and integrated figure production in a streamlined fashion for reproducibility.
- Heavily annotate the Makefile to describe what each python script does, and what the input and output of each python script is.

---

## 11.5. Categorical Variables and Colors in AnnData (Scanpy Convention)

### Rule
For any categorical column `{name}` in `adata.obs`, use the same color convention as Scanpy for all figures:

- **Storage:** Colors are stored in `adata.uns[f'{name}_colors']` as a **list of hex strings** (e.g. `['#FFFFFF', '#66c2a5', ...]`) with **length equal to the number of categories**, in **category order** (same order as `adata.obs[name].cat.categories`).
- **Usage:** When producing a figure that colors by a categorical obs column:
  - If `adata.uns[f'{name}_colors']` **already exists** and has the correct length, **use those colors** and do not overwrite.
  - If it **does not exist** (or length does not match), **create** a default palette (e.g. qualitative colormap), assign one color per category in order, **save** the hex list to `adata.uns[f'{name}_colors']`, and use those colors for the figure.
- **Consistency:** This ensures the same categorical variable is colored identically across all plots and scripts. Prefer setting `adata.uns[f'{name}_colors']` explicitly in project code when a specific palette is desired (e.g. tumor type: Normal / Non-Solid / Solid).

### Implementation
- In `sc_tools.pl.heatmaps`, use `get_obs_category_colors(adata, obs_col, store_if_missing=True)` to obtain a mapping from category value to RGB; it reads or creates `adata.uns[f'{obs_col}_colors']` as above.
- Apply this pattern in any plotting code that colors by a categorical obs column (heatmaps, scatter, legends).

---

## 12. Metadata and Data Sharing

### Core Skills
- Maintain complete, machine-readable metadata in a structured format (e.g., CSV or YAML).
- Share processed objects (AnnData/SpatialData) alongside raw data for peer review.
- Include all spatial coordinate transformations and image references.

---

## 13. Validation and Interpretation

### Best Practices
- Cross-validate computational findings across different platforms or orthogonal modalities (IHC/IF).
- Interpret results strictly within the context of known biology and experimental controls.
- Avoid over-interpreting clustering or spatial domains without independent validation.

---

## 14. Continuous Improvement

- Track the evolution of community standards (e.g., Scverse updates).
- Revisit and update pipelines as new benchmarks and probabilistic models emerge.
- Treat all analysis pipelines as evolving, versioned software products.

---

## 15. Fail-proof Testing and CI
- Try to generate fail-proof test samples associated with code that is generated. Try to perform unit and integration test where applicable.
- Generate a wide coverage of inputs to generated code, including with empty, sub, and full data.
- Always test generated code so that it **compiles, passes lint, and runs** without error.
- **CI:** Use **GitHub Actions** to run the test suite on push/PR; the workflow should run lint (e.g. `ruff check`), the project test suite (e.g. `pytest`), and optionally Snakemake validation (e.g. `snakemake -n` dry-run or run with fixtures). Fail if any step fails. Document the workflow location (e.g. `.github/workflows/tests.yml`) and how to run the same commands locally.

---

## Appendix: Common Libraries

### Data Containers & Multi-Modal
- **AnnData** (single-modality), **MuData** (multi-modal, multi-omics; mudata), **SpatialData** (spatial multi-modal with images, shapes, points)

### Modeling & Integration
- `scvi-tools`, `moscot`, `cell2location`, `harmonypy`

### Spatial Analysis & Frameworks
- `spatialdata`, `squidpy`, `tangram`, `scvelo`

### GPU-Accelerated & Scalability
- **rapids-singlecell** (RAPIDS for single-cell: GPU-accelerated preprocessing, PCA, clustering; scanpy-like API where applicable)
- `cuml`, `cucim` (RAPIDS ecosystem for ML and imaging)

### Visualization & Statistics
- `scanpy`, `matplotlib`, `seaborn`, `statannotations`, `napari`, `vitessce`, `pinguoin`

### Workflow & Infrastructure
- **Workflow engine:** `snakemake` (all environments). **Container:** `apptainer`/Singularity (Linux/HPC, primary); `docker` (macOS/Windows, fallback).
- Also: `zarr`, `dask`

### CI/CD, Linting & Docs
- `ruff` (lint + format), `pytest`, `sphinx`, `snakemake`, GitHub Actions
