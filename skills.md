# Skills for Single-Cell and Spatial Omics Analysis

Author: Junbum Kim

This document defines actionable, reproducible skills and best-practice guidelines for agentic scientific programming in computational oncology. It serves as a standard for AI-assisted analysis, code reviews, and repository maintenance, bridging legacy Scanpy/Squidpy workflows with modern probabilistic modeling and high-resolution spatial frameworks.

---

# Part I: Analysis Pipeline

Core computational steps in execution order (Phases 1-6 of the pipeline).

---

## 1. Data Ingestion, Integrity, and Multi-Modal Containers

### Phase 0 → AnnData / SpatialData Flow

Phase 0 is split into two sub-steps:
- **Phase 0a:** Run platform tools (Space Ranger, Xenium Ranger, IMC pipeline, CosMx export) → `data/{sample_id}/outs/`
- **Phase 0b:** Load per-sample into portable format → `data/{sample_id}/adata.p0.h5ad` (required) and/or `data/{sample_id}/spatialdata.zarr` (optional, for Visium HD / Xenium)

Use `sc_tools.ingest.loaders` for Phase 0b loading. Each loader sets `obs['sample']`, `obs['library_id']`, `obs['raw_data_dir']`, `obsm['spatial']`.

| Modality | Phase 0b loader | Notes |
|----------|----------------|-------|
| Visium | `load_visium_sample()` | H&E in `uns['spatial']` |
| Visium HD | `load_visium_hd_sample()` | 8um bins default; parquet positions |
| Visium HD Cell | `load_visium_hd_cell_sample()` | SpaceRanger 4 cell segmentation |
| Xenium | `load_xenium_sample()` | spatialdata-io preferred; scanpy fallback |
| IMC | `load_imc_sample()` | Reads segmented h5ad |
| CosMx | `load_cosmx_sample()` | Flat CSV/Parquet or RDS (rpy2+anndata2ri); centroid coords |

Phase 1 then loads all per-sample `adata.p0.h5ad` files, applies per-sample QC, and concatenates via `concat_samples()` → `results/adata.raw.p1.h5ad`.

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

### Modality-Specific Normalization (Phase 3)

| Modality | Normalization | Integration | Feature Selection | Special Notes |
|----------|--------------|-------------|-------------------|---------------|
| Visium / Visium HD | Raw counts to scVI (no manual norm for VAE path); `normalize_total` + `log1p` for non-VAE | scVI (default) | HVG with `batch_key`; optionally intersect with SVG | Backup `adata.raw`; filter MT/RP/HB genes |
| Xenium | `normalize_total` + `log1p` | Standard PCA or Harmony/scVI for multi-sample | HVG on log-normalized data | Tune normalization per dataset; n_PCs critical for KNN |
| CosMx | `normalize_total` + `log1p` | Standard or Harmony/scVI | HVG | Phase 0b: flat CSV/Parquet or RDS → `load_cosmx_sample()`; detection algorithm changed 2023-2024 |
| IMC | `arcsinh(X/5)` -- NOT `log1p` | CytoVI (scvi-tools, totalVI-inspired) | All markers (30-50 proteins); no HVG step | Segmentation quality is critical upstream step |

**Implementation:** `sc_tools.pp.preprocess(adata, modality=..., integration=..., ...)` dispatches to modality-specific recipes. Individual steps (`backup_raw`, `normalize_total`, `log_transform`, `arcsinh_transform`, `filter_genes_by_pattern`, `pca`, `cluster`, etc.) are all importable from `sc_tools.pp` for fine-grained control.

### GPU Acceleration
- **rapids-singlecell:** Drop-in scanpy replacement with 5-100x speedups for preprocessing, PCA, neighbors, clustering. `sc_tools.pp` auto-detects GPU availability via `get_backend()` and uses rapids-singlecell when present, scanpy otherwise. No API change.
- Install: `pip install rapids-singlecell`

---

## 4. Feature Selection

### Core Skills
- Identify highly variable genes (HVGs) prior to dimensionality reduction.
- Tune feature counts based on dataset size and biological complexity.
- Utilize spatially aware methods (e.g., `squidpy.gr.spatial_autocorr`) for discovering spatially variable genes.
- For scVI integration: use `seurat_v3` HVG flavor (works on raw counts). For non-VAE pipelines: use `seurat` flavor on log-normalized data.

---

## 5. Dimensionality Reduction and Latent Space

### Core Steps
- **scVI:** Prioritize deep generative modeling for dimensionality reduction and batch integration.
- **Harmony:** Fast PCA-based correction (`harmonypy`); use when scVI is too slow or not applicable.
- **PCA:** Use for initial visualization or as input for traditional graph construction if VAEs are not applicable.
- Choose components based on variance explained (elbow plots) or reconstruction error.

### Auto-detection of Representation
`sc_tools.pp.neighbors()` auto-detects `use_rep` in priority order: `X_scVI` > `X_cytovi` > `X_pca_harmony` > `X_pca`.

### Validation
- Inspect latent space mixing using metrics like ASW (Average Silhouette Width) or iLISI.
- Confirm graph connectivity is biologically reasonable and not driven by technical batches.

---

## 6. Clustering and Annotation

### Core Skills
- Use graph-based methods (Leiden or Louvain).
- **Spatial Clustering (UTAG):** Unsupervised Tissue Architecture with Graphs -- graph-based message passing on spatial adjacency to identify microanatomical domains. Runs *after* standard Leiden as an optional complement. `max_dist`: 10-20 for IMC, 10-100 for transcriptomics. `slide_key` for multi-image batch processing. Install: `pip install git+https://github.com/ElementoLab/utag.git@main`.
- **Automated Validation:** Cross-reference identified marker genes against established cell-type databases (e.g., CellTypist) to assign biological meaning.

---

# Part II: Spatial and Integrative Analysis

Spatial-specific methods, multi-modal integration, and imaging modality considerations.

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

## 9. Imaging-Based Spatial Modalities (IMC, CosMx, Xenium)

### Core Skills
- Validate segmentation masks/quality prior to expression analysis.
- Handle "binless" or subcellular data by linking molecular coordinates to cell identifiers or neighborhood grids.
- Perform biological sanity checks using known tissue structures (e.g., Keratin in epithelium).

---

# Part III: Statistical Analysis and Visualization

Differential expression, statistical testing, color standards, and publication figure production.

---

## 10. Differential Expression and Spatially Variable Features

### Differential Expression (DE)
- Use non-parametric tests or model-based DE (e.g., `model.differential_expression()` in `scvi-tools`).
- **Strict Correction:** Always apply Benjamini-Hochberg (FDR) adjustment for multiple hypothesis testing.
- Report effect sizes (log-fold change) alongside adjusted p-values.

### Statistical Visualization
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

## 11. Colors and Palettes

### 11.1. Categorical Variables in AnnData (Scanpy Convention)

For any categorical column `{name}` in `adata.obs`, use the same color convention as Scanpy for all figures:

- **Storage:** Colors are stored in `adata.uns[f'{name}_colors']` as a **list of hex strings** (e.g. `['#FFFFFF', '#66c2a5', ...]`) with **length equal to the number of categories**, in **category order** (same order as `adata.obs[name].cat.categories`).
- **Usage:** When producing a figure that colors by a categorical obs column:
  - If `adata.uns[f'{name}_colors']` **already exists** and has the correct length, **use those colors** and do not overwrite.
  - If it **does not exist** (or length does not match), **create** a default palette (e.g. qualitative colormap), assign one color per category in order, **save** the hex list to `adata.uns[f'{name}_colors']`, and use those colors for the figure.
- **Consistency:** This ensures the same categorical variable is colored identically across all plots and scripts. Prefer setting `adata.uns[f'{name}_colors']` explicitly in project code when a specific palette is desired (e.g. tumor type: Normal / Non-Solid / Solid).

**Implementation:** In `sc_tools.pl.heatmaps`, use `get_obs_category_colors(adata, obs_col, store_if_missing=True)` to obtain a mapping from category value to RGB; it reads or creates `adata.uns[f'{obs_col}_colors']` as above. Apply this pattern in any plotting code that colors by a categorical obs column (heatmaps, scatter, legends).

### 11.2. Color-Blind Safe Palettes

All figures must use color-blind safe palettes. Never rely on red-green contrast alone.

#### Categorical: Okabe-Ito / Wong palette (default for <= 8 categories)

| Color | Name | Hex |
|-------|------|-----|
| ![](.) | Black | `#000000` |
| ![](.) | Orange | `#E69F00` |
| ![](.) | Sky Blue | `#56B4E9` |
| ![](.) | Bluish Green | `#009E73` |
| ![](.) | Yellow | `#F0E442` |
| ![](.) | Blue | `#0072B2` |
| ![](.) | Vermillion | `#D55E00` |
| ![](.) | Reddish Purple | `#CC79A7` |

```python
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"]
```

#### Categorical: Paul Tol palettes (for more categories)

- **Bright (7):** `#4477AA`, `#EE6677`, `#228833`, `#CCBB44`, `#66CCEE`, `#AA3377`, `#BBBBBB`
- **Vibrant (7):** `#EE7733`, `#0077BB`, `#33BBEE`, `#EE3377`, `#CC3311`, `#009988`, `#BBBBBB`
- **Muted (10):** `#CC6677`, `#332288`, `#DDCC77`, `#117733`, `#88CCEE`, `#882255`, `#44AA99`, `#999933`, `#AA4499`, `#DDDDDD`

### 11.3. Sequential and Diverging Colormaps

| Use case | Recommended | Avoid |
|----------|-------------|-------|
| Sequential (expression, counts) | `viridis`, `mako`, `cividis` | `jet`, `rainbow`, `hot` |
| Diverging (fold change, z-score) | `RdBu_r`, `coolwarm`, `PiYG` | `bwr` (poor perceptual uniformity) |
| Spatial heatmaps | `magma`, `inferno`, `viridis` | `jet` |

**Rule:** Never use `jet` or `rainbow` colormaps. These are not perceptually uniform and are inaccessible to colorblind viewers. Use `viridis` (default), `cividis` (deuteranopia-safe), or `mako`/`magma` for sequential data.

---

## 12. Publication-Quality Figure Production

Guidelines for producing manuscript-ready figures for Nature family, Cell family, and Science/AAAS journals. All projects generating figures for `figures/manuscript/` must follow these standards. For color palette selection, see Section 11.

### 12.1. Dimensions and Resolution

| Property | Nature | Cell | Science |
|----------|--------|------|---------|
| **Single column** | 89 mm (3.5 in) | 85 mm (3.35 in) | 57 mm (2.24 in) |
| **1.5 column** | 120 mm (4.7 in) | — | 121 mm (4.76 in) |
| **Double column (full width)** | 183 mm (7.2 in) | 174 mm (6.85 in) | 184 mm (7.24 in) |
| **Max height** | 247 mm (9.7 in) | 229 mm (9 in) | 229 mm (9 in) |
| **Min DPI (raster/photo)** | 300 | 300 | 300 |
| **Min DPI (line art)** | 600-1000 | 600 | 600 |
| **Preferred format** | EPS, PDF, TIFF | EPS, PDF, TIFF | PDF, EPS |
| **Color mode** | RGB (not CMYK) | RGB | RGB |

**Default in sc_tools:** Save figures at **300 DPI** (raster) or vector (PDF/SVG) at **single-column width (89 mm)**. Use `plt.savefig(path, dpi=300, bbox_inches="tight")` and set figure size explicitly in inches.

### 12.2. Typography

| Property | Nature | Cell | Science |
|----------|--------|------|---------|
| **Font family** | Helvetica, Arial, or sans-serif | Helvetica, Arial | Helvetica, Arial |
| **Axis labels / text** | 5-7 pt | 6-8 pt | 7 pt min |
| **Panel labels** | 8 pt, **bold lowercase** (a, b, c) | 8 pt, **bold UPPERCASE** (A, B, C) | 8 pt, **bold UPPERCASE** (A, B, C) |
| **Title text** | 8 pt max | 8-10 pt | 8 pt max |

**Matplotlib defaults for sc_tools:**

```python
import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7,                # Base font size (7 pt)
    "axes.titlesize": 8,           # Subplot titles
    "axes.labelsize": 7,           # Axis labels
    "xtick.labelsize": 6,          # Tick labels
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "figure.titlesize": 8,
    "pdf.fonttype": 42,            # TrueType (editable text in PDF)
    "ps.fonttype": 42,
    "svg.fonttype": "none",        # Keep text as text in SVG
})
```

**Panel labeling convention:** Default to **bold lowercase** (a, b, c) for Nature-family submissions. Switch to **bold UPPERCASE** (A, B, C) for Cell/Science. Place panel labels in the **top-left corner** outside the plot area. Use `fig.text(x, y, "a", fontsize=8, fontweight="bold")` or equivalent annotation.

### 12.3. Statistical Reporting in Figures

For significance bar definitions and comparison modes, see Section 10.

#### General requirements (all journals)

- **Report exact P-values** when possible (e.g., P = 0.003), not just thresholds.
- **Define error bars** explicitly in figure legends: SD, SEM, 95% CI, or IQR.
- **Report sample sizes** (n) for each group in the figure legend or on the plot.
- **Name the statistical test** in the legend (e.g., two-sided Wilcoxon rank-sum, one-way ANOVA with Tukey HSD).
- **Multiple testing correction:** Always Benjamini-Hochberg FDR for multi-gene/multi-comparison analyses (see Section 10).

#### Journal-specific notes

| Requirement | Nature | Cell | Science |
|-------------|--------|------|---------|
| P-value reporting | Exact values preferred; asterisks acceptable with definition | Exact values in text; asterisks in figures | Exact values preferred |
| Error bars | Must define in legend (SD, SEM, or CI) | Must define; 95% CI preferred | Must define |
| Sample size | n for each group in legend | n in figure or legend | n in legend |
| Effect size | Encouraged (e.g., Cohen d, log2FC) | Required for key comparisons | Encouraged |

When using `statannotations`, configure: `annotator.configure(test="Mann-Whitney", text_format="star", loc="inside")`.

### 12.4. Line Weights, Borders, and Whitespace

| Element | Specification |
|---------|---------------|
| **Axis lines** | 0.5-1.0 pt |
| **Plot borders** | 0.5 pt or remove (borderless preferred for spatial plots) |
| **Tick marks** | 0.5 pt, outward |
| **Error bar caps** | 0.5-1.0 pt |
| **Scale bars** | 1.0 pt, solid black, with text label (e.g., "100 um") |
| **Grid lines** | Remove or use 0.25 pt light gray (`#CCCCCC`) |
| **Whitespace** | Minimize; use `plt.tight_layout()` or `bbox_inches="tight"` |

**Scale bars for spatial plots:** Required for all spatial/tissue images. Place in the **bottom-right** corner. Label with physical units (um or mm). Never rely on axis ticks alone for spatial scale.

### 12.5. Multi-Panel Figure Assembly

#### Using marsilea (recommended for complex composite figures)

**Marsilea** is a declarative, composable visualization library built on matplotlib for assembling multi-panel figures with cross-layout alignment. Recommended for heatmaps with side annotations, grouped comparisons, and composite biological figures.

Key capabilities:
- **Heatmap composites:** `Heatmap`, `CatHeatmap`, `SizedHeatmap` with aligned side plots (dendrograms, bar charts, labels)
- **Cross-layout paradigm:** Central plot with top/bottom/left/right annotation layers via `.add_top()`, `.add_left()`, etc.
- **Multi-dataset alignment:** `ClusterBoard` for side-by-side heatmaps sharing row/column groupings
- **Biological plots:** OncoPrint, UpSet plots, stacked bars, arc diagrams
- **Scanpy integration:** Works with AnnData for expression heatmaps

```python
import marsilea as ma
import marsilea.plotter as mp

# Example: heatmap with dendrogram and grouped annotations
h = ma.Heatmap(data, cmap="viridis", width=4, height=6)
h.add_left(mp.Labels(row_labels))
h.add_top(mp.Dendrogram(linkage))
h.add_right(mp.Bar(scores, color="#0072B2"))
h.render()
```

Install: `pip install marsilea` or `pip install sc-tools[viz]`.

#### Manual assembly with matplotlib

For simpler layouts, use `matplotlib.gridspec` or `fig.subfigures()`:

```python
fig = plt.figure(figsize=(7.2, 9.0))  # Full width, near max height (Nature)
gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
ax_a = fig.add_subplot(gs[0, 0])
# ... add panels
# Panel labels
for ax, label in zip(axes, "abcdef"):
    ax.text(-0.1, 1.1, label, transform=ax.transAxes, fontsize=8, fontweight="bold")
```

### 12.6. Supplementary Figure Standards

Supplementary figures follow the same quality requirements as main figures:
- Same DPI (300+ for raster, vector preferred)
- Same font sizes and families
- Same color-blind safe palettes
- Same statistical annotation standards
- Panel labels continue from main figures or restart with "S" prefix (e.g., S1a, S1b or per journal convention)
- Maximum page size: typically **letter** (8.5 x 11 in) or **A4** (210 x 297 mm)

### 12.7. Figure Checklist (pre-submission)

Before placing any figure in `figures/manuscript/`:

- [ ] Resolution >= 300 DPI (raster) or vector format (PDF/SVG/EPS)
- [ ] Fonts are Helvetica/Arial, 5-8 pt range
- [ ] Panel labels present (bold, correct case for target journal)
- [ ] Color palette is color-blind safe (no red-green only contrast; no jet/rainbow)
- [ ] Scale bars on all spatial/tissue images with physical unit labels
- [ ] Axis labels include units where applicable
- [ ] Error bars defined in legend (SD, SEM, CI)
- [ ] Sample sizes (n) reported per group
- [ ] Statistical test named in legend
- [ ] P-values shown (exact preferred) or significance bars defined
- [ ] White background (no gray figure background)
- [ ] No unnecessary gridlines or chart junk
- [ ] Figure saved with `bbox_inches="tight"` to minimize whitespace

---

# Part IV: Engineering, Reproducibility, and CI

Workflow standards, metadata management, testing, and continuous improvement.

---

## 13. Reproducibility and Workflow Standards

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

## 14. Metadata and Data Sharing

### Core Skills
- Maintain complete, machine-readable metadata in a structured format (e.g., CSV or YAML).
- Share processed objects (AnnData/SpatialData) alongside raw data for peer review.
- Include all spatial coordinate transformations and image references.

---

## 15. Validation and Interpretation

### Best Practices
- Cross-validate computational findings across different platforms or orthogonal modalities (IHC/IF).
- Interpret results strictly within the context of known biology and experimental controls.
- Avoid over-interpreting clustering or spatial domains without independent validation.

---

## 16. Testing and CI

- Try to generate fail-proof test samples associated with code that is generated. Try to perform unit and integration test where applicable.
- Generate a wide coverage of inputs to generated code, including with empty, sub, and full data.
- Always test generated code so that it **compiles, passes lint, and runs** without error.
- **CI:** Use **GitHub Actions** to run the test suite on push/PR; the workflow should run lint (e.g. `ruff check`), the project test suite (e.g. `pytest`), and optionally Snakemake validation (e.g. `snakemake -n` dry-run or run with fixtures). Fail if any step fails. Document the workflow location (e.g. `.github/workflows/tests.yml`) and how to run the same commands locally.

---

## 17. Continuous Improvement

- Track the evolution of community standards (e.g., Scverse updates).
- Revisit and update pipelines as new benchmarks and probabilistic models emerge.
- Treat all analysis pipelines as evolving, versioned software products.

---

# Appendix

## Common Libraries

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
- **marsilea** (composable multi-panel figures: heatmaps with side annotations, cross-layout alignment, OncoPrint, UpSet; built on matplotlib)

### Workflow & Infrastructure
- **Workflow engine:** `snakemake` (all environments). **Container:** `apptainer`/Singularity (Linux/HPC, primary); `docker` (macOS/Windows, fallback).
- Also: `zarr`, `dask`

### CI/CD, Linting & Docs
- `ruff` (lint + format), `pytest`, `sphinx`, `snakemake`, GitHub Actions

## Key References

- [sc-best-practices: Integration](https://www.sc-best-practices.org/cellular_structure/integration.html)
- [rapids-singlecell docs](https://rapids-singlecell.readthedocs.io/en/latest/)
- [scvi-tools models](https://scvi-tools.org/)
- [UTAG: Nature Methods (2022)](https://www.nature.com/articles/s41592-022-01657-2)
