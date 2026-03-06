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

### 12.8. Figure Intent, Insight, and Readability

Every figure must communicate a specific insight. A figure without a clear takeaway is not a meaningful figure. These principles apply to ALL figure production — original work and reproduction alike.

#### 12.8.1 Figure Intent and Insight (General)

- **Every figure has a claim.** Before coding, state in one sentence what the reader should conclude from each panel. If you cannot articulate the insight, the figure is not ready to produce.
- **Direction of effect matters.** When a panel shows a statistical comparison (box plot, bar plot, violin), verify: which group is higher/lower? Is the difference significant? The direction and significance must match the stated insight. Example: if the insight is "Cytotoxic LME has worse survival", the KM curve for Cytotoxic must be below others, and the log-rank p-value must be significant. If the bars go the wrong way, something is wrong with the data, the grouping, or the claim — investigate before publishing.
- **Write the caption first.** The caption is the contract between the code and the narrative. Draft it before coding. Each panel gets a sentence describing what it shows and what the reader should see.
- **Validate insight after generation.** After producing a figure, read it as a reviewer would: Does the visual tell the story the caption claims? Do the statistics support the direction of effect? If a box plot shows "Group A is higher" but the medians are equal, the insight is wrong or the figure needs fixing.

#### 12.8.2 Data Presentation for Readability (General)

- **Heatmaps**: Z-score normalize (row or column depending on comparison axis). Raw intensity values are almost never publication-readable. Row-normalize when comparing features across groups; column-normalize when comparing groups across features.
- **Color scales**: Diverging palettes (RdBu_r, coolwarm) for z-scored data centered at 0. Sequential palettes (viridis, magma, YlOrRd) for proportions, counts, and unsigned metrics. Never use jet or rainbow.
- **Row/column ordering**: Use hierarchical clustering OR biological category grouping (T cells together, myeloid together, stromal together). Never use arbitrary or alphabetical ordering for biological data.
- **Consistent colors**: Define a project-level color dictionary for major categories (e.g., cell types, patient groups, conditions). Store in a shared `figure_config.py` and import in every figure script. Reuse across ALL figures for visual coherence.
- **Violins/boxes**: Show individual data points when n < 50. Add BH-corrected significance annotations (Wilcoxon rank-sum for 2-group, Kruskal-Wallis + Dunn for multi-group).
- **KM curves**: Number-at-risk table below the plot. Exact log-rank p-value on the plot face. n per group in the legend (e.g., "Cold (n=115)").
- **Bar plots**: Define error bars explicitly (SD, SEM, or 95% CI). Label n per group.
- **UMAP/embeddings**: Rasterize points for PDF size. Axis labels (UMAP1, UMAP2). Color-blind safe palettes: Okabe-Ito for <=8 categories, Paul Tol for >8.
- **Forest plots**: Point estimate (HR) with 95% CI. Reference group label. Sort by effect size or biological grouping.

#### 12.8.3 Figure Reproduction (Manuscript-Specific)

When reproducing figures from a published manuscript or aligning with a collaborator's existing figures:

- **Align with manuscript claims.** Read the manuscript text AND figure captions. Map every panel in the caption to a subplot in the script. Missing panels = incomplete figure.
- **Cross-reference Methods.** Use exact methodology from the paper (clustering parameters, statistical tests, normalization). Do not substitute methods without documenting the deviation.
- **Quantitative validation.** If the manuscript says "5 LME classes including Cold (35%)", the output must show all 5 classes with proportions matching within 5%. If the manuscript reports p=0.0064 for OS, the reproduced p-value should be in the same order of magnitude.
- **Direction-of-effect validation.** Compare every statistical panel against the manuscript text. If the paper says "Cytotoxic LME is enriched for M1 macrophages", the violin/bar for M1 in Cytotoxic must be the highest. If it is not, trace the discrepancy (wrong grouping variable, label mismatch, data subset difference) before proceeding.

#### 12.8.4 Iterative Figure QA

- After generating, visually compare with the original or caption description.
- Per-figure checklist: Does it tell the same story? All panels present? Axes readable at print size? Statistical claims match direction and significance?
- Fix and regenerate until quality matches or exceeds the original.
- Document comparison: what matched, what diverged, and why (e.g., different cohort size, updated method).

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

## 18. HPC Cluster Workflow (brb / cayuga)

All production analysis runs on SLURM-managed clusters. This section defines how to structure jobs, exploit parallelism, use GPU acceleration, and orchestrate complex workflows — including multi-agent patterns.

---

### 18.1 Cluster Reference

| Cluster | SLURM | Scratch | CPU partitions | GPU partitions |
|---------|-------|---------|----------------|----------------|
| **brb** | UP | `/athena/elementolab/scratch/juk4007/` | `scu-cpu` (c7–c42) | `scu-gpu` (L40S / A40 / RTX6000) |
| **cayuga** | UP | `/athena/elementolab/scratch/juk4007/` | `scu-cpu` | `scu-gpu` (A100 / A40) |

**Critical rules:**
- **Never write to `~/`** (login node home). All outputs go to `/athena/elementolab/scratch/juk4007/` or project-relative paths under scratch.
- **Never run compute on the login node** — use `srun` (interactive) or `sbatch` (batch). Even `pip install` should be done in an `srun --pty bash` session.
- Scratch is shared Lustre (`/athena/elementolab/scratch/`). Read-heavy workflows (e.g. opening large `.h5ad` mid-job) are slow on Lustre — copy to `/tmp` on the compute node first.
- The `/athena/` Lustre mount is identical on brb and cayuga; files written by one cluster are immediately visible on the other.

**Conda on brb:**
```bash
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate sc_tools    # py3.11, scanpy 1.11.5, scvi 1.4.2, harmonypy
```

**GPU flags on brb** (NOT `gpu:a100:1` — brb has L40S/A40/RTX6000):
```bash
#SBATCH --gres=gpu:1
```

---

### 18.2 SLURM Job Patterns

#### Single job (preprocessing, integration, scoring)

```bash
#!/usr/bin/env bash
#SBATCH --job-name=sc_tools_preprocess
#SBATCH --partition=scu-cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=4:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

mkdir -p logs
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate sc_tools
cd /athena/elementolab/scratch/juk4007/sc_tools
python projects/visium_hd/robin/scripts/preprocess.py
```

#### Array jobs — embarrassingly parallel per sample (SpaceRanger, per-sample QC, cell segmentation)

```bash
#!/usr/bin/env bash
#SBATCH --job-name=spaceranger
#SBATCH --partition=scu-cpu
#SBATCH --array=0-15%4          # 16 samples, 4 concurrent
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/spaceranger_%A_%a.out

SAMPLES=($(cut -f1 metadata/phase0/all_samples.tsv | tail -n +2))
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
spaceranger count --id="$SAMPLE" ...
```

Key flags:
- `--array=0-N%K` — N+1 tasks, at most K concurrent (throttle for I/O-heavy jobs)
- `$SLURM_ARRAY_TASK_ID` — 0-based index into sample list
- `%A_%a` in output path — parent job ID + array index

#### Dependency chain (pipeline phases as SLURM jobs)

```bash
JOB1=$(sbatch --parsable spaceranger_array.sbatch)
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 ingest_concat.sbatch)
JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 preprocess.sbatch)
```

`afterok` — only start if all array tasks of the parent succeeded. Use `afterany` to proceed regardless of failure (e.g. for cleanup jobs).

#### Interactive GPU session (debugging, notebook exploration)

```bash
srun --partition=scu-gpu --gres=gpu:1 --cpus-per-task=8 --mem=64G --time=2:00:00 --pty bash
```

---

### 18.3 Parallelization by Workload Type

#### SpaceRanger / Xenium Ranger — array per sample

- One SLURM array task per sample. No inter-sample dependencies.
- Each task: 32 CPUs, 220–240G RAM, up to 48h wall time.
- Write outputs directly to scratch: `data/{sample_id}/outs/`.
- Use `%4`–`%8` concurrency limit to avoid Lustre saturation.
- Validate completion by checking for `outs/filtered_feature_bc_matrix.h5` before downstream rules.

#### Cell segmentation (SpaceRanger 4) — array per sample, GPU optional

- SR4 segmentation is CPU-only in standard mode; runs within the SpaceRanger array job.
- For custom Cellpose/StarDist segmentation: separate GPU array job.
- `--gres=gpu:1`, 4–8 CPUs, 32G RAM per task.

#### Per-sample QC and loading (Phase 0b → Phase 1) — array per sample

- Load per-sample adata, apply `filter_spots()`, save `adata.p0.h5ad`.
- Each task: 4–8 CPUs, 32–64G RAM (depending on modality; Visium HD needs more).
- Concatenation (`concat_samples()`) runs as a single downstream job after all per-sample jobs complete.

#### scVI integration — single GPU job

- scVI training is single-process; does not benefit from multi-GPU (ELBO computation is sequential across batches).
- Use 1 A40 or L40S GPU, 8–16 CPUs (data loading), 64–128G RAM.
- For large datasets (>500K cells): backed AnnData, increase `batch_size` in scVI, use `num_workers=4` for the data loader.
- `sc_tools.pp.run_scvi()` auto-detects GPU (`use_gpu=True` by default when CUDA available).

```bash
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
```

#### Harmony — CPU only, fast

- Harmony runs on PCA embeddings, not raw counts — no GPU needed.
- Single job, 8–16 CPUs, 32–64G RAM. Runtime: minutes for <1M cells.
- `sc_tools.pp.run_harmony()` handles the harmonypy call.

#### rapids-singlecell — **default** GPU-accelerated clustering and UMAP (post-integration)

**Use rapids-singlecell by default on any GPU node.** It is a drop-in replacement for scanpy `pp.neighbors`, `tl.umap`, `tl.leiden` and is faster even for moderate datasets (>50K cells); the speedup is decisive above 200K cells.

- A40 (48GB VRAM) handles ~1M cells comfortably for neighbors + UMAP + Leiden.
- `sc_tools.pp` auto-detects GPU; rapids path is taken whenever `rapids_singlecell` is importable.
- Workflow: integrate on GPU (scVI/Harmony) → stay in GPU memory → cluster + UMAP → transfer back before saving.
- Do not fall back to scanpy for these steps on GPU nodes — always use rapids.

```bash
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
```

```python
import rapids_singlecell as rsc
rsc.get.anndata_to_GPU(adata)          # move X and obsm to GPU
rsc.pp.neighbors(adata, use_rep="X_scVI")
rsc.tl.umap(adata)
rsc.tl.leiden(adata, resolution=0.5)
rsc.get.anndata_to_CPU(adata)          # move back before saving
```

#### Deconvolution (cell2location / Tangram) — GPU, per library

- **Always batch by `library_id`** to avoid OOM. One job per library or one array task per library.
- cell2location: 1 GPU, 16–32 CPUs, 128G RAM per library. Runtime: 1–4h per library.
- Tangram: 1 GPU, 8–16 CPUs, 64G RAM per library.
- Merge results after all libraries complete.

```bash
#SBATCH --array=0-7%2    # 8 libraries, 2 concurrent (GPU-limited)
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
```

---

### 18.4 Resource Guidelines

| Task | CPUs | RAM | GPU | Wall time |
|------|------|-----|-----|-----------|
| SpaceRanger count | 32 | 240G | — | 48h |
| Cell segmentation (SR4) | 32 | 220G | — | 24h |
| Per-sample QC / load | 4–8 | 32–64G | — | 1h |
| Concat + global QC | 8 | 64–128G | — | 2h |
| scVI training | 8–16 | 64–128G | 1 | 8–12h |
| Harmony | 8 | 32G | — | 30min |
| rapids neighbors+UMAP+leiden | 4 | 64G | 1 | 30min–2h |
| Deconvolution (per library) | 8–16 | 64–128G | 1 | 1–4h |
| Signature scoring | 8 | 32G | — | 30min |
| Phase 5 figures | 4 | 32G | — | 1h |

Always add 20–30% headroom to RAM estimates (Lustre buffering, Python overhead).

---

### 18.5 Monitoring and Debugging

```bash
# Live queue status
squeue -u juk4007
squeue -u juk4007 -o "%.10i %.12j %.8T %.10M %.6C %.8m %R"  # verbose

# Job accounting (finished/failed)
sacct -j JOBID --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS

# Live resource usage of a running job
sstat -j JOBID --format=MaxRSS,MaxVMSize,AveCPU

# Check GPU utilization interactively (ssh to compute node)
ssh NODENAME "nvidia-smi"
watch -n5 "squeue -u juk4007"

# Cancel all jobs
scancel -u juk4007

# Cancel a job array but keep completed tasks
scancel JOBID_ARRAYID
```

**Lustre log-checking pattern** (avoid slow NFS read in a tight loop):
```bash
ssh brb "cp /athena/elementolab/scratch/juk4007/logs/job_12345.out /tmp/ && cat /tmp/job_12345.out | tail -50"
```

---

### 18.6 Package Installation on HPC: Handling Slow or Hanging Conda

`conda install` on HPC Lustre filesystems frequently hangs or takes 10–30+ minutes due to (1) the SAT solver scanning thousands of packages, (2) Lustre locking when conda tries to lock the package cache, and (3) NFS latency on `~/.conda`. **Never wait — use the alternatives below.**

#### Decision tree

```
Need a Python package?
  ├── Already in sc_tools pyproject.toml extras?  →  pip install -e ".[extra]"  (fastest)
  ├── Pure Python / has a wheel on PyPI?           →  pip install PKG  or  uv pip install PKG
  ├── Needs compiled CUDA extension (rapids)?      →  pip + NVIDIA index (see below)
  ├── Needs a compiled system library?             →  Spack  (see below)
  └── Nothing else works?                         →  mamba/micromamba  (still better than conda)
```

#### Option 1 — pip / uv (preferred for Python packages)

pip resolves and installs in seconds; uv is even faster (Rust-based resolver).

```bash
# Plain pip — use inside the activated conda env
pip install rapids-singlecell

# uv — fastest resolver; compatible with conda envs
uv pip install rapids-singlecell
```

#### Option 2 — rapids-singlecell specifically: NVIDIA PyPI index

conda install of RAPIDS reliably hangs or produces solver conflicts. **Use pip with NVIDIA index instead:**

```bash
# CUDA 12 (A40 / L40S / A100 on brb/cayuga)
pip install "rapids-singlecell[rapids12]" --extra-index-url=https://pypi.nvidia.com

# CUDA 11 fallback
pip install "rapids-singlecell[rapids11]" --extra-index-url=https://pypi.nvidia.com

# Or install pre-releases (bundled CUDA kernels, no toolkit needed at install time)
pip install --pre rapids-singlecell-cu12 --extra-index-url=https://pypi.nvidia.com
```

Check CUDA version first: `nvcc --version` or `nvidia-smi | grep "CUDA Version"`.

#### Option 3 — mamba / micromamba (drop-in conda replacement, ~10x faster resolver)

When a package truly requires conda (e.g. compiled bioinformatics tools with complex native deps):

```bash
# mamba — install once into the base conda env
conda install -n base -c conda-forge mamba

# Then use mamba everywhere instead of conda
mamba install -c conda-forge scanpy

# micromamba — statically linked, no base env required
# Install to scratch (avoid home dir):
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)   # follow prompt; set prefix to scratch
micromamba install -c conda-forge -n sc_tools scanpy
```

**Lustre lockfile hang fix** — if mamba/micromamba still hangs on Lustre, disable lockfiles:
```bash
micromamba config set use_lockfiles false
# or for mamba:
conda config --set use_locks false
```

Move the conda package cache off Lustre to local node scratch:
```bash
conda config --add pkgs_dirs /tmp/conda_pkgs_$USER
```

#### Option 4 — Spack (system-level compiled libraries only)

Use Spack when you need a compiled system library (libhdf5, OpenMPI, a compiled bioinformatics tool) and neither conda nor pip can provide it. **Do not use Spack for Python packages — pip/conda handles those better.**

```bash
spack find hdf5                    # list available versions
spack load hdf5@1.12               # add to PATH / LD_LIBRARY_PATH for current shell
spack load --sh hdf5@1.12          # emit shell commands (useful in sbatch)
spack unload hdf5

# In an sbatch script:
. $(spack location -i hdf5@1.12)/share/spack/setup-env.sh
spack load hdf5@1.12
```

Useful Spack packages on brb/cayuga: `hdf5`, `openmpi`, `cuda`, `libffi`, `zlib`, `star`, `bwa`, `samtools`.

#### Option 5 — Apptainer container (nuclear option)

If a package is impossible to install in the conda env (version conflicts, missing compilers), pull an official container image:

```bash
# rapids-singlecell official container (CUDA + full stack)
apptainer pull rsc.sif ghcr.io/scverse/rapids_singlecell:latest

# Run a script inside it
apptainer exec --nv rsc.sif python my_script.py
```

The `--nv` flag passes through the host GPU.

#### Summary: what to try and in what order

| Package type | First try | Second try | Last resort |
|---|---|---|---|
| Pure Python | `pip install` / `uv pip install` | mamba | conda |
| RAPIDS/GPU | `pip --extra-index-url pypi.nvidia.com` | apptainer pull rsc.sif | mamba (slow) |
| Bioinformatics tool | `conda install -c bioconda` via mamba | Spack | build from source |
| System library | Spack | conda | contact sysadmin |

---

### 18.7 Multi-Agent Orchestration

Complex HPC workflows benefit from splitting responsibilities across multiple agents (or Claude Code sessions). Each agent handles an isolated concern with a well-defined input/output contract.

#### Agent decomposition pattern

```
Orchestrator agent
  ├── [Agent A] Job submission — reads manifest, writes sbatch, submits array
  ├── [Agent B] Status monitor — polls squeue/sacct, reports completions and failures
  ├── [Agent C] Per-phase analysis — runs preprocessing/integration for one phase
  └── [Agent D] Figure/report generation — reads checkpoints, produces outputs
```

#### Rules for multi-agent HPC work

1. **Shared state via files, not memory.** Agents communicate through checkpoint files (`adata.p0.h5ad`, `results/integration_method.txt`, sentinel files `.phase1.done`). No shared in-memory state.
2. **One agent per phase.** Each agent handles one pipeline phase (or one cluster of related jobs). Hand off via checkpoint file existence.
3. **Idempotent operations.** Every script and sbatch must be safe to re-run: check for existing output before recomputing. Use Snakemake rules — they enforce input/output contracts automatically.
4. **Explicit status files.** Write `results/.phase1.done` (via `touch`) on success; write `results/.phase1.failed` with error message on failure. Downstream agents gate on these sentinels.
5. **Parallel branches are independent.** Phase 3.5 (Demographics) and 3.5b (Gene scoring) can run as parallel agents from the same Phase 3 checkpoint. Neither blocks the other.
6. **Job submission agent is separate from analysis agent.** One agent writes and submits sbatch files; a separate agent (or the same agent in a later turn) monitors and acts on results. Do not poll in a tight loop — poll with `gh run list` or `sacct` at 60–120 second intervals.

#### Example: launching parallel per-sample ingestion in Claude Code

```python
# In a Task tool call, launch one subagent per batch
# Each subagent submits one sbatch array job for its batch
task_ids = []
for batch in ["batch1", "batch2"]:
    task_ids.append(
        launch_task(
            agent="general-purpose",
            prompt=f"Submit SLURM array job for {batch} samples. "
                   f"Read metadata/phase0/{batch}_samples.tsv. "
                   f"Write sbatch to scripts/ingest_{batch}.sbatch. Submit with sbatch."
        )
    )
# Wait for both, then submit downstream concat job
```

#### Recommended split for sc_tools pipeline on brb

| Stage | Agent role | Parallelism |
|-------|-----------|-------------|
| Phase 0a (SpaceRanger) | Job submission agent | Array: 1 task/sample |
| Phase 0b (load) | Job submission agent | Array: 1 task/sample |
| Phase 1 (QC + concat) | Single analysis agent | Sequential (depends on all 0b) |
| Phase 3 (integration benchmark) | Benchmark agent | Array: 1 task/method (9 methods) |
| Phase 3.5 + 3.5b | Two parallel agents | Independent branches |
| Phase 5 (figures) | Figure agent per figure set | Parallel by figure group |

#### Integration benchmark as parallel array (9 methods)

Rather than running 9 integration methods sequentially, submit as a SLURM array. Each task runs one method from `sc_tools.bm`:
```bash
#SBATCH --array=0-8
METHODS=(scvi harmony cytovi scvi+harmony raw pca bbknn scanorama desc)
METHOD="${METHODS[$SLURM_ARRAY_TASK_ID]}"
python -c "
import sc_tools.bm as bm, anndata as ad
adata = ad.read_h5ad('results/adata.normalized.p3.h5ad')
bm.run_single_method(adata, method='$METHOD', output_dir='results/tmp/integration_test')
"
```

---

### 18.8 Common Pitfalls

| Pitfall | Fix |
|---------|-----|
| `conda install` hangs on Lustre | Use `pip` / `uv pip` instead; for RAPIDS use NVIDIA PyPI index; for conda-only packages use mamba with `use_locks: false` |
| OOM during scVI on large datasets | Use backed AnnData; reduce `batch_size`; increase `--mem` |
| Lustre lock contention (many jobs writing same dir) | Stagger array job start with `%K` throttle; write to per-task subdirs |
| Job array silently skips tasks | Check `sacct -j JOBID` — tasks may be in `PENDING` due to resource limits |
| `conda activate` fails in sbatch | Always include `eval "$(conda shell.bash hook 2>/dev/null)"` before `conda activate` |
| GPU idle but VRAM full | Inspect with `nvidia-smi`; prior job may have leaked a process — kill with `fuser -k /dev/nvidia0` |
| Slow h5ad save to Lustre | Save to `/tmp` then `rsync` to Lustre at end of job |
| Array task IDs off-by-one | Use `tail -n +2` when building sample arrays from TSV to skip header |
| `rapids-singlecell` not available | Install with `pip install rapids-singlecell` in sc_tools conda env on GPU node; not in base env |

---

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
