# Architecture Research

**Domain:** Spatial/UMAP plot integration into sc_tools HTML report pipeline (v2.0 milestone)
**Researched:** 2026-03-27
**Confidence:** HIGH (all findings from direct codebase inspection)

## Standard Architecture

### System Overview

The v2.0 milestone adds two integration points to the existing layered architecture:
(1) plot generation feeding into existing report templates, and
(2) a new `sct concat` CLI command registered as a pipeline phase.

```
┌─────────────────────────────────────────────────────────────────────┐
│                     CLI Surface (sct)                                │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌─────────────────────┐ │
│  │ sct qc   │  │ sct pp   │  │ sct bm   │  │ sct concat  [NEW]   │ │
│  │ sct      │  │          │  │          │  │ sct ingest  [NEW]   │ │
│  │ report   │  │          │  │          │  │                     │ │
│  └────┬─────┘  └────┬─────┘  └────┬─────┘  └──────────┬──────────┘ │
│       │              │             │                   │            │
├───────┴──────────────┴─────────────┴───────────────────┴────────────┤
│                   Command Layer (thin adapters)                      │
│  ┌───────────────────────────────────────────────────────────────┐  │
│  │  CLIResult {status, command, data, artifacts, provenance}     │  │
│  └───────────────────────────────────────────────────────────────┘  │
├─────────────────────────────────────────────────────────────────────┤
│         Report Pipeline  (sc_tools/qc/report.py)  [MODIFIED]        │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │  generate_pre_filter_report()   — add spatial QC 4-panel     │   │
│  │  generate_post_filter_report()  — spatial plots exist [DONE] │   │
│  │  generate_post_integration_report() — UMAP grid exists [DONE]│   │
│  │  generate_post_celltyping_report() — add spatial celltype     │   │
│  └────────────────────────────┬─────────────────────────────────┘   │
│                               │ calls                               │
│  ┌────────────────────────────▼─────────────────────────────────┐   │
│  │  Plot Library  (sc_tools/pl/)  [EXTENDED]                    │   │
│  │  ┌───────────────┐  ┌──────────────────┐  ┌───────────────┐  │   │
│  │  │ pl/spatial.py │  │ pl/qc_plots.py   │  │ qc/plots.py   │  │   │
│  │  │ plot_spatial_ │  │ qc_umap_grid()   │  │ qc_2x2_grid() │  │   │
│  │  │ continuous()  │  │ qc_celltype_     │  │ qc_spatial_   │  │   │
│  │  │ [NEW wrapper] │  │ abundance() etc. │  │ multipage()   │  │   │
│  │  └───────────────┘  └──────────────────┘  └───────────────┘  │   │
│  └──────────────────────────────────────────────────────────────┘   │
│                               │ fig_to_base64()                     │
│  ┌────────────────────────────▼─────────────────────────────────┐   │
│  │  Jinja2 Templates  (sc_tools/assets/*.html)  [MODIFIED]      │   │
│  │  pre_filter_qc_template.html  post_celltyping_qc_template.html│  │
│  └──────────────────────────────────────────────────────────────┘   │
├─────────────────────────────────────────────────────────────────────┤
│                   sc_tools Library (unchanged except additions)      │
│  ┌───────────┐  ┌──────┐  ┌──────┐  ┌──────┐  ┌───────────────┐   │
│  │  ingest/  │  │  qc/ │  │  pp/ │  │  pl/ │  │  pipeline.py  │   │
│  │ loaders.py│  │      │  │      │  │      │  │  [MODIFIED]   │   │
│  │ concat_   │  │      │  │      │  │      │  │  add concat   │   │
│  │ samples() │  │      │  │      │  │      │  │  phase        │   │
│  └───────────┘  └──────┘  └──────┘  └──────┘  └───────────────┘   │
├─────────────────────────────────────────────────────────────────────┤
│                   Data Layer                                         │
│  ┌────────────────────┐  ┌─────────────────────┐                    │
│  │ AnnData checkpoints│  │ Self-contained HTML  │                    │
│  │ (h5ad files)       │  │ reports with base64  │                    │
│  │ obsm['spatial']    │  │ PNG plots embedded   │                    │
│  │ uns['spatial']     │  │                      │                    │
│  └────────────────────┘  └─────────────────────┘                    │
└─────────────────────────────────────────────────────────────────────┘
```

### Component Responsibilities

| Component | Responsibility | Status |
|-----------|----------------|--------|
| `sc_tools/cli/qc.py` | `sct report generate` command — thin adapter calling `report.py` functions | EXISTS — needs no changes for plots |
| `sc_tools/cli/ingest.py` | `sct concat` command — new thin adapter wrapping `ingest.concat_samples()` | NEW FILE |
| `sc_tools/qc/report.py` | Report orchestration — builds plot dict, renders Jinja2 template | MODIFIED — add spatial 4-panel to pre_filter and celltype spatial to post_celltyping |
| `sc_tools/pl/spatial.py` | Low-level spatial plot primitives (`plot_spatial_continuous`, `plot_spatial_categorical`) | EXISTS — add multi-metric 4-panel wrapper |
| `sc_tools/pl/qc_plots.py` | Post-integration UMAP grids, cluster distribution, celltype abundance | EXISTS — used by post_integration/post_celltyping already |
| `sc_tools/qc/plots.py` | Pre/post-filter QC plots (`qc_2x2_grid`, `qc_spatial_multipage`) | EXISTS — re-exported via `sc_tools.pl` |
| `sc_tools/qc/report_utils.py` | `fig_to_base64()`, `render_template()`, `auto_detect_embeddings()` | EXISTS — used as-is |
| `sc_tools/assets/*.html` | Jinja2 templates for self-contained HTML output | MODIFIED — add spatial plot sections |
| `sc_tools/pipeline.py` | Phase DAG with `PhaseSpec` objects including checkpoint paths | MODIFIED — register `concat` phase |

## Recommended Project Structure

```
sc_tools/
├── cli/
│   ├── ingest.py               # NEW: sct concat command (wraps ingest.concat_samples)
│   └── qc.py                   # EXISTS: sct report generate (no changes needed)
├── qc/
│   └── report.py               # MODIFIED: add spatial 4-panel + celltype spatial sections
├── pl/
│   ├── spatial.py              # MODIFIED: add qc_spatial_4panel() wrapper function
│   └── qc_plots.py             # EXISTS: unchanged (UMAP grid already works)
├── assets/
│   ├── pre_filter_qc_template.html      # MODIFIED: add spatial-plots section
│   └── post_celltyping_qc_template.html # MODIFIED: add celltype-spatial section
└── pipeline.py                 # MODIFIED: add concat PhaseSpec
```

### Structure Rationale

- **`cli/ingest.py` is a new file** (not `cli/qc.py`) because concat is an ingestion-layer operation — it merges per-sample AnnDatas before QC runs. The naming mirrors the `sc_tools/ingest/` backend it wraps.
- **Plot code stays in `sc_tools/pl/`** — all new plot functions belong in `pl/spatial.py` (for spatial overlays) or `pl/qc_plots.py` (for UMAP/abundance). The report.py orchestrator calls these functions and converts figures to base64 via `fig_to_base64()`. Do not put plot-generation code directly in report.py.
- **Templates only receive pre-rendered base64 strings** — the Jinja2 templates remain logic-free. All conditionals about whether a section exists happen in report.py before rendering, not inside the template.
- **No new `sc_tools/pl/report.py` module** — the existing `sc_tools/qc/report.py` is the report orchestrator. Adding a second report module would create confusion. Extend the existing one.

## Architectural Patterns

### Pattern 1: Base64 PNG Embedding (existing, authoritative)

**What:** Matplotlib figures are rendered to PNG in memory via `fig_to_base64(fig, dpi=150)` from `sc_tools/qc/report_utils.py`, which returns a base64-encoded string. The string is placed directly into the `plots` dict passed to Jinja2 as `<img src="data:image/png;base64,{...}">`. No external files are written for plots.

**Why base64 PNG over alternatives:**
- Plotly JSON would require the Plotly JS bundle in every report. Current reports are self-contained (no CDN dependency). Plotly is useful for interactive charts but adds ~3MB per report.
- Inline SVG works but matplotlib SVG output for spatial scatter plots with 50K+ points is very large (one `<circle>` element per point). PNG with `rasterized=True` is necessary for spatial data.
- The existing pattern is consistent across all four report types. Maintaining consistency avoids a dual-rendering path.

**Confirmed:** `generate_post_filter_report()` already generates spatial plots this way (lines 302-334 in report.py). Follow the same pattern for pre_filter and post_celltyping.

### Pattern 2: Per-Sample Spatial Loop

**What:** For any spatial plot section, iterate over unique sample IDs, subset the AnnData, and generate one figure per sample. Append `{"sample": lib_id, "img": base64_string}` to a list. Pass the list to Jinja2 as a template variable. The template iterates and renders one plot per sample row.

**Guard conditions before attempting the loop:**
```python
_has_spatial = (
    "spatial" in adata.obsm
    and sample_col in adata.obs.columns
    and adata.uns.get("spatial")  # has per-sample spatial metadata
)
```

**Subset pattern (from existing post_filter code):**
```python
for lib in sorted(adata.obs[sample_col].dropna().unique()):
    if str(lib) not in adata.uns.get("spatial", {}):
        continue
    sub = adata[adata.obs[sample_col] == lib].copy()
    try:
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        sc.pl.spatial(sub, color=color_key, library_id=str(lib),
                      show=False, ax=ax, frameon=False)
        plots_list.append({"sample": str(lib), "img": fig_to_base64(fig)})
        plt.close(fig)
    except Exception:
        logger.debug("Spatial plot failed for %s", lib, exc_info=True)
        plt.close("all")
```

**squidpy.pl.spatial_scatter interaction:** squidpy's `spatial_scatter` is an alternative to `sc.pl.spatial`. Both read the same AnnData structure: `adata.obsm["spatial"]` for coordinates and `adata.uns["spatial"][library_id]` for the H&E image. The existing codebase uses `sc.pl.spatial` (scanpy), not `squidpy.pl.spatial_scatter`. Use `sc.pl.spatial` for consistency — it already handles the `library_id` argument and per-library subsetting. squidpy's plotting API is higher-level but adds a dependency. Only adopt `squidpy.pl` if a specific feature (e.g., multi-panel spatial layout) is not achievable with scanpy.

### Pattern 3: New Spatial 4-Panel Function in pl/spatial.py

**What:** Add `qc_spatial_4panel(adata_sub, library_id, metric_keys)` to `sc_tools/pl/spatial.py`. This function produces a 1x4 (or 2x2) grid of spatial scatter plots for a single sample, one panel per QC metric. `report.py` calls this function per sample inside the per-sample loop.

**Signature design:**
```python
def qc_spatial_4panel(
    adata: AnnData,
    library_id: str,
    metric_keys: list[str] | None = None,  # defaults: log1p_total_counts, log1p_n_genes_by_counts, pct_counts_mt, qc_pass
    figsize: tuple[float, float] = (20, 5),
    dpi: int = 120,
) -> plt.Figure:
```

**Metric key defaults for QC report:** `["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "qc_pass"]`. These need to exist in `adata.obs` before the function is called. `report.py` must ensure they are computed (see data flow section).

**obs column validation:** The function should silently skip panels for missing columns rather than raising — spatial data may not have MT reads (some platforms). Return a partial figure with blank panels for missing metrics.

### Pattern 4: cli_handler Decorator for sct concat

**What:** `sct concat` must follow the exact same `@cli_handler` + `CLIResult` pattern as all other CLI commands. The decorator handles provenance, exit codes, and JSON/human output.

```python
# sc_tools/cli/ingest.py
from sc_tools.cli import _check_deps, cli_handler
from sc_tools.models.result import CLIResult, Provenance, Status
import typer

ingest_app = typer.Typer(help="Data ingestion commands")

@ingest_app.command("concat")
@cli_handler
def concat(
    input_files: list[str] = typer.Argument(..., help="Paths to per-sample h5ad files"),
    output: str = typer.Option(..., "--output", "-o", help="Output concatenated h5ad path"),
    sample_col: str = typer.Option("sample", "--sample-col", help="Sample identifier column"),
    modality: str = typer.Option("visium", "--modality", "-m", help="Data modality"),
    calculate_qc: bool = typer.Option(True, "--qc/--no-qc", help="Run QC metrics after concat"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate inputs without executing"),
) -> None:
    """Concatenate per-sample AnnData files into a single multi-sample h5ad (CMD-concat)."""
    _check_deps(["scanpy", "anndata"])
    from pathlib import Path
    from sc_tools.ingest import concat_samples
    ...
    return CLIResult(
        status=Status.success,
        command="concat",
        data={"n_samples": n_samples, "n_cells": n_cells, "output": str(output_path)},
        artifacts=[str(output_path)],
        provenance=Provenance(command="concat"),
        message=f"Concatenated {n_samples} samples: {n_cells} cells -> {output_path}",
    )
```

**Registration in `__main__.py` or `cli/__init__.py`:** The new `ingest_app` must be added to the root `app` the same way `qc_app` and `report_app` are registered. Check `sc_tools/__main__.py` for the `app.add_typer()` calls.

### Pattern 5: Pipeline Phase Registration for concat

**What:** Add a `concat` PhaseSpec to `STANDARD_PHASES` in `sc_tools/pipeline.py`. This makes `sct status` report concat completion and allows DAG-aware dependency checks.

**Correct placement in DAG:** concat sits between `ingest_load` (per-sample h5ads) and `qc_filter` (multi-sample filtered h5ad). It replaces the implicit concatenation that previously happened inside `qc_filter`.

```python
"concat": PhaseSpec(
    label="Multi-Sample Concatenation",
    depends_on=[_dp("ingest_load")],
    branch="ingestion",
    checkpoint="results/adata.concatenated.h5ad",
    phase_group=_DP,
    required_obs=["sample", "library_id", "raw_data_dir"],
    required_obsm=["spatial"],
    x_format="raw counts, concatenated",
    old_code="p0c",
),
```

**Update `qc_filter` dependency:** Change `qc_filter.depends_on` from `[_dp("ingest_load")]` to `[_dp("concat")]` only if concat becomes mandatory. If concat is optional (projects may still ingest pre-concatenated data), keep it `optional=True` and leave `qc_filter` depending on `ingest_load`.

**Recommendation:** Make concat `optional=True`. Some projects ingest a single-sample or pre-merged h5ad directly. The phase DAG should not force concat for those cases.

## Data Flow

### QC Report Plot Data Flow

```
Checkpoint: results/adata.filtered.h5ad
    |
    | sct report generate post_filter --adata ...
    v
report.py: generate_post_filter_report()
    |
    | 1. Ensure obs columns exist:
    |    - log1p_total_counts (np.log1p of total_counts)
    |    - n_genes_by_counts (from sc.pp.calculate_qc_metrics)
    |    - pct_counts_mt (from sc.pp.calculate_qc_metrics with mt vars)
    |    - qc_pass (from classify_samples())
    |
    | 2. Check _has_spatial:
    |    - "spatial" in adata.obsm
    |    - sample_col in adata.obs.columns
    |    - adata.uns.get("spatial") has per-library entries
    |
    | 3. For each sample: call pl.spatial.qc_spatial_4panel()
    |    -> returns plt.Figure
    | 4. fig_to_base64(fig) -> base64 PNG string
    | 5. Append {"sample": lib, "img": b64} to spatial_qc_plots list
    |
    v
Jinja2 context: {"spatial_qc_plots": [...], "plots": {...}}
    |
    v
render_template("post_filter_qc_template.html", context)
    |
    v
Self-contained HTML with embedded base64 PNGs
```

### Integration Report UMAP Data Flow

```
Checkpoint: results/adata.normalized.h5ad
    Required obs: leiden (cluster labels)
    Required obsm: X_umap, X_scvi (or other integration embeddings)
    |
    v
generate_post_integration_report()
    |
    | Already implemented: qc_umap_grid(adata, color_keys=[sample_col, batch_key, "leiden"])
    | Already implemented: qc_embedding_umap_grid(adata, embedding_keys, color_key=sample_col)
    | Already implemented: qc_cluster_distribution(adata, cluster_key, sample_col)
    |
    | New for v2.0: explicit bio obs columns in color_keys
    |   - Add patient/subject_id if present
    |   - Add any column listed in PhaseSpec.required_obs for preprocess phase
    v
No new code needed — extend color_keys list in existing function call
```

### Celltype Report Spatial Data Flow

```
Checkpoint: results/adata.celltyped.h5ad
    Required obs: celltype, celltype_broad
    Required obsm: spatial (for spatial plots), X_umap (for UMAP)
    |
    v
generate_post_celltyping_report()
    |
    | NEW: Add spatial celltype section
    | For each sample:
    |   sub = adata[adata.obs[sample_col] == lib].copy()
    |   sc.pl.spatial(sub, color="celltype", library_id=lib, ...)
    |   fig_to_base64(fig) -> base64 PNG
    |
    v
Jinja2 context: {"celltype_spatial_plots": [...]}
    |
    v
post_celltyping_qc_template.html with new celltype-spatial section
```

### sct concat Data Flow

```
Per-sample h5ad files: data/s1/adata.ingested.h5ad, data/s2/adata.ingested.h5ad, ...
    |
    | sct concat s1.h5ad s2.h5ad --output results/adata.concatenated.h5ad
    v
cli/ingest.py: concat()
    |
    | Validate: all input files exist
    | Validate: same modality (check adata.uns["spatial"] keys)
    | Load each h5ad (sequential, not parallel -- memory constraint)
    | Call sc_tools.ingest.concat_samples(adatas, sample_col=sample_col, calculate_qc=True)
    | Write output h5ad
    | Write .provenance.json sidecar
    |
    v
CLIResult: {status, data: {n_samples, n_cells, output}, artifacts: [output_path]}
```

### Key Data Flows

1. **Pre-filter spatial plots:** `adata.filtered.h5ad` (post-QC, pre-concat or pre-integration) -> per-sample spatial scatter of log1p_counts, log1p_genes, pct_mt, pass/fail -> 4-panel base64 PNGs -> pre_filter template.

2. **Integration UMAPs:** `adata.normalized.h5ad` with `X_umap` and integration embeddings in `obsm` -> `qc_umap_grid()` colored by leiden, batch, sample, patient -> base64 PNGs -> post_integration template. **Already implemented.**

3. **Celltype spatial:** `adata.celltyped.h5ad` with `celltype` in `obs` and `spatial` in `obsm` -> per-sample spatial scatter colored by cell type -> base64 PNGs -> post_celltyping template.

4. **Concat ingestion:** N per-sample h5ads -> `concat_samples()` -> single multi-sample h5ad with QC metrics -> registered in pipeline DAG as `concat` phase.

## Obs Column / obsm Key Validation

Each report type has a minimum data contract. The CLI command must validate these exist before calling the report generator. Failure should return `CLIResult(status=Status.error, ...)` with a clear message, not a Python traceback.

### Pre-filter / Post-filter Report

| Key | Location | Required For | How to Compute if Missing |
|-----|----------|--------------|---------------------------|
| `total_counts` | `obs` | Spatial color axis | `sc.pp.calculate_qc_metrics()` |
| `log1p_total_counts` | `obs` | Spatial 4-panel | `np.log1p(adata.obs["total_counts"])` inline |
| `n_genes_by_counts` | `obs` | Spatial 4-panel | `sc.pp.calculate_qc_metrics()` |
| `log1p_n_genes_by_counts` | `obs` | Spatial 4-panel | `np.log1p(adata.obs["n_genes_by_counts"])` inline |
| `pct_counts_mt` | `obs` | Spatial 4-panel | `sc.pp.calculate_qc_metrics(qc_vars=[...])` — skip panel if absent |
| `qc_pass` | `obs` | Spatial 4-panel pass/fail | `classify_samples()` output |
| `spatial` | `obsm` | All spatial plots | Must exist; skip spatial section if absent |
| `spatial` | `uns` per library | All spatial plots | Must have per-sample entries; skip sample if absent |
| `library_id` or `sample_col` | `obs` | Per-sample loop | Must exist; raise SCToolsDataError if absent |

### Post-integration Report

| Key | Location | Required For | Notes |
|-----|----------|--------------|-------|
| `X_umap` | `obsm` | UMAP grid | Skip UMAP section if absent |
| `leiden` | `obs` | UMAP color, cluster distribution | Required for integration report |
| `batch_key` (auto-detected) | `obs` | UMAP color, integration metrics | Falls back to `library_id` |
| Any `X_scvi`, `X_harmony`, etc. | `obsm` | Per-embedding UMAP, benchmark | Auto-detected by `auto_detect_embeddings()` |

### Post-celltyping Report

| Key | Location | Required For | Notes |
|-----|----------|--------------|-------|
| `celltype` | `obs` | All sections | Hard required; raise ValueError if absent |
| `celltype_broad` | `obs` | Abundance chart | Optional |
| `X_umap` | `obsm` | UMAP colored by celltype | Skip if absent |
| `spatial` | `obsm` | Spatial celltype plots | Skip spatial section if absent |

### sct concat Pre-execution Validation

```
For each input file:
  - File exists and is readable
  - adata.obs contains sample_col
  - adata.obsm["spatial"] exists (warn if absent for spatial modalities)
  - adata.X is raw counts (check adata.uns["sc_tools"]["x_format"] if present)

Cross-file:
  - All files have consistent var_names (same genes)
  - Warn if uns["spatial"] keys overlap (duplicate library IDs)
```

## Anti-Patterns

### Anti-Pattern 1: Plot Code in report.py

**What people do:** Write matplotlib figure generation directly inside `generate_pre_filter_report()` or similar functions.
**Why it's wrong:** report.py becomes a 1000+ line file mixing orchestration with figure styling. Plot functions cannot be reused outside reports (e.g., in notebooks, MCP tools). The existing separation (`qc/plots.py`, `pl/spatial.py`, `pl/qc_plots.py`) is deliberate.
**Do this instead:** Add new plot functions to the appropriate `pl/` submodule. Call them from report.py and convert with `fig_to_base64()`. The report orchestrator should contain no matplotlib calls directly.

### Anti-Pattern 2: Plotly for Spatial Scatter Plots

**What people do:** Use Plotly for interactive spatial plots to enable zoom/pan.
**Why it's wrong:** Spatial datasets have 50K-2.5M spots. Plotly renders each point as an SVG element; the HTML becomes multi-megabyte and unusable. Rasterized PNG at 120-150 DPI is the correct choice for large spatial datasets.
**Do this instead:** Use `sc.pl.spatial()` with `rasterized=True` (scanpy sets this automatically for spatial plots). Use `fig_to_base64(fig, dpi=120)`. Reserve Plotly for small data: metrics tables (already used in benchmark reports), sample-level summary charts.

### Anti-Pattern 3: Subsetting adata Before Spatial Plot Without uns Copy

**What people do:** `sub = adata[mask]` (view, not copy) and pass to `sc.pl.spatial()`.
**Why it's wrong:** `adata.uns["spatial"]` contains per-library image data. A plain slice view may not have the correct `uns` structure for the selected library. `sc.pl.spatial()` accesses `uns["spatial"][library_id]` and will fail or show wrong image.
**Do this instead:** `sub = adata[adata.obs[sample_col] == lib].copy()` — always `.copy()` when subsetting for spatial plots. This is the pattern already used in `generate_post_filter_report()` at line 313 in report.py.

### Anti-Pattern 4: Loading Full adata for concat When Only obs/var Needed for Validation

**What people do:** Load all N h5ad files fully into memory to validate compatibility, then concatenate.
**Why it's wrong:** For 10-sample Visium HD projects, validation would load 10 × 25G = 250G before a single byte is written.
**Do this instead:** For the validation pass in `sct concat`, use h5py to read only `obs` and `var_names` from each file:
```python
import h5py
with h5py.File(path, "r") as f:
    var_names = f["var"]["_index"][:]
    n_obs = f["obs"]["_index"].shape[0]
```
Load full AnnDatas only for the actual concat operation, sequentially (not all at once), streaming into `anndata.concat()` if memory is tight.

### Anti-Pattern 5: Registering concat as a Mandatory Upstream Dependency

**What people do:** Add concat as a hard dependency for qc_filter, forcing all projects through the concat phase.
**Why it's wrong:** Some projects ingest a single pre-merged h5ad. Others already have a concatenated checkpoint from a previous run. Making concat mandatory breaks backward compat with existing projects.
**Do this instead:** Mark the concat phase `optional=True` in PhaseSpec. qc_filter continues to depend on `ingest_load` (the lowest-common-denominator checkpoint). The pipeline DAG already supports optional phases.

## Integration Points

### New vs Modified Files

| File | Change Type | What Changes |
|------|-------------|--------------|
| `sc_tools/cli/ingest.py` | NEW | `sct concat` command, `ingest_app` Typer group |
| `sc_tools/cli/__init__.py` | MODIFIED | Register `ingest_app` via `app.add_typer(ingest_app, name="ingest")` |
| `sc_tools/pipeline.py` | MODIFIED | Add `concat` PhaseSpec to `STANDARD_PHASES` |
| `sc_tools/pl/spatial.py` | MODIFIED | Add `qc_spatial_4panel()` function |
| `sc_tools/qc/report.py` | MODIFIED | Add spatial 4-panel loop to `generate_pre_filter_report()`; add celltype spatial loop to `generate_post_celltyping_report()` |
| `sc_tools/assets/pre_filter_qc_template.html` | MODIFIED | Add spatial-qc-plots section |
| `sc_tools/assets/post_celltyping_qc_template.html` | MODIFIED | Add celltype-spatial section |
| `sc_tools/pl/__init__.py` | MODIFIED | Export `qc_spatial_4panel` |

### Unchanged Files (confirmed)

| File | Why Unchanged |
|------|---------------|
| `sc_tools/ingest/loaders.py` | `concat_samples()` already exists and is correct |
| `sc_tools/cli/qc.py` | `sct report generate` already calls the right functions; report.py changes propagate automatically |
| `sc_tools/qc/report_utils.py` | `fig_to_base64()`, `render_template()`, `auto_detect_embeddings()` all exist and work |
| `sc_tools/pl/qc_plots.py` | `qc_umap_grid()`, `qc_cluster_distribution()`, `qc_celltype_abundance()` all exist and are already used by post_integration/post_celltyping reports |
| `sc_tools/qc/spatial.py` | Spatial autocorrelation (SVG computation); not relevant to plot generation |

### Internal Boundaries

| Boundary | Communication | Notes |
|----------|---------------|-------|
| `report.py` -> `pl/spatial.py` | Direct function call, returns `plt.Figure` | New `qc_spatial_4panel()` follows same signature pattern as existing `plot_spatial_continuous()` |
| `report.py` -> `pl/qc_plots.py` | Direct function call, returns `plt.Figure` | Already established; no changes |
| `report.py` -> `report_utils.fig_to_base64` | Function call, returns `str` | Existing utility; use as-is |
| `cli/ingest.py` -> `ingest/loaders.concat_samples` | Direct import + call | Existing function; CLI is a thin wrapper |
| `cli/ingest.py` -> `pipeline.py` | Not at call time — pipeline registration is static at import | PhaseSpec added to STANDARD_PHASES; `sct status` reads it |
| `report.py` -> Jinja2 templates | `render_template(template_name, context_dict)` | context keys are the only contract; template changes must match new keys |

## Build Order

Build in this order to respect checkpoint dependencies and minimize rework:

```
Phase A: sct concat + pipeline registration (no report dependencies)
    1. sc_tools/pipeline.py — add concat PhaseSpec (optional=True)
    2. sc_tools/cli/ingest.py — sct concat command
    3. sc_tools/cli/__init__.py — register ingest_app
    VALIDATES: sct concat runs, produces adata.concatenated.h5ad
               sct status shows concat phase

Phase B: Spatial 4-panel plot function (no template changes yet)
    4. sc_tools/pl/spatial.py — add qc_spatial_4panel()
    5. sc_tools/pl/__init__.py — export qc_spatial_4panel
    VALIDATES: function importable, produces correct 4-panel figure on test data

Phase C: Pre-filter report spatial section
    6. sc_tools/qc/report.py — add spatial_qc_plots loop to generate_pre_filter_report()
    7. sc_tools/assets/pre_filter_qc_template.html — add spatial-qc section
    VALIDATES: sct report generate pre_filter produces HTML with spatial panels
    CHECKPOINT: results/adata.filtered.h5ad with obsm["spatial"]

Phase D: Post-celltyping report spatial section
    8. sc_tools/qc/report.py — add celltype_spatial_plots loop to generate_post_celltyping_report()
    9. sc_tools/assets/post_celltyping_qc_template.html — add celltype-spatial section
    VALIDATES: sct report generate post_celltyping produces HTML with celltype spatial panels
    CHECKPOINT: results/adata.celltyped.h5ad with obs["celltype"] and obsm["spatial"]

Phase E: Integration report bio obs column extension (lowest risk, existing infrastructure)
    10. sc_tools/qc/report.py — extend color_keys in generate_post_integration_report()
        to include patient/subject_id obs columns
    VALIDATES: post_integration report shows patient UMAP panel
    CHECKPOINT: results/adata.normalized.h5ad with obs["leiden"] and obsm["X_umap"]
```

**Rationale for this order:**
- Phase A (concat) has no dependencies on the plot pipeline. It can be built, tested, and merged independently before any report work starts.
- Phase B (4-panel function) must precede Phase C (template changes) because the template assumes the function exists.
- Phase C (pre-filter) before Phase D (celltyping) because pre-filter checkpoints are earlier in the pipeline and easier to test with minimal data.
- Phase E (integration report extension) is last because it is the lowest risk change — extending an existing working section — and does not require new plot functions.

## squidpy.pl.spatial_scatter Interaction

squidpy's `spatial_scatter` and scanpy's `sc.pl.spatial` both read from the same AnnData structure:
- `adata.obsm["spatial"]` — N x 2 coordinate array
- `adata.uns["spatial"][library_id]["images"]["hires"]` — H&E image array
- `adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]` — coordinate scale factor

The existing codebase uses `sc.pl.spatial` exclusively for report generation (confirmed in `generate_post_filter_report()` at lines 320-326). `sc_tools/qc/spatial.py` uses squidpy for graph-based analysis (`sq.gr.spatial_neighbors`, `sq.gr.spatial_autocorr`), not for plotting.

**Decision: use `sc.pl.spatial` for all report spatial plots.** Reasons:
1. Consistency with existing code (no dual-path maintenance)
2. `sc.pl.spatial` handles `library_id` subsetting and H&E overlay natively
3. squidpy is a soft dependency (used lazily in `sc_tools/qc/spatial.py`); making it required for reports would break installs without squidpy
4. squidpy's `spatial_scatter` is more powerful for multi-sample plots, but the per-sample subsetting pattern already achieves the same result

If squidpy's `spatial_scatter` is needed in the future (e.g., side-by-side multi-sample layout), add it as an optional code path behind a try/import guard in `pl/spatial.py`, consistent with how squidpy is used elsewhere in the codebase.

## Sources

All findings are HIGH confidence from direct codebase inspection (2026-03-27):

- `sc_tools/qc/report.py` — existing report generators and base64 embedding pattern
- `sc_tools/pl/spatial.py` — existing spatial plot primitives
- `sc_tools/pl/qc_plots.py` — existing UMAP grid and abundance functions
- `sc_tools/pl/__init__.py` — existing pl module exports
- `sc_tools/qc/plots.py` — existing QC plot functions (re-exported via pl)
- `sc_tools/pipeline.py` — PhaseSpec DAG structure and registration API
- `sc_tools/cli/__init__.py` — cli_handler decorator, CLIResult, app structure
- `sc_tools/cli/qc.py` — existing report generate command pattern
- `sc_tools/ingest/__init__.py` — concat_samples function availability
- `sc_tools/ingest/loaders.py` — concat_samples signature (adatas, sample_col, calculate_qc)
- `.planning/PROJECT.md` — v2.0 milestone target features and constraints

---
*Architecture research for: v2.0 spatial/UMAP plot integration into sc_tools report pipeline*
*Researched: 2026-03-27*
