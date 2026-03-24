# Phase 6: Scientific Gaps - Research

**Researched:** 2026-03-24
**Domain:** Pseudobulk DE, marker validation, subject metadata, panel-aware cell typing
**Confidence:** HIGH

## Summary

Phase 6 adds four scientific capabilities to sc_tools: (1) pseudobulk differential expression via PyDESeq2, (2) marker validation reporting integrated into the existing post-celltyping report, (3) subject_id metadata enforcement at ingestion, and (4) panel-aware cell typing dispatch that restricts methods for targeted panels. All four requirements build on well-established patterns already in the codebase -- the CLI handler/CLIResult envelope, Jinja2 HTML report system, cell typing backend registry, and lazy dependency checking.

The primary technical risks are low: PyDESeq2 has a clean pandas-based API, the dotplot/marker validation infrastructure already exists in the post-celltyping report, and the cell typing dispatch function is a simple 77-line file with a clear extension point. The pseudobulk aggregation (sum counts by subject_id + celltype) is straightforward numpy/pandas without needing decoupler.

**Primary recommendation:** Implement in dependency order: SCI-03 (subject_id) first (needed by SCI-01), then SCI-01 (pseudobulk DE), then SCI-02 (marker validation report), then SCI-04 (panel dispatch) -- each is largely independent after SCI-03.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- D-01: Auto-infer design formula with override. Auto-build `~ condition + batch` from obs metadata. Detect batch covariates from `library_id`/`batch` columns. User can override with `--formula` flag for custom covariates (e.g., `~ condition + sex + age`).
- D-02: Aggregation by `subject_id + celltype`. Minimum thresholds: >=3 subjects per condition group, >=10 cells per subject+celltype combination for aggregation. Combinations below threshold are excluded with a warning.
- D-03: Per-celltype CSV output + summary JSON. One CSV per cell type with columns: `gene`, `log2FC`, `pvalue`, `padj`, `baseMean`. Output directory: `{project_dir}/results/de/`. CLIResult.data has summary stats, CLIResult.artifacts lists all CSV paths.
- D-04: PyDESeq2 as the DE engine. No alternatives.
- D-05: HTML report extending existing report system (Phase 3 templates in `sc_tools/assets/`). Dotplot of top 5 markers per assigned cell type. Flag table for types with low canonical marker expression (mean expression below configurable threshold, default 0.1).
- D-06: Integrated into `sct report post-celltyping`. Summary includes: n_types tested, n_flagged, total cells. Canonical marker source is from existing gene set signatures (`sc_tools/tl/gene_sets.py`) or user-provided marker file.
- D-07: Flagging is informational -- does not fail the command. Agent decides whether to act on flags.
- D-08: Warn + validate, don't block. Single-sample projects: warn if `subject_id` missing but don't fail. Multi-sample projects: require `subject_id`, validate it's distinct from `library_id`.
- D-09: Batch-condition confounding check at registration/QC. When `subject_id` and a condition column exist, check if batch perfectly confounds condition. Warn if so -- don't block.
- D-10: `subject_id` is an obs column convention, not a schema migration. Enforced via validation functions, not DB constraints.
- D-11: When `n_vars < 1000`, auto-restrict to panel-validated methods: `sctype` (with custom markers) and `custom_gates`. Warn if user requests whole-transcriptome models.
- D-12: `--force-method` flag allows override with a warning logged in provenance.
- D-13: Panel detection logged in CLIResult provenance. Decision recorded: `{panel_detected: true, n_vars: 300, restricted_methods: [...]}`.

### Claude's Discretion
- PyDESeq2 wrapper implementation details (anndata to pandas conversion, count matrix extraction)
- Dotplot rendering library (scanpy.pl.dotplot, custom matplotlib, or Plotly for HTML)
- Exact confounding detection algorithm (perfect confounding vs partial)
- How `annotate_celltypes()` dispatch changes for panel mode (decorator, guard clause, or config)
- Whether `sct de run` is a new command group or subcommand under existing group
- Test fixture design for pseudobulk (synthetic multi-subject data)

### Deferred Ideas (OUT OF SCOPE)
- Cell-typing benchmarking (ADV-06 in v2) -- compare sctype, celltypist, scArches systematically
- Trajectory / RNA velocity (ADV-03 in v2)
- Bootstrap uncertainty for DE results
- Volcano plot generation for DE results -- could be a future report enhancement
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| SCI-01 | Pseudobulk DE module (`sc_tools.tl.de`) -- aggregate counts by subject_id + celltype, run PyDESeq2 with batch covariates | PyDESeq2 0.5.2 API documented; DeseqDataSet takes pandas count matrix + metadata DataFrame + Wilkinson formula string; results in `stat_res.results_df` |
| SCI-02 | Marker validation report -- after cell typing, compute top N marker genes per assigned type, generate dotplot, flag types with low canonical marker expression | Existing `generate_post_celltyping_report()` already accepts `marker_genes` dict and renders dotplot via `sc.pl.dotplot`; extend with flag table and summary stats |
| SCI-03 | Subject-level metadata model -- `subject_id` distinct from `library_id`, enforced at ingestion for multi-sample projects, batch-condition confounding validation | Validation function in `sc_tools/qc/` or `sc_tools/ingest/`; no schema migration needed (D-10) |
| SCI-04 | Panel-aware cell typing dispatch -- when `n_vars < 1000`, restrict to panel-validated methods | Guard clause in `annotate_celltypes()` before `get_backend()` call; backend registry already categorizes methods (Tier 1-2 vs Tier 3) |
</phase_requirements>

## Project Constraints (from CLAUDE.md)

- Output paths: active projects at `~/Documents/projects/active/<project>/` -- never repo root
- Phase transitions: call `set_phase_status` on completion, dispatch documentor for Plan.md/Journal.md
- CLI commands use `cli_handler` decorator + CLIResult envelope
- Lazy imports for heavy deps; `_check_deps()` fast-fail pattern
- HTML reports use Jinja2 templates in `sc_tools/assets/`

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| pydeseq2 | 0.5.2 | Pseudobulk differential expression | Only pure-Python DESeq2 implementation; scverse ecosystem; accepts pandas DataFrames directly |
| scanpy | 1.11.5 (installed) | Dotplot rendering, rank_genes_groups | Already used throughout sc_tools; `sc.pl.dotplot` produces publication-quality marker plots |
| anndata | 0.11.4 (installed) | Data container | Already the data format throughout sc_tools |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pandas | (installed) | Count matrix aggregation, DE results tables | Pseudobulk aggregation and PyDESeq2 I/O |
| numpy | (installed) | Sparse-to-dense conversion, thresholding | Count extraction from adata.X |
| scipy.sparse | (installed) | Sparse matrix handling | Efficient sum aggregation of counts |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| PyDESeq2 | decoupler.get_pseudobulk + edgeR | decoupler adds aggregation helper but requires R backend; PyDESeq2 is pure Python |
| Manual pseudobulk aggregation | decoupler.pp.pseudobulk() | Adds dependency for a ~20 line groupby+sum operation; not worth it |
| sc.pl.dotplot | Custom Plotly dotplot | Plotly is interactive but sc.pl.dotplot is already used in the post-celltyping report |

**Installation:**
```bash
pip install pydeseq2>=0.5.0
```

**Note:** PyDESeq2 is an optional dependency. Use `_check_deps(['pydeseq2'])` fast-fail pattern before loading data, matching the established pattern in `sc_tools/cli/__init__.py`.

## Architecture Patterns

### New Module: `sc_tools/tl/de.py`
```
sc_tools/tl/
    de.py                    # NEW: pseudobulk DE (SCI-01)
    celltype/
        annotate.py          # MODIFY: add panel guard (SCI-04)
sc_tools/qc/
    sample_qc.py             # MODIFY: add subject_id validation (SCI-03)
    report.py                # MODIFY: extend post-celltyping report (SCI-02)
sc_tools/cli/
    de.py                    # NEW: `sct de run` command (SCI-01)
    __init__.py              # MODIFY: register de_app
sc_tools/assets/
    post_celltyping_qc_template.html  # MODIFY: add marker validation section (SCI-02)
sc_tools/tests/
    test_de.py               # NEW: pseudobulk DE tests
    test_subject_metadata.py # NEW: subject_id validation tests
    test_panel_dispatch.py   # NEW: panel-aware dispatch tests
```

### Pattern 1: PyDESeq2 Wrapper (SCI-01)
**What:** Convert AnnData to PyDESeq2 inputs, run DE per celltype, collect results
**When to use:** `sct de run` command
**Example:**
```python
# Source: PyDESeq2 docs (https://pydeseq2.readthedocs.io/en/latest/)
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def run_pseudobulk_de(
    adata: AnnData,
    condition_key: str,
    subject_key: str = "subject_id",
    celltype_key: str = "celltype",
    formula: str | None = None,
    min_subjects_per_group: int = 3,
    min_cells_per_combo: int = 10,
) -> dict[str, pd.DataFrame]:
    """Run pseudobulk DE per celltype. Returns {celltype: results_df}."""
    results = {}
    for ct in adata.obs[celltype_key].unique():
        # 1. Subset to celltype
        mask = adata.obs[celltype_key] == ct
        sub = adata[mask]

        # 2. Aggregate: sum raw counts by subject_id
        # Group by subject_id, sum counts
        groups = sub.obs.groupby(subject_key)
        # Filter groups with < min_cells_per_combo
        valid_subjects = [s for s, g in groups if len(g) >= min_cells_per_combo]

        # 3. Build count matrix (subjects x genes) and metadata
        count_matrix = pd.DataFrame(...)  # subjects x genes, integer counts
        metadata = pd.DataFrame(...)      # subjects x covariates

        # 4. Check minimum subjects per condition group
        # Skip celltype if < min_subjects_per_group in any condition

        # 5. Build formula
        design = formula or f"~ {condition_key}"

        # 6. Run PyDESeq2
        dds = DeseqDataSet(counts=count_matrix, metadata=metadata, design=design)
        dds.deseq2()

        stat_res = DeseqStats(dds, contrast=[condition_key, "treatment", "control"])
        stat_res.summary()

        results[ct] = stat_res.results_df

    return results
```

### Pattern 2: Panel Guard in annotate_celltypes (SCI-04)
**What:** Guard clause at top of `annotate_celltypes()` checking `adata.n_vars`
**When to use:** Before dispatching to backend
**Example:**
```python
# In sc_tools/tl/celltype/annotate.py
PANEL_VALIDATED_METHODS = {"sctype", "custom_gates"}
WHOLE_TRANSCRIPTOME_METHODS = {"celltypist", "scgpt", "geneformer", "scarches", "singler"}

def annotate_celltypes(adata, method, *, force_method=False, **kwargs):
    panel_detected = adata.n_vars < 1000

    if panel_detected and method in WHOLE_TRANSCRIPTOME_METHODS:
        if force_method:
            import logging
            logging.getLogger(__name__).warning(
                "Panel detected (n_vars=%d) but --force-method used for '%s'. "
                "Results may be unreliable.", adata.n_vars, method
            )
        else:
            raise SCToolsDataError(
                f"Method '{method}' requires whole-transcriptome data but panel detected "
                f"(n_vars={adata.n_vars}). Use sctype or custom_gates, or pass force_method=True.",
                suggestion=f"Use --method sctype or --force-method"
            )

    # ... existing dispatch logic
```

### Pattern 3: Subject ID Validation (SCI-03)
**What:** Validation function checking subject_id presence and distinctness from library_id
**When to use:** At ingestion/QC for multi-sample projects
**Example:**
```python
def validate_subject_metadata(
    adata: AnnData,
    *,
    multi_sample: bool | None = None,
    subject_key: str = "subject_id",
    library_key: str = "library_id",
    condition_key: str | None = None,
    batch_key: str | None = None,
) -> list[str]:
    """Return list of warning messages. Empty = all good."""
    warnings = []

    # Auto-detect multi-sample
    if multi_sample is None:
        multi_sample = (library_key in adata.obs.columns and
                       adata.obs[library_key].nunique() > 1)

    if multi_sample:
        if subject_key not in adata.obs.columns:
            warnings.append(f"Multi-sample project but '{subject_key}' not found in obs.")
        elif library_key in adata.obs.columns:
            # Check distinctness
            if adata.obs[[subject_key, library_key]].drop_duplicates().shape[0] == \
               adata.obs[library_key].nunique():
                warnings.append(
                    f"'{subject_key}' maps 1:1 to '{library_key}' -- "
                    "these may be the same. subject_id should represent biological replicates."
                )
    else:
        if subject_key not in adata.obs.columns:
            warnings.append(f"'{subject_key}' not found in obs (single-sample: informational only).")

    # Batch-condition confounding
    if condition_key and batch_key and condition_key in adata.obs.columns and batch_key in adata.obs.columns:
        crosstab = pd.crosstab(adata.obs[batch_key], adata.obs[condition_key])
        if (crosstab > 0).sum(axis=1).max() == 1:
            warnings.append("Perfect batch-condition confounding detected.")

    return warnings
```

### Pattern 4: Marker Validation Report Extension (SCI-02)
**What:** Extend existing `generate_post_celltyping_report()` with marker validation summary
**When to use:** `sct report post-celltyping`
**Example:**
```python
def compute_marker_validation(
    adata: AnnData,
    celltype_key: str,
    marker_genes: dict[str, list[str]],
    threshold: float = 0.1,
) -> pd.DataFrame:
    """Compute mean expression of canonical markers per celltype.

    Returns DataFrame with columns: celltype, marker_gene, mean_expr, flagged.
    """
    # For each celltype, compute mean expression of its canonical markers
    # Flag celltypes where mean expression of ALL markers is below threshold
    ...
```

### Anti-Patterns to Avoid
- **Using normalized counts for pseudobulk:** PyDESeq2 requires raw integer counts. Extract from `adata.raw.X` or `adata.layers['counts']`, never from `adata.X` after normalization.
- **Pseudobulk by library_id instead of subject_id:** Biological replicates must be at the subject level. Technical replicates (libraries from same subject) should be summed together.
- **Blocking on warnings:** D-07, D-08, D-09 all specify warn-don't-block. Validation functions return warning lists, never raise.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Differential expression | Custom Wilcoxon/t-test pseudobulk | PyDESeq2 `DeseqDataSet` + `DeseqStats` | DESeq2 handles size factor estimation, dispersion shrinkage, Cook's outlier detection, BH correction |
| Dotplot rendering | Custom matplotlib dotplot | `scanpy.pl.dotplot` | Already used in post-celltyping report; handles groupby, gene ordering, dendrogram |
| Count aggregation | Complex sparse matrix ops | `pandas.DataFrame.groupby().sum()` after `.toarray()` | For pseudobulk, the per-subject count matrix is small enough for dense pandas ops |
| Multiple testing correction | Manual BH correction | PyDESeq2 does this internally | `padj` column in `stat_res.results_df` |

**Key insight:** PyDESeq2 handles the entire DESeq2 statistical pipeline (size factors, dispersion, shrinkage, testing, correction). The wrapper only needs to handle AnnData-to-pandas conversion and per-celltype iteration.

## Common Pitfalls

### Pitfall 1: Using Normalized Counts for PyDESeq2
**What goes wrong:** PyDESeq2 expects raw integer counts. Passing log-normalized or scaled data produces garbage results silently.
**Why it happens:** `adata.X` is often normalized; raw counts may be in `adata.raw.X` or `adata.layers['counts']`.
**How to avoid:** Explicitly check for raw counts: verify values are non-negative integers. Try `adata.layers.get('counts')` first, then `adata.raw.X`, then `adata.X` with integer validation.
**Warning signs:** Very low dispersion estimates, all genes called significant or none called significant.

### Pitfall 2: Subject ID == Library ID
**What goes wrong:** If subject_id is just a copy of library_id, pseudobulk treats each library as an independent biological replicate, inflating statistical power (pseudoreplication).
**Why it happens:** Users may not understand the distinction, or single-library-per-subject designs genuinely have 1:1 mapping.
**How to avoid:** Validate at registration that subject_id is not identical to library_id when there are multiple libraries. Warn but don't block for legitimate 1:1 designs.
**Warning signs:** n_subjects == n_libraries AND multiple libraries per condition group.

### Pitfall 3: Panel Detection False Positives
**What goes wrong:** A heavily filtered whole-transcriptome dataset might have `n_vars < 1000` after HVG selection, triggering panel mode incorrectly.
**Why it happens:** HVG filtering can reduce gene count below the 1000 threshold.
**How to avoid:** Check `n_vars` on the FULL adata (or `adata.raw.n_vars` if available), not on the HVG-filtered subset. Document that panel detection should happen before HVG selection.
**Warning signs:** Dataset from Visium/scRNA platform but n_vars < 1000.

### Pitfall 4: Sparse Matrix Memory Spike During Aggregation
**What goes wrong:** Calling `.toarray()` on the full sparse count matrix for aggregation can spike memory.
**Why it happens:** Large datasets (>500K cells, >30K genes) have huge dense matrices.
**How to avoid:** Aggregate per-celltype subsets (already smaller). For each celltype subset, convert to dense and group-sum. The per-subject result matrix is always small.
**Warning signs:** MemoryError during aggregation step.

### Pitfall 5: PyDESeq2 Contrast Specification
**What goes wrong:** The `contrast` parameter requires exact level names from the metadata column, and order matters (tested vs reference).
**Why it happens:** PyDESeq2 uses `[variable, tested_level, reference_level]` format; getting the order wrong reverses log2FC sign.
**How to avoid:** Let the user specify condition column and reference level. Default to alphabetically first level as reference (matching DESeq2 convention).
**Warning signs:** log2FC signs are opposite of expected.

## Code Examples

### Pseudobulk Count Aggregation
```python
# Source: Standard pandas groupby pattern for pseudobulk
import scipy.sparse as sp

def aggregate_pseudobulk(
    adata: AnnData,
    subject_key: str,
    celltype_key: str,
    min_cells: int = 10,
    layer: str | None = None,
) -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    """Aggregate counts by subject+celltype.

    Returns {celltype: (count_df, metadata_df)} where count_df is subjects x genes.
    """
    # Get raw counts
    if layer and layer in adata.layers:
        X = adata.layers[layer]
    elif adata.raw is not None:
        X = adata.raw.X
    else:
        X = adata.X

    results = {}
    for ct in adata.obs[celltype_key].cat.categories:
        ct_mask = adata.obs[celltype_key] == ct
        ct_adata = adata[ct_mask]
        ct_X = X[ct_mask] if not hasattr(X, 'shape') else X[ct_mask.values]

        # Group by subject_id and sum
        subjects = ct_adata.obs[subject_key]
        unique_subjects = subjects.unique()

        counts_list = []
        valid_subjects = []
        for subj in unique_subjects:
            subj_mask = subjects == subj
            n_cells = subj_mask.sum()
            if n_cells < min_cells:
                continue
            subj_counts = ct_X[subj_mask.values]
            if sp.issparse(subj_counts):
                subj_counts = subj_counts.toarray()
            counts_list.append(subj_counts.sum(axis=0).flatten())
            valid_subjects.append(subj)

        if len(valid_subjects) < 2:
            continue

        count_df = pd.DataFrame(
            np.array(counts_list),
            index=valid_subjects,
            columns=adata.var_names if adata.raw is None else adata.raw.var_names,
        ).astype(int)

        # Build metadata for valid subjects (one row per subject)
        meta_df = (ct_adata.obs.loc[ct_adata.obs[subject_key].isin(valid_subjects)]
                   .drop_duplicates(subset=[subject_key])
                   .set_index(subject_key)
                   .loc[valid_subjects])

        results[ct] = (count_df, meta_df)

    return results
```

### PyDESeq2 Invocation Per Celltype
```python
# Source: PyDESeq2 0.5.x API (https://pydeseq2.readthedocs.io/en/latest/)
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def run_deseq2_for_celltype(
    count_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    design: str,
    contrast: list[str],
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Run PyDESeq2 on a single celltype's pseudobulk counts."""
    dds = DeseqDataSet(
        counts=count_df,
        metadata=metadata_df,
        design=design,
        refit_cooks=True,
    )
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=contrast, alpha=alpha)
    stat_res.summary()

    # results_df has: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
    return stat_res.results_df
```

### Batch-Condition Confounding Detection
```python
# Source: Standard epidemiological confounding check
def check_confounding(
    obs: pd.DataFrame,
    batch_key: str,
    condition_key: str,
) -> bool:
    """Return True if batch perfectly confounds condition (every batch has only one condition level)."""
    crosstab = pd.crosstab(obs[batch_key], obs[condition_key])
    # Perfect confounding: each batch maps to exactly one condition
    return bool((crosstab > 0).sum(axis=1).eq(1).all())
```

### CLI Command Pattern for DE
```python
# Following established pattern from sc_tools/cli/qc.py
de_app = typer.Typer(help="Differential expression commands")

@de_app.command("run")
@cli_handler
def de_run(
    file: str = typer.Argument(..., help="Path to cell-typed h5ad"),
    condition: str = typer.Option(..., "--condition", "-c", help="Condition column in obs"),
    project_dir: str = typer.Option(".", "--project-dir", help="Project directory for output"),
    subject_col: str = typer.Option("subject_id", "--subject-col", help="Subject ID column"),
    celltype_col: str = typer.Option("celltype", "--celltype-col", help="Cell type column"),
    reference: str = typer.Option(None, "--reference", help="Reference level for contrast"),
    formula: str = typer.Option(None, "--formula", help="Custom design formula override"),
    min_subjects: int = typer.Option(3, "--min-subjects", help="Min subjects per condition"),
    min_cells: int = typer.Option(10, "--min-cells", help="Min cells per subject+celltype"),
) -> None:
    """Run pseudobulk DE per celltype (SCI-01)."""
    _check_deps(["pydeseq2", "scanpy", "anndata"])
    # ... load, validate, aggregate, run, save CSVs, return CLIResult
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| R DESeq2 via rpy2 | PyDESeq2 (pure Python) | 2023 (v0.1), scverse adoption Dec 2025 | No R dependency; native anndata integration |
| Cell-level DE (Wilcoxon) | Pseudobulk DE | 2020+ (Squair et al., 2021) | Correct type-I error control; treats subjects as statistical units |
| decoupler.get_pseudobulk | Manual pandas groupby | Ongoing | decoupler adds a dependency for ~20 lines of code; manual is simpler |
| PyDESeq2 owkin/PyDESeq2 | PyDESeq2 scverse/PyDESeq2 | Dec 2025 | Maintenance transferred to scverse community |

**Deprecated/outdated:**
- Cell-level Wilcoxon tests for multi-sample DE: inflates sample size, does not account for subject-level variation
- rpy2 + R DESeq2: unnecessary complexity now that PyDESeq2 exists

## Open Questions

1. **Raw count layer name convention**
   - What we know: Some pipelines store raw counts in `adata.layers['counts']`, others in `adata.raw.X`, others leave raw in `adata.X`
   - What's unclear: Which convention sc_tools uses across modalities
   - Recommendation: Try `layers['counts']` > `raw.X` > `X` with integer validation. Add `--layer` flag for explicit override.

2. **Condition column name**
   - What we know: There is no standard obs column for experimental condition in existing sc_tools data
   - What's unclear: Whether to enforce a convention or always require `--condition` flag
   - Recommendation: Require `--condition` flag (no default). This is project-specific and cannot be auto-detected.

3. **Marker gene source for validation**
   - What we know: `sc_tools/tl/gene_sets.py` has MSigDB Hallmark sets, GMT loader, custom JSON loader. `marker_genes` dict is already accepted by `generate_post_celltyping_report()`.
   - What's unclear: Whether canonical celltype markers are already bundled or need to be provided per-project
   - Recommendation: Accept user-provided marker file (JSON: `{celltype: [genes]}`). Fall back to the marker_db used during cell typing if available from `adata.uns`.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| pydeseq2 | SCI-01 (pseudobulk DE) | Not installed | -- | Install as optional dep (`pip install pydeseq2>=0.5.0`); `_check_deps` fast-fail |
| scanpy | SCI-02 (dotplot) | Installed | 1.11.5 | -- |
| anndata | All | Installed | 0.11.4 | -- |
| pandas | All | Installed | -- | -- |
| numpy | All | Installed | -- | -- |
| scipy | SCI-01 (sparse ops) | Installed | -- | -- |

**Missing dependencies with no fallback:**
- None (pydeseq2 is optional and handled by `_check_deps`)

**Missing dependencies with fallback:**
- pydeseq2: Not installed. `_check_deps(['pydeseq2'])` will fast-fail with install instructions. Add to `[project.optional-dependencies]` in pyproject.toml under a `de` extra.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest (configured in pyproject.toml) |
| Config file | `pyproject.toml` `[tool.pytest.ini_options]` |
| Quick run command | `pytest sc_tools/tests/test_de.py -x -v` |
| Full suite command | `pytest sc_tools/tests/ -v --tb=short` |

### Phase Requirements to Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| SCI-01 | Pseudobulk aggregation + PyDESeq2 pipeline | unit | `pytest sc_tools/tests/test_de.py -x` | Wave 0 |
| SCI-01 | CLI `sct de run` end-to-end | integration | `pytest sc_tools/tests/test_de.py::TestDECLI -x` | Wave 0 |
| SCI-02 | Marker validation compute + flag table | unit | `pytest sc_tools/tests/test_marker_validation.py -x` | Wave 0 |
| SCI-02 | Post-celltyping report includes validation section | integration | `pytest sc_tools/tests/test_marker_validation.py::TestReportIntegration -x` | Wave 0 |
| SCI-03 | Subject ID validation (multi-sample vs single) | unit | `pytest sc_tools/tests/test_subject_metadata.py -x` | Wave 0 |
| SCI-03 | Batch-condition confounding detection | unit | `pytest sc_tools/tests/test_subject_metadata.py::TestConfounding -x` | Wave 0 |
| SCI-04 | Panel detection + method restriction | unit | `pytest sc_tools/tests/test_panel_dispatch.py -x` | Wave 0 |
| SCI-04 | Force-method override with warning | unit | `pytest sc_tools/tests/test_panel_dispatch.py::TestForceMethod -x` | Wave 0 |

### Sampling Rate
- **Per task commit:** `pytest sc_tools/tests/test_de.py sc_tools/tests/test_subject_metadata.py sc_tools/tests/test_panel_dispatch.py sc_tools/tests/test_marker_validation.py -x -v`
- **Per wave merge:** `pytest sc_tools/tests/ -v --tb=short`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `sc_tools/tests/test_de.py` -- covers SCI-01 (pseudobulk aggregation, PyDESeq2 wrapper, CLI command)
- [ ] `sc_tools/tests/test_subject_metadata.py` -- covers SCI-03 (subject_id validation, confounding check)
- [ ] `sc_tools/tests/test_panel_dispatch.py` -- covers SCI-04 (panel detection, method restriction, force override)
- [ ] `sc_tools/tests/test_marker_validation.py` -- covers SCI-02 (marker validation compute, flag table, report integration)
- [ ] `sc_tools/tests/conftest.py` -- shared fixtures: synthetic multi-subject AnnData with subject_id, library_id, condition, batch columns; panel-sized AnnData (n_vars=40)

### Test Fixture Design (Claude's Discretion)
```python
# Synthetic multi-subject AnnData for pseudobulk tests
def make_multi_subject_adata(
    n_subjects=6, n_cells_per_subject=50, n_genes=200, n_celltypes=3, seed=42
):
    """Create AnnData with subject_id, library_id, condition, batch, celltype."""
    rng = np.random.default_rng(seed)
    n_obs = n_subjects * n_cells_per_subject
    X = rng.negative_binomial(5, 0.3, (n_obs, n_genes))  # Raw integer counts

    subjects = np.repeat([f"subject_{i}" for i in range(n_subjects)], n_cells_per_subject)
    conditions = np.where(np.repeat(np.arange(n_subjects), n_cells_per_subject) < n_subjects // 2,
                         "control", "treatment")
    celltypes = rng.choice([f"type_{i}" for i in range(n_celltypes)], n_obs)
    libraries = [f"lib_{s}" for s in subjects]  # 1:1 for simplicity

    obs = pd.DataFrame({
        "subject_id": pd.Categorical(subjects),
        "library_id": pd.Categorical(libraries),
        "condition": pd.Categorical(conditions),
        "batch": pd.Categorical([f"batch_{int(i >= n_subjects//2)}" for i in
                                  np.repeat(range(n_subjects), n_cells_per_subject)]),
        "celltype": pd.Categorical(celltypes),
    }, index=[f"cell_{i}" for i in range(n_obs)])

    return ad.AnnData(X.astype(np.float32), obs=obs,
                      var=pd.DataFrame(index=[f"Gene{i}" for i in range(n_genes)]))
```

## Sources

### Primary (HIGH confidence)
- PyDESeq2 official docs (v0.5.4): https://pydeseq2.readthedocs.io/en/latest/auto_examples/plot_minimal_pydeseq2_pipeline.html -- API usage, constructor arguments, results format
- PyDESeq2 PyPI: https://pypi.org/project/pydeseq2/ -- verified version 0.5.2 is latest stable
- sc_tools codebase: `sc_tools/tl/celltype/annotate.py`, `sc_tools/qc/report.py`, `sc_tools/cli/__init__.py` -- existing patterns, extension points, report infrastructure

### Secondary (MEDIUM confidence)
- PyDESeq2 GitHub (scverse/PyDESeq2): https://github.com/scverse/PyDESeq2 -- maintenance transferred to scverse Dec 2025
- decoupler pseudobulk docs: https://decoupler.readthedocs.io/en/latest/notebooks/scell/rna_psbk.html -- alternative aggregation approach (not recommended, adds dependency)

### Tertiary (LOW confidence)
- None

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- PyDESeq2 is the only pure-Python DESeq2 implementation; API verified from official docs
- Architecture: HIGH -- all extension points exist in the codebase and have been read; patterns are established
- Pitfalls: HIGH -- raw counts requirement, pseudoreplication, panel false positives are well-documented in literature

**Research date:** 2026-03-24
**Valid until:** 2026-04-24 (stable domain; PyDESeq2 API unlikely to break)
