# Project Research Summary

**Project:** sc_tools v2.0 — Report Plots and Sample Concat
**Domain:** Spatial transcriptomics HTML report pipeline + agent-native CLI
**Researched:** 2026-03-27
**Confidence:** HIGH

## Executive Summary

The v2.0 milestone extends sc_tools' existing report pipeline with spatial/UMAP plot coverage that is currently absent or incomplete, and adds a `sct concat` CLI command to formalize multi-sample assembly as a tracked pipeline phase. This is fundamentally an integration and wiring problem, not a new-library problem: every capability needed (scanpy spatial plots, Plotly interactive charts, anndata concat, base64 PNG embedding) is already present in the codebase. The work is adding the missing connections — spatial plot loops in two report generators, a 4-panel wrapper function in `pl/spatial.py`, two HTML template section additions, and a new CLI command that wraps an already-implemented backend function (`ingest/loaders.py:concat_samples`).

The recommended approach is to build in strict dependency order: `sct concat` first (no report dependencies, easily validated in isolation), then the `qc_spatial_4panel()` plot function, then the pre-filter report spatial section, then the celltype report spatial section, and finally the lowest-risk change — extending the integration report's UMAP color key list. The two most important correctness invariants are: (1) `uns['spatial']` must be preserved through concat using `uns_merge="unique"` with a post-concat verification check, and (2) the `qc_pass_spot` obs column must be written by `sct qc run` before filtering rather than re-derived in the report. Both are cross-component contracts that, if violated, cause silent downstream failures.

The primary risk is HTML report size when spatial plots scale with sample count. A 10-sample cohort with 15 cell types would produce 150 plots at approximately 1.5 MB each, crashing browsers. The MVP mitigation is strict: aggregate categorical overlay only (one figure per sample showing all cell types), with per-cell-type breakdowns deferred to v1.x as an accordion limited to top-8 by abundance. No new dependencies are required except promoting `plotly` from the `[pipeline]` optional extra to base dependencies and updating a stale CDN version tag.

## Key Findings

### Recommended Stack

All visualization, serialization, and CLI primitives needed for v2.0 are already available in the codebase. Two maintenance items are required before report work begins: promote `plotly>=5.18` from `[pipeline]` optional to base `dependencies` in `pyproject.toml` (prevents ImportError when users invoke `sct qc report` or `sct integration report` without the pipeline extras), and update the hardcoded CDN tag in `report_utils.py` line 517 from `plotly-2.27.0.min.js` to `plotly-latest.min.js`.

**Core technologies:**
- `scanpy.pl.spatial`: spatial overlay plots — already used in post-filter report; extend to pre-filter and celltype reports with the identical `sc.pl.spatial(sub, color=col, library_id=lib, show=False, ax=ax)` call pattern
- `matplotlib / fig_to_base64`: rasterized PNG embedding — authoritative pattern for all report plots; SVG is unusable at 50K–2.5M spots; base64 PNG keeps reports self-contained (no CDN for images)
- `plotly / plotly_to_html`: interactive charts (metrics, radar) — already used for benchmark sections; CDN-based, one `<script>` tag per report
- `anndata.concat`: multi-sample assembly — implemented in `ingest/loaders.py:concat_samples()`; v2.0 adds only the CLI wrapper and phase registration
- `Typer + Pydantic CLIResult`: CLI command pattern — the `@cli_handler` + `CLIResult` pattern established in v1.0 applies directly to `sct concat`
- `h5py` (stdlib): lightweight validation before full load — use for per-file obs count and var_names check in `sct concat` pre-flight without loading full AnnDatas into memory

### Expected Features

**Must have (table stakes for MVP):**
- Spatial plots (log1p_counts, log1p_genes, pct_mt) per sample in pre-filter QC report — populates the dead `spatial_multipage` template slot
- `qc_pass_spot` obs column written by `sct qc run` before the filter step (required for pass/fail spatial overlay)
- Pass/fail spatial overlay in pre-filter report reading `qc_pass_spot` from obs
- `sct concat` command: `--input` list, `--output`, `--batch-key`, obs dedup (`index_unique="-"`), uns spatial merge (`uns_merge="unique"`), CLIResult JSON stdout, provenance sidecar, PhaseSpec registration as `optional=True`
- Celltype spatial overlay per sample in post-celltyping report — one aggregate categorical figure per sample, consistent palette across samples
- `subject_id` UMAP coloring auto-detected in integration report when column is present in `adata.obs`
- Modality guard on all new spatial sections (`"spatial" in adata.uns` + Visium/VisiumHD modality check)

**Should have (differentiators, v1.x after MVP validation):**
- Common color scale across samples (cohort-level 2nd–98th percentile vmin/vmax for continuous spatial metrics)
- 4-panel per-sample spatial grid (2x2: counts/genes/mt/pass-fail) replacing the three-separate-plots approach
- Split batch-vs-bio UMAP layout (two-column figure) in integration report
- `celltype_broad` spatial overlay in celltype report
- Per-celltype spatial overlays — accordion, top-8 by abundance hard cap
- `sct concat --from-dir` with modality consistency validation and lexicographic sort

**Defer (v2+ or not at all for MVP):**
- Per-embedding UMAP colored by Leiden (doubles grid compute time; validate demand first)
- Memory-efficient backed-mode concat for VisiumHD (HIGH complexity; document RAM requirement instead)
- Spatial plots for IMC/Xenium (different coordinate conventions and image backends; separate future feature)

### Architecture Approach

The v2.0 architecture respects the existing layered separation: plot functions belong in `sc_tools/pl/`, report orchestration in `sc_tools/qc/report.py`, CLI commands in `sc_tools/cli/`, and templates in `sc_tools/assets/`. Report orchestrators call plot functions, convert figures to base64 via `fig_to_base64()`, and pass pre-rendered strings to Jinja2 — templates remain logic-free. The new `sct concat` command belongs in a new `sc_tools/cli/ingest.py` file (mirrors the `sc_tools/ingest/` backend it wraps) and is registered as an optional phase in the pipeline DAG to preserve backward compatibility with projects that ingest pre-concatenated h5ads.

**Major components and their changes:**

1. `sc_tools/pl/spatial.py` (MODIFIED) — add `qc_spatial_4panel(adata, library_id, metric_keys, figsize, dpi)` returning a single `plt.Figure`; silently skip missing metric columns (blank panels, no raise)
2. `sc_tools/qc/report.py` (MODIFIED) — add spatial 4-panel loop to `generate_pre_filter_report()`; add celltype spatial loop to `generate_post_celltyping_report()`; extend `color_keys` in `generate_post_integration_report()` for `subject_id`
3. `sc_tools/cli/ingest.py` (NEW) — `sct concat` Typer command following `@cli_handler` + `CLIResult` pattern
4. `sc_tools/cli/__init__.py` (MODIFIED) — register `ingest_app` via `app.add_typer()`
5. `sc_tools/pipeline.py` (MODIFIED) — add `concat` PhaseSpec (`optional=True`, depends on `ingest_load`)
6. `sc_tools/assets/pre_filter_qc_template.html` (MODIFIED) — add spatial-qc-plots section
7. `sc_tools/assets/post_celltyping_qc_template.html` (MODIFIED) — add celltype-spatial section

**Confirmed unchanged:** `sc_tools/ingest/loaders.py` (concat_samples complete), `sc_tools/cli/qc.py` (no changes needed), `sc_tools/qc/report_utils.py` (fig_to_base64 and render_template work as-is), `sc_tools/pl/qc_plots.py` (UMAP grid functions work as-is).

### Critical Pitfalls

1. **`uns['spatial']` loss through concat** — `uns_merge="unique"` is required in `anndata.concat`. After concat, verify `set(input_library_ids) == set(output.uns.get('spatial', {}).keys())` and raise `SCToolsDataError` with a message listing missing samples. Silent loss here breaks all downstream spatial plots with no obvious error.

2. **Re-deriving `qc_pass_spot` in the report** — re-implementing filter thresholds inside report generation creates divergence (report shows different pass/fail than what was actually applied). The `qc_pass_spot` column must be written by `sct qc run` before the filter step; the report reads it from obs.

3. **HTML size explosion with spatial plots** — N_samples × N_celltypes plots at ~1 MB each exceeds browser limits. MVP enforces aggregate categorical overlay only (one plot per sample, all cell types). Per-cell-type breakdown is deferred with a hard cap of top-8 by abundance.

4. **Loading full AnnData for pre-flight validation in `sct concat`** — for 10-sample VisiumHD at 25G each, full loads for validation require 250G RAM before writing a byte. Use h5py to read only `obs/_index` and `var/_index` for compatibility checks; load full AnnDatas only for the actual concat call.

5. **Plotly for spatial scatter plots** — at 50K–2.5M spots, Plotly renders each point as an SVG element and the HTML becomes unusable. Use `sc.pl.spatial` + `fig_to_base64(fig, dpi=120)` for all spatial plots. Reserve Plotly for small-data interactive sections (metrics tables, radar charts).

6. **Subsetting AnnData without `.copy()` for spatial plots** — a plain view slice may not carry `uns['spatial'][library_id]` correctly. Always use `adata[adata.obs[sample_col] == lib].copy()` before passing to any spatial plot function. This is the pattern already used in `generate_post_filter_report()` at report.py line 313.

## Implications for Roadmap

Based on research, the milestone decomposes into five ordered build phases:

### Phase A: sct concat + Pipeline Registration

**Rationale:** No dependencies on the report plot pipeline. Can be built, tested, and merged entirely independently. Produces the `adata.concatenated.h5ad` checkpoint that enables multi-sample QC testing. Validates the `uns['spatial']` preservation invariant before any report work begins. The plotly promotion and CDN tag fix should also be done here as housekeeping before report work starts.

**Delivers:** `sct concat` CLI command, `concat` PhaseSpec (`optional=True`), provenance sidecar, per-file h5py pre-flight validation, uns spatial merge with post-concat correctness check, `plotly` promoted to base deps, CDN tag updated.

**Addresses:** All `sct concat` table stakes from FEATURES.md — obs dedup, batch key injection, uns spatial preservation, CLIResult JSON output, pipeline phase registration.

**Avoids:** Loading full AnnData for validation (h5py only for pre-flight); `uns['spatial']` silent loss (`uns_merge="unique"` + verification).

**Research flag:** Standard patterns — Typer + CLIResult pattern from v1.0 applies directly. No additional research needed.

### Phase B: Spatial Plot Function (pl/spatial.py)

**Rationale:** The `qc_spatial_4panel()` function must exist and be independently validated before any report template changes are made. Isolating the function development prevents rework if the plot API needs revision.

**Delivers:** `qc_spatial_4panel(adata, library_id, metric_keys)` in `pl/spatial.py`, exported from `pl/__init__.py`. Handles missing metric columns gracefully with blank panels.

**Uses:** `sc.pl.spatial` (existing), `plt.subplots` (existing). Follows the `plot_spatial_continuous()` signature pattern already in the module.

**Research flag:** Standard patterns — established within the existing pl/spatial.py module. No additional research needed.

### Phase C: Pre-filter Report Spatial Section

**Rationale:** Pre-filter checkpoints are the earliest in the pipeline and simplest to test with minimal data. This phase also introduces `qc_pass_spot` in `sct qc run` — the most behavior-changing addition in the milestone (modifying a core QC command) — which should be tested in isolation before building further report sections on top of it.

**Delivers:** Spatial QC 4-panel per sample in pre-filter HTML report; `qc_pass_spot` obs column added in `sct qc run`; sample ordering by QC score (failed first) in spatial section; dead `spatial_multipage` template slot activated.

**Implements:** Per-sample spatial loop pattern (Pattern 2 from ARCHITECTURE.md), base64 PNG embedding (Pattern 1), both directly modeled on `generate_post_filter_report()` lines 302–334.

**Avoids:** Recomputing pass/fail in report (read from obs); Plotly for spatial data (use fig_to_base64); subsetting without `.copy()`.

**Research flag:** Standard patterns — directly cloned from existing post-filter spatial loop. No additional research needed.

### Phase D: Post-celltyping Report Spatial Section

**Rationale:** Same per-sample spatial loop pattern validated in Phase C. Placed later because it depends on a later pipeline checkpoint (`adata.celltyped.h5ad`) and adds the categorical color handling and modality guard that are less familiar than the continuous metrics in Phase C.

**Delivers:** Aggregate categorical per-sample celltype spatial overlay in post-celltyping HTML report; consistent palette pre-computed from all unique celltype values; Visium/VisiumHD modality guard via `get_modality_terms()`.

**Addresses:** Celltype spatial table stakes from FEATURES.md — one figure per sample, consistent palette, modality guard.

**Avoids:** Per-cell-type breakdown at MVP (HTML size explosion); IMC spatial path (different code path, future feature).

**Research flag:** Standard patterns — same loop as Phase C; `plot_spatial_categorical` already exists in `pl/spatial.py`. No additional research needed.

### Phase E: Integration Report UMAP Extension

**Rationale:** Lowest-risk change in the entire milestone — adding `subject_id` to the `color_keys` list in `generate_post_integration_report()` with an `if "subject_id" in adata.obs.columns` guard. No new functions, no template changes, no new patterns. Placed last to minimize regression risk.

**Delivers:** `subject_id` / patient-level UMAP coloring auto-detected in integration report when column is present in `adata.obs`.

**Avoids:** Computing UMAP fresh during report generation (require `X_umap` in obsm; skip section if absent); 5-method × 4-color-key UMAP panel explosion.

**Research flag:** Standard patterns — single guard condition extending an existing list. No research needed.

### Phase Ordering Rationale

- Phase A (concat) is fully independent and produces a multi-sample checkpoint enabling realistic test data for all subsequent phases.
- Phase B (plot function) must precede Phases C and D because both templates depend on the function.
- Phase C before Phase D: pre-filter checkpoint is earlier in the pipeline and simpler; the `qc_pass_spot` addition needs isolated testing before celltype work adds complexity.
- Phase E last: zero new plot infrastructure; extending an existing list carries the lowest regression risk.

### Research Flags

All five phases follow standard patterns present in the existing codebase. No phase requires a `/gsd:research-phase` call.

- **Phases A–E:** Standard patterns confirmed by direct codebase inspection. The v1.0 CLI pattern (Typer + CLIResult + provenance sidecar) is directly reusable for Phase A. The spatial loop and base64 embedding patterns from `generate_post_filter_report()` are directly reusable for Phases B–D. Phase E is a one-line change.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | All findings from direct codebase inspection; no new libraries; plotly promotion and CDN fix confirmed against pyproject.toml and report_utils.py source |
| Features | HIGH | Based on direct inspection of all four report generators and HTML templates; gap analysis against existing slots (dead spatial_multipage, missing celltype spatial) confirmed |
| Architecture | HIGH | All component boundaries, unchanged files, and build order confirmed by direct codebase inspection; every claim grounded in a specific file and line range |
| Pitfalls | HIGH | Grounded in both domain research (scanpy memory issues, AnnData concat semantics, HTML size limits) and codebase evidence (existing spatial loop pattern, migration debt) |

**Overall confidence:** HIGH

### Gaps to Address

- **`qc_pass_spot` schema documentation:** The column name and semantics (`True` = spot passes QC) should be documented as a standard obs key in the sc_tools data contract when introduced in Phase C.
- **`sct concat` optional vs. mandatory:** The PhaseSpec registers as `optional=True` and `qc_filter` continues depending on `ingest_load`. If a future project always requires concat, the DAG wiring may need revisiting — treat as a known design decision, revisit after adoption.
- **Plotly CDN on HPC without internet:** The integration report uses a CDN-based Plotly script tag. On HPC nodes without internet access, interactive chart sections fail to render. This is an existing design choice; document as a known limitation in the report footer. Not blocking MVP.
- **VisiumHD concat RAM requirement:** Loading N × 25G files fully is feasible only if RAM >= 4 × N × 25G. The MVP documents this requirement and relies on `sct estimate` to warn before execution. True backed-mode concat deferred to v2.1.

## Sources

### Primary (HIGH confidence)

- `sc_tools/qc/report.py` — existing report generators, spatial loop pattern (lines 302–334), base64 embedding
- `sc_tools/pl/spatial.py` — `plot_spatial_continuous`, `plot_spatial_categorical` signatures
- `sc_tools/pl/qc_plots.py` — `qc_umap_grid`, `qc_cluster_distribution`, `qc_celltype_abundance`
- `sc_tools/qc/report_utils.py` — `fig_to_base64`, `render_template`, `auto_detect_embeddings`, CDN tag (line 517)
- `sc_tools/ingest/loaders.py` — `concat_samples` implementation (lines 903–947)
- `sc_tools/pipeline.py` — PhaseSpec DAG structure and registration API
- `sc_tools/cli/__init__.py` — `cli_handler` decorator, `CLIResult`, app structure
- `sc_tools/assets/*.html` — template slot structure for all four report types
- `pyproject.toml` — `requires-python = ">=3.10"`, dep versions, extras
- `.planning/PROJECT.md` — v2.0 milestone targets, Phase 6 subject metadata status
- [scanpy.pl.spatial docs](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pl.spatial.html) — `ax=`, `vmin`/`vmax`, `library_id` params
- [plotly.io.to_html docs](https://plotly.com/python-api-reference/generated/plotly.io.to_html.html) — `include_plotlyjs` options

### Secondary (MEDIUM confidence)

- [Anthropic: Writing effective tools for AI agents](https://www.anthropic.com/engineering/writing-tools-for-agents) — agent CLI design principles
- [scanpy Issue #2365](https://github.com/scverse/scanpy/issues/2365) — backed mode memory behavior with large h5ad files
- [CLI framework comparison (2025)](https://dasroot.net/posts/2025/12/building-cli-tools-python-click-typer-argparse/) — Typer vs Click vs argparse tradeoffs

### Tertiary (LOW confidence)

- HTML size estimates (~1 MB per 150-DPI spatial PNG, browser crash limits) — estimates based on typical matplotlib output; validate against actual DLBCL project data in Phase C

---
*Research completed: 2026-03-27*
*Ready for roadmap: yes*
