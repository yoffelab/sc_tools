# Phase 1: Benchmark Fixes - Context

**Gathered:** 2026-03-20
**Status:** Ready for planning

<domain>
## Phase Boundary

Fix integration benchmark bugs (BM-01 through BM-07) and add unit tests (TST-01 through TST-03) so that benchmarking produces correct, reproducible, memory-efficient results on real project data. No new benchmark capabilities — this is a fix-and-harden phase.

</domain>

<decisions>
## Implementation Decisions

### h5py Loading API (BM-01)
- Add `embedding_files` parameter to `compare_integrations()`: dict mapping method name to h5ad file path
- Backwards compatible — existing `adata` + `embeddings` dict API remains functional
- When `embedding_files` is provided, load only `obsm` embeddings + required `obs` columns via h5py (no full AnnData load)
- Extract `batch_key` plus a single optional `bio_key` (replaces the celltype-only assumption)
- `bio_key` defaults to `celltype_key` but can be any clinically relevant variable (e.g., `condition`, `disease_status`) for bio conservation metrics
- Cell type is often missing in early pipeline stages — the API must work without it (batch-only metrics)
- Peak memory target: <2GB for 2.5M-cell dataset when using h5py path

### Claude's Discretion
- NaN handling strategy for resolVI embeddings (BM-02) — per-embedding valid-cell masking is the natural approach
- `_stratified_subsample` fix (BM-04) — replace `sorted(indices)[:n]` with proportional group downsampling
- Runtime tracking implementation (BM-05) — add timing around each method in `run_integration_benchmark`
- Parameter provenance storage (BM-06) — store batch_weight, bio_weight, seed, resolution in DataFrame attrs or as columns
- `_recipe_targeted_panel` scVI fix (BM-07) — skip normalization when scVI is selected so raw counts are preserved
- Test fixture design (TST-01/02/03) — synthetic data with controlled properties for unit tests

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Benchmark code
- `sc_tools/bm/integration.py` — All benchmark functions: `compare_integrations`, `compute_integration_metrics`, `run_integration_benchmark`, `_stratified_subsample`, `run_full_integration_workflow`
- `sc_tools/bm/report.py` — HTML report generation for benchmark results (consumer of benchmark DataFrame)
- `sc_tools/pl/benchmarking.py` — Plotting functions for integration comparison (consumer of benchmark DataFrame)

### Preprocessing (targeted panel bug)
- `sc_tools/pp/recipes.py` — `_recipe_targeted_panel` (line 226): normalizes before scVI raw count check
- `sc_tools/pp/strategy.py` — `SmallStrategy` class handling normalization pipeline

### Requirements
- `.planning/REQUIREMENTS.md` — BM-01 through BM-07, TST-01 through TST-03 define exact acceptance criteria

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `_asw_batch_sklearn`, `_pcr_sklearn`, `_ari_sklearn`, `_nmi_sklearn`: sklearn fallback metrics already implemented
- `_compute_scib_metrics` / `_compute_scib_metrics_batch_only`: scib-metrics integration working
- `compute_composite_score`: batch/bio weighting with configurable weights
- `SmallStrategy` in `sc_tools/pp/strategy.py`: modular normalization/integration pipeline

### Established Patterns
- scib-metrics availability checked at import time (`_HAS_SCIB` flag)
- Metrics clamped to [0, 1] range after computation
- `run_integration_benchmark` saves per-method h5ad files to `output/tmp/integration_test/` — this is the natural input for the new h5py loading path
- `DataFrame.attrs` used to annotate metadata (e.g., `scib_fallback`)

### Integration Points
- `compare_integrations()` is called by `run_integration_benchmark()` and the MCP tool `run_full_phase`
- `_recipe_targeted_panel()` called from `preprocess()` dispatcher for xenium/cosmx/visium_hd_cell modalities
- Benchmark DataFrame consumed by `sc_tools/bm/report.py` and `sc_tools/pl/benchmarking.py`
- Existing tests: `test_bm_report.py` (report), `test_pp.py` (preprocessing) — no `test_bm_integration.py` exists

</code_context>

<specifics>
## Specific Ideas

- celltype_key is often absent in early pipeline stages — bio conservation should accept any clinically relevant variable (condition, disease_status, etc.) via a `bio_key` parameter
- The h5py loading path should mirror the existing `save_intermediates` output structure (per-method h5ad files in a directory)

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 01-benchmark-fixes*
*Context gathered: 2026-03-20*
