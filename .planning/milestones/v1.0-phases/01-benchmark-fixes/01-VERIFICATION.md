---
phase: 01-benchmark-fixes
verified: 2026-03-21T14:00:00Z
status: passed
score: 11/11 must-haves verified
re_verification: false
---

# Phase 01: Benchmark Fixes Verification Report

**Phase Goal:** Fix sc_tools internals — integration benchmark bugs (NaN handling, subsampling bias, runtime tracking, parameter provenance), h5py-based embedding loading, and targeted panel recipe fix.
**Verified:** 2026-03-21T14:00:00Z
**Status:** passed
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | NaN embedding rows are filtered before metric computation and logged as warnings | VERIFIED | `_mask_nan_rows` at line 45; per-embedding NaN masking in `compare_integrations` loop at lines 524-531; `TestNaNHandling.test_compare_integrations_with_nan_embedding` passes |
| 2 | `compare_integrations` accepts `subsample_n` parameter and subsamples before metrics | VERIFIED | Signature line 449: `subsample_n: int | None = None`; wired at lines 511-512; `TestCompareIntegrationsExtended.test_subsample_n_parameter` passes |
| 3 | `_stratified_subsample` preserves group proportions without index truncation bias | VERIFIED | Old `sorted(indices)[:n]` line is GONE; replaced with `rng.choice(indices, size=n, replace=False)` at line 934; `TestStratifiedSubsampleFix.test_proportional_allocation` and `test_exact_count` pass |
| 4 | `run_integration_benchmark` output DataFrame includes `runtime_s` column per method | VERIFIED | `runtimes` dict at line 689; `time.perf_counter()` at line 692; `finally` block at lines 776-781; `runtime_s` column added at line 805 |
| 5 | `compare_integrations` stores `benchmark_params` in `DataFrame.attrs` with batch_weight, bio_weight, subsample_n, seed, resolution, and use_scib keys | VERIFIED | Lines 603-610 set all 6 keys; `test_benchmark_params_in_attrs` asserts all 6 keys with correct values |
| 6 | `compute_integration_metrics` returns valid scores for synthetic single-batch and no-celltype data | VERIFIED | `resolution: float = 1.0` added at line 233; `test_single_batch_metrics` (asw_batch=0.0) and `test_no_celltype_metrics` pass |
| 7 | `_stratified_subsample` handles edge cases: n > n_obs returns full copy, single group works | VERIFIED | `n >= adata.n_obs` returns `adata.copy()` at lines 920-921; `test_n_greater_than_nobs` and `test_single_group` pass |
| 8 | `compare_integrations` loads embeddings from h5ad files via h5py without full AnnData load | VERIFIED | `_load_embedding_h5py` at lines 63-120; `embedding_files` parameter at line 452; file-based path at lines 544-590; `TestEmbeddingFilesParameter` passes |
| 9 | `compare_integrations` accepts `bio_key` parameter that defaults to `celltype_key` for bio conservation metrics | VERIFIED | `bio_key: str | None = None` at line 444; resolved at line 505; `TestBioKeyParameter` passes |
| 10 | h5py loading path reconstructs categorical obs columns from codes+categories | VERIFIED | `grp["categories"]` and `grp["codes"]` reconstruction at lines 100-120; `TestLoadEmbeddingH5py.test_loads_categorical_obs` passes |
| 11 | `_recipe_targeted_panel` skips normalization when `integration=='scvi'` so raw counts are preserved | VERIFIED | `if integration == "scvi":` at line 244 precedes any `normalize_total` call; old warning string absent; `TestTargetedPanelScVI` all 3 tests pass |

**Score:** 11/11 truths verified

---

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `sc_tools/bm/integration.py` | Fixed `_stratified_subsample`, NaN masking, `subsample_n`, `runtime_s`, param provenance with resolution, `_load_embedding_h5py`, `embedding_files`, `bio_key` | VERIFIED | All functions present and substantive; wired throughout |
| `sc_tools/tests/test_integration_benchmark.py` | Unit tests for all BM-02..06, TST-01..03, and BM-01 h5py path | VERIFIED | 7 test classes present: `TestStratifiedSubsampleFix`, `TestNaNHandling`, `TestCompareIntegrationsExtended`, `TestLoadEmbeddingH5py`, `TestEmbeddingFilesParameter`, `TestBioKeyParameter` |
| `sc_tools/pp/recipes.py` | Fixed `_recipe_targeted_panel` with scVI branch before normalization | VERIFIED | `if integration == "scvi":` at line 244; `flavor="seurat_v3"` at line 246 |
| `sc_tools/tests/test_pp.py` | Tests for targeted panel scVI normalization skip | VERIFIED | `class TestTargetedPanelScVI` at line 633 with 3 substantive tests |

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `compare_integrations` | `_mask_nan_rows` | Called per-embedding before `compute_integration_metrics` | WIRED | NaN masking applied inline (lines 524-531); `_mask_nan_rows` also exposed as standalone helper |
| `compare_integrations` | `_stratified_subsample` | `subsample_n` parameter triggers subsampling before metric loop | WIRED | Lines 511-512: `if adata is not None and subsample_n is not None and subsample_n < adata.n_obs` |
| `run_integration_benchmark` | `time.perf_counter` | Timing each method integration | WIRED | `import time` at line 12; `t0 = time.perf_counter()` at line 692; `finally` block at line 776-781 |
| `compare_integrations` | `_load_embedding_h5py` | Called when `embedding_files` parameter is provided | WIRED | Lines 561: `_load_embedding_h5py(fpath, obsm_key, obs_keys_needed)` inside embedding_files loop |
| `_recipe_targeted_panel` | scVI branch | `if integration == "scvi":` before `normalize_total` | WIRED | Line 244 branches before line 261 (`normalize_total`) |

---

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| BM-01 | 01-02 | h5py loading of pre-computed embeddings (<2GB for 2.5M cells) | SATISFIED | `_load_embedding_h5py` + `embedding_files` param; memory arithmetic test passes (300MB < 2GB) |
| BM-02 | 01-01 | Filter NaN rows per-embedding before metric computation | SATISFIED | Per-embedding NaN masking in `compare_integrations` loop; `_mask_nan_rows` helper |
| BM-03 | 01-01 | Configurable `subsample_n` parameter in `compare_integrations` | SATISFIED | `subsample_n: int | None = None` parameter wired to `_stratified_subsample` |
| BM-04 | 01-01 | Fix `_stratified_subsample` truncation bias | SATISFIED | `sorted(indices)[:n]` eliminated; `rng.choice(indices, size=n, replace=False)` used |
| BM-05 | 01-01 | Add `runtime_s` column to `run_integration_benchmark` output | SATISFIED | `runtimes` dict + `time.perf_counter()` + `df["runtime_s"] = df["method"].map(runtimes)` |
| BM-06 | 01-01 | Store benchmark parameters alongside results | SATISFIED | `df.attrs["benchmark_params"]` with all 6 keys: batch_weight, bio_weight, subsample_n, seed, resolution, use_scib |
| BM-07 | 01-02 | Fix `_recipe_targeted_panel` — skip normalization for scVI | SATISFIED | `if integration == "scvi":` guard at line 244 before any `normalize_total` call |
| TST-01 | 01-01 | Unit tests for `compute_integration_metrics` (synthetic, single batch, no celltype) | SATISFIED | `test_single_batch_metrics` and `test_no_celltype_metrics` in `TestCompareIntegrationsExtended` |
| TST-02 | 01-01 | Unit tests for `compare_integrations` (NaN embeddings, subsampling) | SATISFIED | `TestNaNHandling.test_compare_integrations_with_nan_embedding`; `test_subsample_n_parameter` |
| TST-03 | 01-01 | Unit tests for `_stratified_subsample` (proportionality, n > n_obs, single group) | SATISFIED | `TestStratifiedSubsampleFix`: `test_proportional_allocation`, `test_n_greater_than_nobs`, `test_single_group` |

All 10 requirements satisfied. No orphaned requirements — REQUIREMENTS.md Traceability table lists BM-01..07 and TST-01..03 as Phase 1 / Complete.

---

### Anti-Patterns Found

No blockers or warnings detected.

- No TODO/FIXME/PLACEHOLDER comments in modified files
- No stub implementations (empty returns, hardcoded empty data)
- No disconnected wiring — all helpers are called from their intended call sites
- Old buggy line `sorted(indices)[:n]` confirmed absent
- Old warning string "scVI expects raw counts" confirmed absent

---

### Test Execution Results

```
pytest sc_tools/tests/test_integration_benchmark.py sc_tools/tests/test_pp.py::TestTargetedPanelScVI
59 passed, 2 skipped in 19.63s
```

All 59 collected tests pass. 2 skipped tests are pre-existing (not related to this phase).

---

### Commit Verification

All 6 task commits confirmed in git log:

| Commit | Type | Description |
|--------|------|-------------|
| `ec905c4` | test | Add failing tests for BM-02..06 and TST-01..03 (RED) |
| `773bd96` | feat | Fix integration benchmark bugs BM-02..06 with TDD (GREEN) |
| `3452b56` | test | Add failing tests for h5py loading, embedding_files, bio_key (RED) |
| `4668827` | feat | Add h5py embedding loading, embedding_files, bio_key (GREEN) |
| `c3e6973` | test | Add failing tests for targeted panel scVI normalization bug (RED) |
| `118cf3b` | fix | Fix `_recipe_targeted_panel` to skip normalization for scVI (GREEN) |

---

### Human Verification Required

None. All phase deliverables are verifiable programmatically via test execution, code inspection, and grep.

---

## Summary

Phase 01 fully achieved its goal. All 10 requirements (BM-01 through BM-07, TST-01 through TST-03) are implemented and tested. The five integration benchmark bugs are fixed with proper TDD discipline (RED commits before GREEN commits). The h5py loading path enables 2.5M-cell benchmarking without OOM. The targeted panel recipe now correctly preserves raw counts for scVI. No anti-patterns, no stubs, no orphaned code. 59 tests pass.

---

_Verified: 2026-03-21T14:00:00Z_
_Verifier: Claude (gsd-verifier)_
