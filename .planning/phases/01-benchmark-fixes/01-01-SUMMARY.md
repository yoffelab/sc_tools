---
phase: 01-benchmark-fixes
plan: 01
subsystem: benchmarking
tags: [numpy, pandas, anndata, sklearn, integration-metrics, tdd]

# Dependency graph
requires: []
provides:
  - "_mask_nan_rows helper for per-embedding NaN filtering"
  - "Fixed _stratified_subsample with proportional allocation"
  - "subsample_n, seed, resolution params in compare_integrations"
  - "benchmark_params provenance in DataFrame.attrs"
  - "runtime_s column in run_integration_benchmark output"
affects: [01-benchmark-fixes, cli-benchmark]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Per-embedding NaN masking before metric computation"
    - "DataFrame.attrs for benchmark parameter provenance"
    - "Proportional stratified subsampling with rng.choice trim"

key-files:
  created: []
  modified:
    - sc_tools/bm/integration.py
    - sc_tools/tests/test_integration_benchmark.py

key-decisions:
  - "NaN masking applied per-embedding inside compare_integrations loop (not globally), so each method's valid cell set is independent"
  - "Subsample trim uses rng.choice(indices, n, replace=False) to avoid index-position bias"
  - "resolution threaded through compute_integration_metrics to _ari_sklearn/_nmi_sklearn for configurable Leiden clustering"

patterns-established:
  - "TDD red-green for sc_tools bug fixes: write failing tests first, then implement"
  - "_mask_nan_rows(X, *arrays) pattern for aligned NaN filtering across arrays"

requirements-completed: [BM-02, BM-03, BM-04, BM-05, BM-06, TST-01, TST-02, TST-03]

# Metrics
duration: 13min
completed: 2026-03-21
---

# Phase 01 Plan 01: Integration Benchmark Fixes Summary

**Fixed five integration benchmark bugs (NaN handling, subsampling bias, subsample_n, runtime tracking, param provenance) with 15 new TDD tests covering edge cases**

## Performance

- **Duration:** 13 min
- **Started:** 2026-03-21T13:07:11Z
- **Completed:** 2026-03-21T13:21:00Z
- **Tasks:** 1 (TDD: red + green)
- **Files modified:** 2

## Accomplishments

- NaN embedding rows (resolVI) filtered per-embedding before metric computation, preventing NaN propagation in benchmark scores
- _stratified_subsample truncation bias eliminated: proportional group allocation with random trim instead of sorted()[:n]
- compare_integrations now accepts subsample_n, seed, resolution parameters with full provenance stored in DataFrame.attrs
- run_integration_benchmark output includes runtime_s column per method via time.perf_counter()
- 15 new tests across 3 test classes (TestStratifiedSubsampleFix, TestNaNHandling, TestCompareIntegrationsExtended), all 45 tests pass

## Task Commits

Each task was committed atomically:

1. **Task 1 RED: TDD failing tests** - `ec905c4` (test)
2. **Task 1 GREEN: Implementation** - `773bd96` (feat)

## Files Created/Modified

- `sc_tools/bm/integration.py` - Added _mask_nan_rows, fixed _stratified_subsample, wired subsample_n/seed/resolution, benchmark_params attrs, runtime_s tracking
- `sc_tools/tests/test_integration_benchmark.py` - Added 3 test classes (15 tests) and _make_nan_embedding_adata fixture

## Decisions Made

- NaN masking applied per-embedding inside compare_integrations loop rather than globally, so each method's valid cell set is independent
- Subsample trim uses rng.choice to randomly select n indices rather than sorted truncation
- resolution parameter threaded through compute_integration_metrics to _ari_sklearn and _nmi_sklearn for configurable Leiden clustering resolution
- Runtime recorded via finally block in run_integration_benchmark to capture timing even on failure

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Integration benchmark fixes (BM-02..06) complete and tested
- Plan 01-02 (BM-01 h5py loading, BM-07 recipe fix) can proceed independently
- IBD spatial project can use fixed compare_integrations with subsample_n for large datasets

---
*Phase: 01-benchmark-fixes*
*Completed: 2026-03-21*
