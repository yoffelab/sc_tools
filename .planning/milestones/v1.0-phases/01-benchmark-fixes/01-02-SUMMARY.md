---
phase: 01-benchmark-fixes
plan: 02
subsystem: benchmarking
tags: [h5py, anndata, numpy, pandas, integration-metrics, tdd, recipes]

# Dependency graph
requires:
  - phase: 01-benchmark-fixes plan 01
    provides: "compare_integrations with subsample_n, seed, resolution, NaN masking, benchmark_params"
provides:
  - "_load_embedding_h5py helper for memory-efficient h5ad reads"
  - "embedding_files parameter for file-based benchmark loading"
  - "bio_key parameter for flexible bio conservation metrics"
  - "Fixed _recipe_targeted_panel with scVI branch before normalization"
affects: [cli-benchmark, ibd-spatial-project]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "h5py selective reads: obsm embeddings + obs categorical reconstruction from codes/categories"
    - "Minimal AnnData construction from h5py arrays for metric computation"
    - "scVI-first branching in preprocessing recipes (normalize after branch, not before)"

key-files:
  created: []
  modified:
    - sc_tools/bm/integration.py
    - sc_tools/pp/recipes.py
    - sc_tools/tests/test_integration_benchmark.py
    - sc_tools/tests/test_pp.py

key-decisions:
  - "embedding_files auto-discovers obsm key from h5ad file (reads first/only obsm key)"
  - "bio_key defaults to celltype_key for backwards compatibility -- both params coexist"
  - "compare_integrations adata param now optional (None when embedding_files provides all methods)"
  - "_recipe_targeted_panel scVI branch mirrors _recipe_visium pattern exactly"

patterns-established:
  - "_load_embedding_h5py pattern for memory-efficient h5ad access (300MB for 2.5M cells vs full AnnData)"
  - "Minimal AnnData from h5py arrays: zeros X + obsm + categorical obs for metric computation"
  - "Mock strategy.select_features/reduce_and_integrate/embed_and_cluster for recipe unit tests"

requirements-completed: [BM-01, BM-07]

# Metrics
duration: 16min
completed: 2026-03-21
---

# Phase 01 Plan 02: h5py Loading Path and Targeted Panel scVI Fix Summary

**h5py-based embedding loading for 2.5M-cell benchmarks without OOM, bio_key parameter for flexible bio conservation, and fixed _recipe_targeted_panel to preserve raw counts for scVI**

## Performance

- **Duration:** 16 min
- **Started:** 2026-03-21T13:23:15Z
- **Completed:** 2026-03-21T13:39:00Z
- **Tasks:** 2 (both TDD: red + green)
- **Files modified:** 4

## Accomplishments

- _load_embedding_h5py reads obsm embeddings and obs columns via h5py with categorical codes/categories reconstruction, enabling benchmarking 2.5M-cell datasets at ~300MB (vs full AnnData OOM)
- compare_integrations now accepts embedding_files dict for file-based loading and bio_key for any obs column as bio conservation metric target
- _recipe_targeted_panel correctly branches on scVI before normalization, preserving raw counts needed by the VAE model
- 14 new tests across 6 test classes, all 100 tests pass (14 new + 86 existing across both test files)

## Task Commits

Each task was committed atomically:

1. **Task 1 RED: TDD failing tests for h5py/embedding_files/bio_key** - `3452b56` (test)
2. **Task 1 GREEN: implement _load_embedding_h5py and updated compare_integrations** - `4668827` (feat)
3. **Task 2 RED: TDD failing tests for targeted panel scVI normalization** - `c3e6973` (test)
4. **Task 2 GREEN: fix _recipe_targeted_panel scVI branch** - `118cf3b` (fix)

## Files Created/Modified

- `sc_tools/bm/integration.py` - Added _load_embedding_h5py, embedding_files param, bio_key param, file-based loading in compare_integrations
- `sc_tools/pp/recipes.py` - Restructured _recipe_targeted_panel to branch on scVI before normalization
- `sc_tools/tests/test_integration_benchmark.py` - Added TestLoadEmbeddingH5py, TestEmbeddingFilesParameter, TestBioKeyParameter (11 tests) + _write_test_h5ad fixture
- `sc_tools/tests/test_pp.py` - Added TestTargetedPanelScVI (3 tests)

## Decisions Made

- embedding_files auto-discovers the obsm key from each h5ad file (reads the first available obsm key) rather than requiring caller to specify -- simpler API since benchmark h5ad files typically contain a single embedding
- bio_key is an additional parameter (not a replacement for celltype_key) -- defaults to celltype_key when not provided, maintaining full backwards compatibility
- compare_integrations adata parameter is now optional (can be None when all methods come from embedding_files) with input validation raising ValueError if neither is provided
- _recipe_targeted_panel fix mirrors _recipe_visium pattern exactly -- scVI branch calls select_features with flavor="seurat_v3" before any normalization

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Fixed test fixture directory creation**
- **Found during:** Task 1 GREEN (TestEmbeddingFilesParameter.test_from_files)
- **Issue:** `_write_test_h5ad(tmp_path / "a")` failed because subdirectory didn't exist
- **Fix:** Added `tmp_path.mkdir(parents=True, exist_ok=True)` in _write_test_h5ad fixture
- **Files modified:** sc_tools/tests/test_integration_benchmark.py
- **Committed in:** `4668827` (Task 1 GREEN commit)

**2. [Rule 3 - Blocking] Mocked select_features in targeted panel tests**
- **Found during:** Task 2 GREEN (TestTargetedPanelScVI)
- **Issue:** `select_features(flavor="seurat_v3")` requires scikit-misc which is not installed in test env
- **Fix:** Added `patch.object(strategy, "select_features")` to all 3 targeted panel tests
- **Files modified:** sc_tools/tests/test_pp.py
- **Committed in:** `118cf3b` (Task 2 GREEN commit)

---

**Total deviations:** 2 auto-fixed (2 blocking)
**Impact on plan:** Both fixes necessary for test execution. No scope creep.

## Issues Encountered

None beyond the auto-fixed items above.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All 10 BM/TST requirements for Phase 01 now complete (BM-01..07, TST-01..03)
- Phase 01 (benchmark-fixes) is fully complete
- IBD spatial project can use embedding_files for 2.5M-cell benchmark comparison
- CLI Phase 02 can wrap compare_integrations as `sct benchmark integration --from-dir`

## Self-Check: PASSED

- All 4 modified files exist on disk
- All 4 task commits found in git log
- All 13 acceptance criteria verified

---
*Phase: 01-benchmark-fixes*
*Completed: 2026-03-21*
