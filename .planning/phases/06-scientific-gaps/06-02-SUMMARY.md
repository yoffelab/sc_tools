---
phase: 06-scientific-gaps
plan: 02
subsystem: analysis
tags: [pydeseq2, pseudobulk, differential-expression, cli, typer]

# Dependency graph
requires:
  - phase: 06-01
    provides: "adata_multi_subject fixture, validate_subject_metadata, _check_deps pattern"
provides:
  - "sc_tools.tl.de module with aggregate_pseudobulk and run_pseudobulk_de"
  - "sct de run CLI command for pseudobulk DE"
  - "Auto-inferred design formula with collinearity guard"
affects: [08-multi-omic]

# Tech tracking
tech-stack:
  added: [pydeseq2]
  patterns: [pseudobulk-aggregation, collinearity-guard, per-celltype-iteration]

key-files:
  created:
    - sc_tools/tl/de.py
    - sc_tools/cli/de.py
    - sc_tools/tests/test_de.py
  modified:
    - sc_tools/cli/__init__.py

key-decisions:
  - "Raw count extraction priority: layers[layer] > layers['counts'] > raw.X > X"
  - "Collinear batch covariates (1:1 with subject_id) excluded from auto-formula via _is_collinear_with_subject check"
  - "Reference level defaults to alphabetically first (matching DESeq2 convention)"
  - "PyDESeq2 failures per-celltype logged and skipped (does not abort entire pipeline)"

patterns-established:
  - "Pseudobulk pattern: aggregate per-celltype subsets (sparse-safe), then run DE on small dense matrices"
  - "Formula inference with collinearity guard for batch covariates"

requirements-completed: [SCI-01]

# Metrics
duration: 5min
completed: 2026-03-24
---

# Phase 06 Plan 02: Pseudobulk DE Summary

**Pseudobulk DE module wrapping PyDESeq2 with per-celltype aggregation, collinearity-guarded design formula, and sct de run CLI command**

## Performance

- **Duration:** 5 min
- **Started:** 2026-03-24T14:31:36Z
- **Completed:** 2026-03-24T14:37:18Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- aggregate_pseudobulk correctly sums raw counts by subject+celltype with min_cells threshold filtering
- run_pseudobulk_de runs PyDESeq2 per celltype with auto-inferred or custom design formula
- Collinear batch covariates (1:1 with subject_id) automatically excluded from auto-formula to prevent rank-deficient design matrix
- sct de run CLI command registered, writes per-celltype CSVs to results/de/ with standardized columns
- 20 tests total: 15 pass, 5 skip gracefully (pydeseq2 not installed)

## Task Commits

Each task was committed atomically:

1. **Task 1: Pseudobulk DE module (TDD)** - `d21a335` (test: failing tests) -> `8991d70` (feat: implementation)
2. **Task 2: CLI command sct de run** - `6169514` (feat: CLI + registration + CLI tests)

## Files Created/Modified
- `sc_tools/tl/de.py` - Pseudobulk DE module: aggregate_pseudobulk, run_pseudobulk_de, _infer_design_formula
- `sc_tools/cli/de.py` - CLI command: sct de run with all flags (condition, subject-col, formula, etc.)
- `sc_tools/cli/__init__.py` - Registered de_app, added pydeseq2 to _DEP_INSTALL
- `sc_tools/tests/test_de.py` - 20 tests: aggregation (7), formula inference (4), PyDESeq2 integration (4), CLI (5)

## Decisions Made
- Raw count extraction uses priority chain: explicit layer > layers['counts'] > raw.X > X
- Collinear batch detection checks 1:1 mapping between batch column values and subject index
- Reference level defaults to alphabetically first condition level (DESeq2 convention)
- PyDESeq2 failures per-celltype are logged and skipped, not fatal to the pipeline
- Batch candidate columns limited to ('batch', 'library_id') for auto-detection

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - pydeseq2 is an optional dependency handled by _check_deps fast-fail pattern.

## Next Phase Readiness
- DE module ready for use in downstream analysis workflows
- Per-celltype CSV output format established for future volcano plot or enrichment analysis
- Existing tests unbroken (118 passed in test_qc, test_subject_metadata, test_panel_dispatch)

## Self-Check: PASSED

---
*Phase: 06-scientific-gaps*
*Completed: 2026-03-24*
