---
phase: 03-core-commands
plan: 03
subsystem: cli
tags: [typer, cli, preprocessing, benchmarking, h5py, integration]

requires:
  - phase: 03-core-commands/01
    provides: "CLI framework with cli_handler, _check_deps, app with subcommand groups"
  - phase: 01-benchmark-fixes
    provides: "compare_integrations with embedding_files, subsample_n, seed, resolution params"
provides:
  - "sct preprocess run command wrapping recipes.preprocess with modality auto-detection"
  - "sct benchmark integration command wrapping compare_integrations with --from-dir"
  - "benchmark --report flag generating HTML with subsampling params"
  - "benchmark_params kwarg on generate_benchmark_report"
affects: [04-self-discovery, 05-provenance]

tech-stack:
  added: []
  patterns: ["CLI commands return CLIResult, no return type annotation to avoid Typer eval issue"]

key-files:
  created: []
  modified:
    - sc_tools/cli/preprocess.py
    - sc_tools/cli/benchmark.py
    - sc_tools/bm/report.py

key-decisions:
  - "Removed -> CLIResult return annotation to avoid Typer get_type_hints NameError on Python 3.10"
  - "benchmark_params rendered as HTML table in report header via context dict"
  - "CLI calls generate_benchmark_report with output_path kwarg to avoid positional arg conflict with aggregated"

patterns-established:
  - "CLI command functions omit return type annotation when wrapped by @cli_handler (Typer compat)"

requirements-completed: [CMD-02, CMD-04]

duration: 3min
completed: 2026-03-23
---

# Phase 03 Plan 03: Preprocess Run and Benchmark Integration CLI Commands Summary

**sct preprocess run with modality auto-detection (D-08) and sct benchmark integration with --from-dir embedding discovery and --report HTML generation (D-12)**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-23T02:58:26Z
- **Completed:** 2026-03-23T03:01:46Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Implemented `sct preprocess run` wrapping `recipes.preprocess` with auto-modality detection from adata.uns keys or panel size heuristic (D-08)
- Implemented `sct benchmark integration --from-dir` wrapping `compare_integrations` with h5ad file auto-discovery and subsampling transparency (D-07, D-12)
- Added `--report` flag generating HTML benchmark report with subsampling params rendered in header (D-12 item 3)
- Added `benchmark_params` keyword argument to `generate_benchmark_report` in bm/report.py

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement sct preprocess run command** - `49e0ec9` (feat)
2. **Task 2: Implement sct benchmark integration command with --report flag** - `04bdef7` (feat)

## Files Created/Modified
- `sc_tools/cli/preprocess.py` - Full preprocess run command with _detect_modality, modality/integration validation, _check_deps
- `sc_tools/cli/benchmark.py` - Full benchmark integration command with h5ad discovery, subsampling, --report flag
- `sc_tools/bm/report.py` - Added benchmark_params kwarg and HTML rendering to generate_benchmark_report

## Decisions Made
- Removed `-> "CLIResult"` return type annotation from CLI command functions because Typer's `get_type_hints()` evaluates the forward reference in module scope where CLIResult is not imported (Python 3.10 compat issue). Used `# returns CLIResult` comment instead.
- `benchmark_params` rendered as a Bootstrap-styled HTML card with parameter table, injected into the report template context
- CLI passes `output_path=` as keyword arg to `generate_benchmark_report` to avoid positional conflict with the `aggregated` parameter

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Removed CLIResult return type annotation**
- **Found during:** Task 1 (preprocess_run)
- **Issue:** `-> "CLIResult"` return type annotation caused `NameError: name 'CLIResult' is not defined` when Typer evaluated type hints via `get_type_hints()` in Python 3.10 -- CLIResult was only imported inside the function body
- **Fix:** Removed return type annotation, added `# returns CLIResult` comment
- **Files modified:** sc_tools/cli/preprocess.py, sc_tools/cli/benchmark.py
- **Verification:** `sct preprocess --help` and `sct benchmark --help` both exit 0, all 24 tests pass
- **Committed in:** 49e0ec9 (Task 1), 04bdef7 (Task 2)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Necessary fix for Python 3.10 runtime compatibility. No scope creep.

## Issues Encountered
None beyond the return type annotation fix documented above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- preprocess run and benchmark integration commands complete and registered
- Ready for Plan 04 (sct validate and sct status commands)
- All existing tests (24) pass

---
*Phase: 03-core-commands*
*Completed: 2026-03-23*
