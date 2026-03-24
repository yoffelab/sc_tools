---
phase: 07-memory-safety
plan: 02
subsystem: cli
tags: [typer, dry-run, memory-estimation, tier-dispatch]

# Dependency graph
requires:
  - phase: 07-01
    provides: IOGateway, DataTier, estimate_from_h5, METHOD_MULTIPLIERS
provides:
  - sct estimate command for pre-execution memory/runtime projection
  - --dry-run flag on all data-touching commands (validates without executing)
  - --force flag on all data-touching commands (future memory guard override)
  - Tier-aware @cli_handler decorator with DataTier metadata
affects: [08-multi-omic]

# Tech tracking
tech-stack:
  added: []
  patterns: [cli_handler(tier=DataTier) decorator pattern, dry-run early-exit with Status.skipped]

key-files:
  created:
    - sc_tools/cli/estimate.py
    - sc_tools/tests/test_cli_estimate.py
    - sc_tools/tests/test_cli_dryrun.py
  modified:
    - sc_tools/cli/__init__.py
    - sc_tools/cli/qc.py
    - sc_tools/cli/preprocess.py
    - sc_tools/cli/benchmark.py
    - sc_tools/cli/de.py
    - sc_tools/cli/validate.py

key-decisions:
  - "cli_handler extended to accept tier kwarg via decorator pattern: @cli_handler(tier=DataTier.T3_FULL)"
  - "Dry-run error handling moved inside try block for proper SCToolsUserError capture"
  - "validate_run dry-run handled inline (not via cli_handler) to preserve custom exit code 2 behavior"

patterns-established:
  - "cli_handler(tier=DataTier) pattern: all data-touching commands declare their loading tier"
  - "Dry-run pattern: --dry-run pops kwargs, validates input, returns Status.skipped with estimate"

requirements-completed: [MEM-02, MEM-03]

# Metrics
duration: 6min
completed: 2026-03-24
---

# Phase 7 Plan 2: CLI Estimate and Dry-Run Summary

**sct estimate command with memory/runtime projection and --dry-run/--force flags on all data-touching CLI commands**

## Performance

- **Duration:** 6 min
- **Started:** 2026-03-24T21:34:34Z
- **Completed:** 2026-03-24T21:40:40Z
- **Tasks:** 2
- **Files modified:** 9

## Accomplishments
- Extended cli_handler to support tier keyword argument and dry_run/force kwargs interception
- Created sct estimate command that projects peak memory and runtime from h5 metadata
- Added --dry-run and --force flags to qc run, preprocess run, benchmark integration, de run, and validate run
- All 11 new tests passing, all 78 existing CLI tests unaffected (backwards compatible)

## Task Commits

Each task was committed atomically:

1. **Task 1: Extend cli_handler and create sct estimate command**
   - `8020899` (test): add failing tests for sct estimate CLI command
   - `aee5271` (feat): implement sct estimate command and extend cli_handler with tier/dry-run/force
2. **Task 2: Add --dry-run and --force flags to all data-touching commands**
   - `6a2a11f` (test): add failing tests for --dry-run on data-touching CLI commands
   - `c9f6ed5` (feat): add --dry-run and --force flags to all data-touching CLI commands

## Files Created/Modified
- `sc_tools/cli/__init__.py` - Extended cli_handler with tier/dry-run/force support
- `sc_tools/cli/estimate.py` - New sct estimate command with register_estimate pattern
- `sc_tools/cli/qc.py` - Added --dry-run, --force, tier=T3_FULL
- `sc_tools/cli/preprocess.py` - Added --dry-run, --force, tier=T3_FULL
- `sc_tools/cli/benchmark.py` - Added --dry-run, --force, tier=T2_SUMMARY
- `sc_tools/cli/de.py` - Added --dry-run, --force, tier=T3_FULL
- `sc_tools/cli/validate.py` - Added --dry-run with inline handling (no --force needed for T1)
- `sc_tools/tests/test_cli_estimate.py` - 5 tests for estimate command
- `sc_tools/tests/test_cli_dryrun.py` - 6 tests for dry-run behavior

## Decisions Made
- cli_handler extended with decorator pattern supporting both `@cli_handler` and `@cli_handler(tier=...)` for backwards compatibility
- Dry-run error handling moved inside the main try block so SCToolsUserError on missing files produces proper error JSON (Rule 1 bug fix)
- validate_run keeps inline dry-run handling to preserve its custom exit code 2 behavior (Phase 3 decision D-03)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Dry-run error path not caught by error handler**
- **Found during:** Task 2 (test_dryrun_missing_file)
- **Issue:** SCToolsUserError raised in dry-run block was outside the try/except, producing empty output instead of error JSON
- **Fix:** Moved dry-run logic inside the main try block of cli_handler
- **Files modified:** sc_tools/cli/__init__.py
- **Verification:** test_dryrun_missing_file passes with proper error JSON
- **Committed in:** c9f6ed5

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Essential fix for correctness. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 7 (memory-safety) complete: IO Gateway + tiered loading (Plan 01) + CLI estimate/dry-run (Plan 02)
- Ready for Phase 8 (multi-omic assembly)
- All commands now declare their data tier and support --dry-run for safe pre-execution validation

---
*Phase: 07-memory-safety*
*Completed: 2026-03-24*

## Self-Check: PASSED
