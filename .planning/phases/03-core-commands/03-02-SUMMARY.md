---
phase: 03-core-commands
plan: 02
subsystem: cli
tags: [typer, cli, qc, validation, pipeline-status, reports]

requires:
  - phase: 03-01
    provides: "CLI app skeleton with cli_handler, _check_deps, 7 subcommand groups"
  - phase: 02-cli-foundation
    provides: "CLIResult model, error taxonomy, Typer app framework"
provides:
  - "sct validate run <phase> <file> -- checkpoint validation with CLIResult JSON"
  - "sct status show -- pipeline DAG with graceful registry fallback (D-09)"
  - "sct qc run <file> -- QC metrics with per-sample classification"
  - "sct report generate <type> -- HTML report generation (pre_filter, post_filter, post_integration, post_celltyping)"
affects: [04-self-discovery, 05-provenance]

tech-stack:
  added: []
  patterns: [cli-handler-decorator, lazy-imports, graceful-registry-fallback]

key-files:
  created: []
  modified:
    - sc_tools/cli/validate.py
    - sc_tools/cli/status.py
    - sc_tools/cli/qc.py
    - sc_tools/cli/__init__.py

key-decisions:
  - "validate_run does NOT use @cli_handler -- needs custom exit code logic (exit 2 for data error)"
  - "status_show uses -> None return type annotation to avoid Typer forward-reference resolution error with CLIResult"
  - "report_app defined in qc.py but registered as separate 'report' subcommand in __init__.py"

patterns-established:
  - "Custom exit codes: commands needing non-0/1 exits skip @cli_handler and call _emit + SystemExit directly"
  - "Registry fallback: D-09 pattern wraps db import in try/except, shows DAG structure even without registry"

requirements-completed: [CMD-01, CMD-03, CMD-05, CMD-06]

duration: 3min
completed: 2026-03-23
---

# Phase 03 Plan 02: Validate, Status, QC, and Report CLI Commands Summary

**Four CLI commands wrapping existing backend: validate (checkpoint checks), status (DAG viewer with registry fallback), qc run (per-sample metrics), and report generate (HTML report dispatch)**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-23T02:58:17Z
- **Completed:** 2026-03-23T03:00:57Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- `sct validate run` wraps validate_file with pass/fail JSON output and exit code 2 on validation failure
- `sct status show` displays full pipeline DAG (10 phases) with graceful registry fallback per D-09
- `sct qc run` computes per-sample QC metrics via compute_sample_metrics + classify_samples
- `sct report generate` dispatches to 4 report generators (pre_filter, post_filter, post_integration, post_celltyping)
- All commands support --project-dir (D-03) and use lazy imports (CLI-08)

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement sct validate and sct status commands** - `f1c5b0d` (feat)
2. **Task 2: Implement sct qc run and sct report commands** - `ed4711e` (feat)

## Files Created/Modified
- `sc_tools/cli/validate.py` - Checkpoint validation command with phase slug dispatch
- `sc_tools/cli/status.py` - Pipeline DAG status with registry graceful fallback
- `sc_tools/cli/qc.py` - QC metrics and HTML report generation commands
- `sc_tools/cli/__init__.py` - Registered report_app as 'report' subcommand

## Decisions Made
- validate_run skips @cli_handler to use custom exit code 2 for data errors (validate failure is a data error, not a user error)
- Return type annotation uses `-> None` instead of `-> "CLIResult"` to avoid Typer's forward-reference resolution attempting to evaluate the string annotation
- report_app colocated in qc.py (report commands share QC infrastructure) but registered as separate top-level 'report' subcommand

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed Typer forward-reference resolution error on CLIResult return type**
- **Found during:** Task 1 (status_show)
- **Issue:** `-> "CLIResult"` string annotation caused NameError when Typer tried to resolve type hints via `get_type_hints()`
- **Fix:** Changed return type to `-> None` since @cli_handler wraps the function and Typer never sees the actual return value
- **Files modified:** sc_tools/cli/status.py
- **Verification:** `sct status --help` exits 0, `sct status show` returns valid JSON with 10 phases
- **Committed in:** f1c5b0d (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Minor type annotation fix for Typer compatibility. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All four commands registered and reachable via --help
- 24 existing CLI tests still pass
- Ready for Plan 03 (preprocess, benchmark) and Plan 04 (integration tests)

---
*Phase: 03-core-commands*
*Completed: 2026-03-23*
