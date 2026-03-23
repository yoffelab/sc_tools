---
phase: 03-core-commands
plan: 04
subsystem: testing
tags: [pytest, typer, cli-runner, integration-tests, e2e, mcp, cli-result]

# Dependency graph
requires:
  - phase: 03-02
    provides: "CLI commands (qc, validate, status, preprocess, benchmark)"
  - phase: 03-03
    provides: "benchmark_params support in report generation"
provides:
  - "23 CLI integration tests covering all Phase 3 commands"
  - "E2E test scaffold with skipif guards for HPC real data"
  - "CMD-07 verification: CLI and MCP share identical CLIResult fields"
affects: [04-self-discovery, 05-provenance]

# Tech tracking
tech-stack:
  added: []
  patterns: ["CliRunner-based integration testing", "skipif-guarded E2E tests via env vars"]

key-files:
  created:
    - sc_tools/tests/test_cli_commands.py
    - sc_tools/tests/test_cli_e2e.py
  modified:
    - sc_tools/cli/qc.py

key-decisions:
  - "CliRunner() without mix_stderr (not supported in this Typer version)"
  - "E2E tests guarded by SCT_TEST_DATA_DIR and SCT_TEST_EMBEDDING_DIR env vars"
  - "qc_run uses sc_tools.qc.metrics.calculate_qc_metrics instead of raw scanpy call"

patterns-established:
  - "Integration test pattern: CliRunner + _parse_result + fixture-based h5ad files"
  - "E2E test pattern: skipif-guarded classes with env var data paths"

requirements-completed: [CMD-07, TST-05, TST-06]

# Metrics
duration: 4min
completed: 2026-03-23
---

# Phase 3 Plan 4: CLI Tests & CMD-07 Verification Summary

**23 integration tests covering all CLI commands with CLIResult JSON verification, E2E scaffold with skipif guards, and verified CLI/MCP shared Result type**

## Performance

- **Duration:** 4 min
- **Started:** 2026-03-23T03:03:44Z
- **Completed:** 2026-03-23T03:08:33Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- 23 integration tests across 8 test classes covering CMD-01 through CMD-08
- CMD-07 verified: CLI and MCP tools produce CLIResult with identical field sets
- E2E test scaffold with 3 skipif-guarded tests for HPC real data
- Full test suite: 47 passed, 3 skipped (24 Phase 2 + 23 Phase 3 + 3 E2E skipped)

## Task Commits

Each task was committed atomically:

1. **Task 1: CLI integration tests for all commands** - `0cdb1d6` (test)
2. **Task 2: E2E test scaffold and full suite verification** - `3bcbe87` (test)

## Files Created/Modified
- `sc_tools/tests/test_cli_commands.py` - 23 integration tests for all CLI commands
- `sc_tools/tests/test_cli_e2e.py` - E2E test scaffold with skipif guards
- `sc_tools/cli/qc.py` - Fixed qc_run to use sc_tools.qc.metrics.calculate_qc_metrics

## Decisions Made
- CliRunner() without mix_stderr: Typer's CliRunner does not support this parameter in current version
- E2E tests use SCT_TEST_DATA_DIR and SCT_TEST_EMBEDDING_DIR environment variables for data paths
- Fixed qc_run to use the project's own calculate_qc_metrics which properly creates mt/hb boolean var columns before calling scanpy

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed qc_run calling scanpy with invalid qc_vars**
- **Found during:** Task 1 (CLI integration tests)
- **Issue:** qc_run passed qc_vars=["MT-", "mt-"] to sc.pp.calculate_qc_metrics, but scanpy expects boolean column names in adata.var, not prefix strings
- **Fix:** Replaced with sc_tools.qc.metrics.calculate_qc_metrics which properly creates boolean var["mt"] and var["hb"] columns first
- **Files modified:** sc_tools/cli/qc.py
- **Verification:** test_qc_run_produces_metrics passes with exit code 0
- **Committed in:** 0cdb1d6 (Task 1 commit)

**2. [Rule 1 - Bug] Removed unsupported mix_stderr parameter**
- **Found during:** Task 1 (CLI integration tests)
- **Issue:** Plan specified CliRunner(mix_stderr=False) but Typer's CliRunner does not accept this parameter
- **Fix:** Changed to CliRunner() without parameters
- **Files modified:** sc_tools/tests/test_cli_commands.py
- **Verification:** All 23 tests pass
- **Committed in:** 0cdb1d6 (Task 1 commit)

---

**Total deviations:** 2 auto-fixed (2 bugs)
**Impact on plan:** Both fixes necessary for test execution. No scope creep.

## Issues Encountered
None beyond the auto-fixed deviations above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 03 (core-commands) fully complete: all 4 plans executed
- Full CLI test coverage established for all commands
- Ready for Phase 04 (self-discovery) with stable, tested CLI interface

---
*Phase: 03-core-commands*
*Completed: 2026-03-23*
