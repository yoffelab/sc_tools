---
phase: 03-core-commands
plan: 01
subsystem: cli
tags: [typer, cli-package, dependency-check, pytest-fixtures]

# Dependency graph
requires:
  - phase: 02-cli-foundation
    provides: "Typer app, CLIResult model, error taxonomy, cli.py monolith"
provides:
  - "cli/ package structure with submodule stubs (qc, preprocess, validate, benchmark, status)"
  - "_check_deps() utility for fast-fail on missing optional dependencies"
  - "adata_100 shared test fixture for CLI integration tests"
  - "validate_app and status_app registered as subcommands"
affects: [03-02, 03-03, 03-04, 04-self-discovery]

# Tech tracking
tech-stack:
  added: []
  patterns: ["cli/ package with lazy submodule imports", "conftest.py shared fixtures for CLI tests"]

key-files:
  created:
    - sc_tools/cli/__init__.py
    - sc_tools/cli/qc.py
    - sc_tools/cli/preprocess.py
    - sc_tools/cli/validate.py
    - sc_tools/cli/benchmark.py
    - sc_tools/cli/status.py
    - sc_tools/tests/conftest.py
  modified: []

key-decisions:
  - "integrate_app and celltype_app kept in __init__.py (no Phase 3 commands for them)"
  - "Submodule imports at bottom of __init__.py with noqa E402 for clean separation"

patterns-established:
  - "CLI submodule pattern: each submodule defines its own Typer app with callback"
  - "_check_deps() pattern: import-test with install hint lookup table"
  - "adata_100 fixture: 100 cells, 200 genes, batch/celltype/spatial/embeddings/MT-genes"

requirements-completed: [CMD-08, TST-05]

# Metrics
duration: 2min
completed: 2026-03-23
---

# Phase 03 Plan 01: CLI Package Migration Summary

**cli.py migrated to cli/ package with _check_deps fast-fail utility and adata_100 shared test fixture**

## Performance

- **Duration:** 2 min
- **Started:** 2026-03-23T02:53:49Z
- **Completed:** 2026-03-23T02:56:13Z
- **Tasks:** 2
- **Files modified:** 7

## Accomplishments
- Migrated monolithic cli.py to cli/ package with 5 stub submodules (D-01)
- Added _check_deps() utility for CMD-08 fast-fail on missing optional dependencies
- Registered validate_app and status_app as new subcommands (BLOCKER 1 fix)
- Created conftest.py with adata_100 fixture for CLI integration tests (TST-05)
- All 24 existing Phase 2 tests pass with backward compatibility preserved (D-02)

## Task Commits

Each task was committed atomically:

1. **Task 1: Migrate cli.py to cli/ package with _check_deps** - `3754eb3` (feat)
2. **Task 2: Create conftest.py with shared test fixtures** - `5dcf211` (feat)

## Files Created/Modified
- `sc_tools/cli/__init__.py` - Main CLI app, helpers, _check_deps, subcommand registration
- `sc_tools/cli/qc.py` - QC sub-app stub
- `sc_tools/cli/preprocess.py` - Preprocess sub-app stub
- `sc_tools/cli/validate.py` - Validate sub-app stub
- `sc_tools/cli/benchmark.py` - Benchmark sub-app stub
- `sc_tools/cli/status.py` - Status sub-app stub
- `sc_tools/tests/conftest.py` - Shared test fixtures (adata_100, adata_100_h5ad, adata_100_preprocess_checkpoint)

## Decisions Made
- integrate_app and celltype_app kept in __init__.py since no Phase 3 commands target them
- Submodule imports placed at bottom of __init__.py with noqa E402 for clean separation of helpers from registrations

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- cli/ package structure ready for Plans 02-04 to add commands to submodules
- _check_deps() available for all command implementations needing optional deps
- adata_100 fixture available for integration tests in Plans 02-04
- validate and status subcommands registered and ready for command implementation

---
*Phase: 03-core-commands*
*Completed: 2026-03-23*
