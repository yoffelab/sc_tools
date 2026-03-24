---
phase: 02-cli-foundation
plan: 02
subsystem: cli
tags: [typer, rich, cli-app, error-handler, mcp-result]

# Dependency graph
requires:
  - phase: 02-01
    provides: CLIResult Pydantic envelope, exception hierarchy with exit codes, sct entry point in pyproject.toml
provides:
  - Typer CLI app (sc_tools/cli.py) with global --human callback, error handler, five stub command groups
  - __main__.py thin Typer wrapper for python -m sc_tools
  - CLI argument parsing tests (24 tests, TST-04)
  - MCP validate_checkpoint returning CLIResult JSON (D-11 proof-of-concept)
affects: [03-core-commands, 04-self-discovery]

# Tech tracking
tech-stack:
  added: [typer CLI app framework]
  patterns: [cli_handler decorator for error-to-exit-code mapping, _emit/_render_rich output helpers, invoke_without_command callbacks for stub groups]

key-files:
  created:
    - sc_tools/cli.py
    - sc_tools/tests/test_cli.py
  modified:
    - sc_tools/__main__.py
    - sc_tools/mcp/tools_server.py
    - pyproject.toml

key-decisions:
  - "is_eager=True on --version callback to run before Typer command resolution"
  - "invoke_without_command=True on all sub-app callbacks so stub groups appear in help"
  - "pyproject.toml [project.scripts] moved after [project].dependencies to fix TOML parsing error"

patterns-established:
  - "cli_handler decorator: wraps command functions with error catching, CLIResult emission, and SystemExit mapping"
  - "_emit() writes JSON to stdout, _render_rich() writes Rich to stderr -- never mixed"
  - "Stub command groups use callback with invoke_without_command=True pattern"

requirements-completed: [CLI-01, CLI-02, CLI-04, CLI-07, CLI-08, TST-04]

# Metrics
duration: 5min
completed: 2026-03-21
---

# Phase 02 Plan 02: CLI App & Tests Summary

**Typer CLI app with global --human flag, error handler mapping exceptions to exit codes 0-3, five stub command groups, 24 argument parsing tests, and MCP CLIResult migration proof-of-concept**

## Performance

- **Duration:** 5 min
- **Started:** 2026-03-21T18:14:52Z
- **Completed:** 2026-03-21T18:19:52Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- Typer CLI app with global --human flag, --version, error handler decorator, and five stub command groups (qc, preprocess, integrate, benchmark, celltype)
- 24 CLI argument parsing tests covering help, version, --human, unknown commands, CLIResult model, error taxonomy, non-interactive constraint, lazy imports
- MCP validate_checkpoint returns CLIResult JSON instead of raw string (D-11 proof-of-concept)
- __main__.py reduced to 4-line thin wrapper

## Task Commits

Each task was committed atomically:

1. **Task 1: Create Typer CLI app with error handler, --human support, and stub command groups** - `6d208ae` (feat)
2. **Task 2: Write CLI argument parsing tests and update MCP validate_checkpoint to return CLIResult** - `2de2ace` (feat)

## Files Created/Modified
- `sc_tools/cli.py` - Typer app with global --human callback, error handler decorator, output helpers, five stub command groups
- `sc_tools/__main__.py` - Thin 4-line wrapper importing and running Typer app
- `sc_tools/tests/test_cli.py` - 24 CLI argument parsing tests (TST-04)
- `sc_tools/mcp/tools_server.py` - validate_checkpoint returns CLIResult JSON
- `pyproject.toml` - Fixed [project.scripts] placement (moved after dependencies to fix TOML parsing)

## Decisions Made
- Used `is_eager=True` on `--version` option callback so it executes before Typer's command resolution (avoids "Missing command" error when no subcommand provided)
- All five stub command groups use `invoke_without_command=True` callbacks to appear in help without requiring placeholder commands
- Fixed pyproject.toml `[project.scripts]` section placement: was between `[project]` and `dependencies`, causing setuptools to interpret `dependencies` as a scripts entry

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed --version flag conflict with Typer command resolution**
- **Found during:** Task 1
- **Issue:** `--version` in the main callback raised `typer.Exit()` but Typer still expected a command, resulting in exit code 2
- **Fix:** Used `is_eager=True` with a dedicated `_version_callback` function so version runs before command parsing
- **Files modified:** sc_tools/cli.py
- **Verification:** `runner.invoke(app, ['--version'])` returns exit code 0
- **Committed in:** 6d208ae

**2. [Rule 3 - Blocking] Fixed pyproject.toml [project.scripts] placement**
- **Found during:** Task 1 (pip install verification)
- **Issue:** `[project.scripts]` was placed between `[project]` header and `dependencies` key, making setuptools interpret `dependencies` as belonging to `project.scripts` rather than `project`
- **Fix:** Moved `[project.scripts]` after the `dependencies` array, before `[project.optional-dependencies]`
- **Files modified:** pyproject.toml
- **Verification:** TOML parsing succeeds (though pip install still fails due to pre-existing Python >=3.11 constraint vs 3.10 env)
- **Committed in:** 6d208ae

---

**Total deviations:** 2 auto-fixed (1 bug, 1 blocking)
**Impact on plan:** Both fixes necessary for CLI to function. No scope creep.

## Issues Encountered
- `pip install -e ".[cli]"` fails because `requires-python = ">=3.11"` but environment has Python 3.10.15. This is a pre-existing pyproject.toml constraint, not caused by this plan. All three CLI deps (typer, pydantic, rich) are already installed in the environment, so verification succeeded via direct Python imports.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- CLI framework complete: app object, error handler, output helpers, stub groups all in place
- Phase 3 (core commands) can add real commands to the stub groups by importing qc_app, preprocess_app, etc. from sc_tools.cli
- MCP CLIResult pattern proven on validate_checkpoint; remaining MCP tools can migrate in Phase 3
- 24 tests provide regression baseline for Phase 3 command additions

## Self-Check: PASSED

All files verified present. Both task commits found in git log.

---
*Phase: 02-cli-foundation*
*Completed: 2026-03-21*
