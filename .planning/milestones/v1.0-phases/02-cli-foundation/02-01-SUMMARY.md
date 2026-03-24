---
phase: 02-cli-foundation
plan: 01
subsystem: cli
tags: [pydantic, exception-hierarchy, cli-result, error-taxonomy]

# Dependency graph
requires: []
provides:
  - CLIResult Pydantic envelope with dual serialization (CLI JSON + MCP dict)
  - Exception hierarchy with category/exit_code metadata (SCToolsError, SCToolsUserError, SCToolsDataError, SCToolsRuntimeError, SCToolsFatalError)
  - ErrorInfo model with retryable/fixable/fatal taxonomy
  - Provenance model with auto-detected version
  - sct entry point registered in pyproject.toml
  - cli optional dependency group (typer, rich, pydantic)
affects: [02-02, 03-core-commands, 04-self-discovery, 05-provenance]

# Tech tracking
tech-stack:
  added: [pydantic v2 CLIResult model]
  patterns: [Pydantic envelope for structured CLI output, exception classes with taxonomy metadata]

key-files:
  created:
    - sc_tools/errors.py
    - sc_tools/models/__init__.py
    - sc_tools/models/result.py
  modified:
    - pyproject.toml

key-decisions:
  - "_get_version() uses importlib.metadata.version('sci-sc-tools') with fallback to 'unknown'"
  - "ConfigDict(use_enum_values=True) on CLIResult for JSON string serialization of Status enum"
  - "timezone.utc with noqa: UP017 per project convention (Python 3.10 compat)"

patterns-established:
  - "Exception hierarchy: SCToolsError base -> User/Data/Runtime/Fatal subclasses with category + exit_code"
  - "CLIResult dual serialization: model_dump_json() for CLI stdout, model_dump(mode='json') for MCP returns"

requirements-completed: [CLI-03, CLI-05, CLI-06]

# Metrics
duration: 2min
completed: 2026-03-21
---

# Phase 02 Plan 01: Type Contracts Summary

**CLIResult Pydantic envelope with Status/ErrorInfo/Provenance models and 5-class exception hierarchy mapping to semantic exit codes 0/1/2/3**

## Performance

- **Duration:** 2 min
- **Started:** 2026-03-21T18:10:11Z
- **Completed:** 2026-03-21T18:12:34Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- CLIResult Pydantic v2 model with all D-06 fields (status, command, data, artifacts, provenance, message, error, partial_failures)
- Exception hierarchy: SCToolsError base with SCToolsUserError(exit 1), SCToolsDataError(exit 2), SCToolsRuntimeError(exit 3), SCToolsFatalError(exit 3)
- ErrorInfo model with retryable/fixable/fatal taxonomy and actionable suggestion text
- Provenance with auto-detected sc_tools version via importlib.metadata
- pyproject.toml: sct entry point + cli optional dependency group (typer, rich, pydantic)

## Task Commits

Each task was committed atomically:

1. **Task 1: Create exception hierarchy and CLIResult model** - `ea08dff` (feat)
2. **Task 2: Add CLI optional dependency group to pyproject.toml** - `7d59c7a` (chore)

## Files Created/Modified
- `sc_tools/errors.py` - Exception hierarchy with category/suggestion/exit_code class attributes
- `sc_tools/models/result.py` - CLIResult, Status, ErrorInfo, Provenance Pydantic models
- `sc_tools/models/__init__.py` - Package re-exports for clean imports
- `pyproject.toml` - Added [project.scripts] sct entry point and cli optional deps

## Decisions Made
- `_get_version()` helper defined before Provenance class, uses try/except for PackageNotFoundError
- `ConfigDict(use_enum_values=True)` on CLIResult ensures JSON serialization uses string enum values
- `timezone.utc` with `# noqa: UP017` per project convention for Python 3.10 compatibility

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Type contracts ready for Plan 02 (cli.py with Typer app, __main__.py wrapper, tests)
- sc_tools/errors.py and sc_tools/models/result.py are importable without heavy dependencies
- ErrorInfo.category aligns with SCToolsError.category values for taxonomy consistency

## Self-Check: PASSED

All 5 files verified present. Both task commits (ea08dff, 7d59c7a) found in git log.

---
*Phase: 02-cli-foundation*
*Completed: 2026-03-21*
