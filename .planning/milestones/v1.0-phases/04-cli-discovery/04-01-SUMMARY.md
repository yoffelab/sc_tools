---
phase: 04-cli-discovery
plan: 01
subsystem: cli
tags: [typer, click, introspection, json-schema, pydantic, discovery]

# Dependency graph
requires:
  - phase: 02-cli-foundation
    provides: "Typer app, CLIResult envelope, cli_handler decorator, error hierarchy"
  - phase: 03-core-commands
    provides: "6 leaf commands (qc run, preprocess run, etc.) to introspect"
provides:
  - "sct list-commands --json: machine-readable catalog of all CLI commands with params"
  - "sct describe <cmd>: JSON schema for specific command params + CLIResult output schema"
  - "sct schema: full CLI contract as single JSON document with $defs"
  - "_walk_commands helper for Typer/Click tree introspection"
  - "_param_to_schema helper for Click param to JSON Schema conversion"
affects: [05-provenance, 06-scientific-gaps, 07-memory-safety, 08-multi-omic]

# Tech tracking
tech-stack:
  added: []
  patterns: [typer-click-introspection, register-pattern-for-circular-import-avoidance]

key-files:
  created:
    - sc_tools/cli/discovery.py
    - sc_tools/tests/test_cli_discovery.py
  modified:
    - sc_tools/cli/__init__.py

key-decisions:
  - "register_discovery(app) pattern avoids circular import (discovery.py never imports app at module level)"
  - "Commands keyed by space-separated names (e.g. 'qc run') matching CLI invocation syntax"
  - "_CLICK_TYPE_MAP handles both lowercase and uppercase Click type names for robustness"

patterns-established:
  - "Discovery registration: register_discovery(app) called at bottom of cli/__init__.py"
  - "Param schema: Click params converted to JSON Schema objects with type, default, description, enum"

requirements-completed: [DSC-01, DSC-02, DSC-03]

# Metrics
duration: 2min
completed: 2026-03-23
---

# Phase 4 Plan 1: CLI Discovery Summary

**Three discovery commands (list-commands, describe, schema) providing machine-readable CLI introspection via Typer/Click tree walking and Pydantic model_json_schema**

## Performance

- **Duration:** 2 min
- **Started:** 2026-03-23T05:59:32Z
- **Completed:** 2026-03-23T06:01:27Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- `sct list-commands --json` returns catalog of 6 leaf commands with full parameter schemas (types, defaults, required, descriptions)
- `sct describe "qc run"` returns per-command JSON schema including CLIResult output schema with $defs
- `sct schema` returns full CLI contract as single JSON document with commands dict, $defs (CLIResult, Status, ErrorInfo, Provenance), schema_version
- All 16 discovery tests pass; all 24 existing CLI tests unbroken

## Task Commits

Each task was committed atomically:

1. **Task 1: Write tests for discovery commands (RED phase)** - `48d0683` (test)
2. **Task 2: Implement discovery commands and register (GREEN phase)** - `4e2ac1c` (feat)

## Files Created/Modified
- `sc_tools/cli/discovery.py` - Discovery commands: list-commands, describe, schema + _walk_commands, _param_to_schema, _command_to_entry helpers
- `sc_tools/tests/test_cli_discovery.py` - 16 tests across 5 test classes (TestListCommands, TestDescribe, TestSchema, TestOutputFormat, TestLazyImports)
- `sc_tools/cli/__init__.py` - Added register_discovery(app) call at bottom

## Decisions Made
- Used `register_discovery(app)` pattern to avoid circular imports (discovery.py imports app lazily inside function bodies)
- Commands keyed by space-separated names ("qc run", "benchmark integration") matching actual CLI invocation syntax
- _CLICK_TYPE_MAP includes both lowercase and uppercase variants of Click type names for robustness
- Meta commands (version, list-commands, describe, schema) excluded from discovery output via _EXCLUDED_COMMANDS set
- Stub groups (integrate, celltype) with no subcommands automatically skipped by _walk_commands

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Known Stubs

None - all discovery commands fully functional with real data from CLI tree.

## Next Phase Readiness
- Discovery commands complete; agents can now programmatically discover all CLI commands
- Ready for Phase 5 (provenance) which can use discovery to validate command registration

## Self-Check: PASSED

- FOUND: sc_tools/cli/discovery.py
- FOUND: sc_tools/tests/test_cli_discovery.py
- FOUND: .planning/phases/04-cli-discovery/04-01-SUMMARY.md
- FOUND: commit 48d0683 (test RED)
- FOUND: commit 4e2ac1c (feat GREEN)

---
*Phase: 04-cli-discovery*
*Completed: 2026-03-23*
