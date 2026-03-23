---
phase: 05-provenance
plan: 02
subsystem: provenance
tags: [lineage, trace, sidecar, h5py, sha256, typer, cli]

# Dependency graph
requires:
  - phase: 05-provenance-01
    provides: ProvenanceRecord model, sidecar read/write, checksum utility, h5ad uns embedding
  - phase: 04-cli-discovery
    provides: register_discovery(app) pattern, cli_handler decorator
provides:
  - trace_lineage BFS engine with cycle detection and SHA256 relocation
  - _read_uns_provenance h5py-based fallback reader
  - sct provenance show command (sidecar + uns fallback)
  - sct provenance trace command (recursive lineage walk)
  - register_provenance(app) registration function
affects: [06-scientific-gaps, benchmark, preprocessing]

# Tech tracking
tech-stack:
  added: []
  patterns: [bfs-lineage-walk, sha256-file-relocation, register-command-group-pattern]

key-files:
  created:
    - sc_tools/provenance/trace.py
    - sc_tools/cli/provenance.py
    - sc_tools/tests/test_cli_provenance.py
  modified:
    - sc_tools/cli/__init__.py

key-decisions:
  - "BFS with visited set for cycle detection (unlimited depth)"
  - "SHA256 relocation scan capped at 1000 files to avoid huge dir traversals"
  - "Relocation note preserved on origin nodes found via SHA256"
  - "register_provenance pattern follows register_discovery convention"
  - "Lazy imports inside CLI commands for fast startup"

patterns-established:
  - "register_{module}(app) pattern for adding command groups to sct CLI"
  - "BFS lineage walk with canonical path cycle detection"
  - "h5py-only uns reader for provenance without AnnData load"

requirements-completed: [PRV-03, PRV-04]

# Metrics
duration: 7min
completed: 2026-03-23
---

# Phase 5 Plan 2: Provenance CLI Commands Summary

**BFS lineage trace engine with cycle detection, SHA256 relocation, and adata.uns fallback; sct provenance show/trace CLI commands registered via register_provenance(app) pattern**

## Performance

- **Duration:** 7 min
- **Started:** 2026-03-23T19:27:01Z
- **Completed:** 2026-03-23T19:34:00Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Lineage trace engine (trace.py) with BFS walk, visited-set cycle detection, chronological output (oldest first), adata.uns h5py fallback, and SHA256-based file relocation
- `sct provenance show` reads sidecar or falls back to adata.uns for h5ad files; raises SCToolsUserError when no provenance found
- `sct provenance trace` walks lineage chain to raw data origins with --project-dir support
- 15 tests: 9 trace engine unit tests + 6 CLI integration tests

## Task Commits

Each task was committed atomically:

1. **Task 1: Lineage trace engine (trace.py)** (TDD)
   - `792d117` (test: failing tests - RED phase)
   - `2877eb1` (feat: implementation - GREEN phase)
2. **Task 2: Provenance CLI commands (show + trace) and registration** (TDD)
   - `af93bc7` (feat: CLI commands + tests - GREEN phase, RED embedded)

## Files Created/Modified
- `sc_tools/provenance/trace.py` - BFS lineage trace engine with _read_uns_provenance and _find_by_sha256
- `sc_tools/cli/provenance.py` - register_provenance(app) with show and trace commands
- `sc_tools/cli/__init__.py` - Added register_provenance(app) call at bottom
- `sc_tools/tests/test_cli_provenance.py` - 15 tests covering all trace and CLI functionality

## Decisions Made
- BFS with visited set for cycle detection rather than depth limit -- handles arbitrary DAGs
- SHA256 relocation scan capped at 1000 files to prevent scanning enormous project directories
- Relocation note preserved even on origin nodes ("input relocated; origin (no provenance)")
- register_provenance follows same pattern as register_discovery for consistency
- All heavy imports (h5py, provenance modules) are lazy inside command functions

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] SHA256 relocation note not preserved on origin nodes**
- **Found during:** Task 1 (GREEN phase)
- **Issue:** When a relocated file had no provenance, the relocation_note was overwritten with "origin (no provenance)"
- **Fix:** Combined notes: "input relocated; origin (no provenance)" when file found by SHA256 but has no provenance
- **Files modified:** sc_tools/provenance/trace.py
- **Committed in:** 2877eb1

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Fix necessary for accurate relocation tracking. No scope creep.

## Issues Encountered
- Pre-existing test failure in test_cli_commands.py (preprocess missing file exit code 3 vs expected 1) -- unrelated to our changes
- Pre-existing test failure in test_provenance.py (missing sc_tools.pp.strategy module) -- unrelated to our changes

## Known Stubs
None -- all provenance commands and trace engine are fully implemented with real logic.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 5 (Provenance) is now complete: all PRV requirements fulfilled
- Sidecar infrastructure (Plan 01) + CLI commands and trace engine (Plan 02)
- Ready for Phase 6 (Scientific Gaps)

## Self-Check: PASSED

- All 4 created/modified files verified on disk
- All 3 task commits verified in git log (792d117, 2877eb1, af93bc7)

---
*Phase: 05-provenance*
*Completed: 2026-03-23*
