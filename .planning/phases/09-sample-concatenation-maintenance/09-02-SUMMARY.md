---
phase: 09-sample-concatenation-maintenance
plan: 02
subsystem: cli-concat
tags: [concat, spatial, provenance, pipeline-dag]
dependency_graph:
  requires: [cli-foundation, provenance-sidecars, pipeline-dag]
  provides: [sct-concat-command, concat-pipeline-phase]
  affects: [ingest-loaders, cli-init, pipeline-phases]
tech_stack:
  added: []
  patterns: [register_X command pattern, uns_merge=unique for spatial, _input_files provenance convention]
key_files:
  created:
    - sc_tools/cli/concat.py
    - sc_tools/tests/test_concat.py
  modified:
    - sc_tools/cli/__init__.py
    - sc_tools/ingest/loaders.py
    - sc_tools/pipeline.py
    - sc_tools/tests/test_pipeline.py
decisions:
  - "Used register_concat(app) direct command pattern (not sub-app) since concat is a single command"
  - "concat PhaseSpec is optional=True so projects can skip it and feed per-sample h5ads directly to qc_filter"
  - "Fixed existing concat_samples uns_merge='same' bug alongside new CLI implementation"
metrics:
  duration: 6min
  completed: "2026-03-28T02:47:25Z"
  tasks: 3
  files: 6
---

# Phase 09 Plan 02: sct concat CLI Command Summary

Implemented `sct concat` CLI command for merging same-modality h5ad files with spatial data preservation via `uns_merge="unique"`, provenance sidecar via `_input_files` convention, and registered concat as an optional pipeline phase between ingest_load and qc_filter.

## Tasks Completed

### Task 1: Create test scaffolds (TDD RED)
- Created `sc_tools/tests/test_concat.py` with 10 test methods across 4 test classes
- TestConcatCommand (3 tests), TestSpatialPreservation (2), TestConcatProvenance (2), TestConcatPipeline (3)
- RED state confirmed: all tests fail with ModuleNotFoundError
- Commit: `cd34715`

### Task 2: Implement sct concat CLI command (CONCAT-01, CONCAT-02, CONCAT-03)
- Created `sc_tools/cli/concat.py` with `register_concat` following the established pattern
- Validates >= 2 inputs, checks file existence, pre-flight gene overlap check via h5py
- Concatenates with `ad.concat(uns_merge="unique")` preserving spatial metadata
- Post-concat verification ensures all spatial keys survived
- Uses `smart_write_checkpoint` for output, `_input_files` convention for provenance
- Registered in `sc_tools/cli/__init__.py`
- Fixed `sc_tools/ingest/loaders.py` line 930: `uns_merge="same"` -> `uns_merge="unique"`
- All 7 implementation tests passing (GREEN state)
- Commit: `a30d882`

### Task 3: Register concat as pipeline phase (CONCAT-04)
- Added `concat` PhaseSpec to STANDARD_PHASES between ingest_load and qc_filter
- Set `optional=True` so projects can skip concat and go directly to qc_filter
- Added TestConcatPhase class with 6 tests in test_pipeline.py
- Updated existing all_phases_complete tests to include concat
- Fixed lint issues (import ordering, zip strict parameter)
- All 77 pipeline+concat tests passing, lint clean
- Commit: `8748c04`

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed test invocation pattern**
- **Found during:** Task 2
- **Issue:** Tests passed "concat" as CLI argument to CliRunner but register_concat registers the command directly on the test app, causing "unexpected extra argument" error
- **Fix:** Removed "concat" from CliRunner invoke args since the command is the only command on the test app
- **Files modified:** sc_tools/tests/test_concat.py

**2. [Rule 1 - Bug] Fixed return type annotation**
- **Found during:** Task 2
- **Issue:** `-> CLIResult` annotation on concat function caused NameError with `from __future__ import annotations` and Typer's inspect.signature
- **Fix:** Changed to `-> None` following the established pattern in assemble.py
- **Files modified:** sc_tools/cli/concat.py

**3. [Rule 1 - Bug] Fixed provenance sidecar path in tests**
- **Found during:** Task 2
- **Issue:** Test used `.with_suffix(".provenance.json")` which replaces `.h5ad`, but actual convention is append `.provenance.json`
- **Fix:** Changed to `Path(str(output) + ".provenance.json")`
- **Files modified:** sc_tools/tests/test_concat.py

**4. [Rule 3 - Blocking] Updated existing pipeline tests**
- **Found during:** Task 3
- **Issue:** `test_all_phases_complete` tests listed only 10 phases; adding concat (optional, non-iterative) caused assertion failure
- **Fix:** Added "concat" to all_phases lists in both flat-slug and tuple-based tests
- **Files modified:** sc_tools/tests/test_pipeline.py

## Verification Results

- `pytest sc_tools/tests/test_concat.py`: 10 passed
- `pytest sc_tools/tests/test_pipeline.py -k concat`: 6 passed
- `grep 'uns_merge="same"' sc_tools/ingest/loaders.py`: no matches (bug fixed)
- `ruff check` on modified files: all checks passed

## Known Stubs

None - all functionality is fully wired.

## Self-Check: PASSED

- [x] sc_tools/cli/concat.py exists
- [x] sc_tools/tests/test_concat.py exists
- [x] Commits cd34715, a30d882, 8748c04 exist
- [x] All tests pass, lint clean
