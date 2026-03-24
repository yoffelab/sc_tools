---
phase: 07-memory-safety
plan: 01
subsystem: io
tags: [h5py, anndata, psutil, memory-guard, tiered-loading, lazy-import]

# Dependency graph
requires:
  - phase: 02-cli-foundation
    provides: "CLI lazy import pattern (CLI-08)"
  - phase: 05-provenance
    provides: "h5py append pattern for lightweight file access"
provides:
  - "IOGateway class with T1/T2/T3 tiered read dispatch"
  - "DataTier enum (T1_METADATA, T2_SUMMARY, T3_FULL)"
  - "read_h5ad_metadata -- h5py-only metadata extraction with fallback chain"
  - "estimate_from_h5 -- pre-load memory estimation for dense and sparse X"
  - "SCToolsRuntimeError for memory guard violations"
  - "METHOD_MULTIPLIERS dict for command-specific estimation"
affects: [07-02, 08-multi-omic]

# Tech tracking
tech-stack:
  added: []
  patterns: ["tiered IO dispatch via enum", "h5py lazy import inside methods", "pre-load memory guard with force override"]

key-files:
  created:
    - sc_tools/io/__init__.py
    - sc_tools/io/gateway.py
    - sc_tools/io/metadata.py
    - sc_tools/io/estimate.py
    - sc_tools/io/errors.py
    - sc_tools/tests/test_io_gateway.py
    - sc_tools/tests/test_estimate.py
  modified: []

key-decisions:
  - "SCToolsRuntimeError created in sc_tools/io/errors.py (sc_tools/errors.py did not exist)"
  - "Memory guard uses 80% of available RAM threshold with psutil.virtual_memory"
  - "estimate_from_h5 uses 2.0x overhead multiplier for peak estimation"
  - "All heavy imports (h5py, anndata, psutil) lazy inside method bodies per CLI-08"

patterns-established:
  - "Tiered IO dispatch: T1 returns dict (h5py only), T2 returns backed AnnData, T3 returns full AnnData"
  - "Pre-load memory guard: estimate before load, block if exceeds threshold, force=True bypasses"
  - "n_obs/n_vars fallback chain: attrs._index_length -> _index dataset length -> first key length"

requirements-completed: [MEM-01, MEM-02]

# Metrics
duration: 5min
completed: 2026-03-24
---

# Phase 7 Plan 1: IO Gateway Summary

**Tiered IO Gateway with h5py metadata reader, memory estimator, and 80%-RAM guard for safe large-file loading**

## Performance

- **Duration:** 5 min
- **Started:** 2026-03-24T21:25:18Z
- **Completed:** 2026-03-24T21:30:25Z
- **Tasks:** 1 (TDD: RED + GREEN)
- **Files created:** 7

## Accomplishments
- IOGateway class dispatches reads by DataTier enum: T1 returns metadata dict via h5py, T2 returns backed AnnData, T3 returns full AnnData with memory guard
- Pre-load memory guard blocks T3 loads when estimated peak memory exceeds 80% of available RAM; force=True bypasses
- estimate_from_h5 computes accurate memory estimates for both dense and sparse h5ad files from h5py metadata alone
- All 11 unit tests pass covering metadata extraction, tiered loading, memory guard blocking, and force override

## Task Commits

Each task was committed atomically:

1. **Task 1 (RED): Failing tests for IO gateway** - `d453875` (test)
2. **Task 1 (GREEN): IO Gateway implementation** - `33b839e` (feat)

## Files Created/Modified
- `sc_tools/io/__init__.py` - Public API exports: IOGateway, DataTier, read_h5ad_metadata, estimate_from_h5
- `sc_tools/io/gateway.py` - IOGateway class with tiered read dispatch and memory guard
- `sc_tools/io/metadata.py` - h5py metadata extraction with n_obs/n_vars fallback chain
- `sc_tools/io/estimate.py` - Memory estimation from h5 metadata (dense + sparse X)
- `sc_tools/io/errors.py` - SCToolsRuntimeError for memory guard violations
- `sc_tools/tests/test_io_gateway.py` - 7 tests: T1/T2/T3 reads, memory guard, force override
- `sc_tools/tests/test_estimate.py` - 4 tests: dense/sparse estimation, required keys, zero-gene edge case

## Decisions Made
- Created `sc_tools/io/errors.py` with `SCToolsRuntimeError` since `sc_tools/errors.py` referenced in plan did not exist (Rule 3 -- blocking issue)
- Memory guard uses `psutil.virtual_memory().available` for available RAM check (consistent with existing `sc_tools/memory/profiling.py` pattern)
- Overhead multiplier of 2.0x applied to base memory estimate for peak estimation (per RESEARCH.md guidance)
- METHOD_MULTIPLIERS dict included in estimate.py for Plan 02 estimate command to use

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Created sc_tools/io/errors.py**
- **Found during:** Task 1 (implementation)
- **Issue:** Plan referenced `from sc_tools.errors import SCToolsRuntimeError` but `sc_tools/errors.py` did not exist in the codebase
- **Fix:** Created `sc_tools/io/errors.py` with `SCToolsRuntimeError` class matching the spec (category="fatal", exit_code=3)
- **Files created:** `sc_tools/io/errors.py`
- **Verification:** Tests import and catch the error correctly
- **Committed in:** 33b839e (Task 1 GREEN commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Error class co-located in io package rather than top-level. No scope creep.

## Issues Encountered
- Test mocking required patching `psutil.virtual_memory` at the psutil module level (not gateway module level) because psutil is lazily imported inside the method. Resolved by mocking `estimate_from_h5` return value to produce large estimates that trigger the guard.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- IO Gateway module ready for Plan 02 (sct estimate command, dry-run mode)
- IOGateway.read() provides the tiered loading API that Plan 02 will wire into CLI commands
- estimate_from_h5 and METHOD_MULTIPLIERS ready for the estimate command

---
*Phase: 07-memory-safety*
*Completed: 2026-03-24*
