---
phase: 07-memory-safety
verified: 2026-03-24T22:00:00Z
status: gaps_found
score: 5/6 must-haves verified
re_verification: false
gaps:
  - truth: "T3 load refuses when estimated memory exceeds 80% of available RAM"
    status: partial
    reason: "Memory guard raises sc_tools.io.errors.SCToolsRuntimeError (base Exception) instead of sc_tools.errors.SCToolsRuntimeError (SCToolsError subclass). CLI handler catches the wrong type -- memory guard errors fall through to the generic except Exception handler, lose the 'Use --force to override' suggestion, and appear as SCToolsFatalError with a generic message."
    artifacts:
      - path: "sc_tools/io/errors.py"
        issue: "SCToolsRuntimeError inherits from Exception directly, not from sc_tools.errors.SCToolsError. No suggestion attribute."
      - path: "sc_tools/io/gateway.py"
        issue: "Imports and raises sc_tools.io.errors.SCToolsRuntimeError, which is not in the CLI's catch chain for SCToolsRuntimeError."
    missing:
      - "sc_tools/io/errors.py: Make SCToolsRuntimeError inherit from sc_tools.errors.SCToolsRuntimeError, or remove the file and import from sc_tools.errors directly in gateway.py"
      - "Verify that the CLI --force suggestion appears in error output when memory guard fires via sct qc run or sct preprocess run"
---

# Phase 7: Memory Safety Verification Report

**Phase Goal:** Large datasets (2.5M cells, 25G h5ad) can be processed without OOM through tiered loading, pre-execution estimation, and dry-run validation
**Verified:** 2026-03-24T22:00:00Z
**Status:** gaps_found
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | IOGateway.read() with T1 returns metadata dict without loading array data | VERIFIED | gateway.py:64-65, test_t1_returns_dict_not_anndata passes |
| 2 | IOGateway.read() with T2 returns backed AnnData (obs/var accessible, X on disk) | VERIFIED | gateway.py:67-70, test_t2_backed_read passes |
| 3 | IOGateway.read() with T3 loads full AnnData into memory | VERIFIED | gateway.py:72-78, test_t3_full_read passes |
| 4 | T3 load refuses when estimated memory exceeds 80% of available RAM | PARTIAL | _check_memory_guard raises correct error in isolation (test passes). However sc_tools.io.errors.SCToolsRuntimeError does not inherit from sc_tools.errors.SCToolsRuntimeError, so when triggered via the CLI the exception falls through to the generic handler and the actionable "Use --force to override" suggestion is lost. |
| 5 | T3 load with force=True bypasses memory guard | VERIFIED | gateway.py:73-74, test_force_override passes |
| 6 | estimate_from_h5 returns estimated_peak_mb for both dense and sparse X | VERIFIED | estimate.py:24-77, all 4 estimate tests pass |
| 7 | sct estimate preprocess run returns JSON with estimated_peak_mb and estimated_runtime_s | VERIFIED | cli/estimate.py, all 5 CLI estimate tests pass |
| 8 | --dry-run on sct qc run validates inputs and exits 0 without modifying data | VERIFIED | cli/__init__.py dry_run block, test_dryrun_qc_run and test_dryrun_no_side_effects pass |
| 9 | --dry-run output is a valid CLIResult JSON with status=skipped | VERIFIED | cli/__init__.py:267 Status.skipped, all dry-run tests confirm |
| 10 | All data-touching commands accept --dry-run and --force flags | VERIFIED | qc.py, preprocess.py, benchmark.py, de.py, validate.py all contain --dry-run |

**Score:** 5/6 core truths fully verified (truth #4 partial due to error class inheritance gap)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `sc_tools/io/__init__.py` | Public API: IOGateway, DataTier, read_h5ad_metadata, estimate_from_h5 | VERIFIED | All 4 exports present, __all__ declared |
| `sc_tools/io/gateway.py` | IOGateway class with tiered read dispatch and memory guard | VERIFIED | DataTier enum + IOGateway.read() + _check_memory_guard() all present |
| `sc_tools/io/metadata.py` | h5py metadata extraction with n_obs/n_vars fallback chain | VERIFIED | read_h5ad_metadata with _index_length -> _index -> first-key fallback |
| `sc_tools/io/estimate.py` | Memory and runtime estimation from h5 metadata | VERIFIED | estimate_from_h5 + METHOD_MULTIPLIERS dict |
| `sc_tools/io/errors.py` | SCToolsRuntimeError for memory guard (deviation from plan) | PARTIAL | Created as standalone class inheriting Exception, not SCToolsError hierarchy |
| `sc_tools/cli/estimate.py` | sct estimate command with memory/runtime projection | VERIFIED | register_estimate + estimate_command + estimate_from_h5 + METHOD_MULTIPLIERS |
| `sc_tools/cli/__init__.py` | Extended cli_handler with tier, dry_run, force support | VERIFIED | cli_handler(func=None, *, tier=None) with dry_run/force pop and Status.skipped |
| `sc_tools/tests/test_io_gateway.py` | Unit tests for T1/T2/T3 reads and memory guard | VERIFIED | 7 tests all pass |
| `sc_tools/tests/test_estimate.py` | Unit tests for memory estimation (dense and sparse) | VERIFIED | 4 tests all pass |
| `sc_tools/tests/test_cli_estimate.py` | Integration tests for sct estimate CLI | VERIFIED | 5 tests all pass |
| `sc_tools/tests/test_cli_dryrun.py` | Integration tests for --dry-run on data-touching commands | VERIFIED | 6 tests all pass |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| sc_tools/io/gateway.py | sc_tools/io/metadata.py | import read_h5ad_metadata | VERIFIED | gateway.py:17 direct import |
| sc_tools/io/gateway.py | sc_tools/io/estimate.py | import estimate_from_h5 | VERIFIED | gateway.py:18 direct import |
| sc_tools/io/gateway.py | psutil | virtual_memory().available for memory guard | VERIFIED | gateway.py:87 psutil lazy-imported inside _check_memory_guard |
| sc_tools/cli/__init__.py | sc_tools/io (estimate) | import estimate_from_h5 for dry-run | VERIFIED | cli/__init__.py:254 lazy import inside wrapper |
| sc_tools/cli/estimate.py | sc_tools/io/estimate.py | import estimate_from_h5 for memory projection | VERIFIED | cli/estimate.py:28 |
| sc_tools/cli/__init__.py | sc_tools/models/result.py | CLIResult with status=skipped for dry-run | VERIFIED | cli/__init__.py:267 Status.skipped |
| sc_tools/io/errors.py | sc_tools/errors.SCToolsRuntimeError | inheritance for CLI error catch chain | NOT_WIRED | sc_tools.io.errors.SCToolsRuntimeError inherits Exception directly, not SCToolsError. CLI handler cannot catch it as SCToolsRuntimeError. |

### Data-Flow Trace (Level 4)

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|---------------|--------|--------------------|--------|
| sc_tools/io/estimate.py | estimated_peak_mb | h5py metadata read via read_h5ad_metadata | Yes -- n_obs, n_vars, nnz from actual h5 file | FLOWING |
| sc_tools/cli/estimate.py | est dict | estimate_from_h5(file_path) | Yes -- reads file metadata, applies multiplier | FLOWING |
| sc_tools/cli/__init__.py (dry_run) | dry_data["estimate"] | estimate_from_h5(file_arg) | Yes -- real h5 metadata when file exists | FLOWING |

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| Import public API | `python -c "from sc_tools.io import IOGateway, DataTier, read_h5ad_metadata, estimate_from_h5; print('OK')"` | OK | PASS |
| No top-level heavy imports in sc_tools/io/ | grep `^import h5py\|^import anndata\|^import psutil` in sc_tools/io/*.py | No matches | PASS |
| All 22 phase 07 tests | `python -m pytest test_io_gateway.py test_estimate.py test_cli_estimate.py test_cli_dryrun.py -v` | 22 passed in 1.64s | PASS |
| Memory guard cross-class compatibility | `issubclass(sc_tools.io.errors.SCToolsRuntimeError, sc_tools.errors.SCToolsRuntimeError)` | False | FAIL |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| MEM-01 | 07-01 | IO Gateway with tiered loading strategy: h5py for metadata/embeddings, backed mode for summaries, full load only for compute | SATISFIED | IOGateway with T1/T2/T3, all tests pass |
| MEM-02 | 07-01, 07-02 | sct estimate -- pre-execution estimation of peak memory and runtime | SATISFIED | sc_tools/cli/estimate.py, test_cli_estimate.py all pass |
| MEM-03 | 07-02 | --dry-run flag for all data-touching commands -- validate inputs, report what would happen, without executing | SATISFIED | All 5 data-touching commands have --dry-run, test_cli_dryrun.py all pass |

No orphaned requirements: all Phase 7 requirement IDs (MEM-01, MEM-02, MEM-03) appear in plan frontmatter. REQUIREMENTS.md traceability table marks all three as Complete.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| sc_tools/io/errors.py | 6 | `SCToolsRuntimeError(Exception)` -- not in SCToolsError hierarchy | Warning | Memory guard exception is not caught by CLI handler's `except SCToolsRuntimeError` branch; falls to generic handler. Exit code 3 still produced but error message loses "Use --force to override" suggestion. |

### Human Verification Required

#### 1. Memory Guard CLI Surface Behavior

**Test:** Run `sct qc run <large_file.h5ad>` on a file large enough to trigger the 80% RAM guard (without --force), then inspect the JSON error output.
**Expected:** Error JSON should contain suggestion "Use --force to override or reduce dataset size." and status="error" with exit code 3.
**Why human:** Requires an actual large h5ad file to trigger the guard in a real CLI invocation. Programmatic check confirmed the exception class mismatch -- the suggestion text will be the generic fallback, not the memory-specific one.

### Gaps Summary

One gap blocking full goal achievement at the CLI surface level:

**Memory guard error class isolation:** The SUMMARY noted that `sc_tools/errors.py` "did not exist" when Plan 01 ran, so `sc_tools/io/errors.py` was created as a workaround. However, `sc_tools/errors.py` does exist (confirmed) and contains `SCToolsRuntimeError`. The io module's `SCToolsRuntimeError` was created as a standalone class inheriting `Exception` rather than re-using or extending the existing one. As a result, when the memory guard fires during a CLI command (e.g., `sct qc run`), the CLI handler catches it via `except Exception` (the generic fallback), not via `except SCToolsRuntimeError` (the intended path). The exit code is still 3, but the actionable suggestion "Use --force to override" is replaced by "This is an unexpected error. Please report it."

**Fix:** In `sc_tools/io/errors.py`, change `class SCToolsRuntimeError(Exception)` to `from sc_tools.errors import SCToolsRuntimeError` and remove the duplicate definition, OR change it to `class SCToolsRuntimeError(SCToolsRuntimeError)` via proper import. The simplest fix is to remove `sc_tools/io/errors.py` entirely and update `gateway.py` to `from sc_tools.errors import SCToolsRuntimeError`.

All other truths are fully satisfied. All 22 tests pass. No stubs or placeholder implementations found. Lazy import invariant (CLI-08) is maintained throughout.

---

_Verified: 2026-03-24T22:00:00Z_
_Verifier: Claude (gsd-verifier)_
