---
phase: 05-provenance
plan: 01
subsystem: provenance
tags: [pydantic, sha256, sidecar, h5py, leiden, reproducibility]

# Dependency graph
requires:
  - phase: 02-cli-foundation
    provides: CLIResult model, cli_handler decorator, _emit() function
  - phase: 03-core-commands
    provides: Command implementations populating artifacts field
provides:
  - InputFile and ProvenanceRecord Pydantic models
  - sc_tools/provenance/ package (checksum, sidecar utilities)
  - Automatic sidecar writing in cli_handler
  - h5ad uns embedding for provenance portability
  - random_state threading through Leiden clustering paths
affects: [05-02-provenance-commands, benchmark, preprocessing]

# Tech tracking
tech-stack:
  added: []
  patterns: [atomic-sidecar-write, _input_files-convention-key, h5py-append-mode-uns]

key-files:
  created:
    - sc_tools/provenance/__init__.py
    - sc_tools/provenance/checksum.py
    - sc_tools/provenance/sidecar.py
    - sc_tools/tests/test_provenance.py
  modified:
    - sc_tools/models/result.py
    - sc_tools/cli/__init__.py
    - sc_tools/bm/integration.py
    - sc_tools/pp/reduce.py

key-decisions:
  - "ProvenanceRecord separate from Provenance class for backwards compatibility"
  - "_input_files convention key in CLIResult.data popped before emission"
  - "h5py append mode for uns embedding avoids full AnnData load"
  - "Atomic sidecar write via tempfile + os.replace"
  - "random_state=0 default for deterministic Leiden clustering"

patterns-established:
  - "Sidecar convention: {artifact}.provenance.json alongside every CLI output"
  - "_input_files list in result.data for provenance input tracking"
  - "h5ad outputs get provenance in both sidecar and adata.uns['sct_provenance']"

requirements-completed: [PRV-01, PRV-02, PRV-05]

# Metrics
duration: 12min
completed: 2026-03-23
---

# Phase 5 Plan 1: Provenance Infrastructure Summary

**ProvenanceRecord model with SHA256 checksums, automatic .provenance.json sidecar writing in cli_handler, h5ad uns embedding via h5py, and random_state threading through Leiden clustering**

## Performance

- **Duration:** 12 min
- **Started:** 2026-03-23T19:11:55Z
- **Completed:** 2026-03-23T19:24:00Z
- **Tasks:** 2
- **Files modified:** 8

## Accomplishments
- InputFile and ProvenanceRecord Pydantic models with full provenance fields (command, params, inputs with SHA256, version, timestamp, runtime_s, peak_memory_mb)
- Provenance utility package: streaming SHA256 hasher, atomic sidecar write/read, peak memory measurement (macOS/Linux aware), h5ad uns embedding
- cli_handler automatically writes .provenance.json sidecars for every successful command with artifacts, and embeds provenance in adata.uns for h5ad outputs
- random_state=0 default threaded through _leiden_cluster, _ari/_nmi, compute_integration_metrics, compare_integrations, and pp.reduce.cluster
- 20 unit tests covering all provenance functionality

## Task Commits

Each task was committed atomically:

1. **Task 1: Provenance model, utilities, sidecar hook, and h5ad uns embedding** (TDD)
   - `a0e265c` (test: failing tests - RED phase)
   - `8c2dbbe` (feat: implementation - GREEN phase)
2. **Task 2: Thread random_state through Leiden clustering** - `49a9e9c` (feat)

## Files Created/Modified
- `sc_tools/provenance/__init__.py` - Package init
- `sc_tools/provenance/checksum.py` - SHA256 streaming file hasher
- `sc_tools/provenance/sidecar.py` - Sidecar write/read/path, peak memory, build_provenance_record, embed_provenance_in_adata
- `sc_tools/models/result.py` - Added InputFile and ProvenanceRecord models
- `sc_tools/cli/__init__.py` - Added _write_provenance_sidecars hook, time.monotonic timing, _input_files handling
- `sc_tools/bm/integration.py` - random_state parameter through _leiden_cluster, _ari, _nmi, compute_integration_metrics, compare_integrations
- `sc_tools/pp/reduce.py` - random_state parameter in cluster()
- `sc_tools/tests/test_provenance.py` - 20 unit tests

## Decisions Made
- ProvenanceRecord is a separate class from the existing minimal Provenance to avoid breaking CLIResult backwards compatibility
- _input_files convention key in result.data is popped before _emit() so it never appears in JSON stdout output
- h5py append mode ("a") for uns embedding avoids loading the full AnnData into memory
- Atomic sidecar write uses tempfile + os.replace to prevent corruption from parallel CLI execution
- random_state defaults to 0 (not 42) per D-14 for consistent reproducibility
- Platform-aware peak memory: resource.getrusage returns bytes on macOS, KB on Linux

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed mock target for h5ad embedding tests**
- **Found during:** Task 1 (GREEN phase)
- **Issue:** Tests patching `sc_tools.cli.embed_provenance_in_adata` failed because the import is lazy (inside _write_provenance_sidecars)
- **Fix:** Changed mock targets to `sc_tools.provenance.sidecar.embed_provenance_in_adata` and `sc_tools.provenance.sidecar.write_sidecar`
- **Files modified:** sc_tools/tests/test_provenance.py
- **Committed in:** 8c2dbbe

**2. [Rule 1 - Bug] Fixed leiden mock in random_state test**
- **Found during:** Task 2
- **Issue:** Mocking sc.tl.leiden prevented creation of 'leiden' obs column, causing KeyError
- **Fix:** Used side_effect mock that creates the leiden column before returning
- **Files modified:** sc_tools/tests/test_provenance.py
- **Committed in:** 49a9e9c

---

**Total deviations:** 2 auto-fixed (2 bugs in tests)
**Impact on plan:** Both fixes necessary for test correctness. No scope creep.

## Issues Encountered
- Pre-existing test failure in test_bm_report.py (Plotly embed assertion) unrelated to our changes -- logged as out-of-scope

## Known Stubs
None -- all provenance functions are fully implemented with real logic.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Provenance infrastructure complete, ready for Plan 02 (`sct provenance show` and `sct provenance trace` commands)
- All sidecar utilities exported and tested
- cli_handler hook operational for all existing commands

---
*Phase: 05-provenance*
*Completed: 2026-03-23*
