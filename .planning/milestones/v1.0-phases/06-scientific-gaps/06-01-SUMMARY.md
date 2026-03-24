---
phase: 06-scientific-gaps
plan: 01
subsystem: qc, celltype
tags: [subject_id, metadata_validation, panel_detection, cell_typing, confounding]

# Dependency graph
requires:
  - phase: 01-sc-tools-fixes
    provides: "Stable sc_tools internals for wrapping"
provides:
  - "validate_subject_metadata() for subject-level metadata enforcement"
  - "check_confounding() for batch-condition confounding detection"
  - "Panel-aware annotate_celltypes() with PANEL_VALIDATED_METHODS guard"
  - "panel_dispatch uns provenance for CLIResult consumption"
affects: [06-scientific-gaps, 08-multi-omic]

# Tech tracking
tech-stack:
  added: []
  patterns: [warn-dont-block validation, panel guard before backend dispatch, uns provenance storage]

key-files:
  created:
    - sc_tools/qc/metadata.py
    - sc_tools/tests/conftest.py
    - sc_tools/tests/test_subject_metadata.py
    - sc_tools/tests/test_panel_dispatch.py
  modified:
    - sc_tools/qc/__init__.py
    - sc_tools/tl/celltype/annotate.py

key-decisions:
  - "validate_subject_metadata returns warning list (never raises) per D-08/D-09"
  - "Panel detection uses adata.raw.n_vars when available to avoid HVG false positives (Pitfall 3)"
  - "panel_dispatch stored in adata.uns as canonical provenance record for future CLI bridging"
  - "1:1 subject-library mapping detected via drop_duplicates comparison"

patterns-established:
  - "Warn-dont-block: validation functions return list[str] warnings, never raise"
  - "Panel guard pattern: constants + guard clause before get_backend() dispatch"
  - "uns provenance: structured dict in adata.uns for machine-readable metadata"

requirements-completed: [SCI-03, SCI-04]

# Metrics
duration: 7min
completed: 2026-03-24
---

# Phase 06 Plan 01: Subject Metadata & Panel Dispatch Summary

**Subject-level metadata validation with confounding detection, and panel-aware cell typing dispatch restricting whole-transcriptome methods for targeted panels (n_vars < 1000)**

## Performance

- **Duration:** 7 min
- **Started:** 2026-03-24T14:20:25Z
- **Completed:** 2026-03-24T14:27:28Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments
- validate_subject_metadata handles all 4 cases: multi-sample missing, 1:1 mapping, single-sample info, batch-condition confounding
- check_confounding uses pd.crosstab to detect perfect confounding
- annotate_celltypes now blocks celltypist/scgpt/geneformer/scarches/singler for panel data (n_vars < 1000)
- force_method=True override with logged warning for advanced users
- Panel dispatch provenance stored in adata.uns['panel_dispatch'] for CLIResult bridging
- 18 new tests all passing, existing 76 QC tests unbroken

## Task Commits

Each task was committed atomically (TDD: RED then GREEN):

1. **Task 1 RED: Subject metadata tests** - `3d99128` (test)
2. **Task 1 GREEN: Subject metadata implementation** - `1335200` (feat)
3. **Task 2 RED: Panel dispatch tests** - `7ff69d4` (test)
4. **Task 2 GREEN: Panel dispatch implementation** - `e333d3a` (feat)

## Files Created/Modified
- `sc_tools/qc/metadata.py` - validate_subject_metadata and check_confounding functions
- `sc_tools/qc/__init__.py` - Added metadata module exports
- `sc_tools/tl/celltype/annotate.py` - Panel guard with PANEL_VALIDATED_METHODS, force_method param, uns provenance
- `sc_tools/tests/conftest.py` - Shared fixtures: adata_multi_subject, adata_panel
- `sc_tools/tests/test_subject_metadata.py` - 10 tests for subject_id validation and confounding
- `sc_tools/tests/test_panel_dispatch.py` - 8 tests for panel detection and method restriction

## Decisions Made
- validate_subject_metadata returns warning list (never raises) per D-08/D-09 warn-dont-block principle
- Panel detection uses adata.raw.n_vars when available to avoid HVG false positives (Pitfall 3 from RESEARCH.md)
- panel_dispatch stored in adata.uns as canonical provenance record; CLIResult bridging deferred to CLI command implementation
- 1:1 subject-library mapping detected by comparing drop_duplicates count vs library nunique

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed test fixture for valid multi-sample scenario**
- **Found during:** Task 1 (GREEN phase)
- **Issue:** adata_multi_subject fixture has 1 lib per subject (1:1 mapping), which correctly triggers distinctness warning. Test expected no warnings for "valid" data.
- **Fix:** Test for valid multi-sample uses inline fixture with 3 subjects, 6 libraries (non-1:1 mapping) instead of the shared fixture.
- **Files modified:** sc_tools/tests/test_subject_metadata.py
- **Verification:** All 10 tests pass
- **Committed in:** 1335200 (Task 1 GREEN commit)

---

**Total deviations:** 1 auto-fixed (1 bug in test expectation)
**Impact on plan:** Minor test adjustment. No scope creep.

## Issues Encountered
None.

## Known Stubs
None - all functions are fully implemented with no placeholder data.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Subject metadata validation ready for pseudobulk DE (Plan 06-02, SCI-01) which requires subject_id
- Panel dispatch ready for CLI celltype command integration when implemented
- Shared test fixtures (adata_multi_subject, adata_panel) available for Plans 06-02 and 06-03

## Self-Check: PASSED

All 4 created files verified on disk. All 4 commit hashes found in git log.

---
*Phase: 06-scientific-gaps*
*Completed: 2026-03-24*
