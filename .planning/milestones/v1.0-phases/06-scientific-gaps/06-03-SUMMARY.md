---
phase: 06-scientific-gaps
plan: 03
subsystem: qc
tags: [marker-validation, celltyping, html-report, scanpy, dotplot]

# Dependency graph
requires:
  - phase: 06-01
    provides: conftest fixtures (adata_multi_subject, adata_panel)
provides:
  - compute_marker_validation function for post-celltyping QC
  - render_marker_dotplot for base64 dotplot rendering
  - Marker Validation section in post-celltyping HTML report
affects: [phase-07, phase-08]

# Tech tracking
tech-stack:
  added: []
  patterns: [threshold-based flagging with informational-only semantics]

key-files:
  created:
    - sc_tools/qc/marker_validation.py
    - sc_tools/tests/test_marker_validation.py
  modified:
    - sc_tools/qc/__init__.py
    - sc_tools/qc/report.py
    - sc_tools/assets/post_celltyping_qc_template.html

key-decisions:
  - "Flagging uses all-markers-below-threshold rule per D-05: celltype flagged only when ALL canonical markers are below threshold"
  - "has_marker_validation stored as True or None (not False) for sidebar section filter compatibility"
  - "NaN genes excluded from all-below decision to avoid false flagging on missing genes"

patterns-established:
  - "Informational flagging pattern: compute returns data, never raises"
  - "Template conditional sections guarded by has_* context vars mapped to None for sidebar filtering"

requirements-completed: [SCI-02]

# Metrics
duration: 7min
completed: 2026-03-24
---

# Phase 06 Plan 03: Marker Validation Report Summary

**Marker validation compute + post-celltyping report integration with threshold-based flagging, dotplot, and flag table**

## Performance

- **Duration:** 7 min
- **Started:** 2026-03-24T14:32:17Z
- **Completed:** 2026-03-24T14:39:27Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- `compute_marker_validation` correctly identifies flagged celltypes where all canonical markers fall below threshold
- `render_marker_dotplot` produces base64 PNG via scanpy.pl.dotplot
- Post-celltyping HTML report conditionally shows Marker Validation section with summary cards, dotplot, and flag table
- Backward compatible: no marker section rendered when marker_genes=None
- 10 tests (7 unit + 3 integration) all passing, existing celltyping tests unaffected

## Task Commits

Each task was committed atomically:

1. **Task 1 RED: Failing tests for marker validation** - `865b1c2` (test)
2. **Task 1 GREEN: Implement marker validation compute and dotplot** - `0647b96` (feat)
3. **Task 2: Integrate marker validation into post-celltyping report** - `000325a` (feat)

_TDD task had separate RED and GREEN commits._

## Files Created/Modified
- `sc_tools/qc/marker_validation.py` - compute_marker_validation and render_marker_dotplot functions
- `sc_tools/qc/__init__.py` - Added marker_validation module and compute_marker_validation export
- `sc_tools/qc/report.py` - Extended generate_post_celltyping_report with marker validation context
- `sc_tools/assets/post_celltyping_qc_template.html` - Added Marker Validation section with cards, dotplot, flag table
- `sc_tools/tests/test_marker_validation.py` - 10 tests covering unit and integration scenarios

## Decisions Made
- Flagging uses all-markers-below-threshold rule: a celltype is flagged only when ALL its canonical markers have mean expression below threshold (NaN genes excluded from this decision)
- `has_marker_validation` stored as `True` or `None` (not `False`) to work with the existing sidebar section filter which checks `is not None`
- Informational note in template explains flagging does not indicate failure

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed has_marker_validation boolean filtering for sidebar nav**
- **Found during:** Task 2 (report integration)
- **Issue:** `has_marker_validation=False` passed `is not None` check, causing "Marker Validation" to appear in sidebar even when disabled
- **Fix:** Changed to `has_marker_validation or None` so False becomes None for sidebar filtering
- **Files modified:** sc_tools/qc/report.py
- **Verification:** test_report_without_marker_genes_no_validation passes
- **Committed in:** 000325a (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Auto-fix necessary for backward compatibility. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Known Stubs
None - all functions are fully wired with real data paths.

## Next Phase Readiness
- Marker validation module available for CLI wrapping in Phase 7+
- Post-celltyping report template fully extended and backward compatible
- All phase 06 plans now complete (06-01, 06-02, 06-03)

---
*Phase: 06-scientific-gaps*
*Completed: 2026-03-24*
