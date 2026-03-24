---
phase: 08-multi-omic-assembly
plan: 01
subsystem: assembly
tags: [mudata, muon, multi-omic, patient-metadata, cross-modal-query]

requires:
  - phase: 06-scientific-gaps
    provides: "validate_subject_metadata, SCToolsDataError error hierarchy"
provides:
  - "MultiOmicAtlas class wrapping MuData with patient-centric access"
  - "join_subject_metadata outer join across modalities"
  - "build_mudata MuData construction with metadata in uns"
  - "celltype_proportions cross-modal query"
  - "aggregate_by_level N-level abstraction stacking"
affects: [08-02-joint-embedding-cli]

tech-stack:
  added: [mudata, muon, mofapy2]
  patterns: [lazy-mudata-import, outer-join-patient-metadata, modality-keyed-mudata]

key-files:
  created:
    - sc_tools/assembly/__init__.py
    - sc_tools/assembly/_metadata.py
    - sc_tools/assembly/_build.py
    - sc_tools/assembly/_atlas.py
    - sc_tools/assembly/_query.py
    - sc_tools/tests/test_assembly.py
  modified:
    - pyproject.toml

key-decisions:
  - "Outer join on subject_id includes all patients even with missing modalities"
  - "Patient and sample metadata stored as DataFrames in mdata.uns"
  - "All mudata imports lazy inside function bodies (CLI-08 pattern)"
  - "Modalities missing celltype_key silently skipped in queries (no error)"
  - "h5mu format for serialization with round-trip preservation"

patterns-established:
  - "Lazy mudata import: all mudata/muon imports inside function bodies"
  - "Metadata in uns: patient_metadata and sample_metadata as DataFrames in mdata.uns"
  - "Modality-keyed MuData: dict keys become mod keys, feature spaces stay separate"

requirements-completed: [MOM-01, MOM-02, MOM-04]

duration: 3min
completed: 2026-03-24
---

# Phase 8 Plan 1: Assembly Module Core Summary

**MultiOmicAtlas class with outer-join patient metadata, MuData construction, and cross-modal celltype proportion queries**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-24T22:36:56Z
- **Completed:** 2026-03-24T22:40:26Z
- **Tasks:** 1 (TDD: RED + GREEN)
- **Files modified:** 7

## Accomplishments
- Assembly module with MultiOmicAtlas wrapping MuData for patient-centric multi-omic access
- Outer join on subject_id across modalities with boolean presence flags per modality
- Cross-modal celltype proportions and N-level abstraction stacking queries
- h5mu serialization round-trip verified
- All 11 tests pass covering metadata join, build, atlas, query, and edge cases

## Task Commits

Each task was committed atomically:

1. **Task 1 (RED): Failing tests for assembly module** - `b23dc35` (test)
2. **Task 1 (GREEN): Implement assembly module** - `3ed8814` (feat)

**Plan metadata:** [pending] (docs: complete plan)

## Files Created/Modified
- `sc_tools/assembly/__init__.py` - Public API with lazy imports for MultiOmicAtlas, build_mudata, query functions
- `sc_tools/assembly/_metadata.py` - join_subject_metadata: outer join with validation and overlap warnings
- `sc_tools/assembly/_build.py` - build_mudata: MuData construction with patient/sample metadata in uns
- `sc_tools/assembly/_atlas.py` - MultiOmicAtlas: wrapper class with patient_view, save/load, celltype_proportions
- `sc_tools/assembly/_query.py` - celltype_proportions and aggregate_by_level cross-modal queries
- `sc_tools/tests/test_assembly.py` - 11 unit tests with 4 modality fixture
- `pyproject.toml` - Added multiomics optional dependency group (mudata, muon, mofapy2)

## Decisions Made
- Outer join on subject_id: all patients included even when present in only one modality (boolean has_mod flags)
- Patient and sample metadata stored as DataFrames in mdata.uns (not separate files)
- All mudata imports lazy inside function bodies per CLI-08 pattern
- Modalities missing celltype_key silently skipped in proportion queries (logged, not raised)
- h5mu format for serialization (mudata.write/read) with round-trip metadata preservation
- validate_subject_metadata called per-modality for consistency checks (imported from sc_tools.qc.metadata)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Worktree was behind main (missing sc_tools/errors.py, sc_tools/qc/metadata.py from phases 6-7). Resolved by merging main into worktree before implementation.

## Known Stubs
None - all functions are fully implemented with no placeholder data or TODO markers.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Assembly data layer complete, ready for Plan 02 (joint embedding and CLI commands)
- MultiOmicAtlas.from_modalities() provides the entry point for embedding backends
- MuData stored in h5mu format can be loaded by Plan 02 CLI commands

## Self-Check: PASSED

- All 7 created/modified files verified on disk
- Both commits (b23dc35, 3ed8814) verified in git log
- All 11 tests pass

---
*Phase: 08-multi-omic-assembly*
*Completed: 2026-03-24*
