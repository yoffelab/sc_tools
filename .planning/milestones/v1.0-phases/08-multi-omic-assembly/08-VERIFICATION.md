---
phase: 08-multi-omic-assembly
verified: 2026-03-24T23:10:00Z
status: passed
score: 12/12 must-haves verified
re_verification: false
gaps: []
human_verification:
  - test: "Run sct assemble embed on a real multi-omic h5mu file with mofapy2 installed"
    expected: "h5mu file gains obsm['X_mofa'], CLIResult status=success with embedding_shape"
    why_human: "mofapy2 is not installed in local dev environment; MOFA+ integration tests skipped (3 skipped)"
---

# Phase 8: Multi-Omic Assembly Verification Report

**Phase Goal:** Independently processed modalities (scRNA, IMC, Visium, Xenium) can be assembled into a unified MuData object for cross-modal patient-level analysis
**Verified:** 2026-03-24T23:10:00Z
**Status:** PASSED
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Patient metadata joins across modalities using outer join on subject_id, including patients with missing modalities | VERIFIED | `join_subject_metadata` builds union of all subject_ids with boolean `has_{mod}` flags; `test_metadata_outer_join` confirms PAT3 present even when absent from imc/xenium |
| 2 | MuData object is created from per-modality AnnData objects with modality keys (rna, imc, visium, xenium) | VERIFIED | `build_mudata` calls `mudata.MuData(modalities)` with dict keys; `test_mudata_build` confirms 4-key MuData with correct shapes |
| 3 | MultiOmicAtlas serializes to h5mu and round-trips correctly | VERIFIED | `atlas.save()` calls `mdata.write(path)`; `atlas.load()` calls `mudata.read(path)`; `test_atlas_roundtrip` passes confirming modalities and patient_metadata preserved |
| 4 | Cell type proportions can be queried per patient across all modalities | VERIFIED | `celltype_proportions` iterates `mdata.mod.items()`, computes within-group proportions; `test_celltype_proportions` confirms columns and proportions sum to 1.0 |
| 5 | N-level abstraction stacking works with arbitrary grouping columns | VERIFIED | `aggregate_by_level(mdata, group_cols=["subject_id"])` and `group_cols=["sample_id"]` both produce correct groupings; `test_nlevel_aggregation` passes |
| 6 | MOFA+ embedding produces obsm['X_mofa'] in the MuData object | VERIFIED (code path confirmed; integration test skipped — mofapy2 not installed) | `MofaBackend.run` calls `mu.tl.mofa(mdata, use_obs="union", ...)` then reads `mdata.obsm["X_mofa"]`; returned array is the embedding |
| 7 | Embedding backend dispatch selects method by name string (mofa, multivi, totalvi) | VERIFIED | `get_embedding_backend(name)` returns from `_BACKENDS` dict; `list_embedding_methods()` returns `['mofa', 'multivi', 'totalvi']` confirmed by runtime check |
| 8 | Invalid method name raises ValueError with list of available methods | VERIFIED | `get_embedding_backend` raises `ValueError` with `f"Available: {known}"`; `test_invalid_method` passes |
| 9 | MultiVI/TotalVI raise clear error when modality constraints are violated | VERIFIED | Both backends check for required mod keys before any scvi import; `SCToolsDataError` raised with "requires RNA and ATAC/protein" message; both tests pass |
| 10 | sct assemble build creates h5mu from per-modality h5ad files | VERIFIED | CLI command builds `MultiOmicAtlas.from_modalities`, saves h5mu; `test_build_command` passes with exit code 0 and valid JSON output |
| 11 | sct assemble embed runs joint embedding on an assembled h5mu file | VERIFIED (code path; integration test skipped — mofapy2 not installed) | CLI loads atlas via `MultiOmicAtlas.load`, calls `atlas.embed(method=method, ...)`, saves updated atlas |
| 12 | sct assemble query returns cell type proportions as JSON | VERIFIED | CLI loads atlas, calls `atlas.celltype_proportions()`, returns `CLIResult` with `data.proportions` as dict records; `test_query_command` passes confirming non-empty proportions |

**Score:** 12/12 truths verified (3 MOFA+ integration tests skipped due to mofapy2 absence — code paths verified by inspection)

---

### Required Artifacts

| Artifact | Expected | Exists | Lines | Status | Details |
|----------|----------|--------|-------|--------|---------|
| `sc_tools/assembly/__init__.py` | Public API exports | Yes | 55 | VERIFIED | Lazy `__getattr__` exports all 7 symbols; `MultiOmicAtlas`, `build_mudata`, `join_subject_metadata`, `celltype_proportions`, `aggregate_by_level`, `get_embedding_backend`, `list_embedding_methods` |
| `sc_tools/assembly/_metadata.py` | Hierarchical metadata join with outer join on subject_id | Yes | 125 | VERIFIED | `def join_subject_metadata` present; outer join logic confirmed; validates `subject_key` presence; calls `validate_subject_metadata` per-modality |
| `sc_tools/assembly/_build.py` | MuData construction from per-modality AnnData dict | Yes | 57 | VERIFIED | `def build_mudata` present; lazy `import mudata`; stores `patient_metadata` and `sample_metadata` in `mdata.uns`; calls `mdata.update()` |
| `sc_tools/assembly/_atlas.py` | MultiOmicAtlas class wrapping MuData | Yes | 240 | VERIFIED | `class MultiOmicAtlas` present; `from_modalities`, `patient_view`, `sample_view`, `embed`, `celltype_proportions`, `save`, `load`, `modalities`, `n_obs` all implemented |
| `sc_tools/assembly/_query.py` | Cross-modal query functions | Yes | 155 | VERIFIED | `def celltype_proportions` and `def aggregate_by_level` both present and substantive |
| `sc_tools/assembly/embed/_base.py` | EmbeddingBackend Protocol and registry | Yes | 76 | VERIFIED | `EmbeddingBackend` Protocol (`runtime_checkable`), `_BACKENDS` dict, `register_embedding_backend`, `get_embedding_backend`, `list_embedding_methods` |
| `sc_tools/assembly/embed/_mofa.py` | MOFA+ embedding backend | Yes | 61 | VERIFIED | `class MofaBackend` with `run` static method; calls `mu.tl.mofa(mdata, use_obs="union", ...)`; `register_embedding_backend("mofa", MofaBackend)` at module level |
| `sc_tools/assembly/embed/_multivi.py` | MultiVI embedding backend | Yes | 75 | VERIFIED | `class MultiviBackend`; validates RNA+ATAC; raises `SCToolsDataError` with clear message |
| `sc_tools/assembly/embed/_totalvi.py` | TotalVI embedding backend | Yes | 75 | VERIFIED | `class TotalviBackend`; validates RNA+protein; raises `SCToolsDataError` with clear message |
| `sc_tools/cli/assemble.py` | sct assemble CLI command group | Yes | 199 | VERIFIED | `register_assemble(app)` function; `assemble_app = typer.Typer(...)`; three subcommands: `build`, `embed`, `query` |
| `sc_tools/tests/test_assembly.py` | Unit tests for assembly module | Yes | 257 | VERIFIED | 11 test functions; all pass |
| `sc_tools/tests/test_assembly_embed.py` | Unit tests for embedding backends | Yes | 151 | VERIFIED | 8 test functions; 6 pass, 2 skip (mofapy2 not installed — expected) |
| `sc_tools/tests/test_cli_assemble.py` | Unit tests for CLI assemble commands | Yes | 201 | VERIFIED | 6 test functions; 5 pass, 1 skip (mofapy2 not installed — expected) |

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `sc_tools/assembly/_metadata.py` | `sc_tools/qc/metadata.py` | `validate_subject_metadata` called per-modality | WIRED | Line 64: `from sc_tools.qc.metadata import validate_subject_metadata`; called at line 67 inside try/except for graceful degradation |
| `sc_tools/assembly/_atlas.py` | `sc_tools/assembly/_build.py` | `from_modalities` classmethod calls `build_mudata` | WIRED | Line 64: `from sc_tools.assembly._build import build_mudata`; line 66: `mdata = build_mudata(modalities, ...)` |
| `sc_tools/assembly/_atlas.py` | `sc_tools/assembly/_query.py` | `celltype_proportions` method delegates to query module | WIRED | Line 207: `from sc_tools.assembly._query import celltype_proportions`; line 209: `return celltype_proportions(self.mdata, ...)` |
| `sc_tools/assembly/embed/_mofa.py` | `sc_tools/assembly/embed/_base.py` | `register_embedding_backend('mofa', MofaBackend)` | WIRED | Line 11: import; line 60: registration at module load time |
| `sc_tools/cli/assemble.py` | `sc_tools/assembly/_atlas.py` | `MultiOmicAtlas.from_modalities` and `.save` for build command | WIRED | Lines 68, 75, 76, 122, 124, 172, 174 — all three subcommands import and use `MultiOmicAtlas` |
| `sc_tools/cli/__init__.py` | `sc_tools/cli/assemble.py` | `register_assemble(app)` pattern | WIRED | Line 403: `from sc_tools.cli.assemble import register_assemble`; line 404: `register_assemble(app)` |

---

### Data-Flow Trace (Level 4)

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|--------------|--------|--------------------|--------|
| `_atlas.py::celltype_proportions` | `props` DataFrame | `_query.celltype_proportions(self.mdata, ...)` iterates live `mdata.mod.items()` | Yes — groups real obs data | FLOWING |
| `_metadata.py::join_subject_metadata` | `patient_metadata` DataFrame | Builds from `adata.obs[subject_key].unique()` across all modalities | Yes — reads actual obs data | FLOWING |
| `_build.py::build_mudata` | `mdata` MuData | `mudata.MuData(modalities)` wraps input AnnData dict | Yes — wraps real data | FLOWING |
| `cli/assemble.py::assemble_query` | `props` records in JSON | `atlas.celltype_proportions()` -> `_query.celltype_proportions(mdata)` | Yes — confirmed by `test_query_command` producing non-empty proportions | FLOWING |

---

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| Module imports cleanly | `python -c "from sc_tools.assembly import MultiOmicAtlas; print('import ok')"` | `import ok` | PASS |
| Embedding registry has 3 backends | `python -c "from sc_tools.assembly.embed import list_embedding_methods; print(list_embedding_methods())"` | `['mofa', 'multivi', 'totalvi']` | PASS |
| sct assemble --help shows subcommands | `runner.invoke(app, ['assemble', '--help'])` | build, embed, query listed | PASS |
| All 22 tests pass (3 skipped for mofapy2) | `python -m pytest sc_tools/tests/test_assembly*.py sc_tools/tests/test_cli_assemble.py -v` | 22 passed, 3 skipped | PASS |
| sct assemble build produces h5mu | `test_build_command` via typer.testing.CliRunner | exit code 0, h5mu file created, status=success | PASS |
| sct assemble query returns proportions | `test_query_command` via typer.testing.CliRunner | exit code 0, non-empty proportions list in JSON | PASS |

---

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| MOM-01 | 08-01-PLAN | Patient/subject metadata join across modalities using validated subject_id | SATISFIED | `join_subject_metadata` outer join; validates subject_key; 3 tests cover normal, zero-overlap, missing-key paths |
| MOM-02 | 08-01-PLAN | MuData assembly from independently processed per-modality AnnData objects | SATISFIED | `build_mudata` and `MultiOmicAtlas.from_modalities`; 2 tests confirm correct shape and separate feature spaces |
| MOM-03 | 08-02-PLAN | Joint embedding via multi-modal integration method (MOFA+, multiVI, or totalVI) | SATISFIED | Three backends with pluggable registry; embed method on atlas; 5 tests cover dispatch, validation, protocol conformance (2 integration tests skipped pending mofapy2) |
| MOM-04 | 08-01-PLAN, 08-02-PLAN | Cross-modal queries at patient/project level (cell type proportions for patient X across all modalities) | SATISFIED | `celltype_proportions` + `aggregate_by_level` in library; `sct assemble query` in CLI with optional `--patient` filter; 4 tests confirm behavior |

No orphaned requirements — all four MOM requirements are claimed by plans and verified.

---

### Anti-Patterns Found

No anti-patterns detected.

- No TODO/FIXME/PLACEHOLDER comments in any assembly or CLI files
- No empty implementations (`return null`, `return {}`, `return []` without logic)
- All mudata/muon imports are inside function bodies (lazy per CLI-08): confirmed in `_build.py:37`, `_atlas.py:115/144/236`, `_mofa.py:41`, `_multivi.py:56`, `_totalvi.py:54`
- No hardcoded empty props passed to rendering functions
- No console.log-only stubs

---

### Human Verification Required

#### 1. MOFA+ Integration Test with mofapy2 Installed

**Test:** In an environment with `mofapy2>=0.7.0` installed, run:
```
python -m pytest sc_tools/tests/test_assembly_embed.py::TestMofaEmbedding::test_mofa_embedding -v
python -m pytest sc_tools/tests/test_cli_assemble.py::TestAssembleEmbed::test_embed_command -v
```
**Expected:** Both tests pass. `mdata.obsm["X_mofa"]` has shape `(n_cells, n_factors)`. CLI returns `status=success` with `embedding_shape` in data.
**Why human:** `mofapy2` is not installed in the local dev environment. The code path is fully implemented and structurally correct, but the MOFA+ end-to-end call via `muon.tl.mofa` cannot be exercised without the package.

---

### Gaps Summary

No gaps. All automated checks pass. The phase goal is fully achieved:

- The `sc_tools.assembly` module provides `MultiOmicAtlas` wrapping `MuData` with patient-centric access.
- Outer join on `subject_id` includes all patients even when absent from some modalities (boolean `has_{mod}` flags in `patient_metadata`).
- `build_mudata` constructs modality-keyed MuData with separate feature spaces and metadata in `uns`.
- Cross-modal queries (`celltype_proportions`, `aggregate_by_level`) iterate live modality obs data.
- Three joint embedding backends (MOFA+, MultiVI, TotalVI) registered via pluggable protocol.
- `sct assemble build | embed | query` CLI commands follow the `CLIResult` envelope pattern with lazy imports.
- 22/25 tests pass; 3 skipped are guarded by `pytest.importorskip("mofapy2")` which is by design.

The only item requiring human follow-up is MOFA+ end-to-end integration testing in an environment with `mofapy2` installed.

---

_Verified: 2026-03-24T23:10:00Z_
_Verifier: Claude (gsd-verifier)_
