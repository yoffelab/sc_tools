---
phase: 05-provenance
verified: 2026-03-23T20:00:00Z
status: passed
score: 9/9 must-haves verified
re_verification: false
---

# Phase 5: Provenance Verification Report

**Phase Goal:** Every CLI output has a traceable lineage back to its inputs, parameters, and environment
**Verified:** 2026-03-23T20:00:00Z
**Status:** PASSED
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Every successful CLI command with artifacts writes a .provenance.json sidecar file | VERIFIED | `_write_provenance_sidecars` called in `cli_handler` on `status==success and artifacts` (cli/__init__.py:239). Test `test_cli_handler_writes_sidecar_on_success` passes. |
| 2 | Sidecar contains command, params, inputs with SHA256, version, timestamp, runtime_s, peak_memory_mb | VERIFIED | `ProvenanceRecord` model has all 7 fields (models/result.py:64-80). `test_write_sidecar_creates_json` validates all fields in the written JSON. |
| 3 | No sidecar is written when command errors or has no artifacts | VERIFIED | Guard `if result.status == Status.success and result.artifacts:` (cli/__init__.py:239). Tests `test_cli_handler_no_sidecar_on_error` and `test_cli_handler_no_sidecar_when_no_artifacts` both pass. |
| 4 | When the output artifact is an h5ad file, provenance is also embedded in adata.uns['sct_provenance'] | VERIFIED | `embed_provenance_in_adata` called via h5py append mode (sidecar.py:160-196). Called from `_write_provenance_sidecars` when `artifact.endswith(".h5ad")` (cli/__init__.py:214-215). Test `test_write_provenance_sidecars_calls_embed_for_h5ad` passes. |
| 5 | Leiden clustering produces identical results for identical random_state + resolution | VERIFIED | `random_state` threaded through `_leiden_cluster`, `_ari`, `_nmi`, `compute_integration_metrics`, `compare_integrations` (integration.py), and `cluster()` (reduce.py:160, 185). Default is 0. |
| 6 | sct provenance show <file> displays the provenance record from sidecar or adata.uns fallback | VERIFIED | `show` command in cli/provenance.py reads sidecar then falls back to `_read_uns_provenance` for h5ad. Tests `test_show_returns_provenance_from_sidecar` and `test_show_uns_fallback` pass. |
| 7 | sct provenance trace <file> walks input references recursively to build full lineage | VERIFIED | BFS engine in `trace_lineage` (provenance/trace.py:30-156). 3-step chain test passes, producing chronological output with correct commands. |
| 8 | Trace stops at raw data origins (no provenance) and handles missing intermediates gracefully | VERIFIED | Origin nodes marked with `"note": "origin (no provenance)"`. Missing file handled via SHA256 relocation or direct origin marking. Tests `test_missing_input_marked_as_origin` and `TestSHA256Relocation` pass. |
| 9 | h5ad adata.uns fallback reads provenance without loading full AnnData | VERIFIED | `_read_uns_provenance` uses `h5py.File(path, "r")` directly (trace.py:159-189), no AnnData import. `_read_uns_provenance` also used inside `embed_provenance_in_adata` write path (sidecar.py:185). |

**Score:** 9/9 truths verified

---

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `sc_tools/models/result.py` | InputFile and ProvenanceRecord Pydantic models | VERIFIED | `class InputFile` at line 51, `class ProvenanceRecord` at line 64. Both have full field sets per D-06 and D-08. |
| `sc_tools/provenance/__init__.py` | Package init | VERIFIED | Exists, one-line docstring. |
| `sc_tools/provenance/checksum.py` | SHA256 streaming file hasher | VERIFIED | `sha256_file` implemented with 64KB streaming reads (line 9-27). |
| `sc_tools/provenance/sidecar.py` | Sidecar utilities and h5ad uns embedding | VERIFIED | `write_sidecar`, `read_sidecar`, `sidecar_path_for`, `get_peak_memory_mb`, `build_provenance_record`, `embed_provenance_in_adata` all present and substantive. |
| `sc_tools/cli/__init__.py` | Automatic sidecar writing via _write_provenance_sidecars | VERIFIED | `_write_provenance_sidecars` defined (line 182) and called (line 244). `time.monotonic` timing at line 227. `_input_files` popped before emit (line 234). |
| `sc_tools/provenance/trace.py` | Lineage trace walker with cycle detection and SHA256 relocation | VERIFIED | `trace_lineage`, `_read_uns_provenance`, `_find_by_sha256` all present. `visited` set for cycle detection. BFS with deque. |
| `sc_tools/cli/provenance.py` | show and trace CLI commands via register_provenance(app) | VERIFIED | `register_provenance` at line 20. Two commands registered: `show` and `trace`. `@cli_handler` on both. |
| `sc_tools/tests/test_provenance.py` | Unit tests for sidecar, checksum, h5ad uns, Leiden reproducibility | VERIFIED | 20 tests — all pass. Covers all behaviors in the plan. |
| `sc_tools/tests/test_cli_provenance.py` | Integration tests for provenance show and trace commands | VERIFIED | 15 tests — all pass. Covers show, show-no-prov, uns-fallback, trace chain, trace project-dir, help text. |

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `sc_tools/cli/__init__.py` | `sc_tools/provenance/sidecar.py` | `_write_provenance_sidecars` called after `_emit` on success | WIRED | Lines 239-246 in cli/__init__.py guard + call confirmed. |
| `sc_tools/cli/__init__.py` | `sc_tools/provenance/sidecar.py` | `embed_provenance_in_adata` called for h5ad artifacts | WIRED | Line 215 in `_write_provenance_sidecars`. Lazy import inside the function. |
| `sc_tools/provenance/sidecar.py` | `sc_tools/provenance/checksum.py` | `write_sidecar` uses `sha256_file` via `build_provenance_record` | WIRED | `build_provenance_record` imports and calls `sha256_file` (sidecar.py:123). |
| `sc_tools/bm/integration.py` | `scanpy.tl.leiden` | `_leiden_cluster` passes `random_state=random_state` | WIRED | Line 232: `sc.tl.leiden(tmp, resolution=resolution, key_added="leiden", random_state=random_state)`. |
| `sc_tools/cli/provenance.py` | `sc_tools/provenance/sidecar.py` | `show` command reads sidecar via `read_sidecar()` | WIRED | Line 38 in cli/provenance.py: `prov = read_sidecar(file)`. |
| `sc_tools/cli/provenance.py` | `sc_tools/provenance/trace.py` | `trace` command calls `trace_lineage()` | WIRED | Line 77 in cli/provenance.py: `steps = trace_lineage(file, project_dir=project_dir)`. |
| `sc_tools/provenance/trace.py` | `sc_tools/provenance/sidecar.py` | trace walker reads sidecars via `read_sidecar` | WIRED | Line 52 in trace.py: `from sc_tools.provenance.sidecar import read_sidecar`. Used in main BFS loop. |
| `sc_tools/cli/__init__.py` | `sc_tools/cli/provenance.py` | `register_provenance(app)` called at bottom | WIRED | Lines 339-340 in cli/__init__.py: `from sc_tools.cli.provenance import register_provenance` + call. |

---

### Data-Flow Trace (Level 4)

Provenance phase delivers infrastructure utilities (no dynamic-data rendering components). The CLI commands return structured provenance data read from disk — data-flow is verified through test execution. Level 4 is satisfied by test passes rather than JSX/API trace.

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|---------------|--------|--------------------|--------|
| `cli/provenance.py show` | `prov` dict | `read_sidecar()` reads real .provenance.json from disk, or h5py reads `uns/sct_provenance` | Yes — 15 integration tests confirm real JSON returned | FLOWING |
| `cli/provenance.py trace` | `steps` list | `trace_lineage()` BFS reads real sidecars from disk | Yes — chain tests confirm real steps returned | FLOWING |
| `cli/__init__.py _write_provenance_sidecars` | sidecar file | `build_provenance_record` computes real SHA256 checksums from input files | Yes — `test_input_files_convention_key` and `test_sidecar_params_include_kwargs` confirm real data written | FLOWING |

---

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| test_provenance.py — 20 unit tests | `python -m pytest sc_tools/tests/test_provenance.py -x -q` | 20 passed in 10.38s | PASS |
| test_cli_provenance.py — 15 integration tests | `python -m pytest sc_tools/tests/test_cli_provenance.py -x -q` | 15 passed in 0.28s | PASS |
| Module imports without errors | `python -c "from sc_tools.cli import app; from sc_tools.provenance.trace import trace_lineage"` | Exit 0 | PASS |

---

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| PRV-01 | 05-01-PLAN.md | JSON sidecar `.provenance.json` written alongside every CLI output file | SATISFIED | `_write_provenance_sidecars` in `cli_handler`. 3 tests confirm sidecar presence/absence logic. |
| PRV-02 | 05-01-PLAN.md | Sidecar includes: command, params, input files with SHA256, sc_tools version, timestamp, runtime_s, peak_memory_mb | SATISFIED | `ProvenanceRecord` model has all 7 fields. `test_write_sidecar_creates_json` validates all fields in written JSON. |
| PRV-03 | 05-02-PLAN.md | `sct provenance show <file>` — display provenance for a single output | SATISFIED | `show` command in cli/provenance.py. 3 CLI tests cover show, no-provenance error, and uns fallback. |
| PRV-04 | 05-02-PLAN.md | `sct provenance trace <file>` — trace full lineage DAG via input file references | SATISFIED | `trace` command calls `trace_lineage`. 9 trace unit tests + 2 CLI trace tests cover all specified behaviors. |
| PRV-05 | 05-01-PLAN.md | Reproducible Leiden clustering — configurable resolution and random_state propagation | SATISFIED | `random_state` parameter threaded through `_leiden_cluster`, `_ari`, `_nmi`, `compute_integration_metrics`, `compare_integrations` (integration.py), and `cluster()` (reduce.py). Default=0. |

No orphaned requirements found. All 5 PRV requirements in REQUIREMENTS.md are marked Complete and covered by plans 05-01 and 05-02.

---

### Anti-Patterns Found

| File | Pattern | Severity | Impact |
|------|---------|----------|--------|
| `sc_tools/tests/test_bm_report.py` | Pre-existing test failure (Plotly assertion) — documented in both SUMMARYs as pre-existing, unrelated to phase 5 | INFO | Zero impact on phase 5 goal. Not introduced by this phase. |

No stubs, placeholders, empty implementations, or hollow props found in any phase 5 files. Both SUMMARYs explicitly state "Known Stubs: None".

---

### Human Verification Required

None. All goal truths are verifiable programmatically and confirmed by passing tests.

The one item that could warrant optional human inspection:

#### 1. End-to-End Sidecar Round-Trip on Real Data

**Test:** Run a real `sct preprocess run` or `sct qc run` command on a real h5ad file with `_input_files` populated, then inspect the resulting `.provenance.json` sidecar and `adata.uns['sct_provenance']`.
**Expected:** Sidecar contains correct SHA256 for the input file, correct runtime_s, and the h5ad embeds the same provenance accessible via `sc.read_h5ad(output)`.uns['sct_provenance'].
**Why human (optional):** Unit tests mock the sidecar writing path in some tests. Real data integration with actual input SHA256 computation would confirm end-to-end. However, `test_cli_handler_writes_sidecar_on_success` and `test_input_files_convention_key` confirm the core path non-mocked.

---

### Commits Verified

All 6 commits documented in SUMMARYs confirmed present in git log:

| Commit | Description |
|--------|-------------|
| `a0e265c` | test(05-01): failing tests RED phase |
| `8c2dbbe` | feat(05-01): implementation GREEN phase |
| `49a9e9c` | feat(05-01): random_state threading |
| `792d117` | test(05-02): trace engine RED phase |
| `2877eb1` | feat(05-02): trace engine GREEN phase |
| `af93bc7` | feat(05-02): CLI commands + registration |

---

### Gaps Summary

No gaps. All 9 observable truths verified. All 9 required artifacts pass all 4 levels (exists, substantive, wired, data-flowing). All 8 key links confirmed wired. All 5 PRV requirements satisfied. 35 tests pass (20 unit + 15 integration).

The single pre-existing test failure in `test_bm_report.py` (Plotly assertion) is unrelated to phase 5 and was present before this phase began.

---

_Verified: 2026-03-23T20:00:00Z_
_Verifier: Claude (gsd-verifier)_
