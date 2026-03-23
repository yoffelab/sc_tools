---
phase: 03-core-commands
verified: 2026-03-22T00:00:00Z
status: passed
score: 13/13 must-haves verified
re_verification: false
---

# Phase 3: Core Commands Verification Report

**Phase Goal:** Agents can run QC, preprocessing, validation, benchmarking, status checks, and report generation through `sct` commands instead of ad-hoc scripts
**Verified:** 2026-03-22
**Status:** PASSED
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | `sct qc run <file>` produces JSON metrics summary with per-sample QC | VERIFIED | `qc_run()` wired to `compute_sample_metrics` + `classify_samples`; tests pass with 2 samples, 100 cells |
| 2 | `sct preprocess run <file>` dispatches modality-aware recipe and writes output h5ad | VERIFIED | `preprocess_run()` wired to `recipes.preprocess`; `_detect_modality()` reads `adata.uns` or panel size |
| 3 | `sct benchmark integration --from-dir <dir>` loads embeddings, computes ranked JSON | VERIFIED | `benchmark_integration()` wired to `compare_integrations(embedding_files=...)`; subsampling logged |
| 4 | `sct validate <phase> <file>` returns pass/fail JSON with issues list | VERIFIED | `validate_run()` wired to `validate_file()`; exit 0 on pass, exit 2 on failure |
| 5 | `sct status` shows pipeline DAG even when registry unavailable | VERIFIED | `status_show()` wired to `get_dag()` + `get_available_next()`; `registry_available: false` path confirmed |
| 6 | `sct report <type> --adata <file>` generates HTML report | VERIFIED | `report_generate()` wired to all four `generate_*_report` functions in `sc_tools/qc/report.py` |
| 7 | CLI commands and MCP tools share the same CLIResult type | VERIFIED | Both return identical field sets: `{status, command, data, artifacts, provenance, message, error, partial_failures}` |
| 8 | `_check_deps()` raises `SCToolsUserError` with install instructions for missing deps | VERIFIED | Test confirmed: `_check_deps(["nonexistent_xyz_pkg"])` raises `SCToolsUserError` with "Missing required dependencies" |
| 9 | `conftest.py` provides 100-cell AnnData fixture with spatial, obs, obsm | VERIFIED | `adata_100` has 100 obs, 200 vars, 3 MT genes, library_id, sample, spatial, X_pca, X_scvi, modality |
| 10 | CLI integration tests pass with 100-cell AnnData fixtures | VERIFIED | `test_cli_commands.py`: 23 passed |
| 11 | E2E test file exists with skipif guard for HPC/real data | VERIFIED | `test_cli_e2e.py` uses `pytest.mark.skipif` on `SCT_TEST_DATA_DIR` and `SCT_TEST_EMBEDDING_DIR`; 3 skipped in CI |
| 12 | All 24 Phase 2 regression tests pass | VERIFIED | `test_cli.py`: 24 passed |
| 13 | All 8 command groups registered and reachable via `sct <cmd> --help` | VERIFIED | qc, preprocess, validate, status, benchmark, report, integrate, celltype all exit 0 on `--help` |

**Score:** 13/13 truths verified

---

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `sc_tools/cli/__init__.py` | Typer app, `cli_handler`, `_emit`, `_render_rich`, `_check_deps`, sub-app registrations | VERIFIED | 273 lines; all helpers defined; all 8 sub-apps registered via `app.add_typer()` |
| `sc_tools/cli/qc.py` | `qc_app`, `qc_run`, `report_app`, `report_generate` | VERIFIED | 201 lines; both apps and both commands implemented |
| `sc_tools/cli/validate.py` | `validate_app`, `validate_run` | VERIFIED | 55 lines; `validate_run` implemented, wired to `validate_file` |
| `sc_tools/cli/status.py` | `status_app`, `status_show` | VERIFIED | 77 lines; `status_show` with `@cli_handler` and `registry_available` fallback |
| `sc_tools/cli/preprocess.py` | `preprocess_app`, `preprocess_run`, `_detect_modality` | VERIFIED | 159 lines; full implementation with D-08 auto-detection |
| `sc_tools/cli/benchmark.py` | `benchmark_app`, `benchmark_integration` | VERIFIED | 184 lines; `--report` flag wired to `generate_benchmark_report(benchmark_params=...)` |
| `sc_tools/tests/conftest.py` | `adata_100`, `adata_100_h5ad`, `adata_100_preprocess_checkpoint` | VERIFIED | All 3 fixtures present with correct schema |
| `sc_tools/tests/test_cli_commands.py` | Integration tests for all commands | VERIFIED | 275 lines; 7 test classes, 23 tests |
| `sc_tools/tests/test_cli_e2e.py` | E2E scaffold with skipif guards | VERIFIED | 69 lines; 2 test classes guarded by env var checks |

Old `sc_tools/cli.py` monolith: ABSENT (correctly replaced by `sc_tools/cli/` package)

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `cli/__init__.py` | `cli/qc.py` | `from sc_tools.cli.qc import qc_app, report_app` | WIRED | Line 259 |
| `cli/__init__.py` | `cli/preprocess.py` | `from sc_tools.cli.preprocess import preprocess_app` | WIRED | Line 260 |
| `cli/__init__.py` | `cli/validate.py` | `from sc_tools.cli.validate import validate_app` | WIRED | Line 261 |
| `cli/__init__.py` | `cli/benchmark.py` | `from sc_tools.cli.benchmark import benchmark_app` | WIRED | Line 262 |
| `cli/__init__.py` | `cli/status.py` | `from sc_tools.cli.status import status_app` | WIRED | Line 263 |
| `cli/validate.py` | `sc_tools/validate.py` | `from sc_tools.validate import validate_file` | WIRED | Line 33 (lazy import) |
| `cli/status.py` | `sc_tools/pipeline.py` | `from sc_tools.pipeline import get_dag, get_available_next, tuple_to_display` | WIRED | Line 27 (lazy import) |
| `cli/qc.py` | `sc_tools/qc/sample_qc.py` | `from sc_tools.qc.sample_qc import compute_sample_metrics, classify_samples` | WIRED | Line 41 (lazy import) |
| `cli/preprocess.py` | `sc_tools/pp/recipes.py` | `from sc_tools.pp.recipes import preprocess` | WIRED | Line 117 (lazy import) |
| `cli/benchmark.py` | `sc_tools/bm/integration.py` | `from sc_tools.bm.integration import compare_integrations` | WIRED | Line 92 (lazy import) |
| `cli/benchmark.py` | `sc_tools/bm/report.py` | `from sc_tools.bm.report import generate_benchmark_report` | WIRED | Line 135 (lazy import, conditional) |
| `test_cli_commands.py` | `cli/__init__.py` | `from sc_tools.cli import app` | WIRED | Line 16 |
| `mcp/tools_server.py` | `sc_tools/validate.py` | `from sc_tools.validate import validate_file` | WIRED | Line 64 |

---

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| CMD-01 | 03-02 | `sct qc run` — QC metrics, JSON output | SATISFIED | `qc_run()` implemented and tested; `test_qc_run_produces_metrics` passes |
| CMD-02 | 03-03 | `sct preprocess run` — modality-aware recipe dispatch | SATISFIED | `preprocess_run()` wired to `recipes.preprocess`; `_detect_modality()` tested |
| CMD-03 | 03-02 | `sct validate <phase> <file>` — pass/fail JSON | SATISFIED | `validate_run()` wired to `validate_file()`; issues list + n_issues in output |
| CMD-04 | 03-03 | `sct benchmark integration --from-dir` — ranked JSON | SATISFIED | `benchmark_integration()` wired to `compare_integrations(embedding_files=...)` |
| CMD-05 | 03-02 | `sct status` — pipeline DAG from registry | SATISFIED | `status_show()` returns 10-phase DAG; graceful when registry unavailable (D-09) |
| CMD-06 | 03-02 | `sct report <type>` — HTML report generation | SATISFIED | `report_generate()` wired to all four `generate_*_report` backends |
| CMD-07 | 03-04 | Shared CLIResult type for CLI + MCP | SATISFIED | MCP `validate_checkpoint` returns same 8-field CLIResult; field equality verified in test |
| CMD-08 | 03-01 | Fast-fail `_check_deps()` with install instructions | SATISFIED | `_check_deps()` in `cli/__init__.py`; raises `SCToolsUserError`; all commands call it |
| TST-05 | 03-01, 03-04 | CLI integration tests with 100-cell fixtures | SATISFIED | `conftest.py` + `test_cli_commands.py`; 23 tests across 7 classes all passing |
| TST-06 | 03-04 | E2E test scaffold with skipif guard | SATISFIED | `test_cli_e2e.py` guards on `SCT_TEST_DATA_DIR` + `SCT_TEST_EMBEDDING_DIR`; 3 skipped in CI |

No orphaned requirements found. All 10 Phase 3 requirements are claimed by plans and verified in code.

---

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `test_cli_commands.py` | 18 | `CliRunner()` instead of `CliRunner(mix_stderr=False)` as specified in plan | Info | Tests still pass because happy-path JSON tests use adata_100 fixtures that produce no stderr; error-path tests only check exit codes |

No blocker or warning-level anti-patterns found. No TODO/FIXME/placeholder comments. No stub return values (all commands have substantive implementations wired to backend functions).

---

### Human Verification Required

None. All Success Criteria from ROADMAP.md Phase 3 are verifiable programmatically and confirmed:

1. `sct qc run <file>` + `sct preprocess run <file>` — confirmed by unit tests
2. `sct benchmark integration --from-dir <dir>` — confirmed by option/error tests
3. `sct validate <phase> <file>` — confirmed; exit 0 on pass, exit 2 on data error
4. `sct status` — confirmed; 10 phases in DAG, 1 available_next, registry_available=false handled
5. CLI + MCP shared Result type — confirmed by `TestSharedResult.test_cli_and_mcp_same_result_type`

---

## Test Run Summary

```
test_cli.py:              24 passed   (Phase 2 regression)
test_cli_commands.py:     23 passed
test_cli_e2e.py:          0 passed, 3 skipped (skipif guards correct)
Total:                    47 passed, 3 skipped, 0 failed
```

---

_Verified: 2026-03-22_
_Verifier: Claude (gsd-verifier)_
