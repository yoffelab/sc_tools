---
phase: 1
slug: benchmark-fixes
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-20
---

# Phase 1 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest >=7.0 |
| **Config file** | pyproject.toml `[tool.pytest.ini_options]` |
| **Quick run command** | `pytest sc_tools/tests/test_integration_benchmark.py -x -v` |
| **Full suite command** | `pytest sc_tools/tests/ -v --tb=short -q` |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run `pytest sc_tools/tests/test_integration_benchmark.py -x -v`
- **After every plan wave:** Run `pytest sc_tools/tests/ -v --tb=short -q`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 01-01-01 | 01 | 1 | BM-01 | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestH5pyLoading -x` | ❌ W0 | ⬜ pending |
| 01-01-02 | 01 | 1 | BM-02 | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestNaNHandling -x` | ❌ W0 | ⬜ pending |
| 01-01-03 | 01 | 1 | BM-03 | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestCompareIntegrations::test_subsample_n -x` | ❌ W0 | ⬜ pending |
| 01-01-04 | 01 | 1 | BM-04 | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestStratifiedSubsample::test_proportional -x` | ❌ W0 | ⬜ pending |
| 01-01-05 | 01 | 1 | BM-05 | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestRunIntegrationBenchmark::test_runtime_column -x` | ❌ W0 | ⬜ pending |
| 01-01-06 | 01 | 1 | BM-06 | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestCompareIntegrations::test_params_in_attrs -x` | ❌ W0 | ⬜ pending |
| 01-01-07 | 01 | 1 | BM-07 | unit | `pytest sc_tools/tests/test_pp.py::TestTargetedPanelScVI -x` | ❌ W0 | ⬜ pending |
| 01-02-01 | 02 | 1 | TST-01 | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestComputeMetricsEdgeCases -x` | ❌ W0 | ⬜ pending |
| 01-02-02 | 02 | 1 | TST-02 | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestCompareIntegrationsEdgeCases -x` | ❌ W0 | ⬜ pending |
| 01-02-03 | 02 | 1 | TST-03 | unit | `pytest sc_tools/tests/test_integration_benchmark.py::TestStratifiedSubsampleEdgeCases -x` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `sc_tools/tests/test_integration_benchmark.py` — new test classes: TestH5pyLoading, TestNaNHandling, TestComputeMetricsEdgeCases, TestCompareIntegrationsEdgeCases, TestStratifiedSubsampleEdgeCases
- [ ] `sc_tools/tests/test_pp.py` — new test class: TestTargetedPanelScVI
- [ ] NaN embedding fixture helper (`_make_nan_embedding_adata`)
- [ ] h5ad fixture helper (write temp h5ad via anndata, return path for h5py loading tests)
- [ ] No new framework install needed — pytest already configured

*Existing infrastructure covers framework and conftest. New fixtures and test classes are needed.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Peak memory <2GB for 2.5M cells | BM-01 | Requires 2.5M-cell h5ad file (not in test suite) | Run `scripts/run_integration_benchmark.py` on real project data, monitor RSS via `resource.getrusage` |

*All other behaviors have automated verification.*

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
