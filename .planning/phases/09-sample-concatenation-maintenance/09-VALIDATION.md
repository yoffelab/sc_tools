---
phase: 9
slug: sample-concatenation-maintenance
status: draft
nyquist_compliant: true
wave_0_complete: false
created: 2026-03-27
---

# Phase 9 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest 7.x |
| **Config file** | `pyproject.toml` |
| **Quick run command** | `pytest tests/ -x -q --tb=short -k "concat or maint"` |
| **Full suite command** | `make test` |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run `pytest tests/ -x -q --tb=short -k "concat or maint"`
- **After every plan wave:** Run `make test`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 9-01-T1 | 01 | 1 | MAINT-01 | unit | `python -c "import plotly"` | ✅ | ⬜ pending |
| 9-01-T2 | 01 | 1 | MAINT-02 | unit | `pytest tests/test_qc.py -x -q -k "plotly"` | ✅ | ⬜ pending |
| 9-02-T1 | 02 | 1 | CONCAT-01..04 | unit | `pytest tests/test_concat.py -x -q` | ❌ W0 | ⬜ pending |
| 9-02-T2 | 02 | 1 | CONCAT-01,02,03 | unit | `pytest tests/test_concat.py -x -q -k "TestConcatCommand or TestSpatialPreservation or TestConcatProvenance"` | ❌ W0 | ⬜ pending |
| 9-02-T3 | 02 | 1 | CONCAT-04 | integration | `pytest tests/test_pipeline.py -x -q -k "concat"` | ✅ | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `tests/test_concat.py` — stubs for CONCAT-01, CONCAT-02, CONCAT-03, CONCAT-04
- [ ] Fixtures in `tests/conftest.py` — minimal multi-sample AnnData with spatial uns
- [ ] RED state confirmed — pytest on test_concat.py exits non-zero before implementation

*Existing infrastructure (pytest, conftest.py) covers maintenance tests — only concat tests need new stubs.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| `sct concat` CLI end-to-end with real .h5ad files | CONCAT-01 | Requires real AnnData files on disk | Run `sct concat --input a.h5ad b.h5ad --output merged.h5ad`; inspect `merged.h5ad` with Python |
| Report HTML renders Plotly charts correctly | MAINT-02 | Browser rendering not testable in CI | Open a generated QC report HTML in browser; verify plots load without console errors |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
