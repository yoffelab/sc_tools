---
phase: 6
slug: scientific-gaps
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-24
---

# Phase 6 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest 7.x |
| **Config file** | pyproject.toml `[tool.pytest]` |
| **Quick run command** | `python -m pytest sc_tools/tests/test_subject_metadata.py sc_tools/tests/test_panel_dispatch.py sc_tools/tests/test_de.py -x -q` |
| **Full suite command** | `python -m pytest sc_tools/tests/test_subject_metadata.py sc_tools/tests/test_panel_dispatch.py sc_tools/tests/test_de.py sc_tools/tests/test_marker_validation.py -v` |
| **Estimated runtime** | ~15 seconds |

---

## Sampling Rate

- **After every task commit:** Run quick command for the relevant test file
- **After every plan wave:** Run full suite
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 20 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 06-01-01 | 01 | 1 | SCI-03 | unit | `pytest sc_tools/tests/test_subject_metadata.py -x` | W0 | pending |
| 06-01-02 | 01 | 1 | SCI-04 | unit | `pytest sc_tools/tests/test_panel_dispatch.py -x` | W0 | pending |
| 06-02-01 | 02 | 2 | SCI-01 | unit | `pytest sc_tools/tests/test_de.py -x` | W0 | pending |
| 06-02-02 | 02 | 2 | SCI-01 | unit | `pytest sc_tools/tests/test_de.py -x` | W0 | pending |
| 06-03-01 | 03 | 2 | SCI-02 | unit | `pytest sc_tools/tests/test_marker_validation.py -x` | W0 | pending |
| 06-03-02 | 03 | 2 | SCI-02 | unit | `pytest sc_tools/tests/test_marker_validation.py -x` | W0 | pending |

*Status: pending / green / red / flaky*

---

## Wave 0 Requirements

- [ ] `sc_tools/tests/test_subject_metadata.py` — test stubs for subject_id validation, confounding checks
- [ ] `sc_tools/tests/test_panel_dispatch.py` — test stubs for panel detection, method restriction
- [ ] `sc_tools/tests/test_de.py` — test stubs for pseudobulk DE aggregation, PyDESeq2 wrapper, CLI command
- [ ] `sc_tools/tests/test_marker_validation.py` — test stubs for marker validation compute and report integration

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| HTML marker report visual quality | SCI-02 | Visual inspection of dotplot rendering | Open report in browser, verify dotplot shows markers per type |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 20s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
