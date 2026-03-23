---
phase: 5
slug: provenance
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-23
---

# Phase 5 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest 7.x |
| **Config file** | pyproject.toml `[tool.pytest]` |
| **Quick run command** | `python -m pytest sc_tools/tests/test_provenance.py -x -q` |
| **Full suite command** | `python -m pytest sc_tools/tests/test_provenance.py sc_tools/tests/test_cli_provenance.py -v` |
| **Estimated runtime** | ~10 seconds |

---

## Sampling Rate

- **After every task commit:** Run `python -m pytest sc_tools/tests/test_provenance.py -x -q`
- **After every plan wave:** Run full suite command
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 15 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 05-01-01 | 01 | 1 | PRV-01, PRV-02 | unit | `pytest sc_tools/tests/test_provenance.py -x` | ❌ W0 | ⬜ pending |
| 05-01-02 | 01 | 1 | PRV-01, PRV-02 | unit | `pytest sc_tools/tests/test_provenance.py -x` | ❌ W0 | ⬜ pending |
| 05-02-01 | 02 | 2 | PRV-03, PRV-04 | unit | `pytest sc_tools/tests/test_cli_provenance.py -x` | ❌ W0 | ⬜ pending |
| 05-02-02 | 02 | 2 | PRV-03, PRV-04 | unit | `pytest sc_tools/tests/test_cli_provenance.py -x` | ❌ W0 | ⬜ pending |
| 05-02-03 | 02 | 2 | PRV-05 | unit | `pytest sc_tools/tests/test_provenance.py::test_leiden_reproducibility -x` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `sc_tools/tests/test_provenance.py` — test stubs for provenance model, sidecar writing, SHA256, adata.uns embedding
- [ ] `sc_tools/tests/test_cli_provenance.py` — test stubs for `sct provenance show` and `sct provenance trace` commands

*Existing pytest infrastructure and conftest.py from Phase 3 cover shared fixtures.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Sidecar survives file move + trace recovers via SHA256 | PRV-01, D-10 | Requires physical file relocation | Move an h5ad file, run `sct provenance trace`, verify SHA256 fallback works |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 15s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
