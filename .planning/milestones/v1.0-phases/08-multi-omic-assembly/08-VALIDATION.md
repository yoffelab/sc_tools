---
phase: 8
slug: multi-omic-assembly
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-24
---

# Phase 8 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest 7.x |
| **Config file** | pyproject.toml `[tool.pytest.ini_options]` |
| **Quick run command** | `python -m pytest sc_tools/tests/test_assembly.py -x -q` |
| **Full suite command** | `python -m pytest sc_tools/tests/test_assembly.py sc_tools/tests/test_cli_assemble.py -x -q` |
| **Estimated runtime** | ~20 seconds |

---

## Sampling Rate

- **After every task commit:** Run `python -m pytest sc_tools/tests/test_assembly.py -x -q`
- **After every plan wave:** Run full suite command
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 20 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 08-01-01 | 01 | 1 | MOM-01, MOM-02 | unit | `pytest sc_tools/tests/test_assembly.py::test_metadata_join -x` | ❌ W0 | ⬜ pending |
| 08-01-02 | 01 | 1 | MOM-02 | unit | `pytest sc_tools/tests/test_assembly.py::test_mudata_build -x` | ❌ W0 | ⬜ pending |
| 08-02-01 | 02 | 2 | MOM-03 | unit | `pytest sc_tools/tests/test_assembly.py::test_joint_embedding -x` | ❌ W0 | ⬜ pending |
| 08-02-02 | 02 | 2 | MOM-04 | unit | `pytest sc_tools/tests/test_cli_assemble.py -x` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `sc_tools/tests/test_assembly.py` — stubs for MOM-01, MOM-02, MOM-03 (metadata join, MuData build, joint embedding)
- [ ] `sc_tools/tests/test_cli_assemble.py` — stubs for MOM-04 (CLI commands)
- [ ] `sc_tools/tests/conftest.py` — multi-modality fixtures (synthetic RNA + IMC AnnData objects with shared subject_id)

*Existing pytest infrastructure covers test framework needs.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| MuData with 4 real modalities round-trips to h5mu | MOM-02 | Requires real multi-modal project data | Build MuData from real project h5ad files, write to h5mu, reload, verify all modalities present |
| MOFA+ convergence on disparate-size modalities | MOM-03 | Requires GPU and real-scale data | Run MOFA+ with RNA (500K cells) + IMC (10K cells) on HPC, check factor loadings |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 20s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
