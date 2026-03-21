---
phase: 02
slug: cli-foundation
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-21
---

# Phase 02 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest 8.x |
| **Config file** | `pyproject.toml` [tool.pytest.ini_options] |
| **Quick run command** | `python -m pytest sc_tools/tests/test_cli.py -q` |
| **Full suite command** | `python -m pytest sc_tools/tests/ -q` |
| **Estimated runtime** | ~15 seconds |

---

## Sampling Rate

- **After every task commit:** Run `python -m pytest sc_tools/tests/test_cli.py -q`
- **After every plan wave:** Run `python -m pytest sc_tools/tests/ -q`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 15 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 02-01-01 | 01 | 1 | CLI-01 | unit | `python -m pytest sc_tools/tests/test_cli.py::test_sct_entry_point -q` | ❌ W0 | ⬜ pending |
| 02-01-02 | 01 | 1 | CLI-02 | unit | `python -m pytest sc_tools/tests/test_cli.py::TestCommandGroups -q` | ❌ W0 | ⬜ pending |
| 02-01-03 | 01 | 1 | CLI-03 | unit | `python -m pytest sc_tools/tests/test_cli.py::TestCLIResult -q` | ❌ W0 | ⬜ pending |
| 02-01-04 | 01 | 1 | CLI-08 | unit | `python -m pytest sc_tools/tests/test_cli.py::test_lazy_imports -q` | ❌ W0 | ⬜ pending |
| 02-02-01 | 02 | 2 | CLI-04 | unit | `python -m pytest sc_tools/tests/test_cli.py::TestHumanOutput -q` | ❌ W0 | ⬜ pending |
| 02-02-02 | 02 | 2 | CLI-05, CLI-06 | unit | `python -m pytest sc_tools/tests/test_cli.py::TestErrorTaxonomy -q` | ❌ W0 | ⬜ pending |
| 02-02-03 | 02 | 2 | CLI-07 | unit | `python -m pytest sc_tools/tests/test_cli.py::test_non_interactive -q` | ❌ W0 | ⬜ pending |
| 02-02-04 | 02 | 2 | TST-04 | unit | `python -m pytest sc_tools/tests/test_cli.py -q` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `sc_tools/tests/test_cli.py` — CLI argument parsing and result envelope tests
- [ ] `sc_tools/tests/conftest.py` — shared fixtures (no conftest.py currently exists)

*Wave 0 creates test infrastructure. Executor plans include these as first tasks.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| `sct help` < 500ms | CLI-08 | Timing depends on system load | `time sct help` — verify < 500ms on dev machine |
| `--human` Rich output rendering | CLI-04 | Visual inspection of formatted output | Run `sct --human status` — verify Rich tables on stderr |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 15s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
