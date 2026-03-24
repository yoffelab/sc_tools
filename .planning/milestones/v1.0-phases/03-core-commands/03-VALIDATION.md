---
phase: 03
slug: core-commands
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-21
---

# Phase 03 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest 8.x |
| **Config file** | `pyproject.toml` [tool.pytest.ini_options] |
| **Quick run command** | `python -m pytest sc_tools/tests/test_cli.py sc_tools/tests/test_cli_commands.py -q` |
| **Full suite command** | `python -m pytest sc_tools/tests/ -q` |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run `python -m pytest sc_tools/tests/test_cli.py sc_tools/tests/test_cli_commands.py -q`
- **After every plan wave:** Run `python -m pytest sc_tools/tests/ -q`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 03-01-01 | 01 | 1 | - | infra | `python -c "from sc_tools.cli import app"` | ✅ | ⬜ pending |
| 03-01-02 | 01 | 1 | TST-05 | unit | `python -m pytest sc_tools/tests/test_cli_commands.py -q` | ❌ W0 | ⬜ pending |
| 03-02-01 | 02 | 2 | CMD-01, CMD-06 | integration | `python -m pytest sc_tools/tests/test_cli_commands.py::TestQC -q` | ❌ W0 | ⬜ pending |
| 03-02-02 | 02 | 2 | CMD-02 | integration | `python -m pytest sc_tools/tests/test_cli_commands.py::TestPreprocess -q` | ❌ W0 | ⬜ pending |
| 03-02-03 | 02 | 2 | CMD-03, CMD-04 | integration | `python -m pytest sc_tools/tests/test_cli_commands.py::TestBenchmarkValidate -q` | ❌ W0 | ⬜ pending |
| 03-03-01 | 03 | 3 | CMD-05 | integration | `python -m pytest sc_tools/tests/test_cli_commands.py::TestStatus -q` | ❌ W0 | ⬜ pending |
| 03-03-02 | 03 | 3 | CMD-07, CMD-08 | integration | `python -m pytest sc_tools/tests/test_cli_commands.py::TestSharedResult -q` | ❌ W0 | ⬜ pending |
| 03-03-03 | 03 | 3 | TST-06 | e2e | `python -m pytest sc_tools/tests/test_cli_e2e.py -q` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `sc_tools/tests/conftest.py` — shared fixtures (100-cell AnnData with batch/celltype/embeddings/spatial)
- [ ] `sc_tools/tests/test_cli_commands.py` — integration tests for CLI commands with small fixtures

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Benchmark HTML report shows subsampling info | D-12 | Visual inspection of report | Run `sct benchmark integration --from-dir <dir>`, open HTML, verify per-sample subsampling note |
| `--human` Rich output formatting | CMD-06 | Visual inspection | Run `sct --human report pre_filter <file>`, verify Rich tables on stderr |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
