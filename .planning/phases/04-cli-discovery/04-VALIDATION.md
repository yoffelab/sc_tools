---
phase: 4
slug: cli-discovery
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-23
---

# Phase 4 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest 7.x |
| **Config file** | pyproject.toml `[tool.pytest]` |
| **Quick run command** | `python -m pytest sc_tools/tests/test_cli_discovery.py -x -q` |
| **Full suite command** | `python -m pytest sc_tools/tests/test_cli_discovery.py -v` |
| **Estimated runtime** | ~5 seconds |

---

## Sampling Rate

- **After every task commit:** Run `python -m pytest sc_tools/tests/test_cli_discovery.py -x -q`
- **After every plan wave:** Run `python -m pytest sc_tools/tests/test_cli_discovery.py -v`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 10 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 04-01-01 | 01 | 1 | DSC-01 | unit | `pytest sc_tools/tests/test_cli_discovery.py::test_list_commands_json -x` | ❌ W0 | ⬜ pending |
| 04-01-02 | 01 | 1 | DSC-02 | unit | `pytest sc_tools/tests/test_cli_discovery.py::test_describe_command -x` | ❌ W0 | ⬜ pending |
| 04-01-03 | 01 | 1 | DSC-03 | unit | `pytest sc_tools/tests/test_cli_discovery.py::test_schema_full_contract -x` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `sc_tools/tests/test_cli_discovery.py` — test stubs for DSC-01, DSC-02, DSC-03
- [ ] Shared fixtures: `runner = CliRunner()` for invoking Typer commands

*Existing pytest infrastructure and conftest.py from Phase 3 cover shared fixtures.*

---

## Manual-Only Verifications

*All phase behaviors have automated verification.*

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 10s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
