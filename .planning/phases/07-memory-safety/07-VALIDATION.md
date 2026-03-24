---
phase: 7
slug: memory-safety
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-24
---

# Phase 7 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest 7.x |
| **Config file** | pyproject.toml `[tool.pytest.ini_options]` |
| **Quick run command** | `python -m pytest sc_tools/tests/test_io_gateway.py -x -q` |
| **Full suite command** | `python -m pytest sc_tools/tests/test_io_gateway.py sc_tools/tests/test_cli_estimate.py sc_tools/tests/test_cli_dryrun.py -x -q` |
| **Estimated runtime** | ~15 seconds |

---

## Sampling Rate

- **After every task commit:** Run `python -m pytest sc_tools/tests/test_io_gateway.py -x -q`
- **After every plan wave:** Run full suite command
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 15 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 07-01-01 | 01 | 1 | MEM-01 | unit | `pytest sc_tools/tests/test_io_gateway.py::test_metadata_read -x` | ❌ W0 | ⬜ pending |
| 07-01-02 | 01 | 1 | MEM-01 | unit | `pytest sc_tools/tests/test_io_gateway.py::test_backed_read -x` | ❌ W0 | ⬜ pending |
| 07-01-03 | 01 | 1 | MEM-01 | unit | `pytest sc_tools/tests/test_io_gateway.py::test_memory_guard -x` | ❌ W0 | ⬜ pending |
| 07-02-01 | 02 | 2 | MEM-02 | unit | `pytest sc_tools/tests/test_cli_estimate.py -x` | ❌ W0 | ⬜ pending |
| 07-02-02 | 02 | 2 | MEM-03 | unit | `pytest sc_tools/tests/test_cli_dryrun.py -x` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `sc_tools/tests/test_io_gateway.py` — stubs for MEM-01 (tiered loading, memory guard)
- [ ] `sc_tools/tests/test_cli_estimate.py` — stubs for MEM-02 (estimate command)
- [ ] `sc_tools/tests/test_cli_dryrun.py` — stubs for MEM-03 (dry-run flag)

*Existing pytest infrastructure and fixtures cover test framework needs.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Peak memory <1GB for metadata queries on 25G file | MEM-01 | Requires 25G test file on HPC | Run `sct validate ingest <25G_file>` on brb, check peak_memory_mb in provenance sidecar |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 15s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
