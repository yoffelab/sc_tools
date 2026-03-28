---
phase: 09-sample-concatenation-maintenance
plan: 01
subsystem: dependencies
tags: [plotly, cdn, dependencies, maintenance]
dependency_graph:
  requires: []
  provides: [plotly-base-dep, plotly-cdn-3.4.0]
  affects: [report-rendering, html-templates]
tech_stack:
  added: []
  patterns: [cdn-version-pinning]
key_files:
  created: []
  modified:
    - pyproject.toml
    - sc_tools/qc/report_utils.py
    - sc_tools/assets/base_report_template.html
    - sc_tools/tests/test_qc.py
    - .claude/skills/k-dense-frontend/SKILL.md
decisions:
  - Pin plotly CDN to explicit 3.4.0 (not plotly-latest which is frozen at 1.58.5)
metrics:
  duration: 4min
  completed: "2026-03-28"
  tasks_completed: 2
  tasks_total: 2
requirements: [MAINT-01, MAINT-02]
---

# Phase 09 Plan 01: Plotly Dependency Promotion & CDN Update Summary

Promoted plotly from optional [pipeline] extra to base dependency and updated CDN pin from 2.27.0 to 3.4.0 matching plotly.py 6.6.0

## Task Summary

| Task | Description | Commit | Key Files |
|------|-------------|--------|-----------|
| 1 | Promote plotly to base dependencies (MAINT-01) | 62d8b9e | pyproject.toml |
| 2 | Update Plotly CDN pin 2.27.0 to 3.4.0 (MAINT-02) | cd9ad01 | report_utils.py, base_report_template.html, test_qc.py, SKILL.md |

## What Changed

### Task 1: Promote plotly to base dependencies
- Moved `plotly>=5.18` from `[pipeline]`, `[benchmark]`, and `[spatial]` optional extras into `[project.dependencies]`
- plotly is now available in any sc_tools install without needing extras

### Task 2: Update Plotly CDN pin
- Updated CDN script tag from `plotly-2.27.0.min.js` to `plotly-3.4.0.min.js` in:
  - `sc_tools/qc/report_utils.py` (render_template wrapper)
  - `sc_tools/assets/base_report_template.html` (Jinja2 base template)
  - `sc_tools/tests/test_qc.py` (TestBaseTemplate assertion)
  - `.claude/skills/k-dense-frontend/SKILL.md` (developer reference docs)

## Verification Results

- `grep -c "plotly-3.4.0.min.js"` returns matches in all 4 source files
- `grep "plotly-2.27.0"` returns zero matches in source files
- `make lint` passes (ruff check + format clean)

## Deviations from Plan

None - plan executed exactly as written.

## Known Stubs

None.

## Self-Check: PASSED
