# Plan: Benchmark and QC Report Upgrade

## Context

sc_tools currently produces 7 isolated HTML report files that must be opened individually through
Finder. There is no shared navigation, no landing page to discover reports, inconsistent Bootstrap
versions across templates, and several known bugs (pinned Plotly CDN version, silent exception
swallowing, broken best-method heuristic). This plan upgrades all 7 templates to use a shared
Bootstrap 5.3 + Bootswatch Flatly visual identity, adds a sticky ScrollSpy sidebar to each report,
introduces a `figures/QC/index.html` landing page, and cleans up the outstanding bugs. Tasks are
ordered by priority: P1 addresses the core navigation pain point, P2 modernizes the template
architecture, P3 resolves bugs.

---

## Tasks

### P1 -- Navigation

- [x] **Task 1: `figures/QC/index.html` generator**
  - Add `generate_report_index(output_dir, project_name)` to `sc_tools/bm/report.py` (or new `sc_tools/reports/index.py`)
  - Scan `output_dir` for `*.html` files; extract title, date, best-method metadata from each
  - Render a Bootstrap 5.3 card grid (one card per report): title, type badge, date, key metric, "Open Report" button
  - Call this function at the end of each phase that produces reports
  - Output path: `figures/QC/index.html`

- [x] **Task 2: Bootstrap ScrollSpy sidebar in all 7 templates**
  - Stack: Bootstrap 5.3.8 + Bootswatch Flatly (CDN) + Bootstrap Bundle JS
  - Sidebar layout: `col-md-2`, `position: sticky; top: 0; height: 100vh; overflow-y: auto`
  - ScrollSpy binding: `data-bs-spy="scroll" data-bs-target="#sidebar" data-bs-offset="60"` on `<body>`
  - Every `<section>` gets a matching `id=` anchor
  - Dark mode toggle via Bootstrap `data-bs-theme` attribute
  - Active link style: `.nav-link.active { color: var(--bs-primary); font-weight: 600; }`
  - Python generators build a `sections` list before calling `render_template`, filtering out entries whose `key` maps to a `None` context value:
    ```python
    context["sections"] = [s for s in _SECTIONS if "key" not in s or context.get(s["key"]) is not None]
    ```
  - Section anchor targets per template:
    - `benchmark_report_template.html`: `exec-summary`, `strategy-comparison`, `method-ranking`, `per-dataset`, `tissue-analysis`*, `generalization`*, `statistics`*, `runtime`*, `failure-gallery`* (*conditional)
    - `integration_report_template.html`: `exec-summary`, `ranking`, `metrics-heatmap`, `radar`*, `batch-bio`*, `umap-grid`*
    - `segmentation_report_template.html`: `exec-summary`, `score-ranking`, `metrics-table`, `morphology`*, `marker-snr`*, `overlay`*

- [x] **Task 3: `_wrap_with_tabs` integration for benchmark reports**
  - Call `_find_latest_report` + `_wrap_with_tabs` from `sc_tools/qc/report_utils.py` at end of each benchmark generator
  - Prerequisite: Task 2 must be complete -- tabs + sidebar requires panel-scoped ScrollSpy
  - Tab panels bind ScrollSpy to the panel div, not `<body>`: `data-bs-spy="scroll" data-bs-target="#sidebar-{{ tab_id }}"` + `height:100vh; overflow-y:auto`
  - Each sidebar uses a unique id: `sidebar-tab-current`, `sidebar-tab-prev-1`, etc.

### P2 -- Template modernization

- [x] **Task 4: Jinja2 Environment upgrade**
  - Current: `Template(template_text)` (bare string loading) in `render_template` in both `bm/report.py` and `qc/report_utils.py`
  - Upgrade to: `jinja2.Environment(loader=FileSystemLoader(assets_dir)).get_template(name)`
  - Update signature: `render_template(template_name, context, assets_dir)` (replaces `template_path`)
  - This is a prerequisite for Task 5 (`{% extends %}` requires a loader)

- [x] **Task 5: Shared base template `base_report_template.html`**
  - New file `sc_tools/assets/base_report_template.html`: HTML shell with `<head>`, sidebar scaffold, Bootstrap 5.3 + Flatly CDN links, footer, dark mode toggle
  - Each child template uses `{% extends "base_report_template.html" %}` + `{% block content %}`
  - Eliminates: duplicated CDN tags, inconsistent header colors, repeated dark-mode toggle code across all 7 templates
  - Prerequisite: Task 4 (Jinja2 Environment loader required for `{% extends %}`)

- [x] **Task 6: Offline Plotly mode**
  - Add `offline: bool = False` parameter to all three `generate_*_report` functions
  - `offline=True`: first figure call uses `include_plotlyjs=True`, all subsequent use `False`
  - `offline=False` (default): CDN tag `<script src="https://cdn.plot.ly/plotly-2.27.0.min.js">` injected in `<head>` via base template
  - Default remains CDN to keep file sizes small (reports are typically read after `scp` to local machine)

### P3 -- Bug fixes

- [x] **Task 7: Fix `plotly-latest` CDN reference in `benchmark_report_template.html`**
  - Replace `plotly-latest.min.js` with `plotly-2.27.0.min.js`
  - One-line change; prevents breakage when Plotly CDN removes `latest` alias

- [x] **Task 8: Add `logger.warning` to silent exception blocks**
  - `bm/report.py`: 5 bare `except Exception: context[...] = None` blocks currently swallow errors silently
  - Add to each: `logger.warning("Failed to generate %s plot: %s", section_name, e, exc_info=True)`
  - Also add a warning when `failure_gallery_img` remains `None` (currently a dead code path)

- [x] **Task 9: Fix best-method heuristic in `generate_benchmark_report`**
  - Current logic: uses `boundary_regularity_mean`, falling back to `n_cells_mean` -- not appropriate for integration or QC reports
  - Fix: add `score_col: str | None = None` parameter; auto-detect from `composite_score`, `overall_score`, `boundary_regularity_mean`, `n_cells_mean` (in preference order)
  - Change affects only the "Best Method" display card -- no scoring logic is modified

- [x] **Task 10: Resolve dead `failure_gallery_img` block**
  - `bm/report.py` line ~386: `context["failure_gallery_img"] = None` is hardcoded and never populated
  - Decision: either wire to a real function from `sc_tools.pl.benchmarking`, or remove the `{% if failure_gallery_img %}` block from the template entirely
  - Recommend removing until a real implementation exists (avoids misleading conditional in template)

---

## Technical Decisions

1. **CDN default, offline opt-in**: `offline=False` keeps generated file sizes small; users on air-gapped HPC clusters pass `offline=True` to produce self-contained files.
2. **Jinja2 Environment required before base template**: `_render_template` must be upgraded to use `FileSystemLoader` before `{% extends %}` is valid -- Task 4 is a hard prerequisite for Task 5.
3. **Sidebar + tabs coexistence**: ScrollSpy must bind to the tab panel div (not `<body>`) with unique sidebar IDs per tab panel; binding to `<body>` with multiple sidebars causes the wrong links to activate.
4. **Dynamic sections list**: The Python generator builds the sections list and filters out entries whose associated context value is `None` before passing to `render_template`. Dead sidebar links break ScrollSpy visual state and confuse users.
5. **Execution order**: P1 Task 2 (sidebar) before P1 Task 3 (tabs); P2 Task 4 (Jinja2 env) before P2 Task 5 (base template). P3 bugs are independent and can be fixed in any order.

---

## File Map

```
sc_tools/assets/
  base_report_template.html             NEW        Task 5
  benchmark_report_template.html        MODIFY     Tasks 2, 6, 7
  integration_report_template.html      MODIFY     Tasks 2, 5
  segmentation_report_template.html     MODIFY     Tasks 2, 5
  post_integration_qc_template.html     MODIFY     Task 2
  pre_filter_qc_template.html           MODIFY     Task 2
  post_filter_qc_template.html          MODIFY     Task 2
  post_celltyping_qc_template.html      MODIFY     Task 2

sc_tools/bm/report.py                   MODIFY     Tasks 1, 3, 6, 8, 9, 10
sc_tools/qc/report_utils.py             MODIFY     Task 4 -- render_template upgrade
sc_tools/qc/report.py                   MODIFY     Task 4 -- call signature update
sc_tools/tests/test_bm_report.py        NEW        Tests for all new behavior
```
