---
name: k-dense-frontend
tier: 3
description: "[Tier 3 — Specialized] Self-contained HTML report design for sc_tools benchmark and QC reports. Covers Bootstrap 5.3 sticky sidebar nav, ScrollSpy, Plotly embedding, Jinja2 template patterns, multi-report index pages, and dark mode. Use when designing, reviewing, or improving any sc_tools HTML report template in sc_tools/assets/."
allowed-tools: Read Write Edit Bash Glob Grep WebFetch
---

# K-Dense Frontend Reference

## Core Principle

sc_tools HTML reports are **self-contained interactive files** that:
- Open at `file://` with no server, no build step
- Embed all Plotly charts inline (no external data calls after load)
- Use Jinja2 for server-side rendering (Python, not a JS framework)
- Prioritize **readability and navigation** over creative aesthetics

The gold standard is ydata-profiling: Bootstrap 5 + Jinja2 + anchor nav + tabs.

---

## Stack

| Layer | Choice | CDN |
|-------|--------|-----|
| CSS framework | Bootstrap 5.3.8 | `https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.8/css/bootstrap.min.css` |
| Theme | Bootswatch Flatly | `https://cdn.jsdelivr.net/npm/bootswatch@5.3.3/dist/flatly/bootstrap.min.css` |
| JS (includes Popper) | Bootstrap Bundle 5.3.8 | `https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.8/js/bootstrap.bundle.min.js` |
| Charts | Plotly 2.27.0 | `https://cdn.plot.ly/plotly-2.27.0.min.js` |
| Templating | Jinja2 (Python) | — |

**NEVER use `plotly-latest`** — the CDN deprecated it and it may break silently.

**Offline mode**: Replace CDN links with inline `<style>` / `<script>` blocks containing the minified source. Use Jinja2 `{{ bootstrap_css | safe }}` pattern (see ydata-profiling for reference).

---

## Left Sidebar Navigation Pattern

The primary UX requirement for all long reports (>3 sections).

```html
<!DOCTYPE html>
<html lang="en" data-bs-theme="light">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{{ title }}</title>
  <!-- Bootswatch Flatly — academic, clean -->
  <link rel="stylesheet"
    href="https://cdn.jsdelivr.net/npm/bootswatch@5.3.3/dist/flatly/bootstrap.min.css">
  <style>
    html { scroll-behavior: smooth; }

    /* Sidebar: full-height sticky column */
    #sidebar {
      position: sticky;
      top: 0;
      height: 100vh;
      overflow-y: auto;
      border-right: 1px solid var(--bs-border-color);
      padding: 1rem 0.5rem;
    }

    /* Offset anchors so sticky header doesn't cover them */
    section[id] {
      scroll-margin-top: 1rem;
    }

    /* Active link highlight driven by ScrollSpy */
    #sidebar .nav-link.active {
      color: var(--bs-primary);
      font-weight: 600;
      background: rgba(var(--bs-primary-rgb), 0.08);
      border-radius: 4px;
    }

    /* Sub-links indent */
    #sidebar .nav-link.sub { font-size: 0.85em; padding-left: 1.5rem; }
  </style>
</head>

<!-- ScrollSpy: highlight active section as user scrolls -->
<body data-bs-spy="scroll" data-bs-target="#sidebar" data-bs-offset="60" tabindex="0">

<div class="container-fluid">
  <div class="row g-0">

    <!-- Sidebar: col-2 on md+, hidden on mobile -->
    <nav id="sidebar" class="col-md-2 d-none d-md-block bg-light">
      <div class="px-2 mb-3">
        <h6 class="text-muted text-uppercase fw-bold small">{{ title }}</h6>
        <small class="text-muted">{{ date }}</small>
      </div>
      <ul class="nav flex-column">
        {% for section in sections %}
        <li class="nav-item">
          <a class="nav-link" href="#{{ section.id }}">{{ section.label }}</a>
        </li>
        {% endfor %}
      </ul>
      <!-- Dark mode toggle -->
      <div class="px-2 mt-3">
        <button class="btn btn-sm btn-outline-secondary w-100"
          onclick="document.documentElement.setAttribute('data-bs-theme',
            document.documentElement.getAttribute('data-bs-theme') === 'dark' ? 'light' : 'dark')">
          Toggle Dark Mode
        </button>
      </div>
    </nav>

    <!-- Main content -->
    <main class="col-md-10 py-4 px-4">
      <!-- sections go here, each with id= -->
    </main>

  </div>
</div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.8/js/bootstrap.bundle.min.js"></script>
</body>
</html>
```

### Section anchor pattern

Every `<section>` must have an `id` for ScrollSpy and deep-linking:

```html
<section id="exec-summary">
  <h2>Executive Summary</h2>
  ...
</section>

<section id="strategy-comparison">
  <h2>Strategy Comparison</h2>
  ...
</section>
```

Pass the section list to Jinja2 context as:
```python
context["sections"] = [
    {"id": "exec-summary",        "label": "Executive Summary"},
    {"id": "strategy-comparison", "label": "Strategy Comparison"},
    {"id": "method-ranking",      "label": "Method Ranking"},
    {"id": "per-dataset",         "label": "Per-Dataset Heatmap"},
    {"id": "tissue-analysis",     "label": "Tissue Analysis"},
    {"id": "generalization",      "label": "Generalization"},
    {"id": "statistics",          "label": "Statistics"},
    {"id": "runtime",             "label": "Runtime"},
]
```

---

## Plotly Embedding in Jinja2

### CDN mode (default — requires internet on first open)

```python
# Python: generate div string for embedding
plotly_div = fig.to_html(full_html=False, include_plotlyjs=False)
context["plot_ranking"] = plotly_div
```

Include Plotly JS once in `<head>`:
```html
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
```

Render in template:
```html
{% if plot_ranking %}
  {{ plot_ranking | safe }}
{% endif %}
```

### Fully offline (self-contained, no internet required)

```python
# First figure: embed plotly.js (~3MB) once
first_plot = fig1.to_html(full_html=False, include_plotlyjs=True)
# All subsequent figures: reference the already-loaded library
other_plots = [fig.to_html(full_html=False, include_plotlyjs=False) for fig in [fig2, fig3]]
```

Use this when the report will be opened on air-gapped HPC systems or shared via email.

---

## Multi-Report Index Page

When multiple benchmark reports exist for a project, generate an `index.html` in `figures/QC/`:

```html
<!-- figures/QC/index.html -->
<section id="benchmarks">
  <h2>Benchmark Reports</h2>
  <div class="row row-cols-1 row-cols-md-3 g-4">
    {% for report in reports %}
    <div class="col">
      <div class="card h-100">
        <div class="card-body">
          <h5 class="card-title">{{ report.title }}</h5>
          <p class="card-text text-muted small">{{ report.date }}</p>
          <p class="card-text">
            Best: <strong>{{ report.best_method }}</strong><br>
            Score: {{ "%.3f" | format(report.best_score) }}
          </p>
        </div>
        <div class="card-footer">
          <a href="{{ report.filename }}" class="btn btn-primary btn-sm">Open Report</a>
        </div>
      </div>
    </div>
    {% endfor %}
  </div>
</section>
```

The Python generator should scan `figures/QC/` for `*_report.html` files and extract metadata from their `<title>` and summary cards.

---

## Tab Pattern (for combining related reports)

sc_tools QC already uses `_wrap_with_tabs` from `sc_tools/qc/report_utils.py`. Apply the same pattern to benchmark reports.

```python
# In bm/report.py — after generating a new report, check for prior reports
from sc_tools.qc.report_utils import _find_latest_report, _wrap_with_tabs

prior = _find_latest_report(output_path.parent, pattern="*_benchmark_report*.html")
if prior:
    combined = _wrap_with_tabs([prior, output_path], tab_labels=["Previous", "Current"])
    output_path.write_text(combined)
```

---

## Design Standards for Scientific Reports

| Principle | Rule |
|-----------|------|
| Theme | Bootswatch Flatly (light) — clean, academic, no purple gradients |
| Typography | Bootstrap defaults (system font stack) — no decorative fonts |
| Color | Bootstrap semantic colors (`text-success`, `text-danger`, `text-warning`) for status |
| Tables | `table table-striped table-hover table-sm` — always |
| Charts | Plotly only (interactive); matplotlib only for spatial overlays (embed as base64 PNG) |
| Metrics | Always show N, metric name, value, and interpretation together |
| Status badges | `<span class="badge bg-success">PASS</span>` etc. |
| Cards | Use Bootstrap `.card` for each major metric block |
| Error states | Never silently hide a failed section — show `<div class="alert alert-warning">` instead |

### DO NOT use the Anthropic `frontend-design` skill for sc_tools reports
That skill targets creative/artistic consumer UIs. sc_tools reports need density, legibility, and navigability — not distinctive typography or animated micro-interactions.

---

## Quick Reference: Template Paths

```
sc_tools/assets/
  benchmark_report_template.html      ← generate_benchmark_report()
  integration_report_template.html    ← generate_integration_report()
  segmentation_report_template.html   ← generate_segmentation_report()
  pre_filter_qc_template.html         ← generate_pre_filter_report()
  post_filter_qc_template.html        ← generate_post_filter_report()
  post_integration_qc_template.html   ← generate_post_integration_report()
  post_celltyping_qc_template.html    ← generate_post_celltyping_report()
  qc_report_template.html             ← generate_qc_report() (legacy)

sc_tools/bm/report.py                 ← benchmark report generators
sc_tools/qc/report.py                 ← QC report generators
sc_tools/qc/report_utils.py           ← _wrap_with_tabs, _find_latest_report, render_template
```
