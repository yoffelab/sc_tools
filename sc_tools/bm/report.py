"""
HTML report generation for segmentation and integration benchmarks.

Uses Plotly for interactive charts and matplotlib for static plots
(embedded as base64 PNGs).
"""

from __future__ import annotations

import base64
import io
import logging
import re
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData

from sc_tools.qc.report_utils import (
    _extract_body_content,
    _extract_head_css,
    _wrap_with_tabs,
    render_template,
)

__all__ = [
    "generate_segmentation_report",
    "generate_integration_report",
    "generate_benchmark_report",
    "generate_report_index",
]

logger = logging.getLogger(__name__)

matplotlib.use("Agg")


def _find_latest_bm_report(
    output_dir: Path,
    pattern: str,
    exclude: str = "",
) -> Path | None:
    """Find the most recently modified HTML report matching *pattern* in *output_dir*.

    Parameters
    ----------
    output_dir
        Directory to search.
    pattern
        Glob pattern, e.g. ``"*integration*report*.html"``.
    exclude
        Filename to exclude from candidates (typically the just-written report).

    Returns
    -------
    Path or None
        The matching file with the latest modification time, or ``None`` if
        no candidates exist after exclusion.
    """
    candidates = [
        p for p in output_dir.glob(pattern) if p.name != exclude and p.name != "index.html"
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


_SEG_TEMPLATE = "segmentation_report_template.html"
_INT_TEMPLATE = "integration_report_template.html"
_BM_TEMPLATE = "benchmark_report_template.html"


def _fig_to_base64(fig: plt.Figure, dpi: int = 150) -> str:
    """Render a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def _plotly_to_html(fig, include_plotlyjs: bool | str = False) -> str:
    """Convert a Plotly figure to an embeddable HTML div."""
    return fig.to_html(include_plotlyjs=include_plotlyjs, full_html=False)


def generate_segmentation_report(
    comparison_df: pd.DataFrame,
    morphology_dict: dict[str, pd.DataFrame] | None = None,
    marker_quality_dict: dict[str, pd.DataFrame] | None = None,
    intensity_image: np.ndarray | None = None,
    masks: dict[str, np.ndarray] | None = None,
    output_path: str | Path = "segmentation_benchmark_report.html",
    title: str = "Segmentation Benchmark Report",
    offline: bool = False,
) -> Path:
    """Generate a self-contained HTML segmentation benchmark report.

    Parameters
    ----------
    comparison_df
        Output from ``compare_segmentations``.
    morphology_dict
        Optional dict of morphology DataFrames per method.
    marker_quality_dict
        Optional dict of marker quality DataFrames per method.
    intensity_image
        Optional intensity image for overlay visualization.
    masks
        Optional dict of masks for overlay visualization.
    output_path
        Where to write the HTML report.
    title
        Report title.

    Returns
    -------
    Path to the generated report.
    """
    from sc_tools.pl.benchmarking import (
        plot_composite_score_bar,
        plot_segmentation_comparison_table,
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    best_row = comparison_df.iloc[0] if len(comparison_df) > 0 else {}
    best_method = best_row.get("method", "N/A")
    best_score = best_row.get("composite_score", 0.0)

    _first_plot = [True]  # mutable container for closure

    def _to_html_offline(fig):
        if _first_plot[0]:
            _first_plot[0] = False
            return _plotly_to_html(fig, include_plotlyjs=True)
        return _plotly_to_html(fig, include_plotlyjs=False)

    _render = _to_html_offline if offline else _plotly_to_html

    context = {
        "title": title,
        "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "n_methods": len(comparison_df),
        "has_gt": "detection_f1" in comparison_df.columns,
        "best_method": best_method,
        "best_score": best_score,
        "offline": offline,
        "plot_composite_bar": _render(plot_composite_score_bar(comparison_df)),
        "plot_comparison_table": _render(plot_segmentation_comparison_table(comparison_df)),
    }

    # Optional morphology plot
    if morphology_dict:
        from sc_tools.pl.benchmarking import plot_morphology_distributions

        context["plot_morphology"] = _render(plot_morphology_distributions(morphology_dict))
    else:
        context["plot_morphology"] = None

    # Optional marker SNR plot
    if marker_quality_dict:
        from sc_tools.pl.benchmarking import plot_marker_snr_heatmap

        context["plot_marker_snr"] = _render(plot_marker_snr_heatmap(marker_quality_dict))
    else:
        context["plot_marker_snr"] = None

    # Optional overlay
    if intensity_image is not None and masks is not None:
        from sc_tools.pl.benchmarking import plot_segmentation_overlay

        fig = plot_segmentation_overlay(intensity_image, masks)
        context["overlay_img"] = _fig_to_base64(fig)
    else:
        context["overlay_img"] = None

    _SEG_SECTIONS = [
        {"id": "exec-summary", "label": "Executive Summary"},
        {"id": "metrics-table", "label": "Metrics Comparison"},
        {"id": "score-ranking", "label": "Score Ranking"},
        {"id": "morphology", "label": "Morphology", "key": "plot_morphology"},
        {"id": "marker-snr", "label": "Marker SNR", "key": "plot_marker_snr"},
        {"id": "overlay", "label": "Segmentation Overlay", "key": "overlay_img"},
    ]
    context["sections"] = [
        s for s in _SEG_SECTIONS if "key" not in s or context.get(s["key"]) is not None
    ]

    html = render_template(_SEG_TEMPLATE, context)
    output_path.write_text(html)
    logger.info("Segmentation report saved to %s", output_path)

    try:
        prior = _find_latest_bm_report(
            output_path.parent,
            pattern="*segmentation*report*.html",
            exclude=output_path.name,
        )
        if prior:
            prior_html = prior.read_text()
            combined = _wrap_with_tabs(
                current_label="Current",
                current_content=_extract_body_content(html),
                previous_tabs=[
                    ("Previous", _extract_body_content(prior_html), _extract_head_css(prior_html))
                ],
                current_css=_extract_head_css(html),
            )
            output_path.write_text(combined)
            logger.info("Tabs combined with previous segmentation report %s", prior.name)
    except Exception as e:
        logger.warning("Could not combine segmentation report with previous: %s", e)

    return output_path


def generate_integration_report(
    comparison_df: pd.DataFrame,
    adata: AnnData | None = None,
    embeddings: dict[str, str] | None = None,
    batch_key: str = "batch",
    celltype_key: str = "celltype",
    color_by: str | None = None,
    output_path: str | Path = "integration_benchmark_report.html",
    title: str = "Integration Benchmark Report",
    offline: bool = False,
) -> Path:
    """Generate a self-contained HTML integration benchmark report.

    Parameters
    ----------
    comparison_df
        Output from ``compare_integrations``.
    adata
        Optional AnnData for UMAP grid.
    embeddings
        Optional dict of embeddings for UMAP grid.
    batch_key
        Batch column name (for display).
    celltype_key
        Cell type column name (for display and UMAP coloring).
    color_by
        Column to color UMAPs by. Defaults to ``celltype_key``.
    output_path
        Where to write the HTML report.
    title
        Report title.

    Returns
    -------
    Path to the generated report.
    """
    from sc_tools.pl.benchmarking import (
        plot_integration_comparison_table,
        plot_integration_radar,
        plot_integration_ranking_bar,
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if color_by is None:
        color_by = celltype_key

    best_row = comparison_df.iloc[0] if len(comparison_df) > 0 else {}
    best_method = best_row.get("method", "N/A")
    best_overall = best_row.get("overall_score", 0.0)
    best_batch = best_row.get("batch_score", 0.0)
    best_bio = best_row.get("bio_score", 0.0)

    _first_plot = [True]  # mutable container for closure

    def _to_html_offline(fig):
        if _first_plot[0]:
            _first_plot[0] = False
            return _plotly_to_html(fig, include_plotlyjs=True)
        return _plotly_to_html(fig, include_plotlyjs=False)

    _render = _to_html_offline if offline else _plotly_to_html

    context = {
        "title": title,
        "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "n_methods": len(comparison_df),
        "batch_key": batch_key,
        "celltype_key": celltype_key,
        "best_method": best_method,
        "best_overall": best_overall,
        "best_batch": best_batch,
        "best_bio": best_bio,
        "color_by": color_by,
        "offline": offline,
        "plot_ranking_bar": _render(plot_integration_ranking_bar(comparison_df)),
        "plot_comparison_table": _render(plot_integration_comparison_table(comparison_df)),
    }

    # Radar
    exclude = {"method", "embedding_key", "batch_score", "bio_score", "overall_score"}
    metric_cols = [c for c in comparison_df.columns if c not in exclude]
    if len(metric_cols) >= 3:
        context["plot_radar"] = _render(plot_integration_radar(comparison_df))
    else:
        context["plot_radar"] = None

    # UMAP grid (static matplotlib)
    if adata is not None and embeddings is not None:
        from sc_tools.pl.benchmarking import plot_embedding_comparison_umap

        fig = plot_embedding_comparison_umap(adata, embeddings, color_by=color_by)
        context["umap_img"] = _fig_to_base64(fig)
    else:
        context["umap_img"] = None

    # Batch × bio conservation table
    _bb_cols = [
        c
        for c in ["method", "batch_score", "bio_score", "overall_score"]
        if c in comparison_df.columns
    ]
    if len(_bb_cols) > 1:
        _bb_df = (
            comparison_df[_bb_cols]
            .copy()
            .sort_values(
                "overall_score" if "overall_score" in _bb_cols else _bb_cols[-1], ascending=False
            )
        )
        for col in _bb_cols[1:]:
            _bb_df[col] = _bb_df[col].map(lambda x: f"{x:.3f}" if isinstance(x, float) else x)
        context["batch_bio_table_html"] = _bb_df.to_html(
            index=False,
            classes="table table-striped table-hover table-sm",
            border=0,
        )
    else:
        context["batch_bio_table_html"] = None

    _INT_SECTIONS = [
        {"id": "exec-summary", "label": "Executive Summary"},
        {"id": "batch-bio-table", "label": "Batch \u00d7 Bio Table", "key": "batch_bio_table_html"},
        {"id": "metrics-heatmap", "label": "Metrics Heatmap"},
        {"id": "umap-grid", "label": "UMAP Comparison", "key": "umap_img"},
        {"id": "ranking", "label": "Integration Ranking"},
        {"id": "radar", "label": "Metrics Radar", "key": "plot_radar"},
    ]
    context["sections"] = [
        s for s in _INT_SECTIONS if "key" not in s or context.get(s["key"]) is not None
    ]

    html = render_template(_INT_TEMPLATE, context)
    output_path.write_text(html)
    logger.info("Integration report saved to %s", output_path)

    try:
        prior = _find_latest_bm_report(
            output_path.parent,
            pattern="*integration*report*.html",
            exclude=output_path.name,
        )
        if prior:
            prior_html = prior.read_text()
            combined = _wrap_with_tabs(
                current_label="Current",
                current_content=_extract_body_content(html),
                previous_tabs=[
                    ("Previous", _extract_body_content(prior_html), _extract_head_css(prior_html))
                ],
                current_css=_extract_head_css(html),
            )
            output_path.write_text(combined)
            logger.info("Tabs combined with previous integration report %s", prior.name)
    except Exception as e:
        logger.warning("Could not combine integration report with previous: %s", e)

    return output_path


def generate_benchmark_report(
    results_df: pd.DataFrame,
    aggregated: dict[str, pd.DataFrame] | None = None,
    output_path: str | Path = "benchmark_report.html",
    title: str = "IMC Segmentation Benchmark Report",
    score_col: str | None = None,
    offline: bool = False,
) -> Path:
    """Generate a comprehensive HTML benchmark report.

    Includes: executive summary, strategy comparison, method ranking,
    per-dataset heatmap, per-tissue analysis, generalization,
    statistical comparisons, runtime comparison, failure gallery.

    Parameters
    ----------
    results_df
        Full benchmark results DataFrame.
    aggregated
        Output from ``aggregate_results()``. Computed if not provided.
    output_path
        Where to write the HTML report.
    title
        Report title.

    Returns
    -------
    Path to the generated report.
    """
    from sc_tools.pl.benchmarking import (
        plot_dataset_heatmap,
        plot_runtime_comparison,
        plot_strategy_comparison,
        plot_tissue_boxplot,
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if aggregated is None:
        from sc_tools.bm.runner import aggregate_results

        aggregated = aggregate_results(results_df)

    # Summary stats
    n_rois = results_df["roi_id"].nunique() if "roi_id" in results_df.columns else len(results_df)
    n_datasets = results_df["dataset"].nunique() if "dataset" in results_df.columns else 0
    n_strategies = results_df["strategy"].nunique() if "strategy" in results_df.columns else 0
    n_methods = results_df["method"].nunique() if "method" in results_df.columns else 0

    # Best method: use caller-supplied score_col, or auto-detect from preferred candidates
    by_method = aggregated.get("by_method", pd.DataFrame())
    best_method = "N/A"
    if len(by_method) > 0:
        _score_col = score_col
        if _score_col is None:
            for candidate in [
                "composite_score",
                "overall_score",
                "boundary_regularity_mean",
                "n_cells_mean",
            ]:
                if candidate in by_method.columns:
                    _score_col = candidate
                    break
        if _score_col:
            best_method = by_method.sort_values(_score_col, ascending=False).iloc[0]["method"]

    _first_plot = [True]  # mutable container for closure

    def _to_html_offline(fig):
        if _first_plot[0]:
            _first_plot[0] = False
            return _plotly_to_html(fig, include_plotlyjs=True)
        return _plotly_to_html(fig, include_plotlyjs=False)

    _render = _to_html_offline if offline else _plotly_to_html

    context = {
        "title": title,
        "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "n_rois": n_rois,
        "n_datasets": n_datasets,
        "n_strategies": n_strategies,
        "n_methods": n_methods,
        "best_method": best_method,
        "offline": offline,
    }

    # Strategy comparison plot
    try:
        fig = plot_strategy_comparison(results_df)
        context["plot_strategy_comparison"] = _render(fig)
    except Exception as e:
        logger.warning("Failed to generate strategy comparison plot: %s", e, exc_info=True)
        context["plot_strategy_comparison"] = None

    # Method ranking
    try:
        fig = plot_dataset_heatmap(results_df)
        context["plot_method_ranking"] = _render(fig)
    except Exception as e:
        logger.warning("Failed to generate method ranking plot: %s", e, exc_info=True)
        context["plot_method_ranking"] = None

    # Dataset heatmap
    try:
        fig = plot_dataset_heatmap(results_df, metric="n_cells")
        context["plot_dataset_heatmap"] = _render(fig)
    except Exception as e:
        logger.warning("Failed to generate dataset heatmap plot: %s", e, exc_info=True)
        context["plot_dataset_heatmap"] = None

    # Tissue boxplot
    try:
        fig = plot_tissue_boxplot(results_df)
        context["plot_tissue_boxplot"] = _render(fig)
    except Exception as e:
        logger.warning("Failed to generate tissue boxplot plot: %s", e, exc_info=True)
        context["plot_tissue_boxplot"] = None

    # Generalization
    try:
        from sc_tools.bm.analysis import compute_generalization_matrix
        from sc_tools.pl.benchmarking import plot_generalization_radar

        gen_matrix = compute_generalization_matrix(results_df)
        fig = plot_generalization_radar(gen_matrix)
        context["plot_generalization"] = _render(fig)
    except Exception as e:
        logger.warning("Failed to generate generalization plot: %s", e, exc_info=True)
        context["plot_generalization"] = None

    # Statistical comparisons
    try:
        from sc_tools.bm.analysis import statistical_comparison

        stats_df = statistical_comparison(results_df)
        if len(stats_df) > 0:
            context["stats_table"] = stats_df.to_html(
                index=False, classes="stats-table", float_format="%.4f"
            )
        else:
            context["stats_table"] = None
    except Exception as e:
        logger.warning("Failed to generate statistical comparisons: %s", e, exc_info=True)
        context["stats_table"] = None

    # Runtime comparison
    if "runtime_s" in results_df.columns:
        try:
            fig = plot_runtime_comparison(results_df)
            context["plot_runtime"] = _render(fig)
        except Exception as e:
            logger.warning("Failed to generate runtime comparison plot: %s", e, exc_info=True)
            context["plot_runtime"] = None
    else:
        context["plot_runtime"] = None

    # Failure gallery (not implemented)
    logger.warning("Failure gallery not implemented; skipping failure-gallery section")
    context["failure_gallery_img"] = None

    _BM_SECTIONS = [
        {"id": "exec-summary", "label": "Executive Summary"},
        {"id": "method-ranking", "label": "Method Ranking", "key": "plot_method_ranking"},
        {"id": "statistics", "label": "Statistical Tests", "key": "stats_table"},
        {"id": "per-dataset", "label": "Per-Dataset Heatmap", "key": "plot_dataset_heatmap"},
        {
            "id": "strategy-comparison",
            "label": "Strategy Comparison",
            "key": "plot_strategy_comparison",
        },
        {"id": "tissue-analysis", "label": "Tissue Analysis", "key": "plot_tissue_boxplot"},
        {"id": "generalization", "label": "Generalization", "key": "plot_generalization"},
        {"id": "runtime", "label": "Runtime", "key": "plot_runtime"},
    ]
    context["sections"] = [
        s for s in _BM_SECTIONS if "key" not in s or context.get(s["key"]) is not None
    ]

    html = render_template(_BM_TEMPLATE, context)
    output_path.write_text(html)
    logger.info("Benchmark report saved to %s", output_path)

    try:
        prior = _find_latest_bm_report(
            output_path.parent,
            pattern="*benchmark*report*.html",
            exclude=output_path.name,
        )
        if prior:
            prior_html = prior.read_text()
            combined = _wrap_with_tabs(
                current_label="Current",
                current_content=_extract_body_content(html),
                previous_tabs=[
                    ("Previous", _extract_body_content(prior_html), _extract_head_css(prior_html))
                ],
                current_css=_extract_head_css(html),
            )
            output_path.write_text(combined)
            logger.info("Tabs combined with previous benchmark report %s", prior.name)
    except Exception as e:
        logger.warning("Could not combine benchmark report with previous: %s", e)

    return output_path


_INDEX_TEMPLATE = """<!DOCTYPE html>
<html lang="en" data-bs-theme="light">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>sc_tools Reports{% if project_name %} \u2014 {{ project_name }}{% endif %}</title>
  <link rel="stylesheet"
    href="https://cdn.jsdelivr.net/npm/bootswatch@5.3.3/dist/flatly/bootstrap.min.css">
</head>
<body class="bg-light">
<div class="container py-5">
  <div class="d-flex justify-content-between align-items-center mb-4">
    <div>
      <h1 class="h3 mb-0">sc_tools Reports</h1>
      {% if project_name %}<p class="text-muted mb-0">{{ project_name }}</p>{% endif %}
    </div>
    <button class="btn btn-sm btn-outline-secondary"
      onclick="var e=document.documentElement;e.setAttribute('data-bs-theme',e.getAttribute('data-bs-theme')==='dark'?'light':'dark')">
      Toggle Dark Mode
    </button>
  </div>
  {% if reports %}
  <div class="row row-cols-1 row-cols-md-3 g-4">
    {% for r in reports %}
    <div class="col">
      <div class="card h-100 shadow-sm">
        <div class="card-body">
          <div class="d-flex justify-content-between align-items-start mb-2">
            <h5 class="card-title mb-0" style="font-size:0.95rem;">{{ r.title }}</h5>
            <span class="badge ms-2 flex-shrink-0
              {%- if r.report_type == 'benchmark' %} bg-primary
              {%- elif r.report_type == 'integration' %} bg-success
              {%- elif r.report_type == 'segmentation' %} bg-warning text-dark
              {%- elif r.report_type == 'celltyping' %} bg-info text-dark
              {%- else %} bg-secondary{%- endif %}">{{ r.report_type }}</span>
          </div>
          <p class="text-muted small mb-2">{{ r.date }}</p>
          {% if r.best_method and r.best_method != '\u2014' %}
          <p class="mb-0 small"><strong>{{ r.best_method }}</strong></p>
          {% endif %}
        </div>
        <div class="card-footer bg-transparent border-top-0">
          <a href="{{ r.filename }}" class="btn btn-primary btn-sm w-100">Open Report</a>
        </div>
      </div>
    </div>
    {% endfor %}
  </div>
  {% else %}
  <div class="alert alert-info">No reports found in this directory.</div>
  {% endif %}
  <footer class="text-center text-muted small mt-5">
    Generated {{ generated_date }} \u2014 sc_tools
  </footer>
</div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.8/js/bootstrap.bundle.min.js"></script>
</body>
</html>"""


def generate_report_index(
    output_dir: str | Path,
    project_name: str = "",
    output_path: str | Path | None = None,
) -> Path:
    """Generate a figures/QC/index.html landing page listing all HTML reports.

    Parameters
    ----------
    output_dir
        Directory to scan for ``*.html`` report files.
    project_name
        Shown in the page title and header.
    output_path
        Where to write ``index.html``; defaults to ``output_dir / "index.html"``.

    Returns
    -------
    Path to the generated index file.
    """
    from jinja2 import Template

    output_dir = Path(output_dir)
    if output_path is None:
        output_path = output_dir / "index.html"
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    reports = []
    for html_file in sorted(output_dir.glob("*.html")):
        if html_file.name == output_path.name:
            continue
        text = html_file.read_text(errors="replace")

        # Title
        m = re.search(r"<title>(.*?)</title>", text, re.DOTALL)
        title = m.group(1).strip() if m else html_file.stem.replace("_", " ").title()

        # Date
        m = re.search(r'class="text-muted"[^>]*>(\d{4}-\d{2}-\d{2}[^<]*)<', text)
        if m:
            date = m.group(1).strip()
        else:
            mtime = html_file.stat().st_mtime
            date = datetime.fromtimestamp(mtime).strftime("%Y-%m-%d")

        # Key metric extraction — try multiple patterns used by different report types.
        # 1. Benchmark/integration reports: .card.best .value (method name or score)
        m = re.search(
            r'class="[^"]*\bcard\b[^"]*\bbest\b[^"]*"[^>]*>.*?'
            r'class="value"[^>]*>(.*?)</div>',
            text,
            re.DOTALL,
        )
        best_method = m.group(1).strip() if m else None

        # 2. QC reports: extract n_cells / n_spots from header span
        if best_method is None:
            m2 = re.search(r"Total (?:cells|spots|observations):\s*([\d,]+)", text)
            if m2:
                best_method = f"{m2.group(1).replace(',', '')} cells"

        if best_method is None:
            best_method = "\u2014"

        # Report type from filename
        name_lower = html_file.name.lower()
        if "benchmark" in name_lower:
            report_type = "benchmark"
        elif "integration" in name_lower:
            report_type = "integration"
        elif "segmentation" in name_lower:
            report_type = "segmentation"
        elif "celltyp" in name_lower or "post_cell" in name_lower:
            report_type = "celltyping"
        else:
            report_type = "qc"

        reports.append(
            {
                "filename": html_file.name,
                "title": title,
                "date": date,
                "best_method": best_method,
                "report_type": report_type,
            }
        )

    reports.sort(key=lambda r: r["date"], reverse=True)

    html = Template(_INDEX_TEMPLATE).render(
        project_name=project_name,
        reports=reports,
        generated_date=datetime.now().strftime("%Y-%m-%d %H:%M"),
    )
    output_path.write_text(html)
    logger.info("Report index saved to %s", output_path)
    return output_path
