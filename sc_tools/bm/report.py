"""
HTML report generation for segmentation and integration benchmarks.

Uses Plotly for interactive charts and matplotlib for static plots
(embedded as base64 PNGs).
"""

from __future__ import annotations

import base64
import io
import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData

__all__ = ["generate_segmentation_report", "generate_integration_report"]

logger = logging.getLogger(__name__)

matplotlib.use("Agg")

_SEG_TEMPLATE_PATH = Path(__file__).parent.parent / "data" / "segmentation_report_template.html"
_INT_TEMPLATE_PATH = Path(__file__).parent.parent / "data" / "integration_report_template.html"


def _fig_to_base64(fig: plt.Figure, dpi: int = 150) -> str:
    """Render a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def _plotly_to_html(fig) -> str:
    """Convert a Plotly figure to an embeddable HTML div."""
    return fig.to_html(include_plotlyjs=False, full_html=False)


def _render_template(template_path: Path, context: dict) -> str:
    """Render a Jinja2 template."""
    try:
        from jinja2 import Template
    except ImportError as e:
        raise ImportError(
            "jinja2 is required for report generation. Install with: pip install jinja2"
        ) from e

    template_text = template_path.read_text()
    template = Template(template_text)
    return template.render(**context)


def generate_segmentation_report(
    comparison_df: pd.DataFrame,
    morphology_dict: dict[str, pd.DataFrame] | None = None,
    marker_quality_dict: dict[str, pd.DataFrame] | None = None,
    intensity_image: np.ndarray | None = None,
    masks: dict[str, np.ndarray] | None = None,
    output_path: str | Path = "segmentation_benchmark_report.html",
    title: str = "Segmentation Benchmark Report",
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

    context = {
        "title": title,
        "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "n_methods": len(comparison_df),
        "has_gt": "detection_f1" in comparison_df.columns,
        "best_method": best_method,
        "best_score": best_score,
        "plot_composite_bar": _plotly_to_html(plot_composite_score_bar(comparison_df)),
        "plot_comparison_table": _plotly_to_html(plot_segmentation_comparison_table(comparison_df)),
    }

    # Optional morphology plot
    if morphology_dict:
        from sc_tools.pl.benchmarking import plot_morphology_distributions

        context["plot_morphology"] = _plotly_to_html(plot_morphology_distributions(morphology_dict))
    else:
        context["plot_morphology"] = None

    # Optional marker SNR plot
    if marker_quality_dict:
        from sc_tools.pl.benchmarking import plot_marker_snr_heatmap

        context["plot_marker_snr"] = _plotly_to_html(plot_marker_snr_heatmap(marker_quality_dict))
    else:
        context["plot_marker_snr"] = None

    # Optional overlay
    if intensity_image is not None and masks is not None:
        from sc_tools.pl.benchmarking import plot_segmentation_overlay

        fig = plot_segmentation_overlay(intensity_image, masks)
        context["overlay_img"] = _fig_to_base64(fig)
    else:
        context["overlay_img"] = None

    html = _render_template(_SEG_TEMPLATE_PATH, context)
    output_path.write_text(html)
    logger.info("Segmentation report saved to %s", output_path)
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
        plot_batch_vs_bio,
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
        "plot_ranking_bar": _plotly_to_html(plot_integration_ranking_bar(comparison_df)),
        "plot_comparison_table": _plotly_to_html(plot_integration_comparison_table(comparison_df)),
    }

    # Radar
    exclude = {"method", "embedding_key", "batch_score", "bio_score", "overall_score"}
    metric_cols = [c for c in comparison_df.columns if c not in exclude]
    if len(metric_cols) >= 3:
        context["plot_radar"] = _plotly_to_html(plot_integration_radar(comparison_df))
    else:
        context["plot_radar"] = None

    # Batch vs Bio scatter
    if "batch_score" in comparison_df.columns and "bio_score" in comparison_df.columns:
        context["plot_batch_vs_bio"] = _plotly_to_html(plot_batch_vs_bio(comparison_df))
    else:
        context["plot_batch_vs_bio"] = None

    # UMAP grid (static matplotlib)
    if adata is not None and embeddings is not None:
        from sc_tools.pl.benchmarking import plot_embedding_comparison_umap

        fig = plot_embedding_comparison_umap(adata, embeddings, color_by=color_by)
        context["umap_img"] = _fig_to_base64(fig)
    else:
        context["umap_img"] = None

    html = _render_template(_INT_TEMPLATE_PATH, context)
    output_path.write_text(html)
    logger.info("Integration report saved to %s", output_path)
    return output_path
