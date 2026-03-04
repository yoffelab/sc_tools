"""
QC HTML report generation.

Produces a self-contained HTML file with per-sample metrics table,
inline-generated QC plots (embedded as base64 PNGs), and pass/fail
summary in a structured 6-row layout.
"""

from __future__ import annotations

import base64
import io
import logging
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from anndata import AnnData  # noqa: E402

__all__ = ["generate_qc_report"]

logger = logging.getLogger(__name__)

_TEMPLATE_PATH = Path(__file__).parent.parent / "data" / "qc_report_template.html"


def _fig_to_base64(fig: plt.Figure, dpi: int = 150) -> str:
    """Render a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def generate_qc_report(
    adata: AnnData,
    metrics: pd.DataFrame,
    classified: pd.DataFrame,
    figures_dir: str | Path,
    output_path: str | Path,
    sample_col: str = "library_id",
    modality: str = "visium",
    title: str = "QC Report",
    adata_post: AnnData | None = None,
) -> Path:
    """
    Generate a self-contained HTML QC report with structured 6-row layout.

    Plots are generated inline (not globbed from disk) and embedded as
    base64 PNGs. The layout has 6 rows:

    1. QC 2x2 pre-filter | QC 2x2 post-filter
    2. QC violin pre-filter | QC violin post-filter
    3. pct_counts_mt per-sample | (empty)
    4. Cross-sample comparison bar (linear) | Same (log10)
    5. QC sample violin grouped (linear) | Same (log10)
    6. QC scatter matrix (full width)

    Parameters
    ----------
    adata : AnnData
        Pre-filter AnnData used for QC.
    metrics : pd.DataFrame
        Output of ``compute_sample_metrics``.
    classified : pd.DataFrame
        Output of ``classify_samples`` (with ``qc_pass``, ``qc_fail_reasons``).
    figures_dir : str or Path
        Base figures directory (unused for plot generation but kept for API compat).
    output_path : str or Path
        Path for the output HTML file.
    sample_col : str
        Sample column name.
    modality : str
        Modality string for display.
    title : str
        Report title.
    adata_post : AnnData or None
        Post-filter AnnData. If provided, row 1 and 2 right panels are populated.

    Returns
    -------
    Path
        Path to the generated HTML file.
    """
    try:
        from jinja2 import Template
    except ImportError as e:
        raise ImportError(
            "jinja2 is required for HTML report generation. Install with: pip install jinja2"
        ) from e

    from .plots import (
        qc_2x2_grid,
        qc_pct_mt_per_sample,
        qc_sample_comparison_bar,
        qc_sample_scatter_matrix,
        qc_sample_violin_grouped,
        qc_violin_metrics,
    )

    figures_dir = Path(figures_dir)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Read template
    template_text = _TEMPLATE_PATH.read_text()
    template = Template(template_text)

    # Build table rows
    table_rows = []
    for sample in classified.index:
        row = {"sample": str(sample)}
        row["qc_pass"] = bool(classified.loc[sample, "qc_pass"])
        row["qc_fail_reasons"] = (
            str(classified.loc[sample, "qc_fail_reasons"])
            if classified.loc[sample, "qc_fail_reasons"]
            else ""
        )
        for col in [
            "n_spots",
            "total_counts_median",
            "n_genes_median",
            "n_genes_detected",
            "pct_mt_median",
            "pct_mt_gt5",
        ]:
            val = classified.loc[sample, col] if col in classified.columns else None
            if val is not None and not (isinstance(val, float) and np.isnan(val)):
                row[col] = val
            else:
                row[col] = None
        table_rows.append(row)

    # Generate plots inline
    named_plots = {}

    # Row 1: QC 2x2 pre/post
    named_plots["qc_2x2_pre"] = _fig_to_base64(qc_2x2_grid(adata))
    if adata_post is not None:
        named_plots["qc_2x2_post"] = _fig_to_base64(qc_2x2_grid(adata_post))

    # Row 2: Violin pre/post
    named_plots["violin_pre"] = _fig_to_base64(qc_violin_metrics(adata))
    if adata_post is not None:
        named_plots["violin_post"] = _fig_to_base64(qc_violin_metrics(adata_post))

    # Row 3: pct_mt per sample
    named_plots["pct_mt"] = _fig_to_base64(
        qc_pct_mt_per_sample(adata, sample_col=sample_col, classified=classified)
    )

    # Row 4: Cross-sample comparison bar (linear + log10)
    named_plots["comparison_bar"] = _fig_to_base64(
        qc_sample_comparison_bar(metrics, classified=classified)
    )
    named_plots["comparison_bar_log"] = _fig_to_base64(
        qc_sample_comparison_bar(metrics, classified=classified, log_scale=True)
    )

    # Row 5: Sample violin grouped (linear + log10)
    named_plots["violin_grouped"] = _fig_to_base64(
        qc_sample_violin_grouped(adata, sample_col=sample_col, classified=classified)
    )
    named_plots["violin_grouped_log"] = _fig_to_base64(
        qc_sample_violin_grouped(
            adata, sample_col=sample_col, classified=classified, log_scale=True
        )
    )

    # Row 6: Scatter matrix
    named_plots["scatter_matrix"] = _fig_to_base64(
        qc_sample_scatter_matrix(metrics, classified=classified)
    )

    # Compute summary stats
    n_pass = int(classified["qc_pass"].sum())
    n_fail = int((~classified["qc_pass"]).sum())
    n_spots_total = int(adata.n_obs)
    median_genes = "N/A"
    if "n_genes_median" in classified.columns:
        median_genes = f"{classified['n_genes_median'].median():.0f}"

    has_mt = "pct_mt_median" in classified.columns

    html = template.render(
        title=title,
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        modality=modality,
        n_samples=len(classified),
        n_spots_total=n_spots_total,
        n_pass=n_pass,
        n_fail=n_fail,
        median_genes=median_genes,
        has_mt=has_mt,
        table_rows=table_rows,
        plots=named_plots,
    )

    output_path.write_text(html)
    logger.info("QC report written: %s", output_path)
    return output_path
