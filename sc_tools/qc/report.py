"""
QC HTML report generation.

Produces a self-contained HTML file with per-sample metrics table,
cross-sample comparison plots (embedded as base64 PNGs), and pass/fail
summary.
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


def _collect_plot_pngs(figures_dir: Path) -> list[dict[str, str]]:
    """Find existing QC PNGs in figures_dir/QC/ and encode as base64."""
    plots = []
    qc_dir = figures_dir / "QC"
    if not qc_dir.exists():
        return plots
    for subdir in ["raw", "post", ""]:
        search_dir = qc_dir / subdir if subdir else qc_dir
        if not search_dir.exists():
            continue
        for png in sorted(search_dir.glob("*.png")):
            data = base64.b64encode(png.read_bytes()).decode("ascii")
            caption = png.stem.replace("_", " ").title()
            plots.append({"data": data, "caption": caption})
    return plots


def generate_qc_report(
    adata: AnnData,
    metrics: pd.DataFrame,
    classified: pd.DataFrame,
    figures_dir: str | Path,
    output_path: str | Path,
    sample_col: str = "library_id",
    modality: str = "visium",
    title: str = "QC Report",
) -> Path:
    """
    Generate a self-contained HTML QC report.

    Parameters
    ----------
    adata : AnnData
        AnnData used for QC (for total spot count).
    metrics : pd.DataFrame
        Output of ``compute_sample_metrics``.
    classified : pd.DataFrame
        Output of ``classify_samples`` (with ``qc_pass``, ``qc_fail_reasons``).
    figures_dir : str or Path
        Base figures directory (contains ``QC/raw/``, ``QC/post/``).
    output_path : str or Path
        Path for the output HTML file.
    sample_col : str
        Sample column name.
    modality : str
        Modality string for display.
    title : str
        Report title.

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

    # Collect embedded plots
    plots = _collect_plot_pngs(figures_dir)

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
        plots=plots,
    )

    output_path.write_text(html)
    logger.info("QC report written: %s", output_path)
    return output_path
