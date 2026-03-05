"""
QC HTML report generation.

Produces self-contained HTML files with per-sample metrics tables,
inline-generated QC plots (embedded as base64 PNGs), and pass/fail
summaries. Three report types align with pipeline phases:

- **Pre-filter QC** (Phase 1 entry): raw data distributions
- **Post-filter QC** (Phase 1-2 exit): pre vs post comparison
- **Post-integration QC** (Phase 3 exit): UMAP, clusters, integration metrics
"""

from __future__ import annotations

import logging
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from anndata import AnnData  # noqa: E402

from .report_utils import (  # noqa: E402
    build_metrics_table_rows,
    compute_integration_section,
    compute_segmentation_section,
    fig_to_base64,
    get_date_stamp,
    get_modality_terms,
    render_template,
)

__all__ = [
    "generate_qc_report",
    "generate_pre_filter_report",
    "generate_post_filter_report",
    "generate_post_integration_report",
    "generate_segmentation_qc_report",
    "generate_all_qc_reports",
]

logger = logging.getLogger(__name__)

_DATA_DIR = Path(__file__).parent.parent / "data"
_TEMPLATE_PATH = _DATA_DIR / "qc_report_template.html"
_PRE_FILTER_TEMPLATE = _DATA_DIR / "pre_filter_qc_template.html"
_POST_FILTER_TEMPLATE = _DATA_DIR / "post_filter_qc_template.html"
_POST_INTEGRATION_TEMPLATE = _DATA_DIR / "post_integration_qc_template.html"


# Keep for backward compat — delegates to fig_to_base64 from report_utils
def _fig_to_base64(fig: plt.Figure, dpi: int = 150) -> str:
    """Render a matplotlib figure to a base64-encoded PNG string."""
    return fig_to_base64(fig, dpi=dpi)


# ---------------------------------------------------------------------------
# Pre-filter QC report
# ---------------------------------------------------------------------------


def generate_pre_filter_report(
    adata: AnnData,
    metrics: pd.DataFrame,
    classified: pd.DataFrame,
    output_dir: str | Path,
    *,
    sample_col: str = "library_id",
    modality: str = "visium",
    title: str = "Pre-filter QC Report",
    date_stamp: str | None = None,
    segmentation_masks_dir: str | Path | None = None,
) -> Path:
    """Generate a pre-filter QC HTML report (Phase 1 entry).

    Parameters
    ----------
    adata
        Pre-filter AnnData (raw counts, concatenated).
    metrics
        Output of ``compute_sample_metrics``.
    classified
        Output of ``classify_samples``.
    output_dir
        Directory for the output HTML file.
    sample_col
        Column in ``adata.obs`` identifying samples.
    modality
        Modality string for display.
    title
        Report title.
    date_stamp
        YYYYMMDD string (default: today).
    segmentation_masks_dir
        Optional path to mask TIFFs for segmentation scoring.

    Returns
    -------
    Path to the generated HTML file.
    """
    from .plots import (
        qc_2x2_grid,
        qc_pct_mt_per_sample,
        qc_sample_comparison_bar,
        qc_sample_scatter_matrix,
        qc_sample_violin_grouped,
        qc_scatter_counts_genes,
        qc_violin_metrics,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ds = get_date_stamp(date_stamp)

    # Build plots
    plots: dict[str, str] = {}
    plots["qc_2x2"] = fig_to_base64(qc_2x2_grid(adata))
    plots["scatter_counts"] = fig_to_base64(qc_scatter_counts_genes(adata))
    plots["violin"] = fig_to_base64(qc_violin_metrics(adata))

    if "pct_counts_mt" in adata.obs.columns and sample_col in adata.obs.columns:
        plots["pct_mt"] = fig_to_base64(
            qc_pct_mt_per_sample(adata, sample_col=sample_col, classified=classified)
        )

    plots["comparison_bar"] = fig_to_base64(
        qc_sample_comparison_bar(metrics, classified=classified)
    )
    plots["comparison_bar_log"] = fig_to_base64(
        qc_sample_comparison_bar(metrics, classified=classified, log_scale=True)
    )
    plots["violin_grouped"] = fig_to_base64(
        qc_sample_violin_grouped(adata, sample_col=sample_col, classified=classified)
    )
    plots["violin_grouped_log"] = fig_to_base64(
        qc_sample_violin_grouped(
            adata, sample_col=sample_col, classified=classified, log_scale=True
        )
    )
    plots["scatter_matrix"] = fig_to_base64(
        qc_sample_scatter_matrix(metrics, classified=classified)
    )

    # Summary stats
    n_pass = int(classified["qc_pass"].sum())
    n_fail = int((~classified["qc_pass"]).sum())
    median_genes = "N/A"
    if "n_genes_median" in classified.columns:
        median_genes = f"{classified['n_genes_median'].median():.0f}"
    has_mt = "pct_mt_median" in classified.columns

    # Segmentation
    seg_data = None
    seg_plots: dict[str, str] | None = None
    if segmentation_masks_dir is not None:
        result = compute_segmentation_section(adata, segmentation_masks_dir, sample_col)
        if result is not None:
            seg_data = result["df"]
            seg_plots = result["plots"]

    terms = get_modality_terms(modality)

    context = {
        "title": title,
        "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "modality": modality,
        "terms": terms,
        "n_samples": len(classified),
        "n_spots_total": int(adata.n_obs),
        "n_pass": n_pass,
        "n_fail": n_fail,
        "median_genes": median_genes,
        "has_mt": has_mt,
        "table_rows": build_metrics_table_rows(classified),
        "plots": plots,
        "segmentation": seg_data,
        "seg_plots": seg_plots,
    }

    html = render_template(_PRE_FILTER_TEMPLATE, context)
    output_path = output_dir / f"pre_filter_qc_{ds}.html"
    output_path.write_text(html)
    logger.info("Pre-filter QC report written: %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# Post-filter QC report
# ---------------------------------------------------------------------------


def generate_post_filter_report(
    adata_pre: AnnData,
    adata_post: AnnData,
    metrics: pd.DataFrame,
    classified: pd.DataFrame,
    output_dir: str | Path,
    *,
    sample_col: str = "library_id",
    modality: str = "visium",
    title: str = "Post-filter QC Report",
    date_stamp: str | None = None,
    segmentation_masks_dir: str | Path | None = None,
) -> Path:
    """Generate a post-filter QC HTML report (Phase 1-2 exit).

    Parameters
    ----------
    adata_pre
        Pre-filter AnnData.
    adata_post
        Post-filter/annotated AnnData.
    metrics
        Output of ``compute_sample_metrics``.
    classified
        Output of ``classify_samples``.
    output_dir
        Directory for the output HTML file.
    sample_col
        Column in ``adata.obs`` identifying samples.
    modality
        Modality string for display.
    title
        Report title.
    date_stamp
        YYYYMMDD string (default: today).
    segmentation_masks_dir
        Optional path to mask TIFFs for segmentation scoring.

    Returns
    -------
    Path to the generated HTML file.
    """
    from .plots import (
        qc_2x4_pre_post,
        qc_sample_comparison_bar,
        qc_violin_metrics,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ds = get_date_stamp(date_stamp)

    plots: dict[str, str] = {}

    # 2x4 pre vs post
    plots["qc_2x4"] = fig_to_base64(qc_2x4_pre_post(adata_pre, adata_post))

    # Violin pre/post
    plots["violin_pre"] = fig_to_base64(qc_violin_metrics(adata_pre))
    plots["violin_post"] = fig_to_base64(qc_violin_metrics(adata_post))

    # Cross-sample comparison on post data
    plots["comparison_bar"] = fig_to_base64(
        qc_sample_comparison_bar(metrics, classified=classified)
    )
    plots["comparison_bar_log"] = fig_to_base64(
        qc_sample_comparison_bar(metrics, classified=classified, log_scale=True)
    )

    # HVG plot (if available in post)
    if "highly_variable" in adata_post.var.columns:
        from .plots import plot_highly_variable_genes

        plots["hvg"] = fig_to_base64(plot_highly_variable_genes(adata_post))

    # SVG plot (if available in post)
    if "spatial_i" in adata_post.var.columns:
        from .plots import plot_spatially_variable_genes

        plots["svg"] = fig_to_base64(plot_spatially_variable_genes(adata_post))

    n_pass = int(classified["qc_pass"].sum())
    n_fail = int((~classified["qc_pass"]).sum())
    has_mt = "pct_mt_median" in classified.columns

    # Segmentation
    seg_data = None
    seg_plots_dict: dict[str, str] | None = None
    if segmentation_masks_dir is not None:
        result = compute_segmentation_section(adata_post, segmentation_masks_dir, sample_col)
        if result is not None:
            seg_data = result["df"]
            seg_plots_dict = result["plots"]

    terms = get_modality_terms(modality)

    context = {
        "title": title,
        "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "modality": modality,
        "terms": terms,
        "n_samples": len(classified),
        "n_spots_pre": int(adata_pre.n_obs),
        "n_spots_post": int(adata_post.n_obs),
        "n_genes_pre": int(adata_pre.n_vars),
        "n_genes_post": int(adata_post.n_vars),
        "n_pass": n_pass,
        "n_fail": n_fail,
        "has_mt": has_mt,
        "table_rows": build_metrics_table_rows(classified),
        "plots": plots,
        "segmentation": seg_data,
        "seg_plots": seg_plots_dict,
    }

    html = render_template(_POST_FILTER_TEMPLATE, context)
    output_path = output_dir / f"post_filter_qc_{ds}.html"
    output_path.write_text(html)
    logger.info("Post-filter QC report written: %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# Post-integration QC report
# ---------------------------------------------------------------------------


def generate_post_integration_report(
    adata: AnnData,
    output_dir: str | Path,
    *,
    embedding_keys: dict[str, str] | None = None,
    batch_key: str | None = None,
    celltype_key: str | None = None,
    sample_col: str = "library_id",
    cluster_key: str = "leiden",
    modality: str = "visium",
    title: str = "Post-integration QC Report",
    date_stamp: str | None = None,
    segmentation_masks_dir: str | Path | None = None,
) -> Path:
    """Generate a post-integration QC HTML report (Phase 3 exit).

    Parameters
    ----------
    adata
        Normalized/integrated AnnData.
    output_dir
        Directory for the output HTML file.
    embedding_keys
        Dict mapping method name to ``obsm`` key (auto-detected if None).
    batch_key
        Batch column in ``obs`` (auto-detected from raw_data_dir/batch/library_id).
    celltype_key
        Cell type column in ``obs`` (optional; skip bio metrics if absent).
    sample_col
        Sample column for cluster distribution plot.
    cluster_key
        Cluster column (default: ``leiden``).
    modality
        Modality string for display.
    title
        Report title.
    date_stamp
        YYYYMMDD string (default: today).
    segmentation_masks_dir
        Optional path to mask TIFFs for segmentation scoring.

    Returns
    -------
    Path to the generated HTML file.
    """
    from ..qc.report_utils import auto_detect_embeddings

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ds = get_date_stamp(date_stamp)

    # Auto-detect batch key
    if batch_key is None:
        for candidate in ["raw_data_dir", "batch", "library_id", "sample"]:
            if candidate in adata.obs.columns:
                batch_key = candidate
                break
    if batch_key is None:
        batch_key = sample_col

    # Auto-detect embeddings
    if embedding_keys is None:
        embedding_keys = auto_detect_embeddings(adata)

    # Resolve celltype
    if celltype_key is not None and celltype_key not in adata.obs.columns:
        logger.warning("celltype_key %r not in adata.obs; skipping bio metrics", celltype_key)
        celltype_key = None

    has_celltype = celltype_key is not None
    has_umap = "X_umap" in adata.obsm

    plots: dict[str, str] = {}

    # UMAP grid
    if has_umap:
        from ..pl.qc_plots import qc_umap_grid

        color_keys = []
        for k in [sample_col, batch_key, cluster_key]:
            if k in adata.obs.columns and k not in color_keys:
                color_keys.append(k)
        if has_celltype and celltype_key not in color_keys:
            color_keys.append(celltype_key)

        if color_keys:
            fig = qc_umap_grid(adata, color_keys=color_keys)
            plots["umap_grid"] = fig_to_base64(fig)

    # Cluster distribution
    if cluster_key in adata.obs.columns and sample_col in adata.obs.columns:
        from ..pl.qc_plots import qc_cluster_distribution

        fig = qc_cluster_distribution(adata, cluster_key=cluster_key, sample_col=sample_col)
        plots["cluster_dist"] = fig_to_base64(fig)

    # Integration metrics
    integration_plots: dict[str, str] | None = None
    best_method: str | None = None
    best_score: float | None = None

    if embedding_keys and batch_key in adata.obs.columns:
        result = compute_integration_section(adata, embedding_keys, batch_key, celltype_key)
        if result is not None:
            integration_plots = result["plots"]
            comp_df = result["comparison_df"]
            if len(comp_df) > 0:
                best_method = str(comp_df.iloc[0]["method"])
                best_score = float(comp_df.iloc[0]["overall_score"])

    # Cluster count
    n_clusters = 0
    if cluster_key in adata.obs.columns:
        n_clusters = int(adata.obs[cluster_key].nunique())

    # Number of samples
    n_samples = 1
    if sample_col in adata.obs.columns:
        n_samples = int(adata.obs[sample_col].nunique())

    # Segmentation
    seg_data = None
    seg_plots_dict: dict[str, str] | None = None
    if segmentation_masks_dir is not None:
        seg_result = compute_segmentation_section(adata, segmentation_masks_dir, sample_col)
        if seg_result is not None:
            seg_data = seg_result["df"]
            seg_plots_dict = seg_result["plots"]

    terms = get_modality_terms(modality)

    context = {
        "title": title,
        "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "modality": modality,
        "terms": terms,
        "n_samples": n_samples,
        "n_spots": int(adata.n_obs),
        "n_clusters": n_clusters,
        "n_embeddings": len(embedding_keys) if embedding_keys else 0,
        "best_method": best_method,
        "best_score": best_score,
        "has_celltype": has_celltype,
        "plots": plots,
        "integration_plots": integration_plots,
        "segmentation": seg_data,
        "seg_plots": seg_plots_dict,
    }

    html = render_template(_POST_INTEGRATION_TEMPLATE, context)
    output_path = output_dir / f"post_integration_qc_{ds}.html"
    output_path.write_text(html)
    logger.info("Post-integration QC report written: %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# Segmentation QC report
# ---------------------------------------------------------------------------


def generate_segmentation_qc_report(
    adata: AnnData,
    masks_dir: str | Path,
    output_dir: str | Path,
    *,
    sample_col: str = "library_id",
    modality: str = "visium",
    title: str = "Segmentation QC Report",
    date_stamp: str | None = None,
) -> Path | None:
    """Generate a standalone date-versioned segmentation QC HTML report.

    Discovers mask files in *masks_dir*, computes per-ROI segmentation
    quality scores via ``sc_tools.bm.segmentation.score_segmentation``,
    and optionally compares multiple mask methods per ROI using
    ``compare_segmentations``.

    Parameters
    ----------
    adata
        AnnData with sample annotations (used for context only).
    masks_dir
        Directory containing mask TIFF / NPY / NPZ files.
    output_dir
        Directory for the output HTML file.
    sample_col
        Column in ``adata.obs`` identifying samples.
    modality
        Modality string for display.
    title
        Report title.
    date_stamp
        YYYYMMDD string (default: today).

    Returns
    -------
    Path or None
        Path to the generated HTML file, or ``None`` if no masks found.
    """
    result = compute_segmentation_section(adata, masks_dir, sample_col)
    if result is None:
        logger.info("No masks found in %s; skipping segmentation QC report", masks_dir)
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ds = get_date_stamp(date_stamp)

    # Attempt to use the full bm.report for richer HTML
    try:
        from sc_tools.bm.report import generate_segmentation_report

        output_path = output_dir / f"segmentation_qc_{ds}.html"
        generate_segmentation_report(
            comparison_df=result["df"],
            output_path=output_path,
            title=title,
        )
        logger.info("Segmentation QC report written: %s", output_path)
        return output_path
    except Exception:
        logger.debug(
            "bm.report.generate_segmentation_report failed; falling back to inline",
            exc_info=True,
        )

    # Fallback: render a minimal inline report using the segmentation section
    # from the pre_filter template (which already handles the segmentation block)
    terms = get_modality_terms(modality)
    from datetime import datetime as _dt

    context = {
        "title": title,
        "date": _dt.now().strftime("%Y-%m-%d %H:%M"),
        "modality": modality,
        "terms": terms,
        "n_samples": int(adata.obs[sample_col].nunique()) if sample_col in adata.obs.columns else 1,
        "n_spots_total": int(adata.n_obs),
        "n_pass": 0,
        "n_fail": 0,
        "median_genes": "N/A",
        "has_mt": False,
        "table_rows": [],
        "plots": {},
        "segmentation": result["df"],
        "seg_plots": result.get("plots"),
    }

    html = render_template(_PRE_FILTER_TEMPLATE, context)
    output_path = output_dir / f"segmentation_qc_{ds}.html"
    output_path.write_text(html)
    logger.info("Segmentation QC report (fallback) written: %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def generate_all_qc_reports(
    adata_pre: AnnData,
    metrics: pd.DataFrame,
    classified: pd.DataFrame,
    output_dir: str | Path,
    *,
    adata_post: AnnData | None = None,
    adata_integrated: AnnData | None = None,
    sample_col: str = "library_id",
    modality: str = "visium",
    date_stamp: str | None = None,
    embedding_keys: dict[str, str] | None = None,
    batch_key: str | None = None,
    celltype_key: str | None = None,
    cluster_key: str = "leiden",
    segmentation_masks_dir: str | Path | None = None,
) -> dict[str, Path]:
    """Generate all three QC reports in one call.

    Parameters
    ----------
    adata_pre
        Pre-filter AnnData (for pre-filter report).
    metrics
        Output of ``compute_sample_metrics``.
    classified
        Output of ``classify_samples``.
    output_dir
        Directory for all HTML reports.
    adata_post
        Post-filter AnnData (for post-filter report). Skipped if None.
    adata_integrated
        Post-integration AnnData (for post-integration report). Skipped if None.
    sample_col, modality, date_stamp, embedding_keys, batch_key, celltype_key,
    cluster_key, segmentation_masks_dir
        Passed through to individual report generators.

    Returns
    -------
    dict[str, Path]
        Mapping of report type to output path.
    """
    ds = get_date_stamp(date_stamp)
    results: dict[str, Path] = {}

    results["pre_filter"] = generate_pre_filter_report(
        adata_pre,
        metrics,
        classified,
        output_dir,
        sample_col=sample_col,
        modality=modality,
        date_stamp=ds,
        segmentation_masks_dir=segmentation_masks_dir,
    )

    if adata_post is not None:
        results["post_filter"] = generate_post_filter_report(
            adata_pre,
            adata_post,
            metrics,
            classified,
            output_dir,
            sample_col=sample_col,
            modality=modality,
            date_stamp=ds,
            segmentation_masks_dir=segmentation_masks_dir,
        )

    if adata_integrated is not None:
        results["post_integration"] = generate_post_integration_report(
            adata_integrated,
            output_dir,
            embedding_keys=embedding_keys,
            batch_key=batch_key,
            celltype_key=celltype_key,
            sample_col=sample_col,
            cluster_key=cluster_key,
            modality=modality,
            date_stamp=ds,
            segmentation_masks_dir=segmentation_masks_dir,
        )

    return results


# ---------------------------------------------------------------------------
# Legacy API (backward compat with deprecation warning)
# ---------------------------------------------------------------------------


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
    """Generate a self-contained HTML QC report (legacy API).

    .. deprecated::
        Use ``generate_pre_filter_report``, ``generate_post_filter_report``,
        or ``generate_post_integration_report`` instead. This function is
        kept for backward compatibility.

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
        Post-filter AnnData.

    Returns
    -------
    Path to the generated HTML file.
    """
    warnings.warn(
        "generate_qc_report() is deprecated. Use generate_pre_filter_report(), "
        "generate_post_filter_report(), or generate_post_integration_report() instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    from .plots import (
        qc_2x2_grid,
        qc_pct_mt_per_sample,
        qc_sample_comparison_bar,
        qc_sample_scatter_matrix,
        qc_sample_violin_grouped,
        qc_violin_metrics,
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build table rows
    table_rows = build_metrics_table_rows(classified)

    # Generate plots inline
    named_plots: dict[str, str] = {}
    named_plots["qc_2x2_pre"] = fig_to_base64(qc_2x2_grid(adata))
    if adata_post is not None:
        named_plots["qc_2x2_post"] = fig_to_base64(qc_2x2_grid(adata_post))

    named_plots["violin_pre"] = fig_to_base64(qc_violin_metrics(adata))
    if adata_post is not None:
        named_plots["violin_post"] = fig_to_base64(qc_violin_metrics(adata_post))

    named_plots["pct_mt"] = fig_to_base64(
        qc_pct_mt_per_sample(adata, sample_col=sample_col, classified=classified)
    )
    named_plots["comparison_bar"] = fig_to_base64(
        qc_sample_comparison_bar(metrics, classified=classified)
    )
    named_plots["comparison_bar_log"] = fig_to_base64(
        qc_sample_comparison_bar(metrics, classified=classified, log_scale=True)
    )
    named_plots["violin_grouped"] = fig_to_base64(
        qc_sample_violin_grouped(adata, sample_col=sample_col, classified=classified)
    )
    named_plots["violin_grouped_log"] = fig_to_base64(
        qc_sample_violin_grouped(
            adata, sample_col=sample_col, classified=classified, log_scale=True
        )
    )
    named_plots["scatter_matrix"] = fig_to_base64(
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

    html = render_template(
        _TEMPLATE_PATH,
        {
            "title": title,
            "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
            "modality": modality,
            "n_samples": len(classified),
            "n_spots_total": n_spots_total,
            "n_pass": n_pass,
            "n_fail": n_fail,
            "median_genes": median_genes,
            "has_mt": has_mt,
            "table_rows": table_rows,
            "plots": named_plots,
        },
    )

    output_path.write_text(html)
    logger.info("QC report written: %s", output_path)
    return output_path
