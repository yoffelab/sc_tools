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
    _extract_body_content,
    _extract_head_css,
    _find_latest_report,
    _wrap_with_tabs,
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
    "generate_post_celltyping_report",
    "generate_segmentation_qc_report",
    "generate_all_qc_reports",
]

logger = logging.getLogger(__name__)

_DATA_DIR = Path(__file__).parent.parent / "data"
_TEMPLATE_PATH = _DATA_DIR / "qc_report_template.html"
_PRE_FILTER_TEMPLATE = _DATA_DIR / "pre_filter_qc_template.html"
_POST_FILTER_TEMPLATE = _DATA_DIR / "post_filter_qc_template.html"
_POST_INTEGRATION_TEMPLATE = _DATA_DIR / "post_integration_qc_template.html"
_POST_CELLTYPING_TEMPLATE = _DATA_DIR / "post_celltyping_qc_template.html"


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

    # Tab wrapping: embed pre_filter report as a tab when available
    _prev_tabs: list[tuple[str, str, str]] = []
    _pre_path = _find_latest_report(output_dir, "pre_filter")
    if _pre_path is not None:
        try:
            _raw = _pre_path.read_text()
            _prev_tabs.append(
                ("Report 1: Pre-filter", _extract_body_content(_raw), _extract_head_css(_raw))
            )
        except Exception:
            logger.debug("Failed to embed pre_filter tab in post_filter report", exc_info=True)
    if _prev_tabs:
        _cur_css = _extract_head_css(html)
        _cur_body = _extract_body_content(html)
        html = _wrap_with_tabs("Report 2: Post-filter", _cur_body, _prev_tabs, _cur_css)

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
    comparison_df: pd.DataFrame | None = None,
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
    comparison_df
        Pre-computed integration benchmark DataFrame. When provided, skips
        recomputing ``compare_integrations()`` (saves significant time on
        large datasets).

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

    # Per-embedding UMAP grid (compute UMAP for each integration method)
    if embedding_keys:
        from ..pl.qc_plots import qc_embedding_umap_grid

        color_by = sample_col if sample_col in adata.obs.columns else None
        if color_by is None:
            for cand in ["library_id", "sample", "batch"]:
                if cand in adata.obs.columns:
                    color_by = cand
                    break
        if color_by is not None:
            try:
                logger.info(
                    "Computing per-embedding UMAPs for %d methods...",
                    len(embedding_keys),
                )
                fig = qc_embedding_umap_grid(adata, embedding_keys, color_key=color_by)
                plots["embedding_umaps"] = fig_to_base64(fig)
                logger.info("Per-embedding UMAP grid generated successfully")
            except Exception:
                logger.warning("Per-embedding UMAP grid failed", exc_info=True)

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
        if comparison_df is not None:
            # Use pre-computed benchmark results (skip expensive recomputation)
            result = compute_integration_section(
                adata,
                embedding_keys,
                batch_key,
                celltype_key,
                comparison_df=comparison_df,
            )
        else:
            result = compute_integration_section(
                adata,
                embedding_keys,
                batch_key,
                celltype_key,
            )
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

    # Build ranking table rows for the static HTML table
    ranking_rows: list[dict] = []
    _rank_df = result["comparison_df"] if result is not None else None
    if _rank_df is not None and len(_rank_df) > 0:
        has_group = "group" in _rank_df.columns
        for rank, (_, row) in enumerate(_rank_df.iterrows(), 1):
            entry: dict = {
                "rank": rank,
                "method": str(row.get("method", "")),
                "overall": float(row["overall_score"]) if "overall_score" in row.index else None,
                "batch": float(row["batch_score"]) if "batch_score" in row.index else None,
                "bio": float(row["bio_score"]) if "bio_score" in row.index else None,
                "asw_batch": float(row["asw_batch"]) if "asw_batch" in row.index else None,
                "pcr": float(row["pcr"]) if "pcr" in row.index else None,
                "graph_conn": float(row["graph_connectivity"])
                if "graph_connectivity" in row.index
                else None,
                "asw_celltype": float(row["asw_celltype"]) if "asw_celltype" in row.index else None,
                "ari": float(row["ari"]) if "ari" in row.index else None,
                "nmi": float(row["nmi"]) if "nmi" in row.index else None,
            }
            if has_group:
                entry["group"] = str(row["group"])
            ranking_rows.append(entry)

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
        "ranking_rows": ranking_rows,
        "segmentation": seg_data,
        "seg_plots": seg_plots_dict,
    }

    html = render_template(_POST_INTEGRATION_TEMPLATE, context)

    # Tab wrapping: embed post_filter and pre_filter reports as tabs when available
    _prev_tabs_integ: list[tuple[str, str, str]] = []
    for _rt, _lbl in [
        ("post_filter", "Report 2: Post-filter"),
        ("pre_filter", "Report 1: Pre-filter"),
    ]:
        _p = _find_latest_report(output_dir, _rt)
        if _p is not None:
            try:
                _raw = _p.read_text()
                _prev_tabs_integ.append(
                    (_lbl, _extract_body_content(_raw), _extract_head_css(_raw))
                )
            except Exception:
                logger.debug(
                    "Failed to embed %s tab in post_integration report", _rt, exc_info=True
                )
    if _prev_tabs_integ:
        _cur_css_integ = _extract_head_css(html)
        _cur_body_integ = _extract_body_content(html)
        html = _wrap_with_tabs(
            "Report 3: Post-integration",
            _cur_body_integ,
            _prev_tabs_integ,
            _cur_css_integ,
            current_head_extras='<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>',
        )

    output_path = output_dir / f"post_integration_qc_{ds}.html"
    output_path.write_text(html)
    logger.info("Post-integration QC report written: %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# Post-celltyping QC report
# ---------------------------------------------------------------------------


def generate_post_celltyping_report(
    adata: AnnData,
    output_dir: str | Path,
    *,
    celltype_key: str = "celltype",
    embedding_keys: dict[str, str] | None = None,
    batch_key: str | None = None,
    sample_col: str = "library_id",
    cluster_key: str = "leiden",
    modality: str = "visium",
    title: str = "Post-celltyping QC Report",
    date_stamp: str | None = None,
    segmentation_masks_dir: str | Path | None = None,
    comparison_df: pd.DataFrame | None = None,
    comparison_df_p3: pd.DataFrame | None = None,
    integration_test_dir: str | Path | None = None,
    marker_genes: dict[str, list[str]] | None = None,
) -> Path:
    """Generate a post-celltyping QC HTML report (Phase 4 exit).

    Re-evaluates integration quality using validated cell type labels,
    making bio conservation metrics (ARI, NMI, ASW celltype) meaningful.
    Optionally re-scores all candidate integration embeddings stored in
    ``integration_test_dir``.

    Parameters
    ----------
    adata
        Cell-typed AnnData (with validated cell type labels).
    output_dir
        Directory for the output HTML file.
    celltype_key
        Column in ``adata.obs`` with validated cell type labels.
        **Required** -- raises ``ValueError`` if missing.
    embedding_keys
        Dict mapping method name to ``obsm`` key (auto-detected if None).
    batch_key
        Batch column in ``obs`` (auto-detected from raw_data_dir/batch/library_id).
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
    comparison_df
        Pre-computed integration benchmark DataFrame (Phase 4, validated
        celltypes). When provided, skips recomputing ``compare_integrations()``.
    comparison_df_p3
        Pre-computed integration benchmark DataFrame from Phase 3 (using
        preliminary Leiden labels). Shown alongside Phase 4 values for
        comparison. Optional.
    integration_test_dir
        Path to ``results/tmp/integration_test/`` directory containing
        per-method ``{method}.h5ad`` files from Phase 3 benchmark. If
        provided, embeddings are loaded and re-scored with validated
        cell type labels.
    marker_genes
        Optional dict mapping cell type name to a list of marker genes.
        When provided, a marker dotplot is included in the report.

    Returns
    -------
    Path to the generated HTML file.

    Raises
    ------
    ValueError
        If *celltype_key* is not found in ``adata.obs``.
    """
    import anndata as ad

    from ..qc.report_utils import auto_detect_embeddings

    if celltype_key not in adata.obs.columns:
        raise ValueError(
            f"celltype_key {celltype_key!r} not found in adata.obs. "
            f"Post-celltyping report requires validated cell type labels."
        )

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

    # Load integration test embeddings if available
    if integration_test_dir is not None:
        integration_test_dir = Path(integration_test_dir)
        if integration_test_dir.exists():
            for h5ad_file in sorted(integration_test_dir.glob("*.h5ad")):
                method_name = h5ad_file.stem
                try:
                    test_adata = ad.read_h5ad(h5ad_file)
                    for obsm_key in test_adata.obsm:
                        if obsm_key.startswith("X_") and obsm_key != "X_pca":
                            if obsm_key not in adata.obsm:
                                common = adata.obs_names.intersection(test_adata.obs_names)
                                if len(common) > 0:
                                    import numpy as np

                                    emb = np.zeros(
                                        (adata.n_obs, test_adata.obsm[obsm_key].shape[1]),
                                        dtype=np.float32,
                                    )
                                    idx_main = [adata.obs_names.get_loc(c) for c in common]
                                    idx_test = [test_adata.obs_names.get_loc(c) for c in common]
                                    emb[idx_main] = test_adata.obsm[obsm_key][idx_test]
                                    adata.obsm[obsm_key] = emb
                                    embedding_keys[method_name] = obsm_key
                                    logger.info(
                                        "Loaded embedding %s from %s (%d common cells)",
                                        obsm_key,
                                        h5ad_file.name,
                                        len(common),
                                    )
                            elif obsm_key not in embedding_keys.values():
                                embedding_keys[method_name] = obsm_key
                            break
                except Exception:
                    logger.warning(
                        "Failed to load integration test %s", h5ad_file.name, exc_info=True
                    )

    # Read selected integration method if available
    selected_method = None
    for candidate_dir in [output_dir.parent, output_dir]:
        method_file = candidate_dir / "integration_method.txt"
        if method_file.exists():
            selected_method = method_file.read_text().strip()
            break
    if integration_test_dir is not None:
        method_file = Path(integration_test_dir).parent / "integration_method.txt"
        if method_file.exists():
            selected_method = method_file.read_text().strip()

    has_umap = "X_umap" in adata.obsm
    plots: dict[str, str] = {}

    # Celltype abundance plot
    if sample_col in adata.obs.columns:
        try:
            from ..pl.qc_plots import qc_celltype_abundance

            fig_abund = qc_celltype_abundance(
                adata, celltype_key=celltype_key, sample_col=sample_col
            )
            plots["celltype_abundance"] = fig_to_base64(fig_abund)
        except Exception:
            logger.debug("Celltype abundance plot failed", exc_info=True)

    # UMAP grid (include celltype and celltype_broad if available)
    if has_umap:
        from ..pl.qc_plots import qc_umap_grid

        color_keys = []
        for k in [celltype_key, "celltype_broad", cluster_key, batch_key, sample_col]:
            if k in adata.obs.columns and k not in color_keys:
                color_keys.append(k)

        if color_keys:
            fig = qc_umap_grid(adata, color_keys=color_keys)
            plots["umap_grid"] = fig_to_base64(fig)

    # Cluster distribution (by celltype)
    if cluster_key in adata.obs.columns and sample_col in adata.obs.columns:
        from ..pl.qc_plots import qc_cluster_distribution

        fig = qc_cluster_distribution(adata, cluster_key=cluster_key, sample_col=sample_col)
        plots["cluster_dist"] = fig_to_base64(fig)

    # Marker dotplot (optional)
    if marker_genes:
        try:
            import scanpy as sc

            all_genes = [g for genes in marker_genes.values() for g in genes]
            valid_genes = [g for g in all_genes if g in adata.var_names]
            if valid_genes and celltype_key in adata.obs.columns:
                fig_dot, _ = sc.pl.dotplot(  # type: ignore[misc]
                    adata,
                    var_names=marker_genes,
                    groupby=celltype_key,
                    return_fig=True,
                    show=False,
                )
                plots["marker_dotplot"] = fig_to_base64(fig_dot)
        except Exception:
            logger.debug("Marker dotplot failed", exc_info=True)

    # Integration metrics WITH validated celltypes (bio metrics now meaningful)
    integration_plots: dict[str, str] | None = None
    best_method: str | None = None
    best_score: float | None = None
    best_bio_method: str | None = None
    best_bio_score: float | None = None
    result = None

    if embedding_keys and batch_key in adata.obs.columns:
        result = compute_integration_section(
            adata,
            embedding_keys,
            batch_key,
            celltype_key,
            comparison_df=comparison_df,
        )
        if result is not None:
            integration_plots = result["plots"]
            comp_df = result["comparison_df"]
            if len(comp_df) > 0:
                best_method = str(comp_df.iloc[0]["method"])
                best_score = float(comp_df.iloc[0]["overall_score"])
                # Bio score is primary for Phase 4
                if "bio_score" in comp_df.columns:
                    bio_sorted = comp_df.sort_values("bio_score", ascending=False)
                    best_bio_method = str(bio_sorted.iloc[0]["method"])
                    best_bio_score = float(bio_sorted.iloc[0]["bio_score"])

    # Build ranking_rows sorted by bio_score (primary for Phase 4)
    ranking_rows: list[dict] = []
    _rank_df = result["comparison_df"] if result is not None else None
    if _rank_df is not None and len(_rank_df) > 0:
        # Sort by bio_score descending; fallback to overall_score
        sort_col = "bio_score" if "bio_score" in _rank_df.columns else "overall_score"
        _rank_df_sorted = _rank_df.sort_values(sort_col, ascending=False)

        # Build p3 bio lookup for comparison
        _p3_bio_lookup: dict[str, str] = {}
        if comparison_df_p3 is not None and "bio_score" in comparison_df_p3.columns:
            for _, r3 in comparison_df_p3.iterrows():
                m = str(r3.get("method", r3.name))
                _p3_bio_lookup[m] = f"{r3['bio_score']:.3f}"

        for rank, (_, row) in enumerate(_rank_df_sorted.iterrows(), 1):
            method_name = str(row.get("method", row.name))
            entry: dict = {
                "rank": rank,
                "method": method_name,
                "bio": f"{row['bio_score']:.3f}" if "bio_score" in row.index else "-",
                "bio_p3": _p3_bio_lookup.get(method_name),
                "batch": f"{row['batch_score']:.3f}" if "batch_score" in row.index else "-",
                "ari": f"{row['ari']:.3f}" if "ari" in row.index else "-",
                "nmi": f"{row['nmi']:.3f}" if "nmi" in row.index else "-",
                "asw_cell": f"{row['asw_celltype']:.3f}" if "asw_celltype" in row.index else "-",
                "asw_batch": f"{row['asw_batch']:.3f}" if "asw_batch" in row.index else "-",
                "pcr": f"{row['pcr']:.3f}" if "pcr" in row.index else "-",
                "graph_conn": (
                    f"{row['graph_connectivity']:.3f}" if "graph_connectivity" in row.index else "-"
                ),
                "is_selected": method_name == selected_method,
            }
            ranking_rows.append(entry)

    # Celltype composition table
    ct_counts = adata.obs[celltype_key].value_counts()
    total_cells = int(adata.n_obs)
    celltype_table = [
        {
            "celltype": str(ct),
            "n_cells": int(n),
            "pct": float(n) / total_cells * 100,
        }
        for ct, n in ct_counts.items()
    ]

    # Cluster count
    n_clusters = 0
    if cluster_key in adata.obs.columns:
        n_clusters = int(adata.obs[cluster_key].nunique())

    n_celltypes = int(adata.obs[celltype_key].nunique())

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
        "n_celltypes": n_celltypes,
        "n_embeddings": len(embedding_keys) if embedding_keys else 0,
        "best_method": best_method,
        "best_score": best_score,
        "best_bio_method": best_bio_method,
        "best_bio_score": best_bio_score,
        "selected_method": selected_method,
        "has_celltype": True,
        "celltype_key": celltype_key,
        "show_p3_bio": bool(_p3_bio_lookup) if "_p3_bio_lookup" in dir() else False,
        "ranking_rows": ranking_rows,
        "celltype_table": celltype_table,
        "plots": plots,
        "integration_plots": integration_plots,
        "segmentation": seg_data,
        "seg_plots": seg_plots_dict,
    }

    html = render_template(_POST_CELLTYPING_TEMPLATE, context)

    # Tab wrapping: embed post_integration, post_filter, pre_filter reports when available
    _prev_tabs_ct: list[tuple[str, str, str]] = []
    for _rt, _lbl in [
        ("post_integration", "Report 3: Post-integration"),
        ("post_filter", "Report 2: Post-filter"),
        ("pre_filter", "Report 1: Pre-filter"),
    ]:
        _p = _find_latest_report(output_dir, _rt)
        if _p is not None:
            try:
                _raw = _p.read_text()
                _prev_tabs_ct.append((_lbl, _extract_body_content(_raw), _extract_head_css(_raw)))
            except Exception:
                logger.debug("Failed to embed %s tab in post_celltyping report", _rt, exc_info=True)
    if _prev_tabs_ct:
        _cur_css_ct = _extract_head_css(html)
        _cur_body_ct = _extract_body_content(html)
        html = _wrap_with_tabs(
            "Report 4: Post-celltyping",
            _cur_body_ct,
            _prev_tabs_ct,
            _cur_css_ct,
            current_head_extras='<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>',
        )

    output_path = output_dir / f"post_celltyping_qc_{ds}.html"
    output_path.write_text(html)
    logger.info("Post-celltyping QC report written: %s", output_path)
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
    **kwargs
        Remaining keyword arguments (``sample_col``, ``modality``, ``date_stamp``,
        ``embedding_keys``, ``batch_key``, ``celltype_key``, ``cluster_key``,
        ``segmentation_masks_dir``) are passed through to individual report
        generators.

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
