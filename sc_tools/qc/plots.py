"""
QC plots: 2x2 metric grid, 2x4 pre vs post, multipage spatial QC, and Scanpy-style
violin, scatter, HVG, and SVG plots.

- 2x2 grid: total_counts, n_genes_by_counts, log1p(total_counts), pct_counts_mt.
- 2x4 pre/post: left 2x2 = pre-filter metrics, right 2x2 = post-filter.
- Multipage spatial: one page per sample (total_count, log1p, % mt); common_scale=True.
- qc_violin_metrics: multi-panel violin for n_genes_by_counts, total_counts, pct_counts_mt.
- qc_scatter_counts_genes: scatter total_counts vs n_genes_by_counts colored by pct_counts_mt.
- plot_highly_variable_genes: mean vs dispersion with HVG highlighted.
- plot_spatially_variable_genes: scatter mean/rank vs Moran's I, colored by spatially_variable/pval.

Save under project figures/QC/raw/ (pre) or figures/QC/post/ (post).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

__all__ = [
    "qc_2x2_grid",
    "qc_2x4_pre_post",
    "qc_spatial_multipage",
    "qc_violin_metrics",
    "qc_scatter_counts_genes",
    "plot_highly_variable_genes",
    "plot_spatially_variable_genes",
    "qc_sample_comparison_bar",
    "qc_sample_violin_grouped",
    "qc_sample_scatter_matrix",
    "qc_pct_mt_per_sample",
]


def qc_2x2_grid(
    adata: AnnData,
    *,
    total_counts_col: str = "total_counts",
    n_genes_col: str = "n_genes_by_counts",
    pct_mt_col: str = "pct_counts_mt",
    output_dir: str | Path | None = None,
    basename: str = "qc_2x2",
    dpi: int = 300,
    figsize: tuple[float, float] = (10, 10),
    modality: str = "visium",
) -> plt.Figure:
    """
    Plot 2x2 QC grid: total_counts, n_genes, log1p(total_counts), pct_counts_mt.

    Panels: (1,1) total_counts histogram, (1,2) n_genes_by_counts histogram,
    (2,1) log1p(total_counts) histogram, (2,2) pct_counts_mt histogram if present,
    else log1p(n_genes_by_counts).

    Parameters
    ----------
    adata : AnnData
        Annotated data with obs containing total_counts and n_genes_by_counts
        (from calculate_qc_metrics). pct_counts_mt optional.
    total_counts_col : str
        Obs column for total counts (default 'total_counts').
    n_genes_col : str
        Obs column for number of genes (default 'n_genes_by_counts').
    pct_mt_col : str
        Obs column for percent mitochondrial (default 'pct_counts_mt').
    output_dir : str or Path or None
        If set, save PDF and PNG here (default None).
    basename : str
        Base name for files (default 'qc_2x2').
    dpi : int
        DPI for PNG (default 300).
    figsize : tuple
        Figure size (default (10, 10)).

    Returns
    -------
    matplotlib.figure.Figure
        The figure (caller may show or save).
    """
    from .report_utils import get_modality_terms

    terms = get_modality_terms(modality)

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    for ax in axes.flat:
        ax.set_visible(True)

    # (1,1) total_counts
    _intensity_title = terms["intensity_label"]
    _feat_per_obs = f"{terms['features']} per {terms['observation_lower']}"
    if total_counts_col in adata.obs.columns:
        x = adata.obs[total_counts_col].values
        x = x[~np.isnan(x) & np.isfinite(x)]
        axes[0, 0].hist(
            x, bins=min(80, max(20, len(x) // 20)), color="steelblue", edgecolor="white"
        )
        axes[0, 0].set_title(_intensity_title, fontsize=12, fontweight="bold")
        axes[0, 0].set_xlabel(total_counts_col)
    else:
        axes[0, 0].text(
            0.5,
            0.5,
            f"{total_counts_col} not in obs",
            ha="center",
            va="center",
            transform=axes[0, 0].transAxes,
        )
        axes[0, 0].set_title(_intensity_title, fontsize=12, fontweight="bold")

    # (1,2) n_genes_by_counts
    if n_genes_col in adata.obs.columns:
        x = adata.obs[n_genes_col].values
        x = x[~np.isnan(x) & np.isfinite(x)]
        axes[0, 1].hist(x, bins=min(80, max(20, len(x) // 20)), color="coral", edgecolor="white")
        axes[0, 1].set_title(_feat_per_obs, fontsize=12, fontweight="bold")
        axes[0, 1].set_xlabel(n_genes_col)
    else:
        axes[0, 1].text(
            0.5,
            0.5,
            f"{n_genes_col} not in obs",
            ha="center",
            va="center",
            transform=axes[0, 1].transAxes,
        )
        axes[0, 1].set_title(_feat_per_obs, fontsize=12, fontweight="bold")

    # (2,1) log1p(total_counts)
    if total_counts_col in adata.obs.columns:
        x = np.log1p(adata.obs[total_counts_col].values.astype(float))
        x = x[~np.isnan(x) & np.isfinite(x)]
        axes[1, 0].hist(
            x, bins=min(80, max(20, len(x) // 20)), color="seagreen", edgecolor="white", alpha=0.8
        )
        axes[1, 0].set_title("log1p(total counts)", fontsize=12, fontweight="bold")
        axes[1, 0].set_xlabel("log1p(total_counts)")
    else:
        axes[1, 0].set_visible(False)

    # (2,2) pct_counts_mt or log1p(n_genes)
    if pct_mt_col in adata.obs.columns:
        x = adata.obs[pct_mt_col].values
        x = x[~np.isnan(x) & np.isfinite(x)]
        axes[1, 1].hist(
            x, bins=min(80, max(20, len(x) // 20)), color="purple", edgecolor="white", alpha=0.7
        )
        axes[1, 1].set_title("% mitochondrial", fontsize=12, fontweight="bold")
        axes[1, 1].set_xlabel(pct_mt_col)
    elif n_genes_col in adata.obs.columns:
        x = np.log1p(adata.obs[n_genes_col].values.astype(float))
        x = x[~np.isnan(x) & np.isfinite(x)]
        axes[1, 1].hist(
            x, bins=min(80, max(20, len(x) // 20)), color="gray", edgecolor="white", alpha=0.7
        )
        axes[1, 1].set_title(
            f"log1p({terms['features_lower']} per {terms['observation_lower']})",
            fontsize=12,
            fontweight="bold",
        )
        axes[1, 1].set_xlabel("log1p(n_genes_by_counts)")
    else:
        axes[1, 1].text(
            0.5,
            0.5,
            "No pct_mt or n_genes",
            ha="center",
            va="center",
            transform=axes[1, 1].transAxes,
        )
        axes[1, 1].set_title("(optional)", fontsize=12, fontweight="bold")

    plt.tight_layout()

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
        fig.savefig(output_dir / f"{basename}.png", bbox_inches="tight", dpi=dpi)

    return fig


def qc_2x4_pre_post(
    adata_pre: AnnData,
    adata_post: AnnData,
    *,
    total_counts_col: str = "total_counts",
    n_genes_col: str = "n_genes_by_counts",
    pct_mt_col: str = "pct_counts_mt",
    output_dir: str | Path | None = None,
    basename: str = "qc_2x4_pre_post",
    dpi: int = 300,
    figsize: tuple[float, float] = (16, 10),
    modality: str = "visium",
) -> plt.Figure:
    """
    Plot pre- vs post-filter QC: 2 rows x 4 columns.
    Left 2x2 = pre-filter (raw) metrics; right 2x2 = post-filter metrics.
    Use this so post-filter distributions (e.g. after filter_cells/filter_genes)
    are directly comparable to pre.

    Panels: row0 = total_counts (pre), n_genes (pre) | total_counts (post), n_genes (post);
    row1 = log1p(total_counts) (pre), pct_mt (pre) | log1p(total_counts) (post), pct_mt (post).

    Parameters
    ----------
    adata_pre : AnnData
        Pre-filter (raw) adata with total_counts, n_genes_by_counts (and optionally pct_counts_mt).
    adata_post : AnnData
        Post-filter (and optionally normalized) adata. Should have been filtered so that
        n_obs and metric distributions differ from pre. Must have same obs column names.
    total_counts_col : str
        Obs column for total counts (default 'total_counts').
    n_genes_col : str
        Obs column for number of genes (default 'n_genes_by_counts').
    pct_mt_col : str
        Obs column for percent mitochondrial (default 'pct_counts_mt').
    output_dir : str or Path or None
        If set, save PDF and PNG here (default None).
    basename : str
        Base name for files (default 'qc_2x4_pre_post').
    dpi : int
        DPI for PNG (default 300).
    figsize : tuple
        Figure size (default (16, 10)).

    Returns
    -------
    matplotlib.figure.Figure
        The figure.
    """
    from .report_utils import get_modality_terms

    terms = get_modality_terms(modality)
    _feat_per_obs = f"{terms['features']} per {terms['observation_lower']}"

    fig, axes = plt.subplots(2, 4, figsize=figsize)

    def _draw_hist(ax, x, color, title, xlabel):
        x = np.asarray(x, dtype=float)
        x = x[~np.isnan(x) & np.isfinite(x)]
        if len(x) == 0:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        else:
            ax.hist(x, bins=min(80, max(20, len(x) // 20)), color=color, edgecolor="white")
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel(xlabel, fontsize=9)

    # Pre: (0,0) total_counts, (0,1) n_genes, (1,0) log1p(total_counts), (1,1) pct_mt
    if total_counts_col in adata_pre.obs.columns:
        _draw_hist(
            axes[0, 0],
            adata_pre.obs[total_counts_col].values,
            "steelblue",
            "Pre: Total counts",
            total_counts_col,
        )
        _draw_hist(
            axes[1, 0],
            np.log1p(adata_pre.obs[total_counts_col].values.astype(float)),
            "seagreen",
            "Pre: log1p(total counts)",
            "log1p(total_counts)",
        )
    else:
        axes[0, 0].set_title("Pre: Total counts (missing)", fontsize=11)
        axes[1, 0].set_title("Pre: log1p (missing)", fontsize=11)
    if n_genes_col in adata_pre.obs.columns:
        _draw_hist(
            axes[0, 1],
            adata_pre.obs[n_genes_col].values,
            "coral",
            f"Pre: {_feat_per_obs}",
            n_genes_col,
        )
    else:
        axes[0, 1].set_title("Pre: n_genes (missing)", fontsize=11)
    if pct_mt_col in adata_pre.obs.columns:
        _draw_hist(
            axes[1, 1],
            adata_pre.obs[pct_mt_col].values,
            "purple",
            "Pre: % mitochondrial",
            pct_mt_col,
        )
    else:
        axes[1, 1].set_title("Pre: % mt (missing)", fontsize=11)

    # Post: (0,2) total_counts, (0,3) n_genes, (1,2) log1p(total_counts), (1,3) pct_mt
    if total_counts_col in adata_post.obs.columns:
        _draw_hist(
            axes[0, 2],
            adata_post.obs[total_counts_col].values,
            "steelblue",
            "Post: Total counts",
            total_counts_col,
        )
        _draw_hist(
            axes[1, 2],
            np.log1p(adata_post.obs[total_counts_col].values.astype(float)),
            "seagreen",
            "Post: log1p(total counts)",
            "log1p(total_counts)",
        )
    else:
        axes[0, 2].set_title("Post: Total counts (missing)", fontsize=11)
        axes[1, 2].set_title("Post: log1p (missing)", fontsize=11)
    if n_genes_col in adata_post.obs.columns:
        _draw_hist(
            axes[0, 3],
            adata_post.obs[n_genes_col].values,
            "coral",
            f"Post: {_feat_per_obs}",
            n_genes_col,
        )
    else:
        axes[0, 3].set_title("Post: n_genes (missing)", fontsize=11)
    if pct_mt_col in adata_post.obs.columns:
        _draw_hist(
            axes[1, 3],
            adata_post.obs[pct_mt_col].values,
            "purple",
            "Post: % mitochondrial",
            pct_mt_col,
        )
    else:
        axes[1, 3].set_title("Post: % mt (missing)", fontsize=11)

    fig.suptitle(
        "QC: Pre-filter (left) vs Post-filter (right)", fontsize=14, fontweight="bold", y=1.02
    )
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
        fig.savefig(output_dir / f"{basename}.png", bbox_inches="tight", dpi=dpi)

    return fig


def qc_spatial_multipage(
    adata: AnnData,
    library_id_col: str,
    output_path: str | Path,
    *,
    total_counts_col: str = "total_counts",
    pct_mt_col: str = "pct_counts_mt",
    figsize: tuple[float, float] = (18, 6),
    dpi: int = 300,
    common_scale: bool = True,
) -> None:
    """
    Multipage spatial QC report: one page per sample with 1x3 panels
    (total_count, log1p(total_count), % mt).

    When common_scale is True (default), the same vmin/vmax is used for each
    metric across all pages so color scales are comparable across samples.

    Requires adata.obs[library_id_col], adata.obsm['spatial'], and
    adata.uns['spatial'] with per-library images for spatial plots.

    Parameters
    ----------
    adata : AnnData
        Annotated data with spatial coords and (optionally) H&E in uns['spatial'].
    library_id_col : str
        Column in adata.obs identifying library/sample.
    output_path : str or Path
        Path to output PDF (e.g. figures/QC/raw/qc_spatial_multipage.pdf).
    total_counts_col : str
        Obs column for total counts (default 'total_counts').
    pct_mt_col : str
        Obs column for percent mitochondrial (default 'pct_counts_mt').
    figsize : tuple
        Figure size per page (default (18, 6)).
    dpi : int
        DPI for saved PDF (default 300).
    common_scale : bool
        If True, use global vmin/vmax (99th percentile) per metric across all
        spots so every page uses the same color scale (default True).
    """
    if total_counts_col not in adata.obs.columns:
        raise ValueError(
            f"adata.obs[{total_counts_col!r}] required. Run calculate_qc_metrics first."
        )

    tc = adata.obs[total_counts_col].values.astype(float)
    tc = tc[~np.isnan(tc) & np.isfinite(tc)]
    log1p_counts = pd.Series(
        np.log1p(adata.obs[total_counts_col].values.astype(float)),
        index=adata.obs_names,
    )

    if common_scale and len(tc) > 0:
        vmin_tc, vmax_tc = 0.0, float(np.nanpercentile(tc, 99))
        log1p_vals = np.log1p(tc)
        vmin_log, vmax_log = 0.0, float(np.nanpercentile(log1p_vals, 99))
        if pct_mt_col in adata.obs.columns:
            pmt = adata.obs[pct_mt_col].values.astype(float)
            pmt = pmt[~np.isnan(pmt) & np.isfinite(pmt)]
            vmin_mt, vmax_mt = (
                0.0,
                float(np.nanpercentile(pmt, 99)) if len(pmt) > 0 else (0.0, 100.0),
            )
        else:
            vmin_mt, vmax_mt = 0.0, 100.0
    else:
        vmin_tc = vmax_tc = vmin_log = vmax_log = vmin_mt = vmax_mt = None

    cmap = "viridis"
    panels = [
        {
            "type": "continuous",
            "title": "Total counts",
            "obs_col": total_counts_col,
            "cmap": cmap,
            "vmin": vmin_tc,
            "vmax": vmax_tc,
        },
        {
            "type": "continuous",
            "title": "log1p(total counts)",
            "values": log1p_counts,
            "cmap": cmap,
            "vmin": vmin_log,
            "vmax": vmax_log,
        },
    ]

    if pct_mt_col in adata.obs.columns:
        panels.append(
            {
                "type": "continuous",
                "title": "% mitochondrial",
                "obs_col": pct_mt_col,
                "cmap": cmap,
                "vmin": vmin_mt,
                "vmax": vmax_mt,
            }
        )
    else:
        panels.append(
            {
                "type": "continuous",
                "title": "N/A",
                "values": pd.Series(0.0, index=adata.obs_names),
                "obs_col": "",
                "cmap": "gray",
                "vmin": 0,
                "vmax": 1,
            }
        )

    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    from ..pl import spatial as st_spatial

    st_spatial.multipage_spatial_pdf(
        adata,
        library_id_col,
        panels,
        str(out_path),
        figsize=figsize,
        dpi=dpi,
    )


def qc_violin_metrics(
    adata: AnnData,
    *,
    keys: list[str] | None = None,
    groupby: str | None = None,
    output_dir: str | Path | None = None,
    basename: str = "qc_violin",
    dpi: int = 300,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """
    Multi-panel violin plot for QC metrics (n_genes_by_counts, total_counts, pct_counts_mt).

    Uses scanpy's violin with show=False and captures the figure for saving.
    Requires adata.obs columns from calculate_qc_metrics.

    Parameters
    ----------
    adata : AnnData
        Annotated data with obs containing total_counts, n_genes_by_counts, pct_counts_mt.
    keys : list of str or None
        Obs columns to plot (default: n_genes_by_counts, total_counts, pct_counts_mt).
    groupby : str or None
        Optional obs column to stratify violins (e.g. library_id, sample).
    output_dir : str or Path or None
        If set, save PDF and PNG here.
    basename : str
        Base name for files (default 'qc_violin').
    dpi : int
        DPI for PNG (default 300).
    figsize : tuple or None
        Figure size; if None, scanpy default is used.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import scanpy as sc

    if keys is None:
        keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    keys = [k for k in keys if k in adata.obs.columns]
    if not keys:
        fig, ax = plt.subplots(figsize=figsize or (6, 4))
        ax.text(
            0.5,
            0.5,
            "No QC metric columns in obs",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        if output_dir is not None:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_dir / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
            fig.savefig(output_dir / f"{basename}.png", bbox_inches="tight", dpi=dpi)
        return fig

    if figsize is not None:
        plt.figure(figsize=figsize)
    sc.pl.violin(
        adata,
        keys=keys,
        groupby=groupby,
        multi_panel=True,
        show=False,
    )
    fig = plt.gcf()
    plt.tight_layout()

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
        fig.savefig(output_dir / f"{basename}.png", bbox_inches="tight", dpi=dpi)

    return fig


def qc_scatter_counts_genes(
    adata: AnnData,
    *,
    x: str = "total_counts",
    y: str = "n_genes_by_counts",
    color: str = "pct_counts_mt",
    output_dir: str | Path | None = None,
    basename: str = "qc_scatter",
    dpi: int = 300,
    figsize: tuple[float, float] = (6, 5),
) -> plt.Figure:
    """
    Scatter plot: total_counts (x) vs n_genes_by_counts (y), colored by pct_counts_mt.

    Uses scanpy's scatter with show=False. Requires adata.obs from calculate_qc_metrics.

    Parameters
    ----------
    adata : AnnData
        Annotated data with obs columns for x, y, and color.
    x, y, color : str
        Obs column names (defaults: total_counts, n_genes_by_counts, pct_counts_mt).
    output_dir : str or Path or None
        If set, save PDF and PNG here.
    basename : str
        Base name for files (default 'qc_scatter').
    dpi : int
        DPI for PNG (default 300).
    figsize : tuple
        Figure size (default (6, 5)).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import scanpy as sc

    for col in (x, y, color):
        if col not in adata.obs.columns:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(
                0.5,
                0.5,
                f"Missing obs column: {col}",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            if output_dir is not None:
                output_dir_p = Path(output_dir)
                output_dir_p.mkdir(parents=True, exist_ok=True)
                fig.savefig(output_dir_p / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
                fig.savefig(output_dir_p / f"{basename}.png", bbox_inches="tight", dpi=dpi)
            return fig

    plt.figure(figsize=figsize)
    sc.pl.scatter(adata, x=x, y=y, color=color, show=False)
    fig = plt.gcf()
    plt.tight_layout()

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
        fig.savefig(output_dir / f"{basename}.png", bbox_inches="tight", dpi=dpi)

    return fig


def plot_highly_variable_genes(
    adata: AnnData,
    *,
    output_dir: str | Path | None = None,
    basename: str = "hvg",
    dpi: int = 300,
    figsize: tuple[float, float] = (6, 4),
) -> plt.Figure:
    """
    Plot mean vs dispersion (or normalized dispersion) with highly variable genes highlighted.

    Requires adata.var with 'highly_variable' and flavor-specific columns (means, dispersions
    or dispersions_norm). Uses sc.pl.highly_variable_genes(show=False).

    Parameters
    ----------
    adata : AnnData
        Annotated data after highly_variable_genes (e.g. seurat flavor).
    output_dir : str or Path or None
        If set, save PDF and PNG here.
    basename : str
        Base name for files (default 'hvg').
    dpi : int
        DPI for PNG (default 300).
    figsize : tuple
        Figure size (default (6, 4)).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import scanpy as sc

    if "highly_variable" not in adata.var.columns:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(
            0.5,
            0.5,
            "highly_variable not in adata.var",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        if output_dir is not None:
            output_dir_p = Path(output_dir)
            output_dir_p.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_dir_p / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
            fig.savefig(output_dir_p / f"{basename}.png", bbox_inches="tight", dpi=dpi)
        return fig

    plt.figure(figsize=figsize)
    sc.pl.highly_variable_genes(adata, show=False)
    fig = plt.gcf()
    plt.tight_layout()

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
        fig.savefig(output_dir / f"{basename}.png", bbox_inches="tight", dpi=dpi)

    return fig


def _plot_svg_single(
    ax: plt.Axes,
    spatial_i: np.ndarray,
    x_vals: np.ndarray,
    color_vals: np.ndarray | None,
    color_by: str,
    color_is_bool: bool,
    title: str = "",
) -> None:
    """Draw one SVG scatter on ax (mean/rank vs spatial_i)."""
    valid = ~np.isnan(spatial_i) & np.isfinite(spatial_i)
    if not np.any(valid):
        ax.text(0.5, 0.5, "No valid spatial_i", ha="center", va="center", transform=ax.transAxes)
        return
    if color_vals is not None and color_is_bool:
        scatter = ax.scatter(
            x_vals[valid],
            spatial_i[valid],
            c=color_vals[valid].astype(int),
            cmap=plt.cm.tab10,
            vmin=0,
            vmax=1,
            s=8,
            alpha=0.7,
        )
        cbar = plt.colorbar(scatter, ax=ax, ticks=[0.25, 0.75])
        cbar.ax.set_yticklabels(["False", "True"])
    elif color_vals is not None:
        scatter = ax.scatter(
            x_vals[valid], spatial_i[valid], c=color_vals[valid], cmap="viridis", s=8, alpha=0.7
        )
        plt.colorbar(scatter, ax=ax, label=color_by)
    else:
        ax.scatter(x_vals[valid], spatial_i[valid], c="gray", s=8, alpha=0.7)
    ax.set_xlabel("Mean expr." if title else "Mean expression" if x_vals is not None else "Rank")
    ax.set_ylabel("Moran's I")
    if title:
        ax.set_title(title)


def plot_spatially_variable_genes(
    adata: AnnData,
    *,
    x_axis: str = "mean",
    color_by: str = "spatially_variable",
    output_dir: str | Path | None = None,
    basename: str = "svg",
    dpi: int = 300,
    figsize: tuple[float, float] = (6, 5),
) -> plt.Figure:
    """
    Scatter: x = mean expression (or rank), y = Moran's I (spatial_i), colored by spatially_variable or pval.

    If adata.uns['spatial_variable_per_library'] exists (from spatially_variable_genes_per_library),
    one subplot per library is drawn. Otherwise requires adata.var with spatial_i.

    Parameters
    ----------
    adata : AnnData
        Annotated data after spatially_variable_genes or with uns['spatial_variable_per_library'].
    x_axis : str
        'mean' or 'rank': x-axis (default 'mean').
    color_by : str
        'spatially_variable' or 'spatial_pval' (default 'spatially_variable').
    output_dir : str or Path or None
        If set, save PDF and PNG here.
    basename : str
        Base name for files (default 'svg').
    dpi : int
        DPI for PNG (default 300).
    figsize : tuple
        Figure size per panel (default (6, 5)).

    Returns
    -------
    matplotlib.figure.Figure
    """
    per_lib = adata.uns.get("spatial_variable_per_library")
    if isinstance(per_lib, dict) and len(per_lib) > 0:
        nlib = len(per_lib)
        ncol = min(3, nlib)
        nrow = (nlib + ncol - 1) // ncol
        fig, axes = plt.subplots(nrow, ncol, figsize=(figsize[0] * ncol, figsize[1] * nrow))
        axes = np.atleast_1d(axes).flat
        means_global = (
            adata.var["means"].values.astype(float) if "means" in adata.var.columns else None
        )
        if means_global is None and hasattr(adata.X, "toarray"):
            means_global = np.asarray(adata.X.mean(axis=0)).ravel()
        elif means_global is None:
            means_global = np.asarray(adata.X.mean(axis=0)).ravel()
        for idx, (lib_id, df) in enumerate(per_lib.items()):
            ax = axes[idx]
            si = (
                df["spatial_i"].reindex(adata.var_names).values.astype(float)
                if "spatial_i" in df.columns
                else None
            )
            if si is None or not np.any(np.isfinite(si)):
                ax.text(
                    0.5, 0.5, f"{lib_id}: no data", ha="center", va="center", transform=ax.transAxes
                )
                continue
            x_vals = (
                means_global
                if x_axis == "mean"
                else np.argsort(
                    np.nan_to_num(means_global, nan=np.nanmin(means_global) - 1)
                ).astype(float)
            )
            c_vals = None
            c_bool = False
            if color_by in df.columns:
                c_vals = df[color_by].reindex(adata.var_names).values
                c_bool = np.issubdtype(c_vals.dtype, np.bool_)
            _plot_svg_single(ax, si, x_vals, c_vals, color_by, c_bool, title=str(lib_id))
        for j in range(len(per_lib), len(axes)):
            axes[j].set_visible(False)
        plt.tight_layout()
        if output_dir is not None:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_dir / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
            fig.savefig(output_dir / f"{basename}.png", bbox_inches="tight", dpi=dpi)
        return fig

    if "spatial_i" not in adata.var.columns:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(
            0.5, 0.5, "spatial_i not in adata.var", ha="center", va="center", transform=ax.transAxes
        )
        if output_dir is not None:
            output_dir_p = Path(output_dir)
            output_dir_p.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_dir_p / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
            fig.savefig(output_dir_p / f"{basename}.png", bbox_inches="tight", dpi=dpi)
        return fig

    x = adata.var["spatial_i"].values.astype(float)
    valid = ~np.isnan(x) & np.isfinite(x)
    if not np.any(valid):
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No valid spatial_i", ha="center", va="center", transform=ax.transAxes)
        if output_dir is not None:
            output_dir_p = Path(output_dir)
            output_dir_p.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_dir_p / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
            fig.savefig(output_dir_p / f"{basename}.png", bbox_inches="tight", dpi=dpi)
        return fig

    if "means" in adata.var.columns:
        x_vals = adata.var["means"].values.astype(float)
    else:
        if hasattr(adata.X, "toarray"):
            x_vals = np.asarray(adata.X.mean(axis=0)).ravel()
        else:
            x_vals = np.asarray(adata.X.mean(axis=0)).ravel()
        if len(x_vals) != adata.n_vars:
            x_vals = np.full(adata.n_vars, np.nan)

    if x_axis == "rank":
        order = np.argsort(np.nan_to_num(x_vals, nan=np.nanmin(x_vals) - 1))
        x_vals = np.empty_like(order, dtype=float)
        x_vals[order] = np.arange(len(order))

    fig, ax = plt.subplots(figsize=figsize)
    if color_by in adata.var.columns:
        c_vals = adata.var[color_by].values[valid]
        if np.issubdtype(adata.var[color_by].dtype, np.bool_):
            scatter = ax.scatter(
                x_vals[valid],
                x[valid],
                c=c_vals.astype(int),
                cmap=plt.cm.tab10,
                vmin=0,
                vmax=1,
                s=8,
                alpha=0.7,
            )
            cbar = plt.colorbar(scatter, ax=ax, ticks=[0.25, 0.75])
            cbar.ax.set_yticklabels(["False", "True"])
        else:
            scatter = ax.scatter(
                x_vals[valid],
                x[valid],
                c=c_vals,
                cmap="viridis",
                s=8,
                alpha=0.7,
            )
            plt.colorbar(scatter, ax=ax, label=color_by)
    else:
        ax.scatter(x_vals[valid], x[valid], c="gray", s=8, alpha=0.7)
    ax.set_xlabel("Mean expression" if x_axis == "mean" else "Rank (mean)")
    ax.set_ylabel("Moran's I (spatial_i)")
    ax.set_title("Spatially variable genes")
    plt.tight_layout()

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
        fig.savefig(output_dir / f"{basename}.png", bbox_inches="tight", dpi=dpi)

    return fig


# ---------------------------------------------------------------------------
# Cross-sample comparison plots
# ---------------------------------------------------------------------------


def _save_fig(fig: plt.Figure, output_dir: str | Path | None, basename: str, dpi: int) -> None:
    """Save figure as PDF + PNG if output_dir is provided."""
    if output_dir is not None:
        od = Path(output_dir)
        od.mkdir(parents=True, exist_ok=True)
        fig.savefig(od / f"{basename}.pdf", bbox_inches="tight", dpi=dpi)
        fig.savefig(od / f"{basename}.png", bbox_inches="tight", dpi=dpi)


def _format_log10_ticks(ax, max_val: float) -> None:
    """Set custom log10(x+1) y-tick labels: 0 (1), 1 (10), 2 (100), etc."""
    tick_map = {0: "0 (1)", 1: "1 (10)", 2: "2 (100)", 3: "3 (1K)", 4: "4 (10K)", 5: "5 (100K)"}
    max_tick = int(np.ceil(max_val)) if np.isfinite(max_val) else 5
    ticks = list(range(min(max_tick + 1, 6)))
    ax.set_yticks(ticks)
    ax.set_yticklabels([tick_map.get(t, str(t)) for t in ticks])


def qc_sample_comparison_bar(
    metrics: pd.DataFrame,
    metric_cols: list[str] | None = None,
    classified: pd.DataFrame | None = None,
    output_dir: str | Path | None = None,
    basename: str = "qc_sample_comparison",
    dpi: int = 300,
    log_scale: bool = False,
) -> plt.Figure:
    """
    Bar chart per metric, one bar per sample, sorted by value.

    Failed samples (from ``classified``) are highlighted in red.

    Parameters
    ----------
    metrics : pd.DataFrame
        Output of ``compute_sample_metrics`` (indexed by sample).
    metric_cols : list of str or None
        Columns to plot (default: n_genes_median, total_counts_median,
        pct_mt_median, n_spots).
    classified : pd.DataFrame or None
        If provided (output of ``classify_samples``), failed samples shown in red.
    output_dir : str or Path or None
        If set, save PDF and PNG.
    basename : str
        Base filename.
    dpi : int
        DPI for PNG.
    log_scale : bool
        If True, transform values with log10(x + 1) and annotate y-ticks
        with original-scale labels (default False).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if metric_cols is None:
        metric_cols = [
            c
            for c in ["n_genes_median", "total_counts_median", "pct_mt_median", "n_spots"]
            if c in metrics.columns
        ]
    if not metric_cols:
        fig, ax = plt.subplots()
        ax.text(
            0.5,
            0.5,
            "No metric columns available",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        _save_fig(fig, output_dir, basename, dpi)
        return fig

    n_metrics = len(metric_cols)
    fig, axes = plt.subplots(n_metrics, 1, figsize=(max(8, len(metrics) * 0.5), 4 * n_metrics))
    if n_metrics == 1:
        axes = [axes]

    fail_set = set()
    if classified is not None and "qc_pass" in classified.columns:
        fail_set = set(classified.index[~classified["qc_pass"]])

    for ax, col in zip(axes, metric_cols, strict=False):
        sorted_df = metrics[[col]].dropna().sort_values(col)
        colors = ["#d62728" if s in fail_set else "#1f77b4" for s in sorted_df.index]
        values = sorted_df[col].values.astype(float)
        if log_scale:
            values = np.log10(values + 1)
        ax.bar(range(len(sorted_df)), values, color=colors)
        ax.set_xticks(range(len(sorted_df)))
        labels = [str(s) for s in sorted_df.index]
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(col)
        suffix = " (log10 scale)" if log_scale else ""
        ax.set_title(f"{col}{suffix}", fontweight="bold")
        if log_scale and len(values) > 0:
            _format_log10_ticks(ax, float(np.nanmax(values)))

    title_suffix = " (log10 scale)" if log_scale else ""
    fig.suptitle(
        f"Cross-sample QC comparison{title_suffix}", fontsize=14, fontweight="bold", y=1.01
    )
    plt.tight_layout()
    _save_fig(fig, output_dir, basename, dpi)
    return fig


def qc_sample_violin_grouped(
    adata: AnnData,
    sample_col: str = "library_id",
    keys: list[str] | None = None,
    classified: pd.DataFrame | None = None,
    output_dir: str | Path | None = None,
    basename: str = "qc_sample_violin",
    dpi: int = 300,
    log_scale: bool = False,
) -> plt.Figure:
    """
    Violin plots grouped by sample for direct distribution comparison.

    Parameters
    ----------
    adata : AnnData
        Annotated data with QC columns in obs.
    sample_col : str
        Column in obs identifying samples.
    keys : list of str or None
        Obs columns to plot (default: n_genes_by_counts, total_counts, pct_counts_mt).
    classified : pd.DataFrame or None
        If provided, failed sample names are marked with ``(FAIL)`` suffix.
    output_dir : str or Path or None
        If set, save PDF and PNG.
    basename : str
        Base filename.
    dpi : int
        DPI for PNG.
    log_scale : bool
        If True, apply log10(x + 1) to all keys except pct_counts_mt
        (which stays linear). Custom y-tick annotations are added (default False).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if sample_col not in adata.obs.columns:
        fig, ax = plt.subplots()
        ax.text(
            0.5, 0.5, f"{sample_col} not in obs", ha="center", va="center", transform=ax.transAxes
        )
        _save_fig(fig, output_dir, basename, dpi)
        return fig

    if keys is None:
        keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    keys = [k for k in keys if k in adata.obs.columns]
    if not keys:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No QC columns", ha="center", va="center", transform=ax.transAxes)
        _save_fig(fig, output_dir, basename, dpi)
        return fig

    # Keys that should NOT be log-transformed (percentage metrics)
    _pct_keys = {"pct_counts_mt", "pct_counts_hb"}

    fail_set = set()
    if classified is not None and "qc_pass" in classified.columns:
        fail_set = set(classified.index[~classified["qc_pass"]])

    samples = sorted(adata.obs[sample_col].dropna().unique())
    n_keys = len(keys)
    fig, axes = plt.subplots(n_keys, 1, figsize=(max(8, len(samples) * 0.7), 4 * n_keys))
    if n_keys == 1:
        axes = [axes]

    for ax, key in zip(axes, keys, strict=False):
        do_log = log_scale and key not in _pct_keys
        data_per_sample = []
        labels = []
        for s in samples:
            vals = adata.obs.loc[adata.obs[sample_col] == s, key].dropna().values.astype(float)
            if do_log:
                vals = np.log10(vals + 1)
            data_per_sample.append(vals)
            label = f"{s} (FAIL)" if s in fail_set else str(s)
            labels.append(label)

        non_empty = [(i, d) for i, d in enumerate(data_per_sample) if len(d) > 0]
        if non_empty:
            parts = ax.violinplot(
                [d for _, d in non_empty],
                positions=[i for i, _ in non_empty],
                showmeans=True,
                showmedians=True,
            )
            for pc in parts.get("bodies", []):
                pc.set_alpha(0.7)

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(key)
        suffix = " (log10 scale)" if do_log else ""
        ax.set_title(f"{key}{suffix}", fontweight="bold")
        if do_log and non_empty:
            all_vals = np.concatenate([d for _, d in non_empty])
            if len(all_vals) > 0:
                _format_log10_ticks(ax, float(np.nanmax(all_vals)))

    title_suffix = " (log10 scale)" if log_scale else ""
    fig.suptitle(
        f"Per-sample QC distributions{title_suffix}", fontsize=14, fontweight="bold", y=1.01
    )
    plt.tight_layout()
    _save_fig(fig, output_dir, basename, dpi)
    return fig


def qc_sample_scatter_matrix(
    metrics: pd.DataFrame,
    metric_cols: list[str] | None = None,
    classified: pd.DataFrame | None = None,
    output_dir: str | Path | None = None,
    basename: str = "qc_sample_scatter_matrix",
    dpi: int = 300,
) -> plt.Figure:
    """
    Pairwise scatter of sample-level metrics with pass/fail coloring.

    Parameters
    ----------
    metrics : pd.DataFrame
        Output of ``compute_sample_metrics``.
    metric_cols : list of str or None
        Columns for scatter matrix (default: n_spots, n_genes_median,
        total_counts_median, pct_mt_median).
    classified : pd.DataFrame or None
        If provided, color points by pass (blue) / fail (red).
    output_dir : str or Path or None
        If set, save PDF and PNG.
    basename : str
        Base filename.
    dpi : int
        DPI for PNG.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if metric_cols is None:
        metric_cols = [
            c
            for c in ["n_spots", "n_genes_median", "total_counts_median", "pct_mt_median"]
            if c in metrics.columns
        ]
    metric_cols = [c for c in metric_cols if c in metrics.columns]

    if len(metric_cols) < 2:
        fig, ax = plt.subplots()
        ax.text(
            0.5, 0.5, "Need >= 2 metric columns", ha="center", va="center", transform=ax.transAxes
        )
        _save_fig(fig, output_dir, basename, dpi)
        return fig

    n = len(metric_cols)
    fig, axes = plt.subplots(n, n, figsize=(4 * n, 4 * n))

    fail_set = set()
    if classified is not None and "qc_pass" in classified.columns:
        fail_set = set(classified.index[~classified["qc_pass"]])

    colors = ["#d62728" if s in fail_set else "#1f77b4" for s in metrics.index]

    for i in range(n):
        for j in range(n):
            ax = axes[i, j] if n > 1 else axes
            if i == j:
                vals = metrics[metric_cols[i]].dropna().values
                ax.hist(
                    vals, bins=max(5, len(vals) // 3), color="#1f77b4", edgecolor="white", alpha=0.7
                )
                ax.set_xlabel(metric_cols[i])
            else:
                ax.scatter(
                    metrics[metric_cols[j]].values,
                    metrics[metric_cols[i]].values,
                    c=colors,
                    s=40,
                    alpha=0.8,
                    edgecolors="white",
                    linewidths=0.5,
                )
                if j == 0:
                    ax.set_ylabel(metric_cols[i])
                if i == n - 1:
                    ax.set_xlabel(metric_cols[j])

    fig.suptitle("Sample QC scatter matrix", fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    _save_fig(fig, output_dir, basename, dpi)
    return fig


def qc_pct_mt_per_sample(
    adata: AnnData,
    sample_col: str = "library_id",
    pct_mt_col: str = "pct_counts_mt",
    classified: pd.DataFrame | None = None,
    output_dir: str | Path | None = None,
    basename: str = "qc_pct_mt_per_sample",
    dpi: int = 300,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """
    Per-sample %MT distribution as box plots colored by pass/fail.

    Parameters
    ----------
    adata : AnnData
        Annotated data with ``pct_counts_mt`` in obs.
    sample_col : str
        Column in obs identifying samples.
    pct_mt_col : str
        Obs column for percent mitochondrial (default ``pct_counts_mt``).
    classified : pd.DataFrame or None
        If provided, failed samples are colored red, pass samples blue.
    output_dir : str or Path or None
        If set, save PDF and PNG.
    basename : str
        Base filename.
    dpi : int
        DPI for PNG.
    figsize : tuple or None
        Figure size (auto-scaled by sample count if None).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if pct_mt_col not in adata.obs.columns:
        fig, ax = plt.subplots(figsize=figsize or (8, 5))
        ax.text(
            0.5,
            0.5,
            f"{pct_mt_col} not in obs",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        _save_fig(fig, output_dir, basename, dpi)
        return fig

    if sample_col not in adata.obs.columns:
        fig, ax = plt.subplots(figsize=figsize or (8, 5))
        ax.text(
            0.5,
            0.5,
            f"{sample_col} not in obs",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        _save_fig(fig, output_dir, basename, dpi)
        return fig

    fail_set = set()
    if classified is not None and "qc_pass" in classified.columns:
        fail_set = set(classified.index[~classified["qc_pass"]])

    samples = sorted(adata.obs[sample_col].dropna().unique())
    if figsize is None:
        figsize = (max(8, len(samples) * 0.7), 5)
    fig, ax = plt.subplots(figsize=figsize)

    data_per_sample = []
    for s in samples:
        vals = adata.obs.loc[adata.obs[sample_col] == s, pct_mt_col].dropna().values
        data_per_sample.append(vals)

    bp = ax.boxplot(
        data_per_sample,
        positions=range(len(samples)),
        patch_artist=True,
        widths=0.6,
    )
    for patch, s in zip(bp["boxes"], samples, strict=False):
        color = "#d62728" if s in fail_set else "#1f77b4"
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    labels = [f"{s} (FAIL)" if s in fail_set else str(s) for s in samples]
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(pct_mt_col)
    ax.set_title("% Mitochondrial per sample", fontweight="bold")
    plt.tight_layout()
    _save_fig(fig, output_dir, basename, dpi)
    return fig
