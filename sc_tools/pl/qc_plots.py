"""
General-purpose QC plot functions for post-integration analysis.

- qc_umap_grid: grid of UMAPs colored by different obs columns
- qc_cluster_distribution: stacked bar chart of cluster composition per sample
- qc_embedding_umap_grid: per-embedding UMAP comparison grid
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

__all__ = [
    "qc_umap_grid",
    "qc_cluster_distribution",
    "qc_embedding_umap_grid",
    "qc_celltype_abundance",
]

logger = logging.getLogger(__name__)


def qc_umap_grid(
    adata: AnnData,
    color_keys: list[str] | None = None,
    umap_key: str = "X_umap",
    n_cols: int = 3,
    *,
    max_cells: int = 50_000,
    output_dir: str | Path | None = None,
    basename: str = "qc_umap_grid",
    dpi: int = 300,
    figsize_per_panel: tuple[float, float] = (5, 4),
) -> plt.Figure:
    """Grid of UMAP embeddings colored by different obs columns.

    Parameters
    ----------
    adata
        AnnData with ``umap_key`` in ``obsm``.
    color_keys
        List of ``obs`` column names to color by. If ``None``, auto-detects
        from common keys: ``sample``, ``library_id``, ``batch``,
        ``raw_data_dir``, ``leiden``, ``celltype``, ``celltype_broad``.
    umap_key
        Key in ``adata.obsm`` for UMAP coordinates.
    n_cols
        Number of columns in the grid.
    max_cells
        Maximum number of cells to plot. When ``adata.n_obs`` exceeds this,
        cells are randomly subsampled. The plot order is always shuffled so
        minority categories are not hidden behind dominant ones.
    output_dir
        If provided, save figure to ``{output_dir}/{basename}.png``.
    basename
        Filename stem for saved figure.
    dpi
        Resolution for saved figure.
    figsize_per_panel
        Width and height per subplot panel.

    Returns
    -------
    matplotlib.Figure
    """
    if umap_key not in adata.obsm:
        raise KeyError(f"UMAP key {umap_key!r} not found in adata.obsm")

    if color_keys is None:
        candidates = [
            "sample",
            "library_id",
            "batch",
            "raw_data_dir",
            "leiden",
            "celltype",
            "celltype_broad",
        ]
        color_keys = [k for k in candidates if k in adata.obs.columns]
    if not color_keys:
        color_keys = ["leiden"] if "leiden" in adata.obs.columns else []

    # Subsample + shuffle for balanced visualization
    rng = np.random.default_rng(42)
    n_obs = adata.n_obs
    if n_obs > max_cells:
        idx = rng.choice(n_obs, size=max_cells, replace=False)
    else:
        idx = np.arange(n_obs)
    # Shuffle so later-drawn points don't systematically hide earlier ones
    rng.shuffle(idx)

    n_panels = max(len(color_keys), 1)
    n_rows = int(np.ceil(n_panels / n_cols))

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(figsize_per_panel[0] * n_cols, figsize_per_panel[1] * n_rows),
        squeeze=False,
    )

    coords = np.asarray(adata.obsm[umap_key])[idx]

    for panel_idx, key in enumerate(color_keys):
        ax = axes[panel_idx // n_cols, panel_idx % n_cols]
        values = adata.obs[key].iloc[idx]

        if values.dtype.name == "category" or values.dtype == object:
            categories = values.astype("category").cat.categories
            codes = values.astype("category").cat.codes
            cmap = plt.cm.tab20 if len(categories) <= 20 else plt.cm.gist_ncar
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=codes,
                cmap=cmap,
                s=1,
                alpha=0.6,
                rasterized=True,
            )
        else:
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=values,
                cmap="viridis",
                s=1,
                alpha=0.6,
                rasterized=True,
            )
            fig.colorbar(scatter, ax=ax, shrink=0.6)

        ax.set_title(key, fontsize=10)
        ax.set_xlabel("UMAP1", fontsize=8)
        ax.set_ylabel("UMAP2", fontsize=8)
        ax.tick_params(labelsize=7)

    # Hide unused axes
    for panel_idx in range(len(color_keys), n_rows * n_cols):
        axes[panel_idx // n_cols, panel_idx % n_cols].set_visible(False)

    fig.tight_layout()

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        fig.savefig(out / f"{basename}.png", dpi=dpi, bbox_inches="tight")

    return fig


def qc_cluster_distribution(
    adata: AnnData,
    cluster_key: str = "leiden",
    sample_col: str = "library_id",
    *,
    normalize: bool = True,
    output_dir: str | Path | None = None,
    basename: str = "qc_cluster_distribution",
    dpi: int = 300,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Stacked bar chart of cluster composition per sample.

    Parameters
    ----------
    adata
        AnnData with cluster and sample annotations.
    cluster_key
        Column in ``obs`` with cluster labels.
    sample_col
        Column in ``obs`` with sample identifiers.
    normalize
        If ``True``, show proportions (0-1). Otherwise raw counts.
    output_dir
        If provided, save figure.
    basename
        Filename stem for saved figure.
    dpi
        Resolution for saved figure.
    figsize
        Figure size. Auto-scales if ``None``.

    Returns
    -------
    matplotlib.Figure
    """
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"Cluster key {cluster_key!r} not found in adata.obs")
    if sample_col not in adata.obs.columns:
        raise KeyError(f"Sample column {sample_col!r} not found in adata.obs")

    ct = pd.crosstab(adata.obs[sample_col], adata.obs[cluster_key])
    if normalize:
        ct = ct.div(ct.sum(axis=1), axis=0)

    n_samples = len(ct)
    if figsize is None:
        figsize = (max(8, n_samples * 0.6), 5)

    fig, ax = plt.subplots(figsize=figsize)
    ct.plot.bar(stacked=True, ax=ax, cmap="tab20", width=0.85, legend=True)
    ax.set_ylabel("Proportion" if normalize else "Count")
    ax.set_xlabel(sample_col)
    ax.set_title(f"{cluster_key} distribution per {sample_col}")
    ax.legend(
        title=cluster_key,
        loc="upper right",
        fontsize=6,
        title_fontsize=7,
        framealpha=0.8,
        ncol=1 if ct.shape[1] <= 15 else 2,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="center", fontsize=8)
    fig.tight_layout()

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        fig.savefig(out / f"{basename}.png", dpi=dpi, bbox_inches="tight")

    return fig


def qc_embedding_umap_grid(
    adata: AnnData,
    embedding_keys: dict[str, str],
    color_key: str = "sample",
    *,
    max_cells: int = 30_000,
    n_cols: int = 3,
    figsize_per_panel: tuple[float, float] = (5, 4),
) -> plt.Figure:
    """Compute and plot a UMAP for each integration embedding.

    For each embedding in *embedding_keys*, computes neighbors + UMAP on a
    subsampled copy, then plots a panel colored by *color_key*. This shows
    how each integration method structures the data.

    Parameters
    ----------
    adata
        AnnData with embeddings in ``obsm``.
    embedding_keys
        Dict mapping method name to ``obsm`` key.
    color_key
        Column in ``adata.obs`` to color each panel by.
    max_cells
        Subsample to this many cells before computing UMAP (for speed).
    n_cols
        Number of columns in the grid.
    figsize_per_panel
        Width and height per subplot panel.

    Returns
    -------
    matplotlib.Figure
    """
    import scanpy as sc

    rng = np.random.default_rng(42)
    n_obs = adata.n_obs

    # Subsample + shuffle
    if n_obs > max_cells:
        idx = rng.choice(n_obs, size=max_cells, replace=False)
        adata_sub = adata[idx].copy()
    else:
        adata_sub = adata.copy()

    # Shuffle for balanced plotting
    shuffle_idx = rng.permutation(adata_sub.n_obs)

    methods = list(embedding_keys.items())
    n_panels = len(methods)
    n_rows = int(np.ceil(n_panels / n_cols))

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(figsize_per_panel[0] * n_cols, figsize_per_panel[1] * n_rows),
        squeeze=False,
    )

    if color_key in adata_sub.obs.columns:
        values = adata_sub.obs[color_key].iloc[shuffle_idx]
        is_cat = values.dtype.name == "category" or values.dtype == object
    else:
        values = None
        is_cat = False

    for panel_idx, (method_name, obsm_key) in enumerate(methods):
        ax = axes[panel_idx // n_cols, panel_idx % n_cols]

        if obsm_key not in adata_sub.obsm:
            ax.text(
                0.5,
                0.5,
                f"{obsm_key}\nnot found",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=9,
                color="#999",
            )
            ax.set_title(method_name, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        # Compute UMAP for this embedding
        try:
            sc.pp.neighbors(adata_sub, use_rep=obsm_key, key_added=obsm_key)
            sc.tl.umap(adata_sub, neighbors_key=obsm_key)
            coords = np.asarray(adata_sub.obsm["X_umap"])[shuffle_idx]
        except Exception:
            logger.debug("UMAP failed for %s", method_name, exc_info=True)
            ax.text(
                0.5,
                0.5,
                "UMAP failed",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=9,
                color="#c0392b",
            )
            ax.set_title(method_name, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        if values is not None and is_cat:
            categories = values.astype("category").cat.categories
            codes = values.astype("category").cat.codes
            cmap = plt.cm.tab20 if len(categories) <= 20 else plt.cm.gist_ncar
            ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=codes,
                cmap=cmap,
                s=0.5,
                alpha=0.5,
                rasterized=True,
            )
        elif values is not None:
            ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=values,
                cmap="viridis",
                s=0.5,
                alpha=0.5,
                rasterized=True,
            )
        else:
            ax.scatter(
                coords[:, 0],
                coords[:, 1],
                s=0.5,
                alpha=0.5,
                color="#333",
                rasterized=True,
            )

        ax.set_title(method_name, fontsize=9)
        ax.set_xlabel("UMAP1", fontsize=7)
        ax.set_ylabel("UMAP2", fontsize=7)
        ax.tick_params(labelsize=6)

    # Hide unused axes
    for panel_idx in range(n_panels, n_rows * n_cols):
        axes[panel_idx // n_cols, panel_idx % n_cols].set_visible(False)

    fig.suptitle(f"Per-embedding UMAP (colored by {color_key})", fontsize=11, y=1.01)
    fig.tight_layout()
    return fig


def qc_celltype_abundance(
    adata: AnnData,
    celltype_key: str = "celltype",
    sample_col: str = "library_id",
    *,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Stacked bar chart of celltype proportions per sample.

    Colors are read from ``adata.uns[f'{celltype_key}_colors']`` when
    available; otherwise the Okabe-Ito colorblind-safe palette is used.

    Parameters
    ----------
    adata
        AnnData with cell type and sample annotations.
    celltype_key
        Column in ``adata.obs`` with cell type labels.
    sample_col
        Column in ``adata.obs`` with sample identifiers.
    figsize
        Figure size. Auto-scales when ``None``.

    Returns
    -------
    matplotlib.Figure
    """
    # Okabe-Ito palette (colorblind-safe; skills.md §11)
    _OKABE_ITO = [
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#000000",
        "#999999",
    ]

    ct_counts = (
        adata.obs.groupby([sample_col, celltype_key], observed=True).size().unstack(fill_value=0)
    )
    ct_prop = ct_counts.div(ct_counts.sum(axis=1), axis=0)

    n_types = ct_prop.shape[1]
    if figsize is None:
        figsize = (max(8, ct_prop.shape[0] * 0.8), max(4, n_types * 0.3))

    colors_key = f"{celltype_key}_colors"
    if colors_key in adata.uns and len(adata.uns[colors_key]) >= n_types:
        colors = list(adata.uns[colors_key])[:n_types]
    else:
        colors = (_OKABE_ITO * ((n_types // len(_OKABE_ITO)) + 1))[:n_types]

    fig, ax = plt.subplots(figsize=figsize)
    ct_prop.plot(kind="bar", stacked=True, ax=ax, color=colors, legend=True)
    ax.set_xlabel(sample_col)
    ax.set_ylabel("Proportion")
    ax.set_title(f"Celltype abundance per {sample_col}")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7)
    plt.tight_layout()
    return fig
