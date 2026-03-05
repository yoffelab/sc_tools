"""
General-purpose QC plot functions for post-integration analysis.

- qc_umap_grid: grid of UMAPs colored by different obs columns
- qc_cluster_distribution: stacked bar chart of cluster composition per sample
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

__all__ = [
    "qc_umap_grid",
    "qc_cluster_distribution",
]


def qc_umap_grid(
    adata: AnnData,
    color_keys: list[str] | None = None,
    umap_key: str = "X_umap",
    n_cols: int = 3,
    *,
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

    n_panels = max(len(color_keys), 1)
    n_rows = int(np.ceil(n_panels / n_cols))

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(figsize_per_panel[0] * n_cols, figsize_per_panel[1] * n_rows),
        squeeze=False,
    )

    coords = np.asarray(adata.obsm[umap_key])

    for idx, key in enumerate(color_keys):
        ax = axes[idx // n_cols, idx % n_cols]
        values = adata.obs[key]

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
    for idx in range(len(color_keys), n_rows * n_cols):
        axes[idx // n_cols, idx % n_cols].set_visible(False)

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
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=7,
        title_fontsize=8,
    )
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    fig.tight_layout()

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        fig.savefig(out / f"{basename}.png", dpi=dpi, bbox_inches="tight")

    return fig
