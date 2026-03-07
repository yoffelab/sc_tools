"""
Heatmap and clustermap utilities.

Generic helpers for signature score heatmaps with annotation bars,
hierarchical sorting, and clustering within groups.

Categorical colors follow the scanpy convention: adata.uns[f'{obs_col}_colors']
is a list of hex strings (one per category in order). If missing, we create
and store it so all plotting stays consistent.
"""

from __future__ import annotations

from typing import Any

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from sc_tools.utils.signatures import get_signature_df

__all__ = [
    "hex_to_rgb",
    "get_obs_category_colors",
    "cluster_within_groups",
    "annotation_colors_from_categories",
    "signature_score_heatmap",
]


# Defaults
DEFAULT_CLUSTER_METHOD = "average"
DEFAULT_CLUSTER_METRIC = "euclidean"
DEFAULT_HEATMAP_FIGSIZE = (16, 12)
DEFAULT_CLUSTERMAP_FIGSIZE = (18, 14)
DEFAULT_DENDROGRAM_RATIO = 0.1
DEFAULT_VMIN, DEFAULT_VMAX = -3, 3


def hex_to_rgb(hex_color: str) -> tuple[float, float, float]:
    """
    Convert hex color to RGB tuple in (0, 1) range.

    Parameters
    ----------
    hex_color : str
        Hex string (e.g. '#66c2a5' or '66c2a5').

    Returns
    -------
    tuple
        (r, g, b) in [0, 1].
    """
    hex_color = hex_color.lstrip("#")
    if len(hex_color) != 6:
        return (0.8, 0.8, 0.8)
    return tuple(int(hex_color[i : i + 2], 16) / 255.0 for i in (0, 2, 4))


def get_obs_category_colors(
    adata,
    obs_col: str,
    store_if_missing: bool = True,
) -> dict[Any, tuple[float, float, float]] | None:
    """
    Get category -> RGB color mapping for a categorical obs column using the
    scanpy convention: adata.uns[f'{obs_col}_colors'] is a list of hex strings
    (one per category in order). If missing or length mismatch, create a default
    palette and optionally store it in adata.uns.

    Parameters
    ----------
    adata : AnnData
        Object with obs[obs_col] categorical and optionally uns[f'{obs_col}_colors'].
    obs_col : str
        Name of the categorical column in adata.obs.
    store_if_missing : bool
        If True (default), when colors are missing or invalid, create a palette
        and set adata.uns[f'{obs_col}_colors'] to a list of hex strings.

    Returns
    -------
    dict or None
        Map from category value to (r, g, b) in [0, 1]. None if obs_col is not
        present or not categorical.
    """
    if obs_col not in adata.obs.columns:
        return None
    ser = adata.obs[obs_col]
    if not isinstance(ser.dtype, pd.CategoricalDtype):
        return None
    categories = ser.cat.categories.tolist()
    n = len(categories)
    uns_key = f"{obs_col}_colors"
    existing = adata.uns.get(uns_key)
    if existing is not None and len(existing) == n:
        return {cat: hex_to_rgb(h) for cat, h in zip(categories, existing, strict=True)}
    # Create default palette (hex list, category order)
    if n <= 12:
        palette = sns.color_palette("Set3", n)
    else:
        palette = sns.color_palette("husl", n)
    hex_list = [mcolors.to_hex(c) for c in palette]
    if store_if_missing:
        adata.uns[uns_key] = hex_list
    return {cat: hex_to_rgb(h) for cat, h in zip(categories, hex_list, strict=True)}


def cluster_within_groups(
    data_matrix: np.ndarray,
    group_labels: np.ndarray,
    method: str = DEFAULT_CLUSTER_METHOD,
    metric: str = DEFAULT_CLUSTER_METRIC,
) -> np.ndarray:
    """
    Cluster rows within each group; preserve group order.

    Parameters
    ----------
    data_matrix : np.ndarray
        Data matrix (n_samples x n_features).
    group_labels : np.ndarray
        Group label per row (same length as data_matrix).
    method : str
        Linkage method (default 'average').
    metric : str
        Distance metric (default 'euclidean').

    Returns
    -------
    np.ndarray
        Reordered row indices (cluster within each group, groups in order).
    """
    unique_groups = np.unique(group_labels)
    reordered_indices = []

    for group in unique_groups:
        group_mask = group_labels == group
        group_indices = np.where(group_mask)[0]
        group_data = data_matrix[group_indices, :]

        if len(group_indices) > 1:
            distances = pdist(group_data, metric=metric)
            if len(distances) > 0:
                linkage_matrix = linkage(distances, method=method)
                leaves = leaves_list(linkage_matrix)
                group_clustered = group_indices[leaves]
            else:
                group_clustered = group_indices
        else:
            group_clustered = group_indices

        reordered_indices.extend(group_clustered.tolist())

    return np.array(reordered_indices)


def annotation_colors_from_categories(
    annotations: pd.DataFrame,
    column_colors: dict[str, dict[str, tuple[float, float, float]]] | None = None,
    default_hex: dict[str, str] | None = None,
) -> dict[str, list[tuple[float, float, float]]]:
    """
    Build per-column color lists (RGB) for annotation bars.

    Parameters
    ----------
    annotations : pd.DataFrame
        Index = sample ids, columns = annotation names; values = category labels.
    column_colors : dict, optional
        Maps column name -> {category: (r,g,b)}. If None, uses seaborn Set3 for each column.
    default_hex : dict, optional
        Maps column name -> {category: hex}. Converted to RGB; overridden by column_colors.

    Returns
    -------
    dict
        column name -> list of (r,g,b) in row order.
    """
    default_gray = (0.8, 0.8, 0.8)
    out = {}

    for col in annotations.columns:
        unique_vals = annotations[col].unique()
        if column_colors and col in column_colors:
            pal = column_colors[col]
        elif default_hex and col in default_hex:
            pal_hex = default_hex[col]
            pal = {k: hex_to_rgb(v) for k, v in pal_hex.items()}
        else:
            pal_list = sns.color_palette("Set3", len(unique_vals))
            pal = dict(zip(unique_vals, pal_list, strict=True))
        colors = [pal.get(v, default_gray) for v in annotations[col]]
        out[col] = colors

    return out


def signature_score_heatmap(
    adata,
    sig_columns: list[str],
    annotation_cols: dict[str, str],
    sort_by: list[str],
    category_orders: dict[str, list] | None = None,
    cluster: bool = False,
    sig_prefix: str = "sig:",
    sig_suffix: str = "_z",
    vmin: float = DEFAULT_VMIN,
    vmax: float = DEFAULT_VMAX,
    figsize: tuple[float, float] | None = None,
    solidity_colors_hex: dict[str, str] | None = None,
    legend_title: str | None = None,
) -> tuple[plt.Figure, Any | None]:
    """
    Build heatmap or clustermap of signature scores with annotation bars.

    Annotation columns are given as display_name -> obs column name.
    sort_by is list of display names (primary, secondary). category_orders
    maps display name -> ordered list of categories; columns not in it are
    ordered by mean score (descending).

    Parameters
    ----------
    adata : AnnData
        Object with signature scores in obs[sig_columns] and annotation columns.
    sig_columns : list
        Obs column names for signature scores.
    annotation_cols : dict
        Display name -> obs column name, e.g. {'Patient': 'library_id', 'Solidity': 'tumor_type'}.
    sort_by : list
        [primary, secondary] display names for sorting.
    category_orders : dict, optional
        Display name -> list of category order. Missing names: order by mean score.
    cluster : bool
        If True, build clustermap with within-group clustering; else heatmap only.
    sig_prefix, sig_suffix : str
        Stripped from sig_columns for row labels.
    vmin, vmax : float
        Color scale for score matrix.
    figsize : tuple, optional
        (width, height). Default heatmap (16,12), clustermap (18,14).
    solidity_colors_hex : dict, optional
        For backward compatibility: category -> hex for second annotation (e.g. Solidity).
    legend_title : str, optional
        Title for legend (e.g. 'Solidity').

    Returns
    -------
    fig : Figure
        The figure (caller can save with st.pl.save_figure).
    g : seaborn.ClusterGrid or None
        If cluster=True, the ClusterGrid for further tweaks; else None.
    """
    category_orders = category_orders or {}
    figsize = figsize or (DEFAULT_CLUSTERMAP_FIGSIZE if cluster else DEFAULT_HEATMAP_FIGSIZE)

    # Build annotations DataFrame with display names
    annotations = pd.DataFrame(
        {disp: adata.obs[obs_col].values for disp, obs_col in annotation_cols.items()}
    )
    annotations.index = adata.obs_names

    # Signature matrix and row labels (from obsm or obs)
    sig_df = get_signature_df(adata)
    cols = [c for c in sig_columns if c in sig_df.columns]
    sig_df = sig_df[cols]
    [c.replace(sig_prefix, "").replace(sig_suffix, "") for c in sig_df.columns]

    # Order for each sort dimension
    orders = {}
    for disp in sort_by:
        col = annotation_cols.get(disp, disp)
        if disp in category_orders:
            orders[disp] = category_orders[disp]
        else:
            # Order by mean score (descending)
            means = {}
            for cat in annotations[disp].unique():
                mask = annotations[disp] == cat
                means[cat] = sig_df.loc[mask].mean().mean()
            orders[disp] = sorted(means.keys(), key=lambda x: means[x], reverse=True)

    # Group key for sorting
    def key_row(row):
        key = []
        for d in sort_by:
            val = row[d]
            if d in orders and val in orders[d]:
                key.append(orders[d].index(val))
            else:
                key.append(999)
        return tuple(key)

    annotations = annotations.copy()
    annotations["_key"] = annotations.apply(key_row, axis=1)
    annotations = annotations.sort_values("_key")
    annotations = annotations.drop(columns=["_key"])
    sorted_idx = annotations.index
    sig_df = sig_df.loc[sorted_idx]

    # Optional within-group clustering
    if cluster:
        group_cols = sort_by
        grouped = annotations.groupby(group_cols, observed=True)
        clustered_idx = []
        for _, grp in grouped:
            idx = grp.index.tolist()
            if len(idx) > 1:
                data_grp = sig_df.loc[idx].values
                pos = cluster_within_groups(
                    data_grp,
                    np.zeros(len(idx), dtype=int),
                    method=DEFAULT_CLUSTER_METHOD,
                    metric=DEFAULT_CLUSTER_METRIC,
                )
                clustered_idx.extend([idx[i] for i in pos])
            else:
                clustered_idx.extend(idx)
        sig_df = sig_df.loc[clustered_idx]
        annotations = annotations.loc[clustered_idx]

    # Colors for annotation bars (RGB). Prefer adata.uns[f'{obs_col}_colors']
    # (scanpy convention) for categorical columns; then solidity_colors_hex for
    # "Solidity" if not from uns; else Set3 in annotation_colors_from_categories.
    column_colors = {}
    for disp, obs_col in annotation_cols.items():
        if disp not in annotations.columns:
            continue
        obs_colors = get_obs_category_colors(adata, obs_col, store_if_missing=True)
        if obs_colors is not None:
            column_colors[disp] = obs_colors
    default_hex_dict = None
    if (
        solidity_colors_hex
        and "Solidity" in annotations.columns
        and "Solidity" not in column_colors
    ):
        default_hex_dict = {"Solidity": solidity_colors_hex}
    color_lists = annotation_colors_from_categories(
        annotations, column_colors=column_colors, default_hex=default_hex_dict
    )
    arrays = {
        col: np.array([color_lists[col][i] for i in range(len(annotations))]).reshape(1, -1, 3)
        for col in annotations.columns
    }

    if cluster:
        g = sns.clustermap(
            sig_df.T,
            cmap="RdBu_r",
            center=0,
            vmin=vmin,
            vmax=vmax,
            figsize=figsize,
            linewidths=0,
            cbar_kws={"label": "Z-scored Signature Score"},
            dendrogram_ratio=DEFAULT_DENDROGRAM_RATIO,
            cbar_pos=(0.02, 0.85, 0.03, 0.12),
            row_cluster=True,
            col_cluster=False,
            yticklabels=True,
            xticklabels=False,
            method=DEFAULT_CLUSTER_METHOD,
            metric=DEFAULT_CLUSTER_METRIC,
        )
        fig = g.fig
        heatmap_pos = g.ax_heatmap.get_position()
        ann_h = 0.015
        for i, col in enumerate(annotations.columns):
            ax_ann = fig.add_axes(
                [
                    heatmap_pos.x0,
                    heatmap_pos.y1 + (i * ann_h),
                    heatmap_pos.width,
                    ann_h,
                ]
            )
            ax_ann.imshow(arrays[col], aspect="auto", interpolation="nearest", rasterized=True)
            ax_ann.set_xticks([])
            ax_ann.set_yticks([])
            ax_ann.set_ylabel(col, rotation=0, ha="right", va="center", fontsize=9)
            for spine in ax_ann.spines.values():
                spine.set_visible(False)
        g.ax_heatmap.set_xlabel(f"Spots (n={len(annotations)})", fontsize=11)
        g.ax_heatmap.set_ylabel("Gene Signatures", fontsize=11)
        if legend_title and sort_by and solidity_colors_hex:
            second = sort_by[1] if len(sort_by) > 1 else None
            if second:
                legend_elems = [
                    Patch(facecolor=hex_to_rgb(v), label=k) for k, v in solidity_colors_hex.items()
                ]
                leg_ax = fig.add_axes([0.95, 0.95, 0.05, 0.05])
                leg_ax.axis("off")
                leg_ax.legend(handles=legend_elems, title=legend_title or second, fontsize=9)
        return fig, g
    else:
        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(
            sig_df.T,
            cmap="RdBu_r",
            center=0,
            vmin=vmin,
            vmax=vmax,
            cbar_kws={"label": "Z-scored Signature Score"},
            linewidths=0,
            ax=ax,
            xticklabels=False,
            yticklabels=True,
        )
        y1 = ax.get_position().y1
        w = ax.get_position().width
        x0 = ax.get_position().x0
        h_ann = 0.02
        for i, col in enumerate(annotations.columns):
            ax_ann = fig.add_axes([x0, y1 + (i * h_ann), w, h_ann])
            ax_ann.imshow(arrays[col], aspect="auto", interpolation="nearest")
            ax_ann.set_xticks([])
            ax_ann.set_yticks([])
            ax_ann.set_ylabel(col, rotation=0, ha="right", va="center", fontsize=10)
        ax.set_title(
            f"Signature Scores (sorted by {sort_by[0]} then {sort_by[1]})",
            fontsize=14,
            fontweight="bold",
            pad=20,
        )
        ax.set_xlabel(f"Spots (n={len(annotations)})", fontsize=12)
        ax.set_ylabel("Gene Signatures", fontsize=12)
        if legend_title and solidity_colors_hex:
            legend_elems = [
                Patch(facecolor=hex_to_rgb(v), label=k) for k, v in solidity_colors_hex.items()
            ]
            leg_ax = fig.add_axes([0.95, 0.95, 0.05, 0.05])
            leg_ax.axis("off")
            leg_ax.legend(handles=legend_elems, title=legend_title, fontsize=9)
        return fig, None
