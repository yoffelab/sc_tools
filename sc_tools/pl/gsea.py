"""
GSEA / enrichment result plotting utilities.

Functions
---------
plot_gsea_dotplot     Dot plot of enrichment results (size = -log10 p_adj, color = NES).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_gsea_dotplot(
    result_df: pd.DataFrame,
    top_n: int = 10,
    groups: list[str] | None = None,
    figsize: tuple = (8, 6),
    output_path: str | Path | None = None,
) -> plt.Figure:
    """
    Dot plot of gene set enrichment results.

    Dot size encodes statistical significance (-log10 adjusted p-value);
    dot color encodes normalized enrichment score (NES) or odds ratio.
    Compatible with output from both :func:`run_ora` and
    :func:`run_gsea_pseudobulk`.

    Parameters
    ----------
    result_df : pd.DataFrame
        Long-format enrichment results. Must contain at minimum the columns
        ``group``, ``gene_set``, ``p_adj``. Optional: ``NES`` (used for color
        if present, otherwise ``odds_ratio`` is used, then defaults to black).
    top_n : int
        Top N gene sets to display per group (ranked by p_adj ascending).
        Default 10.
    groups : list[str] or None
        Subset of groups to display. Defaults to all groups in result_df.
    figsize : tuple
        Figure size in inches (width, height). Default (8, 6).
    output_path : str, Path, or None
        If provided, save the figure to this path at 300 DPI. Parent
        directories are created if needed.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object.
    """
    required_cols = {"group", "gene_set", "p_adj"}
    missing = required_cols - set(result_df.columns)
    if missing:
        raise ValueError(f"result_df is missing required columns: {missing}")

    df = result_df.copy()

    if groups is not None:
        df = df[df["group"].isin(groups)]
    if df.empty:
        raise ValueError("No rows to plot after filtering groups.")

    # Select top_n per group by p_adj
    top_rows = df.sort_values("p_adj").groupby("group", sort=False).head(top_n)

    # Determine color column
    color_col = None
    if "NES" in top_rows.columns:
        color_col = "NES"
    elif "odds_ratio" in top_rows.columns:
        color_col = "odds_ratio"

    # Pivot: rows = gene_set, cols = group
    all_groups = sorted(top_rows["group"].unique().tolist())
    all_sets = (
        top_rows.groupby("gene_set")["p_adj"]
        .min()
        .sort_values()
        .index.tolist()[: top_n * len(all_groups)]
    )
    # Deduplicate and keep order
    seen: set = set()
    unique_sets = []
    for s in all_sets:
        if s not in seen:
            unique_sets.append(s)
            seen.add(s)

    # Subset to the sets that appear in top_rows
    unique_sets = [s for s in unique_sets if s in top_rows["gene_set"].values]

    # Build size (dot) and color arrays
    size_mat = pd.DataFrame(np.nan, index=unique_sets, columns=all_groups)
    color_mat = pd.DataFrame(np.nan, index=unique_sets, columns=all_groups)

    for _, row in top_rows.iterrows():
        if row["gene_set"] not in size_mat.index:
            continue
        p = row["p_adj"]
        size_val = -np.log10(max(p, 1e-300)) if not np.isnan(p) else 0.0
        size_mat.loc[row["gene_set"], row["group"]] = size_val
        if color_col is not None:
            color_mat.loc[row["gene_set"], row["group"]] = row[color_col]

    fig, ax = plt.subplots(figsize=figsize)

    # Dot plot via scatter
    xticklabels = all_groups
    yticklabels = unique_sets

    # Determine color range
    color_vals = color_mat.values.ravel()
    color_vals_valid = color_vals[~np.isnan(color_vals)]
    if len(color_vals_valid) > 0:
        vmin, vmax = color_vals_valid.min(), color_vals_valid.max()
        if vmin == vmax:
            vmin -= 1
            vmax += 1
    else:
        vmin, vmax = -1, 1

    cmap = "RdBu_r"
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])

    # Scatter
    for xi, grp in enumerate(all_groups):
        for yi, gs in enumerate(unique_sets):
            sz = size_mat.loc[gs, grp]
            cv = color_mat.loc[gs, grp] if color_col is not None else np.nan
            if np.isnan(sz):
                continue
            dot_size = max(sz * 30, 10)
            if np.isnan(cv):
                c = "grey"
            else:
                c = sm.to_rgba(cv)
            ax.scatter(xi, yi, s=dot_size, c=[c], zorder=3)

    ax.set_xticks(range(len(all_groups)))
    ax.set_xticklabels(xticklabels, rotation=45, ha="right")
    ax.set_yticks(range(len(unique_sets)))
    ax.set_yticklabels(yticklabels)
    ax.set_xlim(-0.5, len(all_groups) - 0.5)
    ax.set_ylim(-0.5, len(unique_sets) - 0.5)
    ax.grid(True, linestyle="--", alpha=0.4)

    # Color bar
    if color_col is not None:
        cbar = fig.colorbar(sm, ax=ax, shrink=0.6)
        cbar.set_label(color_col)

    # Size legend (manual)
    legend_sizes = [1, 2, 5]
    legend_handles = [
        plt.scatter([], [], s=s * 30, c="grey", label=f"-log10(p_adj) = {s}") for s in legend_sizes
    ]
    ax.legend(
        handles=legend_handles,
        title="Significance",
        loc="upper left",
        bbox_to_anchor=(1.15, 1.0),
        frameon=False,
    )

    ax.set_xlabel("Group")
    ax.set_ylabel("Gene Set")
    ax.set_title("Enrichment Dot Plot")

    fig.tight_layout()

    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")

    return fig
