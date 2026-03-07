"""
Statistical annotation utilities for plots.

Provides functions for adding statistical annotations to plots:
- Significance bars
- Asterisks (*, **, ***, ****)
- P-value annotations
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def get_asterisk(pval: float) -> str:
    """
    Get asterisk string based on p-value.

    Parameters
    ----------
    pval : float
        P-value

    Returns
    -------
    str
        Asterisk string: '****' (< 0.0001), '***' (< 0.001), '**' (< 0.01), '*' (< 0.05), or 'ns'
    """
    if pval < 0.0001:
        return "****"
    elif pval < 0.001:
        return "***"
    elif pval < 0.01:
        return "**"
    elif pval < 0.05:
        return "*"
    else:
        return "ns"


def add_significance_bars(
    ax: plt.Axes,
    comparisons: list[dict],
    y_positions: list[float] | None = None,
    bar_height: float | None = None,
    bar_spacing: float | None = None,
) -> None:
    """
    Add significance bars to a plot.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes
    comparisons : list of dict
        List of comparison dictionaries with keys:
        - 'group1_idx': int, index of first group
        - 'group2_idx': int, index of second group
        - 'pval': float, p-value
        - 'adj_pval': float, adjusted p-value (optional)
        - 'label': str, label for comparison (optional)
    y_positions : list of float, optional
        Y positions for bars. If None, auto-calculated.
    bar_height : float, optional
        Height of bars. If None, auto-calculated.
    bar_spacing : float, optional
        Spacing between bars. If None, auto-calculated.
    """
    if len(comparisons) == 0:
        return

    # Get current y-axis limits
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min

    # Auto-calculate if not provided
    if bar_height is None:
        bar_height = y_range * 0.04
    if bar_spacing is None:
        bar_spacing = y_range * 0.08
    if y_positions is None:
        # Calculate maximum whisker extent
        max_whisker = y_max
        for line in ax.lines:
            if line.get_ydata().size > 0:
                max_whisker = max(max_whisker, max(line.get_ydata()))

        # Start bars above maximum whisker
        whisker_buffer = y_range * 0.10
        y_start = max_whisker + whisker_buffer

        # Calculate y positions for each bar
        y_positions = [y_start + (i * bar_spacing) for i in range(len(comparisons))]

    # Filter to only significant comparisons
    sig_comparisons = [c for c in comparisons if c.get("adj_pval", c.get("pval", 1.0)) < 0.05]

    # Draw bars
    for idx, comp in enumerate(sig_comparisons):
        if idx >= len(y_positions):
            break

        y_bar = y_positions[idx]
        x_start = comp["group1_idx"] + 1  # +1 because boxplot indices are 1-based
        x_end = comp["group2_idx"] + 1

        # Draw horizontal bar
        ax.plot([x_start, x_end], [y_bar, y_bar], "k-", linewidth=1.5)

        # Draw vertical ticks
        tick_height = bar_height * 0.4
        ax.plot([x_start, x_start], [y_bar - tick_height, y_bar], "k-", linewidth=1.5)
        ax.plot([x_end, x_end], [y_bar - tick_height, y_bar], "k-", linewidth=1.5)

        # Add asterisk
        pval = comp.get("adj_pval", comp.get("pval", 1.0))
        asterisk = get_asterisk(pval)
        x_asterisk = (x_start + x_end) / 2
        ax.text(
            x_asterisk,
            y_bar + bar_height * 0.4,
            asterisk,
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
        )

    # Adjust y-axis to accommodate bars
    if len(sig_comparisons) > 0:
        max_bar_y = max(y_positions[: len(sig_comparisons)]) + bar_height * 0.6
        ax.set_ylim(y_min, max_bar_y)


def plot_boxplot_with_stats(
    adata,
    value_col: str,
    group_col: str,
    stats_df: pd.DataFrame | None = None,
    ax: plt.Axes | None = None,
    colors: list[str] | None = None,
    figsize: tuple[int, int] = (6, 4),
    show_comparisons: bool = True,
) -> tuple[plt.Figure, plt.Axes]:
    """
    Create boxplot with significance annotations.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    value_col : str
        Column name for values
    group_col : str
        Column name for groups
    stats_df : DataFrame, optional
        Statistics dataframe with columns: 'group', 'adj_pval', 'program'
    ax : Axes, optional
        Matplotlib axes. If None, creates new figure.
    colors : list of str, optional
        Colors for groups (default: uses seaborn palette)
    figsize : tuple
        Figure size (default: (6, 4))
    show_comparisons : bool
        Whether to show comparison text (default: True)

    Returns
    -------
    tuple
        (figure, axes)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Prepare data
    groups = (
        adata.obs[group_col].cat.categories
        if hasattr(adata.obs[group_col], "cat")
        else adata.obs[group_col].unique()
    )
    plot_data = []

    for group in groups:
        values = adata.obs.loc[adata.obs[group_col] == group, value_col].dropna()
        plot_data.append(values.values)

    # Create boxplot
    bp = ax.boxplot(plot_data, labels=groups, patch_artist=True, showfliers=False, widths=0.6)

    # Color boxes
    if colors is None:
        colors = sns.color_palette("Set2", len(groups))
    for patch, color in zip(bp["boxes"], colors[: len(groups)], strict=True):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Add significance bars if stats provided
    if stats_df is not None and show_comparisons:
        # Prepare comparisons for 1-vs-all
        comparisons = []
        comparison_texts = []

        for i, group in enumerate(groups):
            if "program" in stats_df.columns:
                program_stats = stats_df[stats_df["program"] == value_col]
            else:
                program_stats = stats_df

            row = program_stats[program_stats["group"] == group]
            if len(row) > 0 and pd.notna(row.iloc[0].get("adj_pval", np.nan)):
                adj_pval = row.iloc[0]["adj_pval"]
                asterisk = get_asterisk(adj_pval)

                # For 1-vs-all, compare group to all others
                comparisons.append(
                    {
                        "group1_idx": 0,  # Start of all groups
                        "group2_idx": len(groups) - 1,  # End of all groups
                        "pval": adj_pval,
                        "adj_pval": adj_pval,
                        "asterisk": asterisk,
                        "x_pos": i + 1,
                    }
                )

                rest_label = "Rest"
                comparison_texts.append(f"{group} - {rest_label}: {asterisk} (p={adj_pval:.2e})")

        # Add significance bars
        if len(comparisons) > 0:
            add_significance_bars(ax, comparisons)

        # Add comparison text at bottom
        if len(comparison_texts) > 0:
            comparison_str = "\n".join(comparison_texts)
            n_lines = len(comparison_texts)
            bottom_margin = max(0.15, 0.05 + (n_lines * 0.025))

            fig.text(
                0.5,
                0.01,
                comparison_str,
                ha="center",
                va="bottom",
                fontsize=9,
                bbox={
                    "boxstyle": "round,pad=0.5",
                    "facecolor": "lightyellow",
                    "alpha": 0.8,
                    "edgecolor": "black",
                    "linewidth": 1,
                },
                family="monospace",
            )

            plt.tight_layout(rect=[0, bottom_margin, 1, 1])
    else:
        plt.tight_layout()

    # Labels
    ax.set_ylabel(value_col.replace("sig:", "").replace("_z", ""), fontsize=11)
    ax.set_xlabel(group_col, fontsize=11)
    ax.set_title(value_col.replace("sig:", "").replace("_z", ""), fontsize=12, fontweight="bold")
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    sns.despine(ax=ax)

    return fig, ax


__all__ = ["get_asterisk", "add_significance_bars", "plot_boxplot_with_stats"]
