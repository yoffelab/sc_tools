# -*- coding: utf-8 -*-
import warnings
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

# -----------------------------
# Config supplied by you
# -----------------------------
ADATA_PATH = Path("results/scvi.leiden.phenotyped.h5ad")
cluster_key: str = "cluster"
library_id_key: str = "library_id"
library_ids: Dict[str, List[str]] = {
    "P1_NAT": ["PT01-1_NAT", "PT01-3_NAT", "PT01-5_NAT"],
    "P1_TUM": ["PT01-2_TUM", "PT01-4_TUM", "PT01-6_TUM"],
    "P3_NAT": ["P3_S1_BL_NT_1_iF10_S7"],
    "P3_TUM": ["P3_S2_BL_T_2_iG10_S8"],
    "P5_NAT": ["Pat5_Samp1_BL_2_iD10_S5", "Pat_5_BL_Tu_Samp2_A1_iD9_S4"],
    "P1_LN": ["Pat_1_LN_PRT_Samp8_D1_iC9_S3"],
    "P3_LN": ["Pat_3_LN_PRT_Samp_8_D1_iB9_S2", "Pat_3_LN_Samp_7_D1_iA9_S1"],
}

# Layout params
FIGSIZE = (18, 12)
GROUP_SPACER = 0.6      # x gap between groups
BAR_WIDTH = 0.8
LABEL_FONTSIZE = 8

# Matrixplot colormap for any ancillary plots you might add
MATRIX_CMAP = "RdBu_r"

# -----------------------------
# Helper functions
# -----------------------------
def _validate_inputs(adata, cluster_key: str, library_id_key: str) -> None:
    assert cluster_key in adata.obs, f"obs lacks {cluster_key}"
    assert library_id_key in adata.obs, f"obs lacks {library_id_key}"
    assert pd.api.types.is_categorical_dtype(adata.obs[cluster_key]), f"obs['{cluster_key}'] must be categorical"
    assert f"{cluster_key}_colors" in adata.uns, f"uns lacks {cluster_key}_colors"
    n_cats = len(adata.obs[cluster_key].cat.categories)
    n_cols = len(adata.uns[f"{cluster_key}_colors"])
    assert n_cats == n_cols, f"length mismatch between {cluster_key} categories ({n_cats}) and colors ({n_cols})"

def _build_ordered_samples(library_ids: Dict[str, List[str]], present_samples: Sequence[str]) -> Tuple[List[str], List[int], List[str]]:
    """
    Returns ordered_samples, spacer_positions, tick_group_labels
    - ordered_samples is the x order
    - spacer_positions are x indices after which a spacer should be inserted
    - tick_group_labels one per group positioned in the middle of each group
    """
    ordered = []
    spacer_after = []
    group_centers = []
    x_cursor = 0
    for g, samples in library_ids.items():
        group_samples = [s for s in samples if s in present_samples]
        if len(group_samples) == 0:
            continue
        ordered.extend(group_samples)
        # center position for group label
        start = x_cursor
        end = x_cursor + len(group_samples) - 1
        group_centers.append((g, (start + end) / 2.0))
        x_cursor += len(group_samples)
        spacer_after.append(x_cursor - 1)  # after last sample of this group
    return ordered, spacer_after, group_centers

def _compartment(label: str) -> str:
    """Map cluster label to epithelial vs stromal vs immune."""
    l = label.lower()
    # epithelial programs
    if any(k in l for k in [
        "enterocyte", "goblet", "enteroendocrine", "crypt stem", "regenerative",
        "neuroendocrine", "metabolic / stress", "weak epithelial", "epithelium"
    ]):
        return "epithelial"
    # stromal programs
    if any(k in l for k in ["fibroblast", "smooth muscle", "pericyte", "ecm mix", "matrix remodeling", "ecm high"]):
        return "stromal"
    # immune programs
    if any(k in l for k in ["b cell", "plasma", "gc", "myeloid", "mast", "nk", "t cell", "dendritic"]):
        return "immune"
    return "unknown"

def _stack_df(adata, ordered_samples: List[str], cluster_key: str) -> pd.DataFrame:
    """Counts per sample x cluster, with all categories present and zero filled."""
    cats = list(adata.obs[cluster_key].cat.categories)
    df = (
        adata.obs[[cluster_key, library_id_key]]
        .groupby([library_id_key, cluster_key])
        .size()
        .rename("count")
        .reset_index()
        .pivot(index=library_id_key, columns=cluster_key, values="count")
        .reindex(index=ordered_samples, columns=cats)
        .fillna(0)
        .astype(int)
    )
    return df

def _subset_by_compartment(df_counts: pd.DataFrame, labels: Sequence[str], wanted: str) -> pd.DataFrame:
    """Filter columns by compartment and keep zero columns if none present."""
    comp_cols = [lab for lab in labels if _compartment(lab) == wanted]
    if len(comp_cols) == 0:
        # If no columns classified, return an all zeros frame with same index
        return pd.DataFrame(0, index=df_counts.index, columns=[])
    return df_counts[comp_cols].copy()

def _rowwise_props(df_counts: pd.DataFrame) -> pd.DataFrame:
    denom = df_counts.sum(axis=1).replace(0, np.nan)
    props = df_counts.div(denom, axis=0).fillna(0.0)
    return props

def _draw_stacked(ax, props: pd.DataFrame, colors: Dict[str, str], show_percent: bool, title: str, bar_width: float) -> None:
    bottoms = np.zeros(props.shape[0], dtype=float)
    x = np.arange(props.shape[0], dtype=float)
    for col in props.columns:
        vals = props[col].values
        c = colors.get(col, "#cccccc")
        ax.bar(x, vals, width=bar_width, bottom=bottoms, color=c, edgecolor="none", linewidth=0)
        bottoms += vals
    # annotate percentages at top of bars if requested
    if show_percent:
        for i, total in enumerate(bottoms):
            ax.text(i, 1.005, f"{100.0:.2f}%", ha="center", va="bottom", fontsize=LABEL_FONTSIZE)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Proportion")
    ax.set_title(title)

def _draw_counts(ax, df_counts: pd.DataFrame, colors: Dict[str, str], bar_width: float) -> None:
    bottoms = np.zeros(df_counts.shape[0], dtype=float)
    x = np.arange(df_counts.shape[0], dtype=float)
    for col in df_counts.columns:
        vals = df_counts[col].values.astype(float)
        c = colors.get(col, "#cccccc")
        ax.bar(x, vals, width=bar_width, bottom=bottoms, color=c, edgecolor="none", linewidth=0)
        bottoms += vals
    # annotate counts at the top of each bar
    for i, total in enumerate(bottoms):
        ax.text(i, total * 1.005 if total > 0 else 0, f"{int(total)}", ha="center", va="bottom", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Cells")
    ax.set_title("Cell counts (stacked)")

def _insert_group_spacers(ax, n_bars: int, spacer_after: List[int], spacer: float) -> None:
    # Shift ticks so bars remain at integer positions; just draw vertical lines as visual spacers
    for pos in spacer_after[:-1]:  # no spacer after last group
        ax.axvline(x=pos + 0.5, color="white", linewidth=spacer * 20, alpha=1.0, zorder=0)

def _apply_xticks(ax, ordered_samples: List[str], rotation: int = 90) -> None:
    ax.set_xticks(np.arange(len(ordered_samples)))
    ax.set_xticklabels(ordered_samples, rotation=rotation, ha="right", fontsize=LABEL_FONTSIZE)

def _legend_from_colors(ax, labels: Sequence[str], color_map: Dict[str, str], ncols: int = 4) -> None:
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=color_map.get(lbl, "#cccccc"), ec="none")
        for lbl in labels
    ]
    ax.legend(handles, labels, ncols=ncols, bbox_to_anchor=(1.0, 1.02), loc="upper left", frameon=False, fontsize=9, title="Cell types")

# -----------------------------
# Main plotting routine
# -----------------------------
def plot_five_row_stacked(
    adata,
    cluster_key: str,
    library_id_key: str,
    library_ids: Dict[str, List[str]],
    figsize: Tuple[int, int] = (11, 8),     # smaller, compact
    group_spacer: float = 0.6,
    bar_width: float = 0.8,
) -> None:
    _validate_inputs(adata, cluster_key, library_id_key)

    # Canonical order and colors
    labels = list(adata.obs[cluster_key].cat.categories)
    palette_list = list(adata.uns[f"{cluster_key}_colors"])
    color_map = {lab: col for lab, col in zip(labels, palette_list)}

    # Sample order from grouping
    present = adata.obs[library_id_key].astype(str).unique().tolist()
    ordered_samples, spacer_after, group_centers = _build_ordered_samples(library_ids, present)
    if len(ordered_samples) == 0:
        raise ValueError("No matching library_id values from the grouping dictionary were found in the AnnData object")

    # Counts
    df_counts_all = _stack_df(adata, ordered_samples, cluster_key)

    # Split by compartment
    df_epi = _subset_by_compartment(df_counts_all, labels, "epithelial")
    df_str = _subset_by_compartment(df_counts_all, labels, "stromal")
    df_imm = _subset_by_compartment(df_counts_all, labels, "immune")

    # Proportions
    props_all = _rowwise_props(df_counts_all)
    props_epi = _rowwise_props(df_epi) if df_epi.shape[1] else df_epi
    props_str = _rowwise_props(df_str) if df_str.shape[1] else df_str
    props_imm = _rowwise_props(df_imm) if df_imm.shape[1] else df_imm

    # ------ figure: 2 columns -> [plots | legend] ------
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=figsize, constrained_layout=False)
    gs = GridSpec(
        nrows=5, ncols=2, figure=fig,
        height_ratios=[1.1, 1, 1, 1, 1],   # counts a bit taller
        width_ratios=[1.0, 0.22],          # narrow legend panel
        wspace=0.25, hspace=0.2
    )

    # Axes for the five stacked plots
    # inside plot_five_row_stacked, after creating axes
    ax_counts = fig.add_subplot(gs[0, 0])
    ax_all    = fig.add_subplot(gs[1, 0], sharex=ax_counts)
    ax_epi    = fig.add_subplot(gs[2, 0], sharex=ax_counts)
    ax_str    = fig.add_subplot(gs[3, 0], sharex=ax_counts)
    ax_imm    = fig.add_subplot(gs[4, 0], sharex=ax_counts)

    # turn off xticks for all but bottom
    for ax in [ax_counts, ax_all, ax_epi, ax_str]:
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(axis="x", which="both", length=0)  # hide ticks too

    # keep bottom xticks
    _apply_xticks(ax_imm, ordered_samples, rotation=90)

    # Legend axis (single vertical column)
    ax_leg = fig.add_subplot(gs[:, 1])
    ax_leg.axis("off")
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_map.get(lbl, "#cccccc"), ec="none") for lbl in labels]
    leg = ax_leg.legend(
        handles, labels, ncols=1, title="Cell types",
        loc="upper left", frameon=False, fontsize=9, title_fontsize=10,
        borderaxespad=0.0, handlelength=1.2, handletextpad=0.4
    )

    # Draw stacks
    _draw_counts(ax_counts, df_counts_all, color_map, bar_width)
    if props_all.shape[1] > 0:
        _draw_stacked(ax_all, props_all, color_map, show_percent=True, title="All cells proportion", bar_width=bar_width)
    else:
        ax_all.set_title("All cells proportion (no data)")
    if props_epi.shape[1] > 0:
        _draw_stacked(ax_epi, props_epi, color_map, show_percent=True, title="Epithelial proportion (within epithelial)", bar_width=bar_width)
    else:
        ax_epi.set_title("Epithelial proportion (no epithelial labels)")
    if props_str.shape[1] > 0:
        _draw_stacked(ax_str, props_str, color_map, show_percent=True, title="Stromal proportion (within stromal)", bar_width=bar_width)
    else:
        ax_str.set_title("Stromal proportion (no stromal labels)")
    if props_imm.shape[1] > 0:
        _draw_stacked(ax_imm, props_imm, color_map, show_percent=True, title="Immune proportion (within immune)", bar_width=bar_width)
    else:
        ax_imm.set_title("Immune proportion (no immune labels)")

    # Group spacers and ticks
    for ax in (ax_counts, ax_all, ax_epi, ax_str, ax_imm):
        _insert_group_spacers(ax, len(ordered_samples), spacer_after, group_spacer)
        ax.set_xlim(-0.5, len(ordered_samples) - 0.5)

    _apply_xticks(ax_imm, ordered_samples, rotation=90)

    # Group labels under the bottom axis
    for grp, xmid in group_centers:
        ax_imm.text(xmid, -0.18, grp, ha="center", va="top",
                    transform=ax_imm.get_xaxis_transform(), fontsize=10)

    # tighten just the left column; legend stays fixed width
    fig.tight_layout(rect=(0, 0, 0.86, 1))  # leave room for legend column
    plt.show()

# -----------------------------
# Run
# -----------------------------
adata = sc.read(ADATA_PATH)

# Ensure cluster is categorical and colors exist
if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
    adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")
if f"{cluster_key}_colors" not in adata.uns:
    warnings.warn(f"uns lacks {cluster_key}_colors. A fallback gray will be used for all bars")
    adata.uns[f"{cluster_key}_colors"] = ["#cccccc"] * len(adata.obs[cluster_key].cat.categories)

plot_five_row_stacked(adata, cluster_key=cluster_key, library_id_key=library_id_key, library_ids=library_ids)
