#!/usr/bin/env python3
"""
Figure 2: 5 LME (Lymphoma Microenvironment) classes.

Manuscript caption (adapted):
  a) Cell type abundance x 5 LME heatmap (z-scored per cell type across LMEs).
     Rows = cell types (grouped by lineage), columns = 5 LME classes.
  b) Stacked bar of cell type fractions per LME class.
  c) T cell / myeloid frequencies per LME (violin + BH significance).
  d) LME class proportions across the cohort (barplot with n and %).

Insight: The 5 LME classes capture distinct immune microenvironment patterns.
Cold (35%) is immune-depleted. Cytotoxic (21%) has M1 macrophages + CD8+GzmB+.
Stromal (21%) has PDPN+ CAFs. T cell Regulated (15%) has Tfh/Treg.
CD206 Enriched (8%) has CD206+ M2-like macrophages (protective).

Usage:
    python scripts/fig2_lme_classes.py

Input:
    data/downloaded/metadata/1.13.21.merged.abundance.csv
    metadata/lme_class_assignments.csv

Output:
    figures/manuscript/fig2/
"""

import logging
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import (
    CELLTYPE_COLORS,
    LME_COLORS,
    LME_ORDER,
    apply_figure_style,
    build_celltype_palette,
    significance_label,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "fig2"
METADATA_DIR = PROJECT_DIR / "metadata"


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def load_abundance_and_lme(config: dict) -> pd.DataFrame | None:
    """Load merged abundance + LME assignments."""
    abundance_path = PROJECT_DIR / config["metadata"].get("abundance", "")
    lme_path = METADATA_DIR / "lme_class_assignments.csv"

    if abundance_path.exists() and lme_path.exists():
        abundance = pd.read_csv(abundance_path, index_col=0)
        lme = pd.read_csv(lme_path)
        lme["sample"] = lme["sample"].astype(str)
        abundance.index = abundance.index.astype(str)

        merged = abundance.join(lme.set_index("sample")[["LME_class"]], how="inner")
        merged = merged.dropna(subset=["LME_class"])
        logger.info(f"Abundance + LME: {merged.shape[0]} samples, {merged.shape[1] - 1} cell types")

        # Validate all 5 LME classes present
        present = set(merged["LME_class"].unique())
        missing = set(LME_ORDER) - present
        if missing:
            logger.warning(f"Missing LME classes: {missing}")

        return merged

    # Fallback: compute from P4 AnnData
    for panel in ["stromal", "immune"]:
        p4_path = RESULTS_DIR / f"adata.{panel}.celltyped.p4.h5ad"
        if p4_path.exists():
            adata = ad.read_h5ad(p4_path)
            if "LME_class" in adata.obs.columns and "sample" in adata.obs.columns:
                ct_col = next(
                    (c for c in ["celltype", "celltype_broad", "labels", "meta"]
                     if c in adata.obs.columns),
                    None,
                )
                if ct_col:
                    ct_counts = adata.obs.groupby(["sample", ct_col]).size().unstack(fill_value=0)
                    ct_props = ct_counts.div(ct_counts.sum(axis=1), axis=0)
                    lme_map = adata.obs.groupby("sample")["LME_class"].first()
                    ct_props["LME_class"] = lme_map
                    ct_props = ct_props.dropna(subset=["LME_class"])
                    return ct_props

    logger.warning("Could not load abundance + LME data")
    return None


def fig2a_heatmap(df: pd.DataFrame):
    """Cell type abundance x LME heatmap (z-scored per cell type).

    Insight: Each LME class has a characteristic cell type composition pattern.
    Z-scoring per cell type (column) reveals relative enrichment/depletion.
    """
    logger.info("  Fig 2a: LME heatmap (z-scored per cell type)")

    feature_cols = [c for c in df.columns if c != "LME_class"]

    # Mean abundance per LME class
    data = df.copy()
    row_sums = data[feature_cols].sum(axis=1)
    if (row_sums > 1.5).any():
        data[feature_cols] = data[feature_cols].div(row_sums, axis=0)

    mean_by_lme = data.groupby("LME_class")[feature_cols].mean()

    # Reorder rows to manuscript order
    lme_present = [l for l in LME_ORDER if l in mean_by_lme.index]
    mean_by_lme = mean_by_lme.loc[lme_present]

    # Z-score per cell type (column) — highlights which LMEs are enriched for each cell type
    z_scored = mean_by_lme.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=0)

    fig, ax = plt.subplots(figsize=(max(8, len(feature_cols) * 0.5), max(4, len(lme_present) * 0.6)))
    sns.heatmap(
        z_scored,
        cmap="RdBu_r",
        center=0,
        vmin=-2,
        vmax=2,
        ax=ax,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "z-score", "shrink": 0.8},
        linewidths=0.5,
        linecolor="white",
    )
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("")
    ax.set_title("Cell Type Abundance by LME Class")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    # Add LME color patches on y-axis
    for i, lme in enumerate(lme_present):
        color = LME_COLORS.get(lme, "#999999")
        ax.add_patch(plt.Rectangle((-0.8, i), 0.3, 1, color=color,
                                   transform=ax.get_yaxis_transform(), clip_on=False))

    out = FIG_DIR / "fig2a_heatmap.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig2b_composition_stacked(df: pd.DataFrame):
    """Stacked barplot of cell type fractions per LME class.

    Insight: The relative proportion of cell types varies dramatically across LMEs.
    """
    logger.info("  Fig 2b: Stacked barplot")

    feature_cols = [c for c in df.columns if c != "LME_class"]
    data = df.copy()
    row_sums = data[feature_cols].sum(axis=1)
    if (row_sums > 1.5).any():
        data[feature_cols] = data[feature_cols].div(row_sums, axis=0)

    mean_props = data.groupby("LME_class")[feature_cols].mean()
    lme_present = [l for l in LME_ORDER if l in mean_props.index]
    mean_props = mean_props.loc[lme_present]

    # Build color list with fuzzy matching (handles varied naming conventions)
    palette = build_celltype_palette(feature_cols)
    colors = [palette[ct] for ct in feature_cols]

    fig, ax = plt.subplots(figsize=(8, 5))
    mean_props.plot(kind="bar", stacked=True, ax=ax, width=0.8, color=colors)
    ax.set_ylabel("Mean Proportion")
    ax.set_xlabel("")
    ax.set_title("Mean Cell Type Composition per LME Class")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=6, ncol=1)
    plt.xticks(rotation=30, ha="right")

    out = FIG_DIR / "fig2b_composition_stacked.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig2c_abundance_violin(df: pd.DataFrame):
    """Cell type abundance per LME class with BH-corrected Wilcoxon tests.

    Insight: Specific cell types are significantly enriched or depleted in
    particular LME classes (e.g., CD206+ macrophages in CD206 Enriched).
    Direction of effect must match the class definition.
    """
    logger.info("  Fig 2c: Abundance violins with statistics")

    feature_cols = [c for c in df.columns if c != "LME_class"]
    data = df.copy()
    row_sums = data[feature_cols].sum(axis=1)
    if (row_sums > 1.5).any():
        data[feature_cols] = data[feature_cols].div(row_sums, axis=0)

    # Select top variable cell types
    variances = data[feature_cols].var().sort_values(ascending=False)
    top_features = variances.head(min(12, len(feature_cols))).index.tolist()

    lme_present = [l for l in LME_ORDER if l in data["LME_class"].unique()]
    palette = {l: LME_COLORS[l] for l in lme_present}

    n_features = len(top_features)
    ncols = 4
    nrows = (n_features + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 3))
    if nrows == 1:
        axes = axes.reshape(1, -1)
    axes_flat = axes.flatten()

    for i, feat in enumerate(top_features):
        ax = axes_flat[i]
        sns.violinplot(
            data=data, x="LME_class", y=feat, order=lme_present, ax=ax,
            cut=0, inner="box", linewidth=0.5, scale="width", palette=palette,
        )
        # Add individual points if n < 50 per group
        max_n = data.groupby("LME_class").size().max()
        if max_n < 50:
            sns.stripplot(
                data=data, x="LME_class", y=feat, order=lme_present, ax=ax,
                size=2, alpha=0.5, color="black", jitter=True,
            )

        ax.set_title(feat, fontsize=8)
        ax.set_xlabel("")
        ax.set_ylabel("Proportion")
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=6)

        # 1-vs-rest significance for top group
        if len(lme_present) >= 2:
            p_vals = []
            for lme in lme_present:
                group_vals = data.loc[data["LME_class"] == lme, feat].dropna()
                rest_vals = data.loc[data["LME_class"] != lme, feat].dropna()
                if len(group_vals) > 2 and len(rest_vals) > 2:
                    _, p = stats.mannwhitneyu(group_vals, rest_vals, alternative="two-sided")
                    p_vals.append((lme, p))

            if p_vals:
                _, adj_p, _, _ = multipletests([p for _, p in p_vals], method="fdr_bh")
                # Annotate the most significant
                min_idx = np.argmin(adj_p)
                best_lme, best_p = p_vals[min_idx]
                best_adj = adj_p[min_idx]
                sig = significance_label(best_adj)
                if sig != "ns":
                    lme_idx = lme_present.index(best_lme)
                    ymax = data[feat].max()
                    ax.text(lme_idx, ymax * 0.95, sig, ha="center", fontsize=7, fontweight="bold")

    for i in range(n_features, len(axes_flat)):
        axes_flat[i].set_visible(False)

    plt.suptitle("Cell Type Abundance by LME Class", fontsize=10, y=1.02)
    plt.tight_layout()

    out = FIG_DIR / "fig2c_abundance_violin.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig2d_lme_proportions(df: pd.DataFrame):
    """LME class proportions across cohort.

    Insight: Cold is the most prevalent class (35%), followed by Stromal (21%)
    and Cytotoxic (21%). CD206 Enriched is the rarest (8%).
    """
    logger.info("  Fig 2d: LME proportions")

    lme_counts = df["LME_class"].value_counts()
    # Reorder to manuscript order
    lme_present = [l for l in LME_ORDER if l in lme_counts.index]
    lme_counts = lme_counts[lme_present]
    lme_pcts = lme_counts / lme_counts.sum() * 100

    colors = [LME_COLORS.get(c, "#999999") for c in lme_present]

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(range(len(lme_pcts)), lme_pcts.values, color=colors)
    ax.set_xticks(range(len(lme_pcts)))
    ax.set_xticklabels(lme_present, rotation=30, ha="right")
    ax.set_ylabel("% of Cohort")
    ax.set_title("LME Class Distribution")

    for bar, pct, count in zip(bars, lme_pcts.values, lme_counts.values, strict=False):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f"{pct:.1f}%\n(n={count})", ha="center", va="bottom", fontsize=7)

    ax.set_ylim(0, max(lme_pcts.values) * 1.25)
    sns.despine()

    out = FIG_DIR / "fig2d_lme_proportions.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def main():
    apply_figure_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()

    df = load_abundance_and_lme(config)
    if df is None:
        logger.error("Cannot generate Figure 2 without abundance + LME data")
        return

    fig2a_heatmap(df)
    fig2b_composition_stacked(df)
    fig2c_abundance_violin(df)
    fig2d_lme_proportions(df)

    logger.info("Figure 2 complete.")


if __name__ == "__main__":
    main()
