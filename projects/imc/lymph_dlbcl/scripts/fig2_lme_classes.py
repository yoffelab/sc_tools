#!/usr/bin/env python3
"""
Figure 2: 5 LME (Lymphoma Microenvironment) classes.

Panels:
  a) Complex heatmap: cell type abundance per LME class (z-scored, clustered)
  b) Barplot of LME class proportions across cohort
  c) Cell type abundance per LME class (violin + BH-Wilcoxon)
  d) Per-LME cell type composition (stacked barplot)

Usage:
    python scripts/fig2_lme_classes.py

Input:
    results/adata.immune.celltyped.p4.h5ad
    results/adata.stromal.celltyped.p4.h5ad
    data/downloaded/metadata/1.13.21.merged.abundance.csv
    metadata/lme_class_assignments.csv

Output:
    figures/manuscript/fig2/
"""

import logging
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

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "fig2"
METADATA_DIR = PROJECT_DIR / "metadata"


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def significance_label(p: float) -> str:
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def load_abundance_and_lme(config: dict) -> pd.DataFrame | None:
    """Load merged abundance + LME assignments."""
    abundance_path = PROJECT_DIR / config["metadata"].get("abundance", "")
    lme_path = METADATA_DIR / "lme_class_assignments.csv"

    if abundance_path.exists() and lme_path.exists():
        abundance = pd.read_csv(abundance_path, index_col=0)
        lme = pd.read_csv(lme_path)
        lme["sample"] = lme["sample"].astype(str)
        abundance.index = abundance.index.astype(str)

        # Join
        merged = abundance.join(lme.set_index("sample")[["LME_display"]], how="inner")
        merged = merged.rename(columns={"LME_display": "LME_class"})
        merged = merged.dropna(subset=["LME_class"])
        logger.info(f"Abundance + LME: {merged.shape[0]} samples, {merged.shape[1] - 1} cell types")
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
                    logger.info(f"Computed abundance from {panel} P4: {ct_props.shape}")
                    return ct_props

    logger.warning("Could not load abundance + LME data")
    return None


def fig2a_heatmap(df: pd.DataFrame):
    """Complex heatmap: z-scored cell type abundance per sample, annotated by LME."""
    logger.info("  Fig 2a: Complex heatmap")

    lme_col = "LME_class"
    feature_cols = [c for c in df.columns if c != lme_col]

    # Normalize to proportions
    data = df[feature_cols].copy()
    row_sums = data.sum(axis=1)
    if (row_sums > 1.5).any():  # Counts, not proportions
        data = data.div(row_sums, axis=0)

    # Z-score per column
    z_data = data.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=0)

    # Sort by LME class
    order = df[lme_col].sort_values()
    z_data = z_data.loc[order.index]

    # LME color annotation
    lme_colors = {
        "Cold": "#4575b4",
        "CD206 Enriched": "#d73027",
        "Cytotoxic": "#fc8d59",
        "Stromal": "#91bfdb",
        "T cell Regulated": "#fee090",
    }
    row_colors = df.loc[order.index, lme_col].map(lme_colors)
    row_colors.name = "LME"

    g = sns.clustermap(
        z_data,
        cmap="RdBu_r",
        center=0,
        vmin=-2,
        vmax=2,
        row_cluster=False,
        col_cluster=True,
        row_colors=row_colors,
        figsize=(14, 10),
        dendrogram_ratio=(0.05, 0.1),
        cbar_kws={"label": "z-score"},
        xticklabels=True,
        yticklabels=False,
    )
    plt.setp(g.ax_heatmap.get_xticklabels(), fontsize=8, rotation=90)

    # Legend
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=c, label=l) for l, c in lme_colors.items()]
    g.ax_heatmap.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.15, 1),
                        title="LME Class", fontsize=8)

    out = FIG_DIR / "fig2a_heatmap.pdf"
    g.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig2b_lme_proportions(df: pd.DataFrame):
    """Barplot of LME class proportions."""
    logger.info("  Fig 2b: LME proportions")

    lme_counts = df["LME_class"].value_counts()
    lme_pcts = lme_counts / lme_counts.sum() * 100

    colors = {
        "Cold": "#4575b4",
        "CD206 Enriched": "#d73027",
        "Cytotoxic": "#fc8d59",
        "Stromal": "#91bfdb",
        "T cell Regulated": "#fee090",
    }

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(
        range(len(lme_pcts)),
        lme_pcts.values,
        color=[colors.get(c, "#999999") for c in lme_pcts.index],
    )
    ax.set_xticks(range(len(lme_pcts)))
    ax.set_xticklabels(lme_pcts.index, rotation=30, ha="right")
    ax.set_ylabel("% of Cohort")
    ax.set_title("LME Class Distribution")

    for bar, pct, count in zip(bars, lme_pcts.values, lme_counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f"{pct:.1f}%\n(n={count})", ha="center", va="bottom", fontsize=9)

    ax.set_ylim(0, max(lme_pcts.values) * 1.2)
    plt.tight_layout()

    out = FIG_DIR / "fig2b_lme_proportions.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig2c_abundance_violin(df: pd.DataFrame):
    """Cell type abundance per LME class with Wilcoxon tests."""
    logger.info("  Fig 2c: Abundance violins")

    lme_col = "LME_class"
    feature_cols = [c for c in df.columns if c != lme_col]

    # Normalize
    data = df[feature_cols].copy()
    row_sums = data.sum(axis=1)
    if (row_sums > 1.5).any():
        data = data.div(row_sums, axis=0)
    data[lme_col] = df[lme_col]

    # Select top variable cell types
    variances = data[feature_cols].var().sort_values(ascending=False)
    top_features = variances.head(min(12, len(feature_cols))).index.tolist()

    n_features = len(top_features)
    ncols = 4
    nrows = (n_features + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5))
    axes = axes.flatten()

    lme_order = ["Cold", "Stromal", "Cytotoxic", "CD206 Enriched", "T cell Regulated"]
    lme_order = [l for l in lme_order if l in data[lme_col].unique()]

    for i, feat in enumerate(top_features):
        ax = axes[i]
        sns.violinplot(data=data, x=lme_col, y=feat, order=lme_order, ax=ax,
                       cut=0, inner="box", linewidth=0.5, scale="width")
        ax.set_title(feat, fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("Proportion")
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)

    # Hide empty subplots
    for i in range(n_features, len(axes)):
        axes[i].set_visible(False)

    plt.suptitle("Cell Type Abundance by LME Class", fontsize=12)
    plt.tight_layout()

    out = FIG_DIR / "fig2c_abundance_violin.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig2d_composition_stacked(df: pd.DataFrame):
    """Stacked barplot of cell type composition per LME class."""
    logger.info("  Fig 2d: LME composition stacked barplot")

    lme_col = "LME_class"
    feature_cols = [c for c in df.columns if c != lme_col]

    data = df[feature_cols].copy()
    row_sums = data.sum(axis=1)
    if (row_sums > 1.5).any():
        data = data.div(row_sums, axis=0)
    data[lme_col] = df[lme_col]

    mean_props = data.groupby(lme_col)[feature_cols].mean()

    fig, ax = plt.subplots(figsize=(10, 6))
    mean_props.plot(kind="bar", stacked=True, ax=ax, width=0.8)
    ax.set_ylabel("Mean Proportion")
    ax.set_xlabel("")
    ax.set_title("Mean Cell Type Composition per LME Class")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()

    out = FIG_DIR / "fig2d_composition_stacked.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()

    df = load_abundance_and_lme(config)
    if df is None:
        logger.error("Cannot generate Figure 2 without abundance + LME data")
        return

    fig2a_heatmap(df)
    fig2b_lme_proportions(df)
    fig2c_abundance_violin(df)
    fig2d_composition_stacked(df)

    logger.info("Figure 2 complete.")


if __name__ == "__main__":
    main()
