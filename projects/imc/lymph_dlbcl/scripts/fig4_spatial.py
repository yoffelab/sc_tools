#!/usr/bin/env python3
"""
Figure 4: Spatial community analysis (20 communities from k=30 kNN).

Manuscript caption (adapted):
  a) Cell type x community composition heatmap (z-scored per community)
  b) Neighbor diversity violins per community
  c) Community x LME enrichment bubble plot (orange=enriched, purple=depleted)
  d) Gradient boosted model AUC + feature ranking for LME prediction

Insight: Spatial communities define local tissue niches that are differentially
distributed across LME classes. Community composition predicts LME classification.

Usage:
    python scripts/fig4_spatial.py

Input:
    results/adata.stromal.spatial.communities.h5ad
    data/downloaded/metadata/7.31.22.k30_per_patient_community_abundance.csv
    metadata/lme_class_assignments.csv

Output:
    figures/manuscript/fig4/
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
    LME_COLORS,
    LME_ORDER,
    apply_figure_style,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "fig4"
METADATA_DIR = PROJECT_DIR / "metadata"


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def fig4a_community_composition(adata: ad.AnnData):
    """Cell type x community composition heatmap (z-scored per community).

    Insight: Each community has a distinct cell type composition.
    Z-scoring per community reveals relative enrichment.
    """
    logger.info("  Fig 4a: Community composition heatmap")

    comm_col = None
    for col in ["community", "cluster", "community_id", "seurat_clusters"]:
        if col in adata.obs.columns:
            comm_col = col
            break

    ct_col = None
    for col in ["celltype", "celltype_broad", "major_group", "labels", "meta"]:
        if col in adata.obs.columns and 2 <= adata.obs[col].nunique() <= 50:
            ct_col = col
            break

    if comm_col is None or ct_col is None:
        logger.warning(f"  Missing columns: community={comm_col}, celltype={ct_col}")
        return

    ct_counts = adata.obs.groupby([comm_col, ct_col]).size().unstack(fill_value=0)
    ct_props = ct_counts.div(ct_counts.sum(axis=1), axis=0)

    # Z-score per community (row) — highlights enriched cell types within each community
    z_scored = ct_props.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=1)

    # Remove NaN/Inf values for clustermap compatibility
    z_scored = z_scored.dropna(how="all", axis=0).dropna(how="all", axis=1).fillna(0)

    g = sns.clustermap(
        z_scored,
        cmap="RdBu_r",
        center=0,
        vmin=-2,
        vmax=2,
        figsize=(max(8, len(z_scored.columns) * 0.5), max(6, len(z_scored) * 0.3)),
        dendrogram_ratio=(0.08, 0.08),
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "z-score", "shrink": 0.6},
        linewidths=0.3,
    )
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
    g.ax_heatmap.set_xlabel("Cell Type")
    g.ax_heatmap.set_ylabel("Community")

    out = FIG_DIR / "fig4a_community_composition.pdf"
    g.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig4b_neighbor_diversity(adata: ad.AnnData):
    """Neighbor diversity (Shannon entropy) per community.

    Insight: Some communities are highly homogeneous (low diversity)
    while others are mixed niches (high diversity).
    """
    logger.info("  Fig 4b: Neighbor diversity")

    comm_col = None
    for col in ["community", "cluster", "community_id", "seurat_clusters"]:
        if col in adata.obs.columns:
            comm_col = col
            break

    ct_col = None
    for col in ["celltype", "celltype_broad", "labels", "meta"]:
        if col in adata.obs.columns and 2 <= adata.obs[col].nunique() <= 50:
            ct_col = col
            break

    if comm_col is None or ct_col is None:
        return

    from scipy.stats import entropy as shannon_entropy

    records = []

    for name, group in adata.obs.groupby(comm_col):
        ct_counts = group[ct_col].value_counts(normalize=True)
        h = shannon_entropy(ct_counts)
        records.append({"community": name, "shannon_entropy": h, "n_cells": len(group)})

    div_df = pd.DataFrame(records).sort_values("shannon_entropy", ascending=False)

    fig, ax = plt.subplots(figsize=(max(6, len(div_df) * 0.3), 4))
    colors = ["#fc8d59" if h > div_df["shannon_entropy"].median() else "#91bfdb"
              for h in div_df["shannon_entropy"]]
    ax.bar(range(len(div_df)), div_df["shannon_entropy"], color=colors)
    ax.set_xticks(range(len(div_df)))
    ax.set_xticklabels(div_df["community"], rotation=45, ha="right")
    ax.set_ylabel("Shannon Entropy")
    ax.set_xlabel("Community")
    ax.set_title("Cell Type Diversity per Community")
    ax.axhline(div_df["shannon_entropy"].median(), color="gray", linestyle="--",
               linewidth=0.8, label=f"Median={div_df['shannon_entropy'].median():.2f}")
    ax.legend(fontsize=6)
    sns.despine()

    out = FIG_DIR / "fig4b_neighbor_diversity.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig4c_community_lme_enrichment(config: dict):
    """Community x LME enrichment bubble plot.

    Insight: Specific communities are enriched (orange) or depleted (purple)
    in particular LME classes. This links spatial organization to TME phenotype.
    """
    logger.info("  Fig 4c: Community x LME enrichment bubble")

    comm_path = PROJECT_DIR / config["metadata"].get("community_abundance", "")
    lme_path = METADATA_DIR / "lme_class_assignments.csv"

    if not (comm_path.exists() and lme_path.exists()):
        logger.warning(f"  Missing: community={comm_path.exists()}, lme={lme_path.exists()}")
        return

    comm = pd.read_csv(comm_path, index_col=0)
    lme = pd.read_csv(lme_path)
    lme["sample"] = lme["sample"].astype(str)
    comm.index = comm.index.astype(str)

    merged = comm.join(lme.set_index("sample")[["LME_class"]], how="inner")
    merged = merged.dropna(subset=["LME_class"])

    if merged.empty:
        logger.warning("  No overlap between community and LME data")
        return

    feature_cols = [c for c in merged.columns if c != "LME_class"]
    lme_present = [l for l in LME_ORDER if l in merged["LME_class"].unique()]

    # Compute enrichment: z-score of mean per LME vs overall mean
    overall_mean = merged[feature_cols].mean()
    overall_std = merged[feature_cols].std()

    # Also compute Fisher-like p-values (above/below median per community)
    records = []
    for lme in lme_present:
        lme_data = merged[merged["LME_class"] == lme][feature_cols]
        for comm_name in feature_cols:
            z = (lme_data[comm_name].mean() - overall_mean[comm_name]) / (overall_std[comm_name] + 1e-10)
            # Wilcoxon test: LME vs rest
            rest_data = merged[merged["LME_class"] != lme][comm_name]
            if len(lme_data) > 2 and len(rest_data) > 2:
                _, p = stats.mannwhitneyu(lme_data[comm_name], rest_data, alternative="two-sided")
            else:
                p = 1.0
            records.append({
                "LME": lme, "Community": comm_name,
                "z_score": z, "p_value": p,
                "mean_abundance": lme_data[comm_name].mean(),
            })

    enrich_df = pd.DataFrame(records)

    # BH correction
    _, adj_p, _, _ = multipletests(enrich_df["p_value"], method="fdr_bh")
    enrich_df["adj_p"] = adj_p
    enrich_df["neg_log_p"] = -np.log10(enrich_df["adj_p"] + 1e-300)

    # Pivot for bubble plot
    z_pivot = enrich_df.pivot(index="LME", columns="Community", values="z_score")
    p_pivot = enrich_df.pivot(index="LME", columns="Community", values="neg_log_p")

    z_pivot = z_pivot.loc[lme_present]
    p_pivot = p_pivot.loc[lme_present]

    fig, ax = plt.subplots(figsize=(max(8, len(feature_cols) * 0.4), max(4, len(lme_present) * 0.6)))

    for i, lme in enumerate(lme_present):
        for j, comm_name in enumerate(feature_cols):
            z = z_pivot.loc[lme, comm_name]
            neg_log_p = p_pivot.loc[lme, comm_name]
            size = min(neg_log_p * 15, 200)  # cap size
            color = "#fc8d59" if z > 0 else "#8856a7"  # orange=enriched, purple=depleted
            alpha = min(0.3 + abs(z) * 0.2, 0.9)
            if neg_log_p > -np.log10(0.05):  # significant
                ax.scatter(j, i, s=size, c=color, alpha=alpha, edgecolors="black", linewidth=0.3)
            else:
                ax.scatter(j, i, s=size * 0.3, c="#cccccc", alpha=0.3, edgecolors="none")

    ax.set_xticks(range(len(feature_cols)))
    ax.set_xticklabels(feature_cols, rotation=90)
    ax.set_yticks(range(len(lme_present)))
    ax.set_yticklabels(lme_present)
    ax.set_title("Community x LME Enrichment")
    ax.set_xlabel("Community")

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#fc8d59",
               markersize=8, label="Enriched"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#8856a7",
               markersize=8, label="Depleted"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#cccccc",
               markersize=6, label="Not significant"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=6)
    sns.despine()

    out = FIG_DIR / "fig4c_community_lme_enrichment.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig4d_ml_prediction(config: dict):
    """Gradient boosted model for LME prediction from community abundance.

    Insight: Community composition accurately predicts LME class,
    confirming the link between spatial organization and TME phenotype.

    Note: Manuscript uses gradient boosted model, not random forest.
    """
    logger.info("  Fig 4d: ML prediction (gradient boosted)")

    comm_path = PROJECT_DIR / config["metadata"].get("community_abundance", "")
    lme_path = METADATA_DIR / "lme_class_assignments.csv"

    if not (comm_path.exists() and lme_path.exists()):
        logger.warning("  Missing community or LME data for ML")
        return

    try:
        from sklearn.ensemble import GradientBoostingClassifier
        from sklearn.metrics import roc_auc_score, roc_curve
        from sklearn.model_selection import StratifiedKFold
        from sklearn.preprocessing import LabelEncoder, StandardScaler
    except ImportError:
        logger.warning("  sklearn not available")
        return

    comm = pd.read_csv(comm_path, index_col=0)
    lme = pd.read_csv(lme_path)
    lme["sample"] = lme["sample"].astype(str)
    comm.index = comm.index.astype(str)

    merged = comm.join(lme.set_index("sample")[["LME_class"]], how="inner")
    merged = merged.dropna(subset=["LME_class"])

    feature_cols = [c for c in merged.columns if c != "LME_class"]
    X = merged[feature_cols].values
    le = LabelEncoder()
    y = le.fit_transform(merged["LME_class"])
    class_names = le.classes_

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Cross-validated ROC
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: ROC curves
    ax = axes[0]
    mean_aucs = []
    for train_idx, test_idx in cv.split(X_scaled, y):
        clf = GradientBoostingClassifier(
            n_estimators=200, max_depth=4, random_state=42,
            learning_rate=0.1, subsample=0.8,
        )
        clf.fit(X_scaled[train_idx], y[train_idx])
        y_proba = clf.predict_proba(X_scaled[test_idx])

        for i, cls_name in enumerate(class_names):
            if i < len(class_names):
                fpr, tpr, _ = roc_curve((y[test_idx] == i).astype(int), y_proba[:, i])
                auc = roc_auc_score((y[test_idx] == i).astype(int), y_proba[:, i])
                mean_aucs.append((cls_name, auc))

    # Average AUC per class
    auc_by_class = {}
    for cls_name, auc in mean_aucs:
        auc_by_class.setdefault(cls_name, []).append(auc)

    # Final model on all data for feature importance
    clf_final = GradientBoostingClassifier(
        n_estimators=200, max_depth=4, random_state=42,
        learning_rate=0.1, subsample=0.8,
    )
    clf_final.fit(X_scaled, y)

    # ROC on full data (for display)
    y_proba_full = clf_final.predict_proba(X_scaled)
    for i, cls_name in enumerate(class_names):
        fpr, tpr, _ = roc_curve((y == i).astype(int), y_proba_full[:, i])
        mean_auc = np.mean(auc_by_class.get(cls_name, [0]))
        color = LME_COLORS.get(cls_name, None)
        ax.plot(fpr, tpr, color=color, linewidth=1.2,
                label=f"{cls_name} (AUC={mean_auc:.2f})")

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Gradient Boosted Model — One-vs-Rest ROC")
    ax.legend(fontsize=6, loc="lower right")
    sns.despine(ax=ax)

    # Panel 2: Feature importance
    ax2 = axes[1]
    importances = clf_final.feature_importances_
    sorted_idx = np.argsort(importances)[::-1][:15]

    ax2.barh(range(len(sorted_idx)), importances[sorted_idx][::-1],
             color="#4575b4", alpha=0.8)
    ax2.set_yticks(range(len(sorted_idx)))
    ax2.set_yticklabels([feature_cols[i] for i in sorted_idx][::-1])
    ax2.set_xlabel("Feature Importance")
    ax2.set_title("Top 15 Features")
    sns.despine(ax=ax2)

    plt.tight_layout()

    out = FIG_DIR / "fig4d_ml_prediction.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def main():
    apply_figure_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()

    # Load spatial community object
    spatial_path = RESULTS_DIR / "adata.stromal.spatial.communities.h5ad"
    if not spatial_path.exists():
        spatial_path = PROJECT_DIR / config["objects"]["spatial_communities"]

    if spatial_path.exists():
        adata = ad.read_h5ad(spatial_path)
        logger.info(f"Spatial object: {adata.shape}")
        fig4a_community_composition(adata)
        fig4b_neighbor_diversity(adata)
    else:
        logger.warning(f"Spatial object not found: {spatial_path}")

    fig4c_community_lme_enrichment(config)
    fig4d_ml_prediction(config)

    logger.info("Figure 4 complete.")


if __name__ == "__main__":
    main()
