#!/usr/bin/env python3
"""Supp Fig 5: TME clustering sensitivity analysis."""

import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "supp_fig5"


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()

    # Load abundance data
    abundance_path = PROJECT_DIR / config["metadata"].get("abundance", "")
    if not abundance_path.exists():
        logger.error(f"Abundance file not found: {abundance_path}")
        return

    abundance = pd.read_csv(abundance_path, index_col=0)
    logger.info(f"Abundance: {abundance.shape}")

    # Normalize + scale
    row_sums = abundance.sum(axis=1)
    if (row_sums > 1.5).any():
        abundance_norm = abundance.div(row_sums, axis=0)
    else:
        abundance_norm = abundance

    scaler = StandardScaler()
    X = scaler.fit_transform(abundance_norm)

    # a) Silhouette analysis for k=2..20
    k_range = range(2, 21)
    silhouettes = []
    inertias = []

    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X)
        silhouettes.append(silhouette_score(X, labels))
        inertias.append(kmeans.inertia_)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(list(k_range), silhouettes, "o-", color="#4575b4")
    ax1.axvline(x=10, color="red", linestyle="--", alpha=0.5, label="k=10 (chosen)")
    ax1.set_xlabel("Number of Clusters (k)")
    ax1.set_ylabel("Silhouette Score")
    ax1.set_title("Silhouette Analysis")
    ax1.legend()

    ax2.plot(list(k_range), inertias, "o-", color="#d73027")
    ax2.axvline(x=10, color="red", linestyle="--", alpha=0.5, label="k=10 (chosen)")
    ax2.set_xlabel("Number of Clusters (k)")
    ax2.set_ylabel("Inertia")
    ax2.set_title("Elbow Plot")
    ax2.legend()

    plt.tight_layout()
    plt.savefig(FIG_DIR / "supp5a_cluster_metrics.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    # b) Cluster stability (bootstrap resampling)
    from sklearn.metrics import adjusted_rand_score

    n_bootstrap = 50
    k_test = [5, 8, 10, 12, 15]
    stability = {k: [] for k in k_test}

    base_labels = {}
    for k in k_test:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        base_labels[k] = kmeans.fit_predict(X)

    rng = np.random.RandomState(42)
    for _ in range(n_bootstrap):
        idx = rng.choice(len(X), size=len(X), replace=True)
        X_boot = X[idx]
        for k in k_test:
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            boot_labels = kmeans.fit_predict(X_boot)
            ari = adjusted_rand_score(base_labels[k][idx], boot_labels)
            stability[k].append(ari)

    stab_df = pd.DataFrame(stability)
    fig, ax = plt.subplots(figsize=(8, 5))
    stab_df.boxplot(ax=ax)
    ax.set_xlabel("Number of Clusters (k)")
    ax.set_ylabel("Adjusted Rand Index (bootstrap)")
    ax.set_title("Clustering Stability (50 bootstrap iterations)")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "supp5b_stability.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    # c) Heatmaps for alternative k values
    for k in [5, 8, 10]:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X)
        abundance_norm["cluster"] = labels

        mean_by_cluster = abundance_norm.groupby("cluster").mean()
        z_scored = mean_by_cluster.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=0)

        g = sns.clustermap(z_scored, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
                           figsize=(12, max(4, k * 0.5)),
                           xticklabels=True, yticklabels=True)
        g.ax_heatmap.set_title(f"k={k} cluster means (z-scored)")
        g.savefig(FIG_DIR / f"supp5c_heatmap_k{k}.pdf", dpi=300, bbox_inches="tight")
        plt.close()

        abundance_norm.drop(columns=["cluster"], inplace=True)

    logger.info("Supp Fig 5 complete.")


if __name__ == "__main__":
    main()
