#!/usr/bin/env python3
"""
Figure 4: Spatial community analysis (20 communities from k=30 kNN).

Panels:
  a) Heatmap of community composition (cell type proportions per community)
  b) Spatial maps of community assignments
  c) Neighborhood enrichment analysis
  d) Community abundance per LME class

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
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

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
    """Heatmap of cell type proportions per spatial community."""
    logger.info("  Fig 4a: Community composition heatmap")

    comm_col = None
    for col in ["community", "cluster", "community_id"]:
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

    g = sns.clustermap(
        ct_props,
        cmap="YlOrRd",
        figsize=(12, max(6, len(ct_props) * 0.35)),
        dendrogram_ratio=(0.1, 0.1),
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "Proportion"},
    )
    plt.setp(g.ax_heatmap.get_xticklabels(), fontsize=8, rotation=90)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=8)
    g.ax_heatmap.set_xlabel("Cell Type")
    g.ax_heatmap.set_ylabel("Community")

    out = FIG_DIR / "fig4a_community_composition.pdf"
    g.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig4b_spatial_maps(adata: ad.AnnData):
    """Spatial maps showing community assignments."""
    logger.info("  Fig 4b: Spatial maps")

    if "spatial" not in adata.obsm:
        # Try to find spatial coords
        spatial_cols = [c for c in adata.obs.columns if "centroid" in c.lower() or "x" == c.lower() or "y" == c.lower()]
        if len(spatial_cols) >= 2:
            adata.obsm["spatial"] = adata.obs[spatial_cols[:2]].values
        else:
            logger.warning("  No spatial coordinates found")
            return

    comm_col = None
    for col in ["community", "cluster", "community_id"]:
        if col in adata.obs.columns:
            comm_col = col
            break

    if comm_col is None:
        return

    # Plot a few representative samples
    if "sample" in adata.obs.columns or "orig.ident" in adata.obs.columns:
        sample_col = "sample" if "sample" in adata.obs.columns else "orig.ident"
        samples = adata.obs[sample_col].value_counts().head(4).index.tolist()

        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()

        for i, sample in enumerate(samples):
            mask = adata.obs[sample_col] == sample
            sub = adata[mask]
            coords = sub.obsm["spatial"]

            ax = axes[i]
            scatter = ax.scatter(
                coords[:, 0], coords[:, 1],
                c=pd.Categorical(sub.obs[comm_col]).codes,
                cmap="tab20",
                s=1, alpha=0.7,
            )
            ax.set_title(f"Sample: {sample}", fontsize=10)
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set_xlabel("X")
            ax.set_ylabel("Y")

        plt.suptitle("Spatial Community Assignments", fontsize=12)
        plt.tight_layout()

        out = FIG_DIR / "fig4b_spatial_maps.pdf"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"    Saved: {out}")
    else:
        logger.warning("  No sample column for spatial maps")


def fig4c_neighborhood_enrichment(adata: ad.AnnData):
    """Neighborhood enrichment analysis."""
    logger.info("  Fig 4c: Neighborhood enrichment")

    comm_col = None
    for col in ["community", "cluster"]:
        if col in adata.obs.columns:
            comm_col = col
            break

    if comm_col is None:
        return

    try:
        import squidpy as sq

        if "spatial" in adata.obsm:
            sq.gr.spatial_neighbors(adata, coord_type="generic")
            sq.gr.nhood_enrichment(adata, cluster_key=comm_col)

            fig, ax = plt.subplots(figsize=(10, 8))
            sq.pl.nhood_enrichment(adata, cluster_key=comm_col, ax=ax)

            out = FIG_DIR / "fig4c_nhood_enrichment.pdf"
            plt.savefig(out, dpi=300, bbox_inches="tight")
            plt.close()
            logger.info(f"    Saved: {out}")
            return
    except ImportError:
        logger.info("  squidpy not available; computing manually")
    except Exception as e:
        logger.warning(f"  squidpy nhood_enrichment failed: {e}")

    # Manual co-occurrence matrix
    if "spatial" in adata.obsm:
        from scipy.spatial import cKDTree

        coords = adata.obsm["spatial"]
        communities = adata.obs[comm_col].values
        unique_comm = sorted(adata.obs[comm_col].unique())

        tree = cKDTree(coords)
        pairs = tree.query_pairs(r=50)  # 50 pixel radius

        cooccurrence = pd.DataFrame(0, index=unique_comm, columns=unique_comm, dtype=float)
        for i, j in pairs:
            c1, c2 = communities[i], communities[j]
            cooccurrence.loc[c1, c2] += 1
            cooccurrence.loc[c2, c1] += 1

        # Normalize
        row_sums = cooccurrence.sum(axis=1)
        cooccurrence_norm = cooccurrence.div(row_sums + 1e-10, axis=0)

        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(cooccurrence_norm, cmap="YlOrRd", ax=ax, xticklabels=True, yticklabels=True)
        ax.set_title("Spatial Co-occurrence (normalized)")
        plt.tight_layout()

        out = FIG_DIR / "fig4c_cooccurrence.pdf"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"    Saved: {out}")


def fig4d_community_by_lme(config: dict):
    """Community abundance per LME class."""
    logger.info("  Fig 4d: Community abundance by LME")

    comm_path = PROJECT_DIR / config["metadata"].get("community_abundance", "")
    lme_path = METADATA_DIR / "lme_class_assignments.csv"

    if not (comm_path.exists() and lme_path.exists()):
        logger.warning(f"  Missing files: community={comm_path.exists()}, lme={lme_path.exists()}")
        return

    comm = pd.read_csv(comm_path, index_col=0)
    lme = pd.read_csv(lme_path)
    lme["sample"] = lme["sample"].astype(str)
    comm.index = comm.index.astype(str)

    merged = comm.join(lme.set_index("sample")[["LME_display"]], how="inner")
    merged = merged.rename(columns={"LME_display": "LME_class"}).dropna(subset=["LME_class"])

    if merged.empty:
        logger.warning("  No overlap between community and LME data")
        return

    feature_cols = [c for c in merged.columns if c != "LME_class"]
    mean_by_lme = merged.groupby("LME_class")[feature_cols].mean()

    # Z-score per community
    z_scored = mean_by_lme.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=0)

    fig, ax = plt.subplots(figsize=(max(10, len(feature_cols) * 0.4), 6))
    sns.heatmap(z_scored, cmap="RdBu_r", center=0, vmin=-2, vmax=2, ax=ax,
                xticklabels=True, yticklabels=True, cbar_kws={"label": "z-score"})
    ax.set_xlabel("Community")
    ax.set_ylabel("LME Class")
    ax.set_title("Community Abundance by LME Class (z-scored)")
    plt.setp(ax.get_xticklabels(), fontsize=8, rotation=90)
    plt.tight_layout()

    out = FIG_DIR / "fig4d_community_by_lme.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def main():
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
        fig4b_spatial_maps(adata)
        fig4c_neighborhood_enrichment(adata)
    else:
        logger.warning(f"Spatial object not found: {spatial_path}")

    fig4d_community_by_lme(config)

    logger.info("Figure 4 complete.")


if __name__ == "__main__":
    main()
