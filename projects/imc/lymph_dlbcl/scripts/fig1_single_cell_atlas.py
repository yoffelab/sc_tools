#!/usr/bin/env python3
"""
Figure 1: Single-cell atlas of DLBCL tumor microenvironment.

Panels:
  a) UMAP colored by cell type (both panels)
  b) Marker expression heatmap (mean per cell type)
  c) Dotplot of key markers per cell type
  d) Stacked barplot of cell type proportions per sample

Usage:
    python scripts/fig1_single_cell_atlas.py

Input:
    results/adata.immune.celltyped.p4.h5ad
    results/adata.stromal.celltyped.p4.h5ad

Output:
    figures/manuscript/fig1/fig1a_umap_immune.pdf
    figures/manuscript/fig1/fig1a_umap_stromal.pdf
    figures/manuscript/fig1/fig1b_heatmap.pdf
    figures/manuscript/fig1/fig1c_dotplot.pdf
    figures/manuscript/fig1/fig1d_proportions.pdf
"""

import logging
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "fig1"


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def load_panel(panel: str) -> ad.AnnData | None:
    """Load celltyped panel data."""
    for suffix in ["celltyped.p4", "annotated.p2", "raw.p1"]:
        path = RESULTS_DIR / f"adata.{panel}.{suffix}.h5ad"
        if path.exists():
            logger.info(f"Loading {panel}: {path}")
            return ad.read_h5ad(path)
    logger.warning(f"No checkpoint found for {panel}")
    return None


def get_celltype_col(adata: ad.AnnData) -> str:
    """Find best cell type column."""
    for col in ["celltype", "celltype_broad", "labels", "meta"]:
        if col in adata.obs.columns and 2 <= adata.obs[col].nunique() <= 50:
            return col
    return "seurat_clusters"


def fig1a_umap(adata: ad.AnnData, panel: str, celltype_col: str):
    """UMAP colored by cell type."""
    logger.info(f"  Fig 1a: UMAP for {panel}")

    # Compute UMAP if not present
    if "X_umap" not in adata.obsm:
        if "X_pca" not in adata.obsm:
            logger.info("    Computing PCA...")
            sc.pp.pca(adata, n_comps=30)
        logger.info("    Computing neighbors + UMAP...")
        sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")
        sc.tl.umap(adata)

    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color=celltype_col,
        ax=ax,
        show=False,
        frameon=False,
        title=f"{panel.title()} Panel - Cell Types",
        legend_loc="on data" if adata.obs[celltype_col].nunique() <= 12 else "right margin",
        size=1,
    )
    out = FIG_DIR / f"fig1a_umap_{panel}.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig1b_heatmap(adatas: dict, celltype_cols: dict):
    """Marker expression heatmap (mean per cell type, both panels)."""
    logger.info("  Fig 1b: Marker heatmap")

    all_mean = []
    for panel, adata in adatas.items():
        col = celltype_cols[panel]
        groups = adata.obs[col].unique()
        for group in groups:
            mask = adata.obs[col] == group
            if hasattr(adata.X, "toarray"):
                mean_vals = np.array(adata.X[mask].toarray().mean(axis=0)).flatten()
            else:
                mean_vals = np.array(adata.X[mask].mean(axis=0)).flatten()
            for j, marker in enumerate(adata.var_names):
                all_mean.append({
                    "celltype": str(group),
                    "marker": marker,
                    "panel": panel,
                    "mean_expr": mean_vals[j],
                })

    df = pd.DataFrame(all_mean)
    pivot = df.pivot_table(index="celltype", columns="marker", values="mean_expr", aggfunc="mean")

    # Z-score per marker
    z_scored = pivot.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=0)

    fig_h = max(8, len(z_scored) * 0.35)
    fig_w = max(12, len(z_scored.columns) * 0.25)

    g = sns.clustermap(
        z_scored,
        cmap="RdBu_r",
        center=0,
        vmin=-2,
        vmax=2,
        figsize=(fig_w, fig_h),
        dendrogram_ratio=(0.1, 0.15),
        cbar_kws={"label": "z-score"},
        xticklabels=True,
        yticklabels=True,
    )
    g.ax_heatmap.set_xlabel("Markers")
    g.ax_heatmap.set_ylabel("Cell Type")
    plt.setp(g.ax_heatmap.get_xticklabels(), fontsize=7, rotation=90)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=8)

    out = FIG_DIR / "fig1b_heatmap.pdf"
    g.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig1c_dotplot(adatas: dict, celltype_cols: dict):
    """Dotplot of key markers per cell type."""
    logger.info("  Fig 1c: Dotplot")

    # Key markers for IMC
    key_markers_immune = [
        "CD20", "PAX5", "CD3", "CD4", "CD8", "PD1", "CD68", "CD163",
        "CD11c", "HLA-DR", "Ki67", "CD31",
    ]
    key_markers_stromal = [
        "CD20", "PDPN", "VIM", "FN", "CD31", "vWF", "CD68", "CD163",
        "CD206", "CXCL12", "CD73",
    ]

    for panel, adata in adatas.items():
        col = celltype_cols[panel]
        markers = key_markers_immune if panel == "immune" else key_markers_stromal

        # Filter to markers present in data
        available = [m for m in markers if m in adata.var_names]
        if not available:
            # Try case-insensitive
            var_lower = {v.lower(): v for v in adata.var_names}
            available = [var_lower[m.lower()] for m in markers if m.lower() in var_lower]

        if not available:
            logger.warning(f"    No key markers found for {panel}")
            continue

        fig, ax = plt.subplots(figsize=(len(available) * 0.6 + 2, adata.obs[col].nunique() * 0.4 + 2))
        sc.pl.dotplot(adata, var_names=available, groupby=col, ax=ax, show=False)

        out = FIG_DIR / f"fig1c_dotplot_{panel}.pdf"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"    Saved: {out}")


def fig1d_proportions(adatas: dict, celltype_cols: dict):
    """Stacked barplot of cell type proportions per sample."""
    logger.info("  Fig 1d: Cell type proportions")

    for panel, adata in adatas.items():
        col = celltype_cols[panel]
        if "sample" not in adata.obs.columns:
            logger.warning(f"    No 'sample' column in {panel}")
            continue

        # Compute proportions
        ct_counts = adata.obs.groupby(["sample", col]).size().unstack(fill_value=0)
        ct_props = ct_counts.div(ct_counts.sum(axis=1), axis=0)

        # Sort samples by dominant cell type for visual clarity
        ct_props = ct_props.loc[ct_props.max(axis=1).sort_values(ascending=False).index]

        fig, ax = plt.subplots(figsize=(max(12, len(ct_props) * 0.15), 6))
        ct_props.plot(kind="bar", stacked=True, ax=ax, width=0.9, legend=True)
        ax.set_ylabel("Proportion")
        ax.set_xlabel("Sample")
        ax.set_title(f"{panel.title()} Panel - Cell Type Proportions")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)

        if len(ct_props) > 50:
            ax.set_xticklabels([])
        else:
            plt.xticks(rotation=90, fontsize=6)

        out = FIG_DIR / f"fig1d_proportions_{panel}.pdf"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"    Saved: {out}")


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()

    adatas = {}
    celltype_cols = {}

    for panel in config["panels"]:
        adata = load_panel(panel)
        if adata is not None:
            adatas[panel] = adata
            celltype_cols[panel] = get_celltype_col(adata)
            logger.info(f"  {panel}: {adata.n_obs:,} cells, celltype_col='{celltype_cols[panel]}'")

    if not adatas:
        logger.error("No panel data loaded")
        return

    # Generate figure panels
    for panel, adata in adatas.items():
        fig1a_umap(adata, panel, celltype_cols[panel])

    if len(adatas) > 0:
        fig1b_heatmap(adatas, celltype_cols)
        fig1c_dotplot(adatas, celltype_cols)
        fig1d_proportions(adatas, celltype_cols)

    logger.info("Figure 1 complete.")


if __name__ == "__main__":
    main()
