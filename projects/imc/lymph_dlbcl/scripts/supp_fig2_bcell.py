#!/usr/bin/env python3
"""Supp Fig 2: B cell subcluster analysis."""

import logging
import sys
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import apply_figure_style

RESULTS_DIR = PROJECT_DIR / "results"
SEURAT_DIR = RESULTS_DIR / "seurat_converted"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "supp_fig2"


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def main():
    apply_figure_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Load B cell subsets from both panels
    bcell_objects = {
        "immune_T2": SEURAT_DIR / "tcell_2_preprocessing/t2_bcell2_seurat.h5ad",
        "immune_merged": SEURAT_DIR / "tcell_merged/SO2_seurat_bcell.h5ad",
        "stromal_S2": SEURAT_DIR / "stroma_2_preprocessing/S1_seurat_bcell.h5ad",
        "stromal_merged": SEURAT_DIR / "stroma_merged/SO2_seurat_bcell.h5ad",
    }

    for name, path in bcell_objects.items():
        if not path.exists():
            logger.warning(f"B cell object not found: {path}")
            continue

        logger.info(f"Processing {name}: {path}")
        adata = ad.read_h5ad(path)
        logger.info(f"  Shape: {adata.shape}")

        # Find cluster/label column
        cluster_col = None
        for col in ["labels", "meta", "recluster", "seurat_clusters"]:
            if col in adata.obs.columns and 2 <= adata.obs[col].nunique() <= 30:
                cluster_col = col
                break

        if cluster_col is None:
            logger.warning(f"  No suitable cluster column for {name}")
            continue

        # UMAP
        if "X_umap" not in adata.obsm:
            if "X_pca" in adata.obsm:
                sc.pp.neighbors(adata, use_rep="X_pca")
                sc.tl.umap(adata)
            else:
                sc.pp.pca(adata, n_comps=min(30, adata.n_vars - 1))
                sc.pp.neighbors(adata)
                sc.tl.umap(adata)

        fig, ax = plt.subplots(figsize=(8, 7))
        sc.pl.umap(adata, color=cluster_col, ax=ax, show=False, size=2, frameon=False,
                    title=f"B cell subclusters ({name})")
        plt.savefig(FIG_DIR / f"supp2a_umap_{name}.pdf", dpi=300, bbox_inches="tight")
        plt.close()

        # Marker heatmap
        groups = adata.obs[cluster_col].unique()
        if len(groups) <= 30:
            mean_expr = pd.DataFrame(index=groups, columns=adata.var_names, dtype=float)
            # Use layers['raw'] if X is empty (p4 may have zeroed X)
            data_matrix = adata.layers.get("raw", adata.X) if adata.layers else adata.X
            X_dense = data_matrix.toarray() if hasattr(data_matrix, "toarray") else np.array(data_matrix)
            for g in groups:
                mask = adata.obs[cluster_col] == g
                mean_expr.loc[g] = X_dense[mask].mean(axis=0)

            z_scored = mean_expr.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=0)

            # Remove NaN/Inf values for clustermap compatibility
            z_scored = z_scored.dropna(how="all", axis=0).dropna(how="all", axis=1).fillna(0)

            g = sns.clustermap(z_scored, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
                               figsize=(12, max(4, len(groups) * 0.35)),
                               xticklabels=True, yticklabels=True)
            plt.setp(g.ax_heatmap.get_xticklabels(), fontsize=7, rotation=90)
            g.savefig(FIG_DIR / f"supp2b_heatmap_{name}.pdf", dpi=300, bbox_inches="tight")
            plt.close()

        logger.info(f"  {name} done")

    logger.info("Supp Fig 2 complete.")


if __name__ == "__main__":
    main()
