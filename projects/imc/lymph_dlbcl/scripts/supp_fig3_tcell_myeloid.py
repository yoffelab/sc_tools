#!/usr/bin/env python3
"""Supp Fig 3: T cell and myeloid characterization."""

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
SEURAT_DIR = PROJECT_DIR / "results" / "seurat_converted"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "supp_fig3"


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    subsets = {
        # T cell subsets
        "tcell_T2": SEURAT_DIR / "tcell_2_preprocessing/1.29.23.tcell_seurat.h5ad",
        "tcell_T2_v2": SEURAT_DIR / "tcell_2_preprocessing/t2_tcell2_seurat.h5ad",
        # Myeloid subsets
        "myeloid_T2": SEURAT_DIR / "tcell_2_preprocessing/t2_myeloid2_seurat.h5ad",
        "myeloid_S1": SEURAT_DIR / "stroma_1_preprocessing/S1_seurat_myeloid2.h5ad",
    }

    for name, path in subsets.items():
        if not path.exists():
            logger.warning(f"Not found: {path}")
            continue

        logger.info(f"Processing {name}: {path}")
        adata = ad.read_h5ad(path)

        cluster_col = None
        for col in ["labels", "meta", "recluster", "seurat_clusters"]:
            if col in adata.obs.columns and 2 <= adata.obs[col].nunique() <= 30:
                cluster_col = col
                break

        if cluster_col is None:
            continue

        # UMAP
        if "X_umap" not in adata.obsm:
            if "X_pca" in adata.obsm:
                sc.pp.neighbors(adata, use_rep="X_pca")
            else:
                sc.pp.pca(adata, n_comps=min(30, adata.n_vars - 1))
                sc.pp.neighbors(adata)
            sc.tl.umap(adata)

        fig, ax = plt.subplots(figsize=(8, 7))
        sc.pl.umap(adata, color=cluster_col, ax=ax, show=False, size=3, frameon=False,
                    title=f"{name} subclusters")
        plt.savefig(FIG_DIR / f"supp3_umap_{name}.pdf", dpi=300, bbox_inches="tight")
        plt.close()

        # Marker heatmap
        groups = adata.obs[cluster_col].unique()
        if len(groups) <= 30:
            X_dense = adata.X.toarray() if hasattr(adata.X, "toarray") else np.array(adata.X)
            mean_expr = pd.DataFrame(index=groups, columns=adata.var_names, dtype=float)
            for g in groups:
                mean_expr.loc[g] = X_dense[adata.obs[cluster_col] == g].mean(axis=0)

            z_scored = mean_expr.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=0)
            g = sns.clustermap(z_scored, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
                               figsize=(12, max(4, len(groups) * 0.35)),
                               xticklabels=True, yticklabels=True)
            g.savefig(FIG_DIR / f"supp3_heatmap_{name}.pdf", dpi=300, bbox_inches="tight")
            plt.close()

        # Proportions per sample
        sample_col = "sample" if "sample" in adata.obs.columns else "orig.ident"
        if sample_col in adata.obs.columns:
            props = adata.obs.groupby([sample_col, cluster_col]).size().unstack(fill_value=0)
            props = props.div(props.sum(axis=1), axis=0)
            fig, ax = plt.subplots(figsize=(max(10, len(props) * 0.1), 5))
            props.plot(kind="bar", stacked=True, ax=ax, width=0.9, legend=True)
            ax.set_ylabel("Proportion")
            ax.set_title(f"{name} proportions per sample")
            if len(props) > 50:
                ax.set_xticklabels([])
            ax.legend(bbox_to_anchor=(1.02, 1), fontsize=7)
            plt.tight_layout()
            plt.savefig(FIG_DIR / f"supp3_proportions_{name}.pdf", dpi=300, bbox_inches="tight")
            plt.close()

        logger.info(f"  {name} done")

    logger.info("Supp Fig 3 complete.")


if __name__ == "__main__":
    main()
