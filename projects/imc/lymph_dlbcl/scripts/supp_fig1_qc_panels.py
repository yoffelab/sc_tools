#!/usr/bin/env python3
"""Supp Fig 1: Panel composition and QC metrics."""

import logging
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import apply_figure_style

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "supp_fig1"


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def main():
    apply_figure_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config()

    for panel in config["panels"]:
        for suffix in ["celltyped.p4", "annotated.p2", "raw.p1"]:
            path = RESULTS_DIR / f"adata.{panel}.{suffix}.h5ad"
            if path.exists():
                logger.info(f"Loading {panel}: {path}")
                adata = ad.read_h5ad(path)
                break
        else:
            logger.warning(f"No checkpoint for {panel}")
            continue

        # Use layers['raw'] if X is empty (p4 may have zeroed X)
        data_matrix = adata.layers.get("raw", adata.X) if adata.layers else adata.X

        # a) Marker list
        fig, ax = plt.subplots(figsize=(8, max(4, adata.n_vars * 0.25)))
        ax.barh(range(adata.n_vars), np.array(data_matrix.mean(axis=0)).flatten()
                if not hasattr(data_matrix, "toarray") else np.array(data_matrix.toarray().mean(axis=0)).flatten())
        ax.set_yticks(range(adata.n_vars))
        ax.set_yticklabels(adata.var_names, fontsize=7)
        ax.set_xlabel("Mean Intensity")
        ax.set_title(f"{panel.title()} Panel Markers (n={adata.n_vars})")
        plt.tight_layout()
        plt.savefig(FIG_DIR / f"supp1a_markers_{panel}.pdf", dpi=300, bbox_inches="tight")
        plt.close()

        # b) Intensity distributions per marker
        n_markers = adata.n_vars
        ncols = 6
        nrows = (n_markers + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3, nrows * 2.5))
        axes = axes.flatten()

        X_dense = data_matrix.toarray() if hasattr(data_matrix, "toarray") else np.array(data_matrix)
        for i in range(n_markers):
            ax = axes[i]
            ax.hist(X_dense[:, i], bins=50, density=True, alpha=0.7, color="#4575b4")
            ax.set_title(adata.var_names[i], fontsize=8)
            ax.tick_params(labelsize=6)

        for i in range(n_markers, len(axes)):
            axes[i].set_visible(False)

        plt.suptitle(f"{panel.title()} Panel - Marker Intensity Distributions", fontsize=12)
        plt.tight_layout()
        plt.savefig(FIG_DIR / f"supp1b_distributions_{panel}.pdf", dpi=300, bbox_inches="tight")
        plt.close()

        # c) Per-sample nCount/nFeature
        if "nCount_Protein" in adata.obs.columns and "nFeature_Protein" in adata.obs.columns:
            sample_col = "sample" if "sample" in adata.obs.columns else "orig.ident"
            if sample_col in adata.obs.columns:
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
                adata.obs.boxplot(column="nCount_Protein", by=sample_col, ax=ax1, rot=90)
                ax1.set_title(f"{panel.title()} - nCount per Sample")
                adata.obs.boxplot(column="nFeature_Protein", by=sample_col, ax=ax2, rot=90)
                ax2.set_title(f"{panel.title()} - nFeature per Sample")
                plt.tight_layout()
                plt.savefig(FIG_DIR / f"supp1c_persample_qc_{panel}.pdf", dpi=300, bbox_inches="tight")
                plt.close()

        logger.info(f"  {panel} done")

    logger.info("Supp Fig 1 complete.")


if __name__ == "__main__":
    main()
