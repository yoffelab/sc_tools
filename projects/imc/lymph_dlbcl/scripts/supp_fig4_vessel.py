#!/usr/bin/env python3
"""Supp Fig 4: Vessel analysis — endothelial subsets, vessel density per LME."""

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
SEURAT_DIR = RESULTS_DIR / "seurat_converted"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "supp_fig4"
METADATA_DIR = PROJECT_DIR / "metadata"


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Load stromal panel (has endothelial markers: CD31, vWF)
    for suffix in ["celltyped.p4", "annotated.p2", "raw.p1"]:
        path = RESULTS_DIR / f"adata.stromal.{suffix}.h5ad"
        if path.exists():
            adata = ad.read_h5ad(path)
            break
    else:
        # Try stroma subset objects
        path = SEURAT_DIR / "stroma_2_preprocessing/S1_seurat_stroma.h5ad"
        if path.exists():
            adata = ad.read_h5ad(path)
        else:
            logger.error("No stromal data found")
            return

    logger.info(f"Loaded: {adata.shape}")

    # Find endothelial markers
    endo_markers = [m for m in adata.var_names if any(
        e in m.lower() for e in ["cd31", "vwf", "pecam", "endoglin"]
    )]
    logger.info(f"Endothelial markers found: {endo_markers}")

    # Find celltype column
    ct_col = None
    for col in ["celltype", "celltype_broad", "labels", "meta"]:
        if col in adata.obs.columns:
            ct_col = col
            break

    # Endothelial cell analysis
    if ct_col and endo_markers:
        # Violin of endothelial markers by cell type
        fig, axes = plt.subplots(1, len(endo_markers), figsize=(len(endo_markers) * 4, 5))
        if len(endo_markers) == 1:
            axes = [axes]

        X_dense = adata.X.toarray() if hasattr(adata.X, "toarray") else np.array(adata.X)
        for i, marker in enumerate(endo_markers):
            idx = list(adata.var_names).index(marker)
            plot_df = pd.DataFrame({
                "expression": X_dense[:, idx],
                "celltype": adata.obs[ct_col].values,
            })
            sns.violinplot(data=plot_df, x="celltype", y="expression", ax=axes[i],
                           cut=0, inner="box", scale="width")
            axes[i].set_title(marker)
            axes[i].set_xlabel("")
            plt.setp(axes[i].get_xticklabels(), rotation=45, ha="right", fontsize=7)

        plt.suptitle("Endothelial Marker Expression by Cell Type")
        plt.tight_layout()
        plt.savefig(FIG_DIR / "supp4a_endo_markers.pdf", dpi=300, bbox_inches="tight")
        plt.close()

    # Vessel density per LME
    lme_path = METADATA_DIR / "lme_class_assignments.csv"
    if "LME_class" in adata.obs.columns or lme_path.exists():
        if "LME_class" not in adata.obs.columns and lme_path.exists():
            lme = pd.read_csv(lme_path)
            if "sample" in adata.obs.columns:
                lme_map = dict(zip(lme["sample"].astype(str), lme["LME_display"]))
                adata.obs["LME_class"] = adata.obs["sample"].astype(str).map(lme_map)

        if "LME_class" in adata.obs.columns and ct_col:
            # Compute endothelial proportion per sample
            sample_col = "sample" if "sample" in adata.obs.columns else "orig.ident"
            if sample_col in adata.obs.columns:
                endo_keywords = ["endothelial", "vessel", "cd31", "vwf"]
                adata.obs["is_endothelial"] = adata.obs[ct_col].astype(str).str.lower().apply(
                    lambda x: any(k in x for k in endo_keywords)
                )

                vessel_density = adata.obs.groupby(sample_col).agg(
                    endo_pct=("is_endothelial", "mean"),
                    LME_class=("LME_class", "first"),
                ).dropna()

                if not vessel_density.empty and vessel_density["LME_class"].nunique() > 1:
                    fig, ax = plt.subplots(figsize=(8, 5))
                    sns.boxplot(data=vessel_density, x="LME_class", y="endo_pct", ax=ax)
                    sns.stripplot(data=vessel_density, x="LME_class", y="endo_pct",
                                  ax=ax, color="black", size=3, alpha=0.5)
                    ax.set_ylabel("Endothelial Cell Fraction")
                    ax.set_title("Vessel Density by LME Class")
                    plt.xticks(rotation=30, ha="right")
                    plt.tight_layout()
                    plt.savefig(FIG_DIR / "supp4b_vessel_density_lme.pdf", dpi=300, bbox_inches="tight")
                    plt.close()

    logger.info("Supp Fig 4 complete.")


if __name__ == "__main__":
    main()
