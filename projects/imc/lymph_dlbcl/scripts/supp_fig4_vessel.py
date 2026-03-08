#!/usr/bin/env python3
"""Supp Fig 4: Vessel analysis — endothelial markers by cell type, vessel density per LME.

Uses immune panel (has CD31 + proper cell type labels).
Stromal panel has mostly 'Unknown' cell types so is not useful here.
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

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import (
    CELLTYPE_COLORS,
    LME_COLORS,
    LME_ORDER,
    apply_figure_style,
    build_celltype_palette,
    filter_protein_markers,
)

RESULTS_DIR = PROJECT_DIR / "results"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "supp_fig4"
METADATA_DIR = PROJECT_DIR / "metadata"


def _find_best_ct_col(obs_columns):
    """Find best cell type column — prefer one with actual labels over all-Unknown."""
    # Priority order of columns to check
    candidates = ["celltype", "celltype_broad", "celltype_from_subset", "labels", "meta"]
    for col in candidates:
        if col in obs_columns:
            return col
    return None


def main():
    apply_figure_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Use immune panel — has CD31 + proper cell type labels (19 types)
    for suffix in ["celltyped.p4", "annotated.p2", "raw.p1"]:
        path = RESULTS_DIR / f"adata.immune.{suffix}.h5ad"
        if path.exists():
            logger.info(f"Loading immune panel: {path}")
            adata = ad.read_h5ad(path)
            break
    else:
        # Fallback to stromal
        for suffix in ["celltyped.p4", "annotated.p2", "raw.p1"]:
            path = RESULTS_DIR / f"adata.stromal.{suffix}.h5ad"
            if path.exists():
                logger.info(f"Loading stromal panel (fallback): {path}")
                adata = ad.read_h5ad(path)
                break
        else:
            logger.error("No panel data found")
            return

    logger.info(f"Loaded: {adata.shape}")

    # Use layers['raw'] if X is empty (p4 may have zeroed X)
    data_matrix = adata.layers.get("raw", adata.X) if adata.layers else adata.X
    X_dense = data_matrix.toarray() if hasattr(data_matrix, "toarray") else np.array(data_matrix)
    logger.info(f"Data range: {X_dense.min():.3f} to {X_dense.max():.3f}")

    # Find endothelial markers
    endo_markers = [m for m in adata.var_names if any(
        e in m.lower() for e in ["cd31", "vwf", "pecam", "endoglin"]
    )]
    logger.info(f"Endothelial markers found: {endo_markers}")

    # Find cell type column — skip if all Unknown
    ct_col = _find_best_ct_col(adata.obs.columns)
    if ct_col:
        n_unique = adata.obs[ct_col].nunique()
        all_unknown = n_unique == 1 and str(adata.obs[ct_col].iloc[0]).lower() == "unknown"
        if all_unknown:
            logger.warning(f"{ct_col} is all 'Unknown', trying alternatives")
            for alt in ["celltype_broad", "celltype_from_subset", "merged_meta", "merged_recluster"]:
                if alt in adata.obs.columns and alt != ct_col:
                    alt_unique = adata.obs[alt].nunique()
                    if alt_unique > 1:
                        ct_col = alt
                        logger.info(f"Using {ct_col} ({alt_unique} categories)")
                        break
            else:
                ct_col = None
                logger.warning("No cell type column with multiple categories found")

    # --- Panel a: Endothelial marker expression by cell type ---
    if ct_col and endo_markers:
        # Get cell types with reasonable counts, skip Unknown/other for clarity
        ct_values = adata.obs[ct_col].value_counts()
        # Keep types with >500 cells, exclude generic labels
        skip_labels = {"unknown", "other", "nan"}
        cell_types = [
            ct for ct, n in ct_values.items()
            if n > 500 and str(ct).lower() not in skip_labels
        ]
        logger.info(f"Cell types for violin (n={len(cell_types)}): {cell_types}")

        if cell_types:
            mask = adata.obs[ct_col].isin(cell_types)
            palette = build_celltype_palette(cell_types, normalize=False)

            fig, axes = plt.subplots(1, len(endo_markers),
                                     figsize=(max(len(cell_types) * 0.8, 6) * len(endo_markers), 5))
            if len(endo_markers) == 1:
                axes = [axes]

            for i, marker in enumerate(endo_markers):
                idx = list(adata.var_names).index(marker)
                plot_df = pd.DataFrame({
                    "expression": X_dense[mask, idx],
                    "celltype": adata.obs.loc[mask, ct_col].astype(str).values,
                })
                # Order by median expression (endothelial types should be highest)
                order = (plot_df.groupby("celltype", observed=True)["expression"]
                         .median().sort_values(ascending=False).index.tolist())

                sns.violinplot(data=plot_df, x="celltype", y="expression", hue="celltype",
                               ax=axes[i], order=order, palette=palette,
                               cut=0, inner="box", scale="width", linewidth=0.5, legend=False)
                axes[i].set_title(marker, fontsize=10, fontweight="bold")
                axes[i].set_xlabel("")
                axes[i].set_ylabel("Expression")
                plt.setp(axes[i].get_xticklabels(), rotation=45, ha="right", fontsize=7)

            plt.suptitle("Endothelial Marker Expression by Cell Type", fontsize=12, fontweight="bold")
            plt.tight_layout()
            plt.savefig(FIG_DIR / "supp4a_endo_markers.pdf", dpi=300, bbox_inches="tight")
            plt.close()
            logger.info("Saved supp4a_endo_markers.pdf")

    # --- Panel b: Vessel density per LME ---
    lme_path = METADATA_DIR / "lme_class_assignments.csv"
    if "LME_class" not in adata.obs.columns and lme_path.exists():
        lme = pd.read_csv(lme_path)
        if "sample" in adata.obs.columns:
            lme_map = dict(zip(lme["sample"].astype(str), lme["LME_display"], strict=False))
            adata.obs["LME_class"] = adata.obs["sample"].astype(str).map(lme_map)

    if "LME_class" in adata.obs.columns and ct_col:
        sample_col = "sample" if "sample" in adata.obs.columns else "orig.ident"
        if sample_col in adata.obs.columns:
            # Identify endothelial-like cells by cell type label or high CD31 expression
            endo_keywords = ["endothelial", "vessel", "bec", "lec", "cd31", "vwf"]
            is_endo = adata.obs[ct_col].astype(str).str.lower().apply(
                lambda x: any(k in x for k in endo_keywords)
            )
            # Also flag cells with high CD31 in "stroma" category
            if "CD31-mem" in adata.var_names:
                cd31_idx = list(adata.var_names).index("CD31-mem")
                cd31_vals = X_dense[:, cd31_idx]
                cd31_high = cd31_vals > np.percentile(cd31_vals[cd31_vals > 0], 90) if (cd31_vals > 0).any() else np.zeros(len(cd31_vals), dtype=bool)
                stroma_mask = adata.obs[ct_col].astype(str).str.lower().str.contains("stroma")
                is_endo = is_endo | (stroma_mask & cd31_high)

            adata.obs["is_endothelial"] = is_endo

            vessel_density = adata.obs.groupby(sample_col).agg(
                endo_pct=("is_endothelial", "mean"),
                LME_class=("LME_class", "first"),
            ).dropna()

            if not vessel_density.empty and vessel_density["LME_class"].nunique() > 1:
                # Order by LME_ORDER
                order = [l for l in LME_ORDER if l in vessel_density["LME_class"].values]
                palette = {k: LME_COLORS.get(k, "#999999") for k in order}

                fig, ax = plt.subplots(figsize=(8, 5))
                sns.boxplot(data=vessel_density, x="LME_class", y="endo_pct",
                            order=order, palette=palette, ax=ax, linewidth=0.8)
                sns.stripplot(data=vessel_density, x="LME_class", y="endo_pct",
                              order=order, ax=ax, color="black", size=3, alpha=0.5)
                ax.set_ylabel("Endothelial Cell Fraction")
                ax.set_xlabel("")
                ax.set_title("Vessel Density by LME Class", fontsize=12, fontweight="bold")
                # Add n per group
                for i, lme in enumerate(order):
                    n = (vessel_density["LME_class"] == lme).sum()
                    ax.text(i, ax.get_ylim()[1] * 0.95, f"n={n}", ha="center", fontsize=7)
                plt.xticks(rotation=30, ha="right")
                plt.tight_layout()
                plt.savefig(FIG_DIR / "supp4b_vessel_density_lme.pdf", dpi=300, bbox_inches="tight")
                plt.close()
                logger.info("Saved supp4b_vessel_density_lme.pdf")

    logger.info("Supp Fig 4 complete.")


if __name__ == "__main__":
    main()
