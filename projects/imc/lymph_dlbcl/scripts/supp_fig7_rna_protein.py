#!/usr/bin/env python3
"""Supp Fig 7: RNA-protein comparison (IMC vs CosMx).

Correlates mean IMC protein expression vs mean CosMx RNA expression
for overlapping markers across matched cell types.

Panels:
  7a: Scatter plots per marker (IMC protein vs CosMx RNA, Pearson r)
  7b: Summary correlation heatmap across all markers

Input:
  - IMC: results/adata.immune.celltyped.p4.h5ad (IMC protein expression)
  - CosMx: projects/cosmx_1k/lymph_dlbcl/results/cosmx_rna_annotated.h5ad
Output: figures/manuscript/supp_fig7/
"""

import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
COSMX_DIR = Path(__file__).resolve().parent.parent.parent.parent / "cosmx_1k" / "lymph_dlbcl"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import apply_figure_style, normalize_celltype

FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "supp_fig7"

# Overlapping markers between IMC protein panel and CosMx 1000-plex RNA
# IMC var_name (Protein-compartment) -> gene symbol in CosMx
IMC_TO_GENE = {
    "CD3-mem": "CD3E",
    "CD4-mem": "CD4",
    "CD8-mem": "CD8A",
    "CD20-mem": "MS4A1",
    "CD68-mem": "CD68",
    "CD163-mem": "CD163",
    "PDPN-mem": "PDPN",
    "CD31-mem": "PECAM1",
    "Ki67-nuc": "MKI67",
    "BCL2-mem": "BCL2",
    "BCL6-nuc": "BCL6",
    "FoxP3-nuc": "FOXP3",
    "GrB-mem": "GZMB",
    "PD1-mem": "PDCD1",
    "PDL1-mem": "CD274",
    "CD45-mem": "PTPRC",
    "HLA-DR-mem": "HLA-DRA",
    "CD11c-mem": "ITGAX",
    "CD56-mem": "NCAM1",
    "ICOS-mem": "ICOS",
    "CD206-mem": "MRC1",
}

# Cell type mapping for aggregation
CELLTYPE_GROUPS = [
    "CD4 T cell",
    "CD8 T cell",
    "B cell",
    "Macrophage",
    "NK",
    "DC",
    "Endothelial",
    "Fibroblast",
]


def load_imc_data():
    """Load IMC protein expression data."""
    p4_path = PROJECT_DIR / "results" / "adata.immune.celltyped.p4.h5ad"
    if not p4_path.exists():
        logger.error("IMC p4 not found: %s", p4_path)
        return None
    logger.info("Loading IMC: %s", p4_path)
    adata = ad.read_h5ad(p4_path)
    # Use layers['raw'] if X is zeros (known issue)
    if adata.layers and "raw" in adata.layers:
        x_sum = np.abs(adata.X[:100].sum()) if hasattr(adata.X, "sum") else 0
        if x_sum == 0:
            logger.info("  Using layers['raw'] (X is zeros)")
            adata.X = adata.layers["raw"].copy()
    return adata


def load_cosmx_data():
    """Load CosMx RNA expression data."""
    for fname in ["cosmx_rna_annotated.h5ad", "cosmx_rna.h5ad"]:
        path = COSMX_DIR / "results" / fname
        if path.exists():
            logger.info("Loading CosMx: %s", path)
            return ad.read_h5ad(path)
    logger.error("CosMx RNA h5ad not found")
    return None


def compute_celltype_means(adata, markers, ct_col):
    """Compute mean expression per cell type for given markers."""
    genes_in_data = set(adata.var_names)
    present = {k: v for k, v in markers.items() if v in genes_in_data}
    if not present:
        return pd.DataFrame()

    celltypes = sorted(adata.obs[ct_col].unique())
    means = {}
    for ct in celltypes:
        mask = adata.obs[ct_col] == ct
        sub = adata[mask]
        ct_means = {}
        for label, gene in present.items():
            expr = sub[:, gene].X
            if hasattr(expr, "toarray"):
                expr = expr.toarray()
            ct_means[label] = float(np.mean(expr))
        means[ct] = ct_means

    return pd.DataFrame(means).T


def main():
    apply_figure_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    imc = load_imc_data()
    cosmx = load_cosmx_data()

    if imc is None or cosmx is None:
        logger.info("Data not available; creating placeholder")
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(
            0.5, 0.5,
            "RNA-Protein Comparison\n\nRequires both IMC and CosMx data",
            ha="center", va="center", transform=ax.transAxes,
            fontsize=12, style="italic",
        )
        ax.axis("off")
        plt.savefig(FIG_DIR / "supp7_placeholder.pdf", dpi=300, bbox_inches="tight")
        plt.close()
        return

    # Normalize CosMx if needed
    if cosmx.X.max() > 50:
        sc.pp.normalize_total(cosmx, target_sum=1e4)
        sc.pp.log1p(cosmx)

    # Find cell type columns
    imc_ct_col = None
    for col in ["celltype", "celltype_broad", "DLC_code", "cluster"]:
        if col in imc.obs.columns:
            imc_ct_col = col
            break
    cosmx_ct_col = "celltype" if "celltype" in cosmx.obs.columns else "cluster"

    logger.info("IMC celltype col: %s, CosMx celltype col: %s", imc_ct_col, cosmx_ct_col)

    # Normalize cell type labels
    if imc_ct_col:
        imc.obs["ct_norm"] = imc.obs[imc_ct_col].map(normalize_celltype)
    cosmx.obs["ct_norm"] = cosmx.obs[cosmx_ct_col].map(normalize_celltype)

    # Find overlapping markers
    imc_vars = set(imc.var_names)
    cosmx_vars = set(cosmx.var_names)
    overlap = {}
    for imc_name, gene in IMC_TO_GENE.items():
        if imc_name in imc_vars and gene in cosmx_vars:
            overlap[imc_name] = gene
    logger.info("Overlapping markers: %d", len(overlap))

    if len(overlap) < 3:
        logger.warning("Too few overlapping markers (%d); creating placeholder", len(overlap))
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(
            0.5, 0.5,
            f"RNA-Protein Comparison\n\nOnly {len(overlap)} overlapping markers found",
            ha="center", va="center", transform=ax.transAxes,
            fontsize=12, style="italic",
        )
        ax.axis("off")
        plt.savefig(FIG_DIR / "supp7_placeholder.pdf", dpi=300, bbox_inches="tight")
        plt.close()
        return

    # --- Panel 7a: Per-marker scatter plots ---
    # Compute mean expression per cell type
    imc_markers = {k: k for k in overlap}  # IMC var_name -> IMC var_name
    cosmx_markers = dict(overlap.items())  # IMC var_name -> gene symbol

    imc_means = compute_celltype_means(imc, imc_markers, "ct_norm")
    cosmx_means = compute_celltype_means(cosmx, cosmx_markers, "ct_norm")

    # Find common cell types
    common_ct = sorted(set(imc_means.index) & set(cosmx_means.index))
    common_ct = [ct for ct in common_ct if ct not in ("Unknown", "Other", "nan")]
    logger.info("Common cell types: %s", common_ct)

    if len(common_ct) < 2:
        logger.warning("Too few common cell types; using all")
        common_ct = sorted(set(imc_means.index) & set(cosmx_means.index))

    # Per-marker scatter
    n_markers = len(overlap)
    ncols = min(4, n_markers)
    nrows = (n_markers + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows))
    if nrows == 1:
        axes = [axes] if ncols == 1 else list(axes)
    else:
        axes = [ax for row in axes for ax in row]

    correlations = {}
    for idx, (imc_name, gene) in enumerate(sorted(overlap.items())):
        ax = axes[idx]
        x_vals, y_vals = [], []
        for ct in common_ct:
            if ct in imc_means.index and ct in cosmx_means.index:
                x_vals.append(imc_means.loc[ct, imc_name])
                y_vals.append(cosmx_means.loc[ct, imc_name])

        x_vals = np.array(x_vals)
        y_vals = np.array(y_vals)

        ax.scatter(x_vals, y_vals, s=20, alpha=0.7)
        # Add cell type labels
        for i, ct in enumerate(common_ct[:len(x_vals)]):
            ax.annotate(ct, (x_vals[i], y_vals[i]), fontsize=4, alpha=0.7)

        if len(x_vals) > 2 and np.std(x_vals) > 0 and np.std(y_vals) > 0:
            r, p = stats.pearsonr(x_vals, y_vals)
            correlations[gene] = r
            ax.set_title(f"{gene}\nr={r:.2f}, p={p:.2e}", fontsize=6)
        else:
            ax.set_title(gene, fontsize=6)

        ax.set_xlabel("IMC protein", fontsize=5)
        ax.set_ylabel("CosMx RNA", fontsize=5)

    # Hide unused axes
    for idx in range(len(overlap), len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle("IMC protein vs CosMx RNA expression per cell type", fontsize=8)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "supp7a_per_marker_scatter.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("  Saved supp7a_per_marker_scatter.pdf (%d markers)", len(overlap))

    # --- Panel 7b: Summary correlation bar ---
    if correlations:
        corr_df = pd.Series(correlations).sort_values(ascending=True)
        fig, ax = plt.subplots(figsize=(4, max(3, len(corr_df) * 0.3)))
        colors = ["#d73027" if v < 0 else "#4575b4" for v in corr_df.values]
        ax.barh(range(len(corr_df)), corr_df.values, color=colors)
        ax.set_yticks(range(len(corr_df)))
        ax.set_yticklabels(corr_df.index, fontsize=6)
        ax.set_xlabel("Pearson r (IMC protein vs CosMx RNA)")
        ax.set_title("Cross-platform marker correlation")
        ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
        plt.savefig(FIG_DIR / "supp7b_correlation_summary.pdf", dpi=300, bbox_inches="tight")
        plt.close()
        logger.info("  Saved supp7b_correlation_summary.pdf")

        # Save correlation table
        corr_df.to_csv(FIG_DIR / "supp7_correlations.csv")

    logger.info("Supp Fig 7 complete.")


if __name__ == "__main__":
    main()
