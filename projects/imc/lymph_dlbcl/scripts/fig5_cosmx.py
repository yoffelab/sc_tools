#!/usr/bin/env python3
"""Fig 5: CosMx spatial transcriptomics — multicellular interactions.

Panels:
  5a: Schematic (placeholder — BioRender original)
  5b: UMAP of ~350K LME cells colored by cell type (B cells omitted)
  5c: Dotplot of cell populations x key marker genes
  5d: Spatial map (x_slide_mm, y_slide_mm) colored by cell type
  5e: Pathway scores per community (heatmap)
  5f: Receptor-ligand interactions (CD8 T cell context)

Input:  projects/cosmx_1k/lymph_dlbcl/results/cosmx_rna_annotated.h5ad
Output: figures/manuscript/fig5/
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

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
COSMX_DIR = Path(__file__).resolve().parent.parent.parent.parent / "cosmx_1k" / "lymph_dlbcl"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import (
    apply_figure_style,
    build_celltype_palette,
)

FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "fig5"

# Key lineage markers for dotplot (panel 5c)
DOTPLOT_MARKERS = {
    "B cell": ["MS4A1", "PAX5", "BACH2", "CD79A"],
    "CD4 T cell": ["CD4", "IL7R", "FOXP3", "TCF7"],
    "CD8 T cell": ["CD8A", "CD8B", "GZMB", "PRF1"],
    "NK": ["NCAM1", "KLRD1", "NKG7", "GNLY"],
    "Macrophage": ["CD68", "CD163", "MRC1", "CSF1R"],
    "DC": ["ITGAX", "FLT3", "IRF8"],
    "Fibroblast": ["PDPN", "FAP", "COL1A1", "ACTA2"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
}

# Hallmark pathway gene sets (curated subset for CosMx 1000-plex)
PATHWAY_GENES = {
    "TNFa_NFkB": [
        "TNFRSF1A", "TNFRSF1B", "NFKB1", "NFKB2", "NFKBIA", "RELA", "RELB",
        "TRAF1", "TRAF2", "BIRC2", "BIRC3", "TNFAIP3", "BCL2A1", "ICAM1",
    ],
    "JAK_STAT": [
        "JAK1", "JAK2", "JAK3", "STAT1", "STAT3", "STAT5A", "STAT5B",
        "SOCS1", "SOCS3", "IRF1", "CISH",
    ],
    "TGFb": [
        "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3",
        "SMAD4", "SMAD7", "ACVR1", "ACVR1B",
    ],
    "VEGF": [
        "VEGFA", "VEGFB", "VEGFC", "FLT1", "KDR", "FLT4", "NRP1", "NRP2",
        "PGF", "HIF1A",
    ],
    "Hypoxia": [
        "HIF1A", "VEGFA", "SLC2A1", "LDHA", "PGK1", "ENO1", "ALDOA",
        "PKM", "CA9", "BNIP3",
    ],
    "Interferon_gamma": [
        "STAT1", "IRF1", "GBP1", "GBP2", "CXCL9", "CXCL10", "CXCL11",
        "IDO1", "TAP1", "PSMB9",
    ],
}

# Receptor-ligand pairs for CD8 T cell context (panel 5f)
CD8_RECEPTOR_LIGAND = {
    "PD1-PDL1": {"receptor": "PDCD1", "ligand": "CD274"},
    "CTLA4-CD80": {"receptor": "CTLA4", "ligand": "CD80"},
    "LAG3-FGL1": {"receptor": "LAG3", "ligand": "FGL1"},
    "TIGIT-PVR": {"receptor": "TIGIT", "ligand": "PVR"},
    "CD28-CD86": {"receptor": "CD28", "ligand": "CD86"},
    "ICOS-ICOSL": {"receptor": "ICOS", "ligand": "ICOSLG"},
    "4-1BB-4-1BBL": {"receptor": "TNFRSF9", "ligand": "TNFSF9"},
    "OX40-OX40L": {"receptor": "TNFRSF4", "ligand": "TNFSF4"},
}


def load_cosmx_rna():
    """Load annotated CosMx RNA h5ad."""
    path = COSMX_DIR / "results" / "cosmx_rna_annotated.h5ad"
    if not path.exists():
        # Fall back to unannotated
        path = COSMX_DIR / "results" / "cosmx_rna.h5ad"
    if not path.exists():
        logger.error("CosMx RNA h5ad not found at %s", path)
        sys.exit(1)
    logger.info("Loading %s", path)
    adata = ad.read_h5ad(path)
    logger.info("Shape: %s", adata.shape)
    return adata


def panel_5b_umap(adata):
    """UMAP of ~350K LME cells colored by cell type (B cells omitted)."""
    logger.info("Panel 5b: UMAP by cell type (excluding B cells)")

    # Determine cell type column
    ct_col = "celltype" if "celltype" in adata.obs.columns else "cluster"

    # Exclude B cells
    if ct_col == "celltype":
        mask = ~adata.obs[ct_col].str.contains("B cell|Plasma", case=False, na=False)
    else:
        mask = pd.Series(True, index=adata.obs.index)

    sub = adata[mask].copy()
    logger.info("Cells after B cell exclusion: %d", sub.n_obs)

    # Build palette
    categories = sub.obs[ct_col].cat.categories if hasattr(sub.obs[ct_col], "cat") else sub.obs[ct_col].unique()
    palette = build_celltype_palette(categories)

    # Check for UMAP coordinates
    if "X_umap" not in sub.obsm:
        logger.info("Computing UMAP (no precomputed embedding found)")
        sc.pp.pca(sub, n_comps=30)
        sc.pp.neighbors(sub, n_pcs=30)
        sc.tl.umap(sub)

    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(
        sub,
        color=ct_col,
        palette=palette,
        ax=ax,
        show=False,
        frameon=False,
        size=0.5,
        legend_loc="right margin",
        legend_fontsize=6,
    )
    ax.set_title(f"CosMx 1000-plex RNA ({sub.n_obs:,} non-B cells)", fontsize=8)
    plt.savefig(FIG_DIR / "fig5b_umap_celltype.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("  Saved fig5b_umap_celltype.pdf")


def panel_5c_dotplot(adata):
    """Dotplot of cell populations x key marker genes."""
    logger.info("Panel 5c: Dotplot markers x cell types")

    ct_col = "celltype" if "celltype" in adata.obs.columns else "cluster"
    genes_in_data = set(adata.var_names)

    # Collect markers present in data
    marker_list = []
    for _celltype, genes in DOTPLOT_MARKERS.items():
        present = [g for g in genes if g in genes_in_data]
        marker_list.extend(present)

    if not marker_list:
        logger.warning("No dotplot markers found in data")
        return

    marker_list = list(dict.fromkeys(marker_list))  # deduplicate preserving order
    logger.info("  %d markers available for dotplot", len(marker_list))

    # Group labels for dotplot
    var_group_labels = []
    var_group_positions = []
    pos = 0
    for celltype, genes in DOTPLOT_MARKERS.items():
        present = [g for g in genes if g in genes_in_data]
        if present:
            var_group_labels.append(celltype)
            var_group_positions.append((pos, pos + len(present) - 1))
            pos += len(present)

    dp = sc.pl.dotplot(
        adata,
        var_names=marker_list,
        groupby=ct_col,
        var_group_labels=var_group_labels,
        var_group_positions=var_group_positions,
        standard_scale="var",
        show=False,
        return_fig=True,
    )
    dp.savefig(FIG_DIR / "fig5c_dotplot.pdf")
    plt.close("all")
    logger.info("  Saved fig5c_dotplot.pdf")


def panel_5d_spatial(adata):
    """Spatial map colored by cell type for one representative sample."""
    logger.info("Panel 5d: Spatial map")

    ct_col = "celltype" if "celltype" in adata.obs.columns else "cluster"

    # Check for spatial coordinates
    if "x_slide_mm" not in adata.obs.columns:
        logger.warning("No x_slide_mm column; skipping spatial panel")
        return

    # Use CTMA121 (smaller, cleaner for visualization)
    sample_col = "sample" if "sample" in adata.obs.columns else "Run_Tissue_name"
    samples = adata.obs[sample_col].unique()
    target = [s for s in samples if "CTMA121" in str(s)]
    if target:
        sub = adata[adata.obs[sample_col] == target[0]].copy()
        logger.info("  Using sample %s (%d cells)", target[0], sub.n_obs)
    else:
        sub = adata.copy()
        logger.info("  Using all cells (%d)", sub.n_obs)

    # Exclude B cells for consistency with 5b
    if ct_col == "celltype":
        mask = ~sub.obs[ct_col].str.contains("B cell|Plasma", case=False, na=False)
        sub = sub[mask].copy()

    # Subsample if too many for PDF
    if sub.n_obs > 200000:
        idx = np.random.RandomState(42).choice(sub.n_obs, 200000, replace=False)
        sub = sub[idx].copy()

    categories = sub.obs[ct_col].cat.categories if hasattr(sub.obs[ct_col], "cat") else sub.obs[ct_col].unique()
    palette = build_celltype_palette(categories)

    fig, ax = plt.subplots(figsize=(8, 6))
    for ct in categories:
        mask = sub.obs[ct_col] == ct
        color = palette.get(str(ct), "#999999")
        ax.scatter(
            sub.obs.loc[mask, "x_slide_mm"],
            sub.obs.loc[mask, "y_slide_mm"],
            c=[color],
            s=0.1,
            alpha=0.5,
            label=ct,
            rasterized=True,
        )
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(f"CosMx spatial map ({sub.n_obs:,} cells)")
    ax.set_aspect("equal")
    ax.legend(
        markerscale=10,
        fontsize=5,
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        frameon=False,
    )
    plt.savefig(FIG_DIR / "fig5d_spatial_map.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("  Saved fig5d_spatial_map.pdf")


def panel_5e_pathways(adata):
    """Pathway scores per cluster/community (heatmap)."""
    logger.info("Panel 5e: Pathway scores")

    ct_col = "celltype" if "celltype" in adata.obs.columns else "cluster"
    genes_in_data = set(adata.var_names)

    # Score each pathway per cluster
    clusters = sorted(adata.obs[ct_col].unique())
    pathway_scores = {}

    for pathway, genes in PATHWAY_GENES.items():
        present = [g for g in genes if g in genes_in_data]
        if len(present) < 3:
            logger.warning("  %s: only %d genes present, skipping", pathway, len(present))
            continue
        logger.info("  %s: %d/%d genes", pathway, len(present), len(genes))

        cluster_means = {}
        for c in clusters:
            mask = adata.obs[ct_col] == c
            subset = adata[mask, present].X
            if hasattr(subset, "toarray"):
                subset = subset.toarray()
            cluster_means[c] = np.mean(subset)

        pathway_scores[pathway] = cluster_means

    if not pathway_scores:
        logger.warning("No pathways scored; skipping panel 5e")
        return

    score_df = pd.DataFrame(pathway_scores, index=clusters)
    # Z-score per pathway
    score_z = (score_df - score_df.mean()) / score_df.std().replace(0, 1)

    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(score_z.T.values, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
    ax.set_xticks(range(len(clusters)))
    ax.set_xticklabels(clusters, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(len(pathway_scores)))
    ax.set_yticklabels(list(pathway_scores.keys()), fontsize=6)
    ax.set_xlabel("Cell type / cluster")
    ax.set_title("Pathway activity (z-scored)")
    plt.colorbar(im, ax=ax, label="Z-score", shrink=0.8)
    plt.savefig(FIG_DIR / "fig5e_pathway_scores.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("  Saved fig5e_pathway_scores.pdf")


def panel_5f_receptor_ligand(adata):
    """Receptor-ligand interactions in CD8 T cell context."""
    logger.info("Panel 5f: Receptor-ligand (CD8 context)")

    ct_col = "celltype" if "celltype" in adata.obs.columns else "cluster"
    genes_in_data = set(adata.var_names)

    # Find CD8 T cells
    if ct_col == "celltype":
        cd8_mask = adata.obs[ct_col].str.contains("CD8", case=False, na=False)
    else:
        # Use all cells if no celltype annotation
        cd8_mask = pd.Series(True, index=adata.obs.index)

    other_mask = ~cd8_mask
    cd8_cells = adata[cd8_mask]
    other_cells = adata[other_mask]

    if cd8_cells.n_obs == 0:
        logger.warning("No CD8 T cells found; skipping panel 5f")
        return

    logger.info("  CD8 cells: %d, Other: %d", cd8_cells.n_obs, other_cells.n_obs)

    # Score receptor-ligand pairs
    pairs_data = []
    for pair_name, genes in CD8_RECEPTOR_LIGAND.items():
        rec = genes["receptor"]
        lig = genes["ligand"]
        if rec not in genes_in_data or lig not in genes_in_data:
            continue

        # Receptor expression in CD8 T cells
        rec_expr = cd8_cells[:, rec].X
        if hasattr(rec_expr, "toarray"):
            rec_expr = rec_expr.toarray()
        rec_mean = float(np.mean(rec_expr))
        rec_pct = float(np.mean(rec_expr > 0) * 100)

        # Ligand expression in other cells (by cell type)
        for ct in sorted(other_cells.obs[ct_col].unique()):
            ct_mask = other_cells.obs[ct_col] == ct
            lig_expr = other_cells[ct_mask, lig].X
            if hasattr(lig_expr, "toarray"):
                lig_expr = lig_expr.toarray()
            lig_mean = float(np.mean(lig_expr))
            lig_pct = float(np.mean(lig_expr > 0) * 100)

            pairs_data.append({
                "pair": pair_name,
                "receptor": rec,
                "ligand": lig,
                "cd8_receptor_mean": rec_mean,
                "cd8_receptor_pct": rec_pct,
                "source_celltype": ct,
                "ligand_mean": lig_mean,
                "ligand_pct": lig_pct,
                "interaction_score": rec_mean * lig_mean,
            })

    if not pairs_data:
        logger.warning("No receptor-ligand pairs found in data")
        return

    rl_df = pd.DataFrame(pairs_data)
    rl_df.to_csv(FIG_DIR / "fig5f_receptor_ligand_scores.csv", index=False)

    # Heatmap: pairs x source cell types (interaction score)
    pivot = rl_df.pivot_table(
        index="pair", columns="source_celltype", values="interaction_score", aggfunc="mean"
    )

    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(pivot.values, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=6)
    ax.set_title("CD8 T cell receptor-ligand interactions")
    ax.set_xlabel("Ligand source cell type")
    ax.set_ylabel("Receptor-ligand pair")
    plt.colorbar(im, ax=ax, label="Interaction score", shrink=0.8)
    plt.savefig(FIG_DIR / "fig5f_receptor_ligand.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("  Saved fig5f_receptor_ligand.pdf")


def main():
    apply_figure_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    adata = load_cosmx_rna()

    # Ensure data is normalized (log1p)
    if adata.X.max() > 50:
        logger.info("Normalizing data (log1p)")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    panel_5b_umap(adata)
    panel_5c_dotplot(adata)
    panel_5d_spatial(adata)
    panel_5e_pathways(adata)
    panel_5f_receptor_ligand(adata)

    logger.info("Fig 5 complete. Panels saved to %s", FIG_DIR)


if __name__ == "__main__":
    main()
