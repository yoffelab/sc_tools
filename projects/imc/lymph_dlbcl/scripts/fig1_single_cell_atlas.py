#!/usr/bin/env python3
"""
Figure 1: Single-cell phenotypic atlas of DLBCL tumor microenvironment.

Manuscript caption (adapted):
  a) UMAP of immune (T2) and stromal (S2) panels colored by cell type
  b) Marker x cell type heatmap (z-scored per marker) — rows=markers, cols=cell types
  c) Cell type prevalence as fraction of non-B cells
  d) B cell marker expression prevalence
  e) COO distribution across cell types (from DLC380_clinical.tsv)

Insight: The TME of DLBCL is composed of 12 major non-B cell populations
(NESC, DC, Mac, Mono, BEC, LEC, CD4 T, CD8 T, Treg, NK, B cell, Other)
identifiable by IMC protein markers across two complementary panels.

Usage:
    python scripts/fig1_single_cell_atlas.py

Input:
    results/adata.immune.celltyped.p4.h5ad (or p2/p1 fallback)
    results/adata.stromal.celltyped.p4.h5ad

Output:
    figures/manuscript/fig1/
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
import scanpy as sc
import seaborn as sns
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_config import (
    COO_COLORS,
    apply_figure_style,
    build_celltype_palette,
    filter_protein_markers,
    is_bcell_label,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
FIG_DIR = PROJECT_DIR / "figures" / "manuscript" / "fig1"


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def load_panel(panel: str) -> ad.AnnData | None:
    """Load best available panel checkpoint."""
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


def resolve_markers(adata: ad.AnnData, markers: list[str]) -> list[str]:
    """Resolve marker names with case-insensitive fallback."""
    available = []
    var_lower = {v.lower(): v for v in adata.var_names}
    for m in markers:
        if m in adata.var_names:
            available.append(m)
        elif m.lower() in var_lower:
            available.append(var_lower[m.lower()])
    return available


def fig1a_umap(adata: ad.AnnData, panel: str, celltype_col: str):
    """UMAP colored by cell type."""
    logger.info(f"  Fig 1a: UMAP for {panel}")

    # Subsample a copy for UMAP (don't modify original adata)
    if "X_umap" not in adata.obsm:
        if adata.n_obs > 200_000:
            n_sub = 100_000
            logger.info(f"    Subsampling {adata.n_obs:,} -> {n_sub:,} for UMAP")
            rng = np.random.default_rng(42)
            idx = rng.choice(adata.n_obs, n_sub, replace=False)
            adata_sub = adata[sorted(idx)].copy()
        else:
            adata_sub = adata.copy()

        # Use layers['raw'] if X is empty (p4 may have zeroed X)
        if adata_sub.layers and "raw" in adata_sub.layers:
            x_check = adata_sub.X[:min(100, adata_sub.n_obs)]
            if hasattr(x_check, "toarray"):
                x_check = x_check.toarray()
            if np.abs(x_check).max() < 1e-10:
                logger.info("    X is empty; copying layers['raw'] to X for UMAP")
                adata_sub.X = adata_sub.layers["raw"].copy()

        if "X_pca" not in adata_sub.obsm:
            sc.pp.pca(adata_sub, n_comps=30)
        sc.pp.neighbors(adata_sub, n_neighbors=15, use_rep="X_pca")
        sc.tl.umap(adata_sub)
    else:
        adata_sub = adata

    # Build color palette — distinct color per category (no fuzzy collapse)
    categories = adata_sub.obs[celltype_col].astype(str).unique()
    palette = build_celltype_palette(categories, normalize=False)
    logger.info(f"    Palette: {len(palette)} distinct categories")

    # Ensure celltype column is string (not numeric/categorical mismatch)
    adata_sub.obs[celltype_col] = adata_sub.obs[celltype_col].astype(str)

    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(
        adata_sub,
        color=celltype_col,
        ax=ax,
        show=False,
        frameon=False,
        title=f"{panel.title()} Panel",
        legend_loc="right margin",
        size=0.5,
        palette=palette,
    )
    # Rasterize scatter points for smaller PDF
    for child in ax.get_children():
        if hasattr(child, "set_rasterized"):
            child.set_rasterized(True)
    out = FIG_DIR / f"fig1a_umap_{panel}.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig1b_heatmap(adatas: dict, celltype_cols: dict):
    """Marker x cell type heatmap (z-scored per marker).

    Insight: Each cell type is defined by a distinct marker expression profile.
    Z-scoring reveals relative enrichment patterns across cell types.
    """
    logger.info("  Fig 1b: Marker heatmap (z-scored)")

    all_mean = []
    for panel, adata in adatas.items():
        col = celltype_cols[panel]
        n_types = adata.obs[col].nunique()

        # Skip panels with <5 unique cell types (broken annotation)
        if n_types < 5:
            logger.info(f"    Skipping {panel} for heatmap ({n_types} types — insufficient)")
            continue

        # Filter to protein markers only (exclude morphological features)
        protein_markers = filter_protein_markers(adata.var_names)
        protein_idx = [list(adata.var_names).index(m) for m in protein_markers]
        logger.info(f"    {panel}: {len(protein_markers)}/{len(adata.var_names)} protein markers, {n_types} cell types")

        for group in adata.obs[col].unique():
            group_str = str(group)
            # Only include fine-grained subtypes (T*, M* prefixed) for heatmap
            # Skip broad categories (bcell, stroma, other, Unknown) that dilute z-scores
            if not (group_str.startswith("T") or group_str.startswith("M")):
                logger.info(f"    Skipping '{group_str}' from heatmap (broad category)")
                continue
            mask = adata.obs[col] == group
            # Use layers['raw'] if X is empty (p4 may have zeroed X)
            data_matrix = adata.layers.get("raw", adata.X) if adata.layers else adata.X
            if hasattr(data_matrix, "toarray"):
                raw_vals = np.array(data_matrix[mask].toarray())
            else:
                raw_vals = np.array(data_matrix[mask])
            # Data is already Seurat-normalized (centered/scaled); take mean directly
            mean_vals = raw_vals.mean(axis=0).flatten()
            for j, marker_idx in enumerate(protein_idx):
                all_mean.append({
                    "celltype": group_str,  # keep raw labels
                    "marker": protein_markers[j],
                    "panel": panel,
                    "mean_expr": mean_vals[marker_idx],
                })

    df = pd.DataFrame(all_mean)
    pivot = df.pivot_table(index="marker", columns="celltype", values="mean_expr", aggfunc="mean")

    # Z-score per marker (row) — highlights which cell types are enriched for each marker
    z_scored = pivot.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=1)

    # Drop rows/cols that are all NaN; fill remaining NaN with 0
    z_scored = z_scored.dropna(how="all", axis=0).dropna(how="all", axis=1)
    z_scored = z_scored.fillna(0)

    # Log marker variation stats for debugging
    marker_std = pivot.std(axis=1)
    logger.info(f"    Marker std range: {marker_std.min():.4f} - {marker_std.max():.4f}")
    logger.info(f"    Markers with std > 0.001: {(marker_std > 0.001).sum()}/{len(marker_std)}")
    logger.info(f"    z_scored shape: {z_scored.shape}")

    fig_h = max(8, len(z_scored) * 0.3)
    fig_w = max(8, len(z_scored.columns) * 0.5)

    # Build column color bar for cell types (distinct colors, no collapse)
    col_palette = build_celltype_palette(z_scored.columns, normalize=False)
    col_colors = pd.Series(
        [col_palette.get(ct, "#999999") for ct in z_scored.columns],
        index=z_scored.columns,
        name="Cell Type",
    )

    g = sns.clustermap(
        z_scored,
        cmap="RdBu_r",
        center=0,
        vmin=-2.5,
        vmax=2.5,
        figsize=(fig_w, fig_h),
        dendrogram_ratio=(0.08, 0.12),
        cbar_kws={"label": "z-score", "shrink": 0.6},
        xticklabels=True,
        yticklabels=True,
        row_cluster=True,
        col_cluster=True,
        col_colors=col_colors,
    )
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=6)
    g.ax_heatmap.set_xlabel("Cell Type")
    g.ax_heatmap.set_ylabel("Marker")

    out = FIG_DIR / "fig1b_heatmap.pdf"
    g.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig1c_prevalence(adatas: dict, celltype_cols: dict):
    """Cell type prevalence as fraction of non-B cells.

    Insight: Non-B cell immune and stromal populations constitute the TME.
    Prevalence shows relative abundance of each population.
    """
    logger.info("  Fig 1c: Cell type prevalence (non-B cell fraction)")

    # Use immune panel only for prevalence (stromal panel has broken annotations)
    # If immune panel has fine-grained labels, use them directly
    all_counts = {}
    for panel, adata in adatas.items():
        col = celltype_cols[panel]
        n_unique = adata.obs[col].nunique()
        # Skip panels with <5 unique types (indicates broken annotation)
        if n_unique < 5:
            logger.info(f"    Skipping {panel} for prevalence ({n_unique} types — insufficient)")
            continue
        ct_counts = adata.obs[col].value_counts()
        for ct, count in ct_counts.items():
            ct_str = str(ct)
            all_counts[ct_str] = all_counts.get(ct_str, 0) + count
        logger.info(f"    Using {panel} for prevalence ({n_unique} types)")

    # Exclude B cells AND Unknown for non-B fraction
    exclude_keys = [
        k for k in all_counts
        if is_bcell_label(k) or k.lower() in ("unknown", "other")
    ]
    non_b_total = sum(c for k, c in all_counts.items() if k not in exclude_keys)
    logger.info(f"    Excluded from prevalence: {exclude_keys}")

    prevalence = {}
    for ct, count in all_counts.items():
        if ct not in exclude_keys:
            prevalence[ct] = count / non_b_total if non_b_total > 0 else 0

    prev_df = pd.DataFrame.from_dict(prevalence, orient="index", columns=["fraction"])
    prev_df = prev_df.sort_values("fraction", ascending=True)

    # Build palette with distinct colors per subtype
    palette = build_celltype_palette(prev_df.index, normalize=False)
    colors = [palette.get(ct, "#999999") for ct in prev_df.index]

    fig, ax = plt.subplots(figsize=(5, max(4, len(prev_df) * 0.3)))
    ax.barh(range(len(prev_df)), prev_df["fraction"], color=colors)
    ax.set_yticks(range(len(prev_df)))
    ax.set_yticklabels(prev_df.index)
    ax.set_xlabel("Fraction of Non-B Cells")
    ax.set_title("Cell Type Prevalence (Non-B, Non-Unknown)")

    for i, (_ct, row) in enumerate(prev_df.iterrows()):
        ax.text(row["fraction"] + 0.005, i, f"{row['fraction']:.1%}", va="center", fontsize=6)

    out = FIG_DIR / "fig1c_prevalence.pdf"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"    Saved: {out}")


def fig1d_bcell_markers(adatas: dict, celltype_cols: dict):
    """B cell marker expression prevalence.

    Insight: B cell markers (CD20, PAX5, Ki67, BCL2, BCL6, MYC) define
    the tumor cell population and its molecular subtypes.
    """
    logger.info("  Fig 1d: B cell marker prevalence")

    b_markers = ["CD20", "PAX5", "Ki67", "BCL2", "BCL6", "MYC", "MUM1", "CD10"]

    for panel, adata in adatas.items():
        col = celltype_cols[panel]
        available = resolve_markers(adata, b_markers)
        if len(available) < 2:
            continue

        # Get B cell subset
        b_mask = adata.obs[col].astype(str).str.lower().str.contains("b cell|b_cell|bcell")
        if b_mask.sum() < 10:
            continue

        b_adata = adata[b_mask]
        # Use layers['raw'] if X is empty (p4 may have zeroed X)
        b_data = b_adata.layers.get("raw", b_adata.X) if b_adata.layers else b_adata.X
        if hasattr(b_data, "toarray"):
            expr = pd.DataFrame(b_data.toarray(), columns=b_adata.var_names)
        else:
            expr = pd.DataFrame(np.array(b_data), columns=b_adata.var_names)

        # Fraction positive (> median of all cells)
        full_data = adata.layers.get("raw", adata.X) if adata.layers else adata.X
        full_medians = {}
        for m in available:
            idx = list(adata.var_names).index(m)
            if hasattr(full_data, "toarray"):
                full_medians[m] = np.median(full_data[:, idx].toarray())
            else:
                full_medians[m] = np.median(full_data[:, idx])

        prevalence = {}
        for m in available:
            prevalence[m] = (expr[m] > full_medians[m]).mean()

        prev_df = pd.Series(prevalence).sort_values(ascending=True)

        fig, ax = plt.subplots(figsize=(5, max(3, len(prev_df) * 0.35)))
        ax.barh(range(len(prev_df)), prev_df.values, color="#0072B2")
        ax.set_yticks(range(len(prev_df)))
        ax.set_yticklabels(prev_df.index)
        ax.set_xlabel("Fraction Positive (B cells)")
        ax.set_title(f"B Cell Marker Prevalence ({panel})")
        ax.set_xlim(0, 1)

        out = FIG_DIR / f"fig1d_bcell_markers_{panel}.pdf"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"    Saved: {out}")
        break  # only need one panel for B cell markers


def fig1e_coo(adatas: dict, celltype_cols: dict):
    """COO distribution if available.

    Insight: Cell of origin (GCB vs ABC) subtypes show distinct
    microenvironment compositions.
    """
    logger.info("  Fig 1e: COO comparison")

    for panel, adata in adatas.items():
        if "COO" not in adata.obs.columns:
            continue

        col = celltype_cols[panel]
        data = adata.obs.dropna(subset=["COO"])
        if len(data) < 10:
            continue

        ct = pd.crosstab(data[col], data["COO"])
        ct_norm = ct.div(ct.sum(axis=1), axis=0)

        # Use COO colors
        coo_cols = [COO_COLORS.get(c, "#999999") for c in ct_norm.columns]

        fig, ax = plt.subplots(figsize=(max(6, ct_norm.shape[0] * 0.5), 5))
        ct_norm.plot(kind="bar", stacked=True, ax=ax, color=coo_cols, width=0.8)
        ax.set_ylabel("Proportion")
        ax.set_xlabel("")
        ax.set_title("COO Distribution by Cell Type")
        ax.legend(title="COO", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.xticks(rotation=45, ha="right")

        out = FIG_DIR / f"fig1e_coo_{panel}.pdf"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"    Saved: {out}")
        break


def main():
    apply_figure_style()
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

    for panel, adata in adatas.items():
        fig1a_umap(adata, panel, celltype_cols[panel])

    fig1b_heatmap(adatas, celltype_cols)
    fig1c_prevalence(adatas, celltype_cols)
    fig1d_bcell_markers(adatas, celltype_cols)
    fig1e_coo(adatas, celltype_cols)

    logger.info("Figure 1 complete.")


if __name__ == "__main__":
    main()
