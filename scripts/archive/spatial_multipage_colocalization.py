"""
Spatial multipage colocalization: multipage PDFs and proportion boxplots for
macrophage (proliferation x M2) and neutrophil–cytotoxic T-cell (SLC16A3+ neut x cytotoxic T-cell).

Reads scores from obsm['signature_score_z']. Outputs:
- figures/manuscript/macrophage_localization/spatial_multipage.pdf + proportion boxplot
- figures/manuscript/neutrophil_cytotoxic_tcell_localization/spatial_multipage.pdf + proportion boxplot

Run from project root: python scripts/spatial_multipage_colocalization.py
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

import sc_tools.pl.spatial as st_spatial
import sc_tools.pl.statistical as st_statistical
import sc_tools.tl.colocalization as st_coloc
from sc_tools.utils.signatures import get_signature_df

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# Signature columns (obsm full-path names)
PROLIFERATION_COL = "Tumor_Cells/Proliferative_Tumor"
MACROPHAGE_COL = "Liron/M2-macrophages"
NEUTROPHIL_COL = "Liron/SLC16A3+ Neutrophil"
CYTOTOXIC_COL = "Liron/T cell cytotoxicity"

TUMOR_TYPE_MAPPING = {
    "Solid": "Solid",
    "Non-Solid": "Non-Solid",
    "Normal": "Normal",
    "Solid Blood Vessel": "Solid",
    "Solid Bronchus": "Solid",
    "Solid Scar Tissue": "Solid",
    "Non-Solid Blood Vessel": "Non-Solid",
    "Non-Solid Bronchus": "Non-Solid",
    "Normal Blood Vessel": "Normal",
    "Normal Bronchus": "Normal",
    "Solid TLS": "Solid",
    "TLS Solid": "Solid",
    "Non-Solid TLS": "Non-Solid",
    "TLS Non-Solid": "Non-Solid",
    "TLS Normal": "Normal",
}

SOLIDITY_COLORS = {
    "Normal": "#66c2a5",
    "Non-Solid": "#fc8d62",
    "Solid": "#8da0cb",
}


def _prepare_solidity(adata):
    """Add solidity from pathologist_annotation."""
    if "solidity_type" in adata.obs.columns:
        adata.obs["solidity"] = adata.obs["solidity_type"]
    else:
        adata.obs["solidity"] = (
            adata.obs["pathologist_annotation"]
            .astype(str)
            .replace(TUMOR_TYPE_MAPPING)
        )
    adata.obs["solidity"] = pd.Categorical(
        adata.obs["solidity"],
        categories=["Normal", "Non-Solid", "Solid"],
    )
    return adata


def _proportion_nonzero_df(adata, truncated_series):
    """Proportion of spots with truncated_similarity > 0 per (library_id, solidity)."""
    df = adata.obs[["library_id", "solidity"]].copy()
    df["trunc"] = truncated_series.reindex(adata.obs_names).values
    results = []
    for (lib_id, solidity), grp in df.groupby(["library_id", "solidity"]):
        n = len(grp)
        n_pos = (grp["trunc"] > 0).sum()
        results.append({
            "library_id": lib_id,
            "solidity": solidity,
            "proportion_nonzero": n_pos / n if n else 0,
            "n_spots": n,
            "n_nonzero": int(n_pos),
        })
    return pd.DataFrame(results)


def _plot_proportion_boxplot(prop_df, output_dir, title, basename="proportion_nonzero_boxplot"):
    """Boxplot of proportion nonzero by solidity with Mann-Whitney + FDR."""
    solidity_order = ["Normal", "Non-Solid", "Solid"]
    prop_df = prop_df[prop_df["solidity"].isin(solidity_order)].copy()
    prop_df["solidity"] = pd.Categorical(prop_df["solidity"], categories=solidity_order)

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.boxplot(
        data=prop_df,
        x="solidity",
        y="proportion_nonzero",
        ax=ax,
        palette=[SOLIDITY_COLORS.get(s, "#808080") for s in solidity_order],
        width=0.6,
    )
    sns.stripplot(
        data=prop_df,
        x="solidity",
        y="proportion_nonzero",
        ax=ax,
        color="black",
        size=5,
        alpha=0.6,
        jitter=True,
    )
    ax.set_xlabel("Solidity", fontsize=12, fontweight="bold")
    ax.set_ylabel("Proportion of Non-Zero Truncated Similarity", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=15)

    solidity_groups = {
        s: prop_df[prop_df["solidity"] == s]["proportion_nonzero"].values
        for s in solidity_order
    }
    comparisons = []
    p_values = []
    for i, g1 in enumerate(solidity_order):
        for j, g2 in enumerate(solidity_order):
            if i >= j:
                continue
            v1, v2 = solidity_groups[g1], solidity_groups[g2]
            if len(v1) < 2 or len(v2) < 2:
                continue
            _, p = mannwhitneyu(v1, v2, alternative="two-sided")
            comparisons.append({"group1": g1, "group2": g2, "group1_idx": i, "group2_idx": j})
            p_values.append(p)

    if p_values:
        _, p_adj, _, _ = multipletests(p_values, method="fdr_bh", alpha=0.05)
        y_max = prop_df["proportion_nonzero"].max()
        y_range = max(prop_df["proportion_nonzero"].max() - prop_df["proportion_nonzero"].min(), 1e-6)
        bar_height = 0.02 * y_range
        bar_y = y_max + 0.05 * y_range
        for (comp, p_adj_val) in zip(comparisons, p_adj, strict=True):
            if p_adj_val >= 0.05:
                continue
            i1, i2 = comp["group1_idx"], comp["group2_idx"]
            ax.plot([i1, i2], [bar_y, bar_y], "k-", linewidth=1.5)
            ax.plot([i1, i1], [bar_y - bar_height / 2, bar_y], "k-", linewidth=1.5)
            ax.plot([i2, i2], [bar_y - bar_height / 2, bar_y], "k-", linewidth=1.5)
            ax.text(
                (i1 + i2) / 2,
                bar_y + bar_height,
                st_statistical.get_asterisk(p_adj_val),
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )
            ax.text(
                (i1 + i2) / 2,
                bar_y + 2.5 * bar_height,
                f"p={p_adj_val:.3e}",
                ha="center",
                va="bottom",
                fontsize=8,
            )
            bar_y += 3 * bar_height
        ax.set_ylim(
            bottom=prop_df["proportion_nonzero"].min() - 0.05 * y_range,
            top=bar_y + 2 * bar_height,
        )

    ax.grid(True, alpha=0.3, axis="y")
    sns.despine(ax=ax)
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, f"{basename}.pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join(output_dir, f"{basename}.png"), bbox_inches="tight", dpi=300)
    plt.close(fig)

    if p_values:
        stats_df = pd.DataFrame({
            "Comparison": [f"{c['group1']} vs {c['group2']}" for c in comparisons],
            "p_value": p_values,
            "p_adjusted": p_adj,
        })
        stats_df.to_csv(os.path.join(output_dir, f"{basename}_stats.csv"), index=False)


def _build_panels(adata, solidity_col, score1_series, score2_series, trunc_series):
    """Build 6-panel spec for multipage_spatial_pdf."""
    return [
        {"type": "he"},
        {"type": "categorical", "obs_col": "pathologist_annotation", "title": "Pathologist Annotation"},
        {"type": "categorical", "obs_col": solidity_col, "title": "Solidity", "palette": SOLIDITY_COLORS},
        {
            "type": "continuous",
            "title": "Score 1",
            "values": score1_series,
            "cmap": "coolwarm",
            "vmin": -2,
            "vmax": 2,
        },
        {
            "type": "continuous",
            "title": "Score 2",
            "values": score2_series,
            "cmap": "coolwarm",
            "vmin": -2,
            "vmax": 2,
        },
        {
            "type": "continuous",
            "title": "Truncated Similarity",
            "values": trunc_series,
            "cmap": "YlOrRd",
            "vmin": 0,
            "vmax": 1,
        },
    ]


def main():
    _scored_new = Path("results/adata.scored.h5ad")
    _scored_old = Path("results/adata.normalized.scored.p35.h5ad")
    adata_path = _scored_new if _scored_new.exists() else _scored_old
    if not adata_path.exists():
        raise FileNotFoundError(f"Not found: {adata_path}. Run score_gene_signatures.py first.")

    print("Loading adata and signature scores...")
    adata = sc.read_h5ad(adata_path)
    sig_df = get_signature_df(adata, use_z=True)
    if sig_df is None or sig_df.empty:
        raise ValueError("No signature scores in obsm['signature_score_z']. Run score_gene_signatures.py.")

    for col in [PROLIFERATION_COL, MACROPHAGE_COL, NEUTROPHIL_COL, CYTOTOXIC_COL]:
        if col not in sig_df.columns:
            print(f"Warning: signature column '{col}' not found. Available: {list(sig_df.columns)[:15]}...")

    adata = _prepare_solidity(adata)

    library_id_col = "library_id"
    if library_id_col not in adata.obs.columns:
        raise KeyError(f"Missing {library_id_col} in adata.obs")

    # --- Macrophage: Proliferation x M2 ---
    out_mac = "figures/manuscript/macrophage_localization"
    if PROLIFERATION_COL in sig_df.columns and MACROPHAGE_COL in sig_df.columns:
        prolif = sig_df[PROLIFERATION_COL]
        macro = sig_df[MACROPHAGE_COL]
        trunc_mac = pd.Series(
            st_coloc.truncated_similarity(prolif.values, macro.values),
            index=sig_df.index,
        )
        panels_mac = _build_panels(
            adata,
            "solidity",
            prolif,
            macro,
            trunc_mac,
        )
        panels_mac[3]["title"] = "Tumor Proliferation"
        panels_mac[4]["title"] = "M2 Macrophage"

        print("Writing macrophage multipage PDF...")
        st_spatial.multipage_spatial_pdf(
            adata,
            library_id_col,
            panels_mac,
            os.path.join(out_mac, "spatial_multipage.pdf"),
        )
        prop_mac = _proportion_nonzero_df(adata, trunc_mac)
        _plot_proportion_boxplot(
            prop_mac,
            out_mac,
            "Proportion of Non-Zero Truncated Similarity (Proliferation x M2)\nby Solidity",
            "proportion_nonzero_boxplot_macrophage",
        )
    else:
        print("Skipping macrophage analysis (missing signature columns).")

    # --- Neutrophil x Cytotoxic T-cell ---
    out_neut = "figures/manuscript/neutrophil_cytotoxic_tcell_localization"
    if NEUTROPHIL_COL in sig_df.columns and CYTOTOXIC_COL in sig_df.columns:
        neut = sig_df[NEUTROPHIL_COL]
        cyto = sig_df[CYTOTOXIC_COL]
        trunc_neut = pd.Series(
            st_coloc.truncated_similarity(neut.values, cyto.values),
            index=sig_df.index,
        )
        panels_neut = _build_panels(
            adata,
            "solidity",
            neut,
            cyto,
            trunc_neut,
        )
        panels_neut[3]["title"] = "SLC16A3+ Neutrophil"
        panels_neut[4]["title"] = "Cytotoxic T-cell"

        print("Writing neutrophil–cytotoxic T-cell multipage PDF...")
        st_spatial.multipage_spatial_pdf(
            adata,
            library_id_col,
            panels_neut,
            os.path.join(out_neut, "spatial_multipage.pdf"),
        )
        prop_neut = _proportion_nonzero_df(adata, trunc_neut)
        _plot_proportion_boxplot(
            prop_neut,
            out_neut,
            "Proportion of Non-Zero Truncated Similarity (SLC16A3+ Neutrophil x Cytotoxic T-cell)\nby Solidity",
            "proportion_nonzero_boxplot_neutrophil_cytotoxic",
        )
    else:
        print("Skipping neutrophil–cytotoxic analysis (missing signature columns).")

    print("Done.")


if __name__ == "__main__":
    main()
