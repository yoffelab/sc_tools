"""
Neutrophil-Cytotoxic T-cell Localization: Measure correlation between SLC16A3+ neutrophils and cytotoxic T-cells.

Tasks:
1. Measure correlation of SLC16A3+ neutrophil score (x) with cytotoxic T-cell score (y)
2. Stratify by Normal, Non-Solid, and Solid tumors
3. Scatterplot with regression lines and Pearson correlation per tumor type
4. Save PDF/PNG and correlation summary CSV

Uses obsm['signature_score_z'] via get_signature_df (obsm-first); fallback to obs if needed.
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr

from sc_tools.utils.signatures import get_signature_df

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# ============================================================================
# Configuration
# ============================================================================

# Signature column names (obsm full-path; fallback obs uses sig:..._z)
NEUTROPHIL_SIG = "Liron/SLC16A3+ Neutrophil"
CYTOTOXIC_TCELL_SIG = "Liron/T cell cytotoxicity"

# Tumor type mapping (same as macrophage_localization)
TUMOR_TYPE_MAPPING = {
    "Solid": "Solid Tumor",
    "Non-Solid": "Non-Solid Tumor",
    "Normal": "Normal Alveolar Cells",
    "Solid Blood Vessel": "Solid Tumor",
    "Solid Bronchus": "Solid Tumor",
    "Solid Scar Tissue": "Solid Tumor",
    "Non-Solid Blood Vessel": "Non-Solid Tumor",
    "Non-Solid Bronchus": "Non-Solid Tumor",
    "Normal Blood Vessel": "Normal Alveolar Cells",
    "Normal Bronchus": "Normal Alveolar Cells",
    "Solid TLS": "Solid Tumor",
    "TLS Solid": "Solid Tumor",
    "Non-Solid TLS": "Non-Solid Tumor",
    "TLS Non-Solid": "Non-Solid Tumor",
    "TLS Normal": "Normal Alveolar Cells",
}

TUMOR_COLORS = {
    "Normal Alveolar Cells": "#66c2a5",
    "Non-Solid Tumor": "#fc8d62",
    "Solid Tumor": "#8da0cb",
}

OUTPUT_DIR = "figures/manuscript/neutrophil_cytotoxic_tcell_localization"
STATS_DIR = "figures/stats"

# ============================================================================
# Helper Functions
# ============================================================================


def keep_first_unique(input_list):
    """Keep first occurrence of each unique item."""
    seen = set()
    unique_list = []
    for item in input_list:
        if item not in seen:
            unique_list.append(item)
            seen.add(item)
    return unique_list


def calculate_correlation(x, y):
    """Calculate Pearson correlation, handling NaN values."""
    mask = ~(np.isnan(x) | np.isnan(y))
    if np.sum(mask) < 3:
        return np.nan, np.nan
    r, p = pearsonr(x[mask], y[mask])
    return r, p


def _resolve_sig_column(adata, sig_name_obsm, sig_name_obs_fallback=None):
    """
    Resolve signature column: obsm (full path) or obs (sig:..._z).
    Returns (values array, column name used) or (None, None) if not found.
    """
    sig_df = get_signature_df(adata, use_z=True)
    if sig_df is not None and not sig_df.empty and sig_name_obsm in sig_df.columns:
        return sig_df[sig_name_obsm].values, sig_name_obsm
    if sig_name_obs_fallback and sig_name_obs_fallback in adata.obs.columns:
        return adata.obs[sig_name_obs_fallback].values, sig_name_obs_fallback
    # Try obs with sig: prefix and _z suffix
    obs_candidate = f"sig:{sig_name_obsm}_z" if not sig_name_obsm.startswith("sig:") else sig_name_obsm
    if obs_candidate in adata.obs.columns:
        return adata.obs[obs_candidate].values, obs_candidate
    return None, None


def plot_neutrophil_cytotoxic_colocalization(
    adata,
    x_sig_name,
    y_sig_name,
    x_data,
    y_data,
    group_col="tumor_type",
):
    """
    Plot scatterplot of SLC16A3+ neutrophil (x) vs cytotoxic T-cell (y), stratified by tumor type.
    """
    valid_types = ["Normal Alveolar Cells", "Non-Solid Tumor", "Solid Tumor"]
    adata_sub = adata[adata.obs[group_col].isin(valid_types)].copy()
    adata_sub.obs[group_col] = adata_sub.obs[group_col].cat.remove_unused_categories()

    if adata_sub.n_obs == 0:
        print("Warning: No data after filtering by tumor type.")
        return None

    # x_data, y_data are aligned with adata (same length and order as adata.obs)
    x_aligned = np.asarray(x_data)
    y_aligned = np.asarray(y_data)
    groups = adata_sub.obs[group_col].values

    filter_mask = (x_aligned >= 0) & (y_aligned >= 0)
    x_filtered = x_aligned[filter_mask]
    y_filtered = y_aligned[filter_mask]
    groups_filtered = groups[filter_mask]

    x_excluded = x_aligned[~filter_mask]
    y_excluded = y_aligned[~filter_mask]

    fig, ax = plt.subplots(figsize=(7, 6))

    if len(x_excluded) > 0:
        ax.scatter(
            x_excluded,
            y_excluded,
            alpha=0.3,
            s=15,
            color="lightgray",
            edgecolors="none",
            label="_nolegend_",
            zorder=0,
        )

    ax.axvline(x=0, color="black", linestyle="--", linewidth=0.8, alpha=0.5, zorder=1)
    ax.axhline(y=0, color="black", linestyle="--", linewidth=0.8, alpha=0.5, zorder=1)

    correlations = {}
    for tumor_type in valid_types:
        if tumor_type not in adata_sub.obs[group_col].cat.categories:
            continue

        mask = groups_filtered == tumor_type
        x_vals = x_filtered[mask]
        y_vals = y_filtered[mask]

        valid_mask = ~(np.isnan(x_vals) | np.isnan(y_vals))
        x_vals = x_vals[valid_mask]
        y_vals = y_vals[valid_mask]

        if len(x_vals) < 3:
            print(f"Warning: Too few points for {tumor_type} ({len(x_vals)}). Skipping...")
            continue

        r, p = calculate_correlation(x_vals, y_vals)
        correlations[tumor_type] = {"r": r, "p": p, "n": len(x_vals)}

        color = TUMOR_COLORS.get(tumor_type, "gray")
        label = tumor_type.replace(" Alveolar Cells", "").replace(" Tumor", "")
        ax.scatter(
            x_vals, y_vals, alpha=0.5, s=20, color=color, label=label, edgecolors="none", zorder=3
        )

        if len(x_vals) >= 3:
            sort_idx = np.argsort(x_vals)
            x_sorted = x_vals[sort_idx]
            y_sorted = y_vals[sort_idx]
            try:
                slope, intercept, _, _, std_err = stats.linregress(x_sorted, y_sorted)
                if np.isfinite(slope) and np.isfinite(intercept) and np.isfinite(std_err):
                    x_min, x_max = x_sorted.min(), x_sorted.max()
                    if x_max > x_min:
                        x_line = np.linspace(x_min, x_max, 100)
                        y_line = intercept + slope * x_line
                        ax.plot(
                            x_line,
                            y_line,
                            color=color,
                            linewidth=2,
                            alpha=0.8,
                            linestyle="-",
                            zorder=2,
                            label="_nolegend_",
                        )
            except Exception as e:
                print(f"Warning: Error fitting regression for {tumor_type}: {e}")

    x_label = x_sig_name.replace("Liron/", "").replace("sig:", "").replace("_z", "")
    y_label = y_sig_name.replace("Liron/", "").replace("sig:", "").replace("_z", "")
    ax.set_xlabel(f"{x_label} Score (z-scored)", fontsize=12)
    ax.set_ylabel(f"{y_label} Score (z-scored)", fontsize=12)
    ax.set_title("SLC16A3+ Neutrophil vs Cytotoxic T-cell", fontsize=13, fontweight="bold")
    ax.grid(alpha=0.3, linestyle="--", zorder=0)
    ax.legend(
        title="Tumor Type", fontsize=10, title_fontsize=11, frameon=True, fancybox=True, shadow=True
    )
    sns.despine(ax=ax)

    correlation_text = []
    for tumor_type in valid_types:
        if tumor_type in correlations:
            r = correlations[tumor_type]["r"]
            p = correlations[tumor_type]["p"]
            n = correlations[tumor_type]["n"]
            if not np.isnan(r):
                p_str = "p<0.001" if p < 0.001 else f"p={p:.3f}" if p < 0.05 else f"p={p:.3f} (n.s.)"
                label = tumor_type.replace(" Alveolar Cells", "").replace(" Tumor", "")
                correlation_text.append(f"{label}: r={r:.3f}, {p_str} (n={n})")

    if correlation_text:
        fig.text(
            0.5,
            0.02,
            "\n".join(correlation_text),
            ha="center",
            va="bottom",
            fontsize=9,
            bbox=dict(
                boxstyle="round,pad=0.5",
                facecolor="lightyellow",
                alpha=0.8,
                edgecolor="black",
                linewidth=1,
            ),
            family="monospace",
        )

    plt.tight_layout(rect=[0, 0.12, 1, 1])
    return correlations


# ============================================================================
# Main
# ============================================================================


def main():
    print("=" * 70)
    print("Neutrophil vs Cytotoxic T-cell Localization Analysis")
    print("=" * 70)

    _scored_new = Path("results/adata.scored.h5ad")
    _scored_old = Path("results/adata.normalized.scored.p35.h5ad")
    adata_path = _scored_new if _scored_new.exists() else _scored_old
    if not adata_path.exists():
        raise FileNotFoundError(f"AnnData file not found: {adata_path}. Run score_gene_signatures.py first.")

    print("\n1. Loading data...")
    adata = sc.read_h5ad(adata_path)
    print(f"   Loaded: {adata.shape} (spots x genes)")

    print("\n2. Creating tumor type annotations...")
    adata.obs["tumor_type"] = pd.Categorical(
        adata.obs["pathologist_annotation"].astype(str).replace(TUMOR_TYPE_MAPPING),
        categories=keep_first_unique(["Normal Alveolar Cells", "Non-Solid Tumor", "Solid Tumor"]),
    )
    adata = adata[
        adata.obs["tumor_type"].isin(["Normal Alveolar Cells", "Non-Solid Tumor", "Solid Tumor"])
    ].copy()
    adata.obs["tumor_type"] = adata.obs["tumor_type"].cat.remove_unused_categories()
    print("   Tumor type distribution:")
    print(adata.obs["tumor_type"].value_counts())

    print("\n3. Resolving signature columns...")
    x_vals, x_col = _resolve_sig_column(
        adata, NEUTROPHIL_SIG, sig_name_obs_fallback="sig:Liron/SLC16A3+ Neutrophil_z"
    )
    y_vals, y_col = _resolve_sig_column(
        adata, CYTOTOXIC_TCELL_SIG, sig_name_obs_fallback="sig:Liron/T cell cytotoxicity_z"
    )

    if x_vals is None:
        raise ValueError(
            f"Neutrophil signature not found. Expected obsm column '{NEUTROPHIL_SIG}' or obs. "
            "Re-run score_gene_signatures.py after adding SLC16A3+ Neutrophil to gene_signatures.json."
        )
    if y_vals is None:
        raise ValueError(
            f"Cytotoxic T-cell signature not found. Expected obsm column '{CYTOTOXIC_TCELL_SIG}'."
        )
    print(f"   X: {x_col}")
    print(f"   Y: {y_col}")

    print("\n4. Generating colocalization plot...")
    correlations = plot_neutrophil_cytotoxic_colocalization(
        adata, x_col, y_col, x_vals, y_vals, group_col="tumor_type"
    )

    if correlations is None:
        print("   No plot generated (no valid data).")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    Path(STATS_DIR).mkdir(parents=True, exist_ok=True)

    out_pdf = os.path.join(OUTPUT_DIR, "SLC16A3_neutrophil_vs_cytotoxic_tcell.pdf")
    out_png = os.path.join(OUTPUT_DIR, "SLC16A3_neutrophil_vs_cytotoxic_tcell.png")
    plt.savefig(out_pdf, bbox_inches="tight", dpi=300)
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"   Saved: {out_pdf}")
    print(f"   Saved: {out_png}")

    print("\n5. Saving correlation summary...")
    summary_data = [
        {
            "Neutrophil_Signature": x_col,
            "Cytotoxic_Tcell_Signature": y_col,
            "Tumor_Type": tumor_type,
            "Pearson_r": st["r"],
            "p_value": st["p"],
            "n_spots": st["n"],
        }
        for tumor_type, st in correlations.items()
    ]
    summary_df = pd.DataFrame(summary_data)
    summary_path = Path(STATS_DIR) / "neutrophil_cytotoxic_tcell_correlations.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"   Saved: {summary_path}")

    print("\n" + "=" * 70)
    print("Neutrophil vs cytotoxic T-cell localization analysis complete.")
    print(f"   Figures: {OUTPUT_DIR}/")
    print(f"   Correlations CSV: {summary_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()
