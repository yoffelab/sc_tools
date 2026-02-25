"""
Differential Program Analysis: 1-vs-all statistical testing across tumor types.

Performs 1-vs-all comparisons for:
- Tumor programs (EMT, Hypoxia, Proliferative)
- Immune programs (Immune_Lymphoid, Innate_Other, Processes:TLS_Formation, Liron)

All results use Benjamini-Hochberg (FDR) correction and include significance bars
with asterisks in boxplots.
"""

import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

sc.settings.set_figure_params(dpi=300, dpi_save=400)

# Load data
print("Loading data...")
adata = sc.read("results/adata.normalized.scored.p35.h5ad")

# Load gene signatures to identify program columns
with open("metadata/gene_signatures.json") as file:
    spatial_signatures = json.load(file)

# Define tumor type column (using pathologist annotation mapping)
pa = {
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


def keep_first_unique(input_list):
    seen = set()
    unique_list = []
    for item in input_list:
        if item not in seen:
            unique_list.append(item)
            seen.add(item)
    return unique_list


adata.obs["tumor_type"] = pd.Categorical(
    adata.obs["pathologist_annotation"].astype(str).replace(pa),
    categories=keep_first_unique(["Normal Alveolar Cells", "Non-Solid Tumor", "Solid Tumor"]),
)

# Filter to only the three main categories
adata = adata[
    adata.obs["tumor_type"].isin(["Normal Alveolar Cells", "Non-Solid Tumor", "Solid Tumor"])
].copy()
adata.obs["tumor_type"] = adata.obs["tumor_type"].cat.remove_unused_categories()

print(f"Tumor type distribution:\n{adata.obs['tumor_type'].value_counts()}")

# Define program signatures to test
tumor_programs = [
    "sig:Tumor_Cells/EMT_Tumor_z",
    "sig:Tumor_Cells/Hypoxia_Tumor_z",
    "sig:Tumor_Cells/Proliferative_Tumor_z",
]

immune_programs = [
    "sig:Immune_Lymphoid/Naive_CD4_T_z",
    "sig:Immune_Lymphoid/Naive_CD8_T_z",
    "sig:Immune_Lymphoid/Effector_CD8_T_z",
    "sig:Immune_Lymphoid/Exhausted_CD8_T_z",
    "sig:Immune_Lymphoid/Memory_CD4_T_z",
    "sig:Immune_Lymphoid/B_Naive_z",
    "sig:Immune_Lymphoid/B_Memory_z",
    "sig:Innate_Other/Mast_Cells_z",
    "sig:Processes/TLS_Formation_z",
    # Liron T cell signatures
    "sig:Liron/T cell exhaustion_z",
    "sig:Liron/T cell cytotoxicity_z",
    "sig:Liron/Tregs_z",
    "sig:Liron/Suppressive Tregs_z",
    "sig:Liron/CD8.4 TRM pre-exhausted_z",
    "sig:Liron/CD8.1 GZMK pre-exhausted_z",
    # Liron myeloid/macrophage signatures
    "sig:Liron/M2-macrophages_z",
    "sig:Liron/MoMacs Merad_z",
    "sig:Liron/Alveolar macrophages Merad_z",
    "sig:Liron/Mac.2 MoMac M2-like_z",
    "sig:Liron/Mac.6 MoMAc M2-like_z",
]

# Filter to only programs that exist in the data
all_programs = tumor_programs + immune_programs
available_programs = [p for p in all_programs if p in adata.obs.columns]
missing = [p for p in all_programs if p not in adata.obs.columns]
if missing:
    print(f"Warning: Missing program columns: {missing}")

print(
    f"Testing {len(available_programs)} programs across {len(adata.obs['tumor_type'].cat.categories)} tumor types"
)


# Function to perform 1-vs-all statistical testing
def test_one_vs_all(adata, program_col, group_col="tumor_type", method="mannwhitneyu"):
    """
    Perform 1-vs-all comparisons for each group.

    Returns DataFrame with columns: group, pval, adj_pval, n_group, n_rest
    """
    results = []
    groups = adata.obs[group_col].cat.categories

    for group in groups:
        group_mask = adata.obs[group_col] == group
        rest_mask = ~group_mask

        group_values = adata.obs.loc[group_mask, program_col].dropna()
        rest_values = adata.obs.loc[rest_mask, program_col].dropna()

        if len(group_values) < 3 or len(rest_values) < 3:
            results.append(
                {
                    "group": group,
                    "pval": np.nan,
                    "adj_pval": np.nan,
                    "n_group": len(group_values),
                    "n_rest": len(rest_values),
                    "mean_group": group_values.mean() if len(group_values) > 0 else np.nan,
                    "mean_rest": rest_values.mean() if len(rest_values) > 0 else np.nan,
                }
            )
            continue

        if method == "mannwhitneyu":
            stat, pval = mannwhitneyu(group_values, rest_values, alternative="two-sided")
        else:
            raise ValueError(f"Unknown method: {method}")

        results.append(
            {
                "group": group,
                "pval": pval,
                "adj_pval": np.nan,  # Will be filled after FDR correction
                "n_group": len(group_values),
                "n_rest": len(rest_values),
                "mean_group": group_values.mean(),
                "mean_rest": rest_values.mean(),
            }
        )

    df = pd.DataFrame(results)

    # Apply Benjamini-Hochberg correction
    valid_pvals = df["pval"].dropna()
    if len(valid_pvals) > 0:
        _, adj_pvals, _, _ = multipletests(valid_pvals, method="fdr_bh", alpha=0.05)
        df.loc[df["pval"].notna(), "adj_pval"] = adj_pvals

    return df


# Perform statistical testing for all programs
print("\nPerforming 1-vs-all statistical testing...")
all_stats = []

for program in available_programs:
    stats_df = test_one_vs_all(adata, program, group_col="tumor_type")
    stats_df["program"] = program
    all_stats.append(stats_df)

stats_combined = pd.concat(all_stats, ignore_index=True)

# Save statistics
os.makedirs("figures/stats", exist_ok=True)
stats_combined.to_csv("figures/stats/tumor_differences_1vsall_stats.csv", index=False)
print("\nSaved statistics to figures/stats/tumor_differences_1vsall_stats.csv")


# Function to plot boxplots with significance annotations
def plot_program_boxplot(
    adata, program_col, group_col="tumor_type", stats_df=None, output_path=None, figsize=(6, 4)
):
    """
    Create boxplot with significance bars and asterisks.
    Bars are properly dodged to avoid overlap, and all comparisons are listed.
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Prepare data
    plot_data = []
    groups = adata.obs[group_col].cat.categories

    for group in groups:
        values = adata.obs.loc[adata.obs[group_col] == group, program_col].dropna()
        plot_data.append(values.values)

    # Create boxplot
    bp = ax.boxplot(plot_data, tick_labels=groups, patch_artist=True, showfliers=False, widths=0.6)

    # Color boxes
    colors = ["#66c2a5", "#fc8d62", "#8da0cb"]  # Normal, Non-Solid, Solid
    for patch, color in zip(bp["boxes"], colors[: len(groups)]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Get initial y-axis limits for positioning (before adding bars)
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min

    # Calculate maximum whisker extent to avoid overlap
    max_whisker = y_max
    for values in plot_data:
        if len(values) > 0:
            q1 = np.percentile(values, 25)
            q3 = np.percentile(values, 75)
            iqr = q3 - q1
            upper_whisker = q3 + 1.5 * iqr
            max_whisker = max(max_whisker, upper_whisker)

    # Add significance bars and asterisks if stats provided
    comparison_texts = []
    if stats_df is not None:
        program_stats = stats_df[stats_df["program"] == program_col]

        # Collect all valid comparisons
        valid_comparisons = []
        for i, group in enumerate(groups):
            row = program_stats[program_stats["group"] == group]
            if len(row) > 0 and pd.notna(row.iloc[0]["adj_pval"]):
                adj_pval = row.iloc[0]["adj_pval"]

                # Determine asterisk level
                if adj_pval < 0.0001:
                    asterisk = "****"
                elif adj_pval < 0.001:
                    asterisk = "***"
                elif adj_pval < 0.01:
                    asterisk = "**"
                elif adj_pval < 0.05:
                    asterisk = "*"
                else:
                    asterisk = "ns"

                # Format comparison text: "Group A - Rest: p*-value"
                # For 1-vs-all, "Rest" means all other groups combined
                rest_label = "Rest"
                comparison_text = f"{group} - {rest_label}: {asterisk} (p={adj_pval:.2e})"
                comparison_texts.append(comparison_text)

                valid_comparisons.append(
                    {
                        "group_idx": i,
                        "group_name": group,
                        "adj_pval": adj_pval,
                        "asterisk": asterisk,
                        "x_pos": i + 1,
                    }
                )

        # Draw significance bars with proper dodging
        if len(valid_comparisons) > 0:
            # Calculate bar spacing to avoid overlap
            bar_spacing = y_range * 0.08  # Space between bars
            bar_height = y_range * 0.04  # Height of each bar

            # Start bars well above the maximum whisker to avoid interference
            # Add buffer of 10% of y-range above max whisker
            whisker_buffer = y_range * 0.10
            y_start = max_whisker + whisker_buffer

            # Filter to only significant comparisons for visualization
            sig_comparisons = [c for c in valid_comparisons if c["asterisk"] != "ns"]

            for idx, comp in enumerate(sig_comparisons):
                # Calculate y position for this bar (dodged)
                y_bar = y_start + (idx * bar_spacing)

                # Draw horizontal bar spanning all groups (group vs rest)
                x_start = 1
                x_end = len(groups)

                # Draw main horizontal bar
                ax.plot([x_start, x_end], [y_bar, y_bar], "k-", linewidth=1.5)

                # Draw vertical ticks at both ends
                tick_height = bar_height * 0.4
                ax.plot([x_start, x_start], [y_bar - tick_height, y_bar], "k-", linewidth=1.5)
                ax.plot([x_end, x_end], [y_bar - tick_height, y_bar], "k-", linewidth=1.5)

                # Add asterisk above bar, positioned over the specific group being compared
                x_asterisk = comp["x_pos"]
                ax.text(
                    x_asterisk,
                    y_bar + bar_height * 0.4,
                    comp["asterisk"],
                    ha="center",
                    va="bottom",
                    fontsize=11,
                    fontweight="bold",
                )

            # Adjust y-axis to accommodate bars (extend upward)
            if len(sig_comparisons) > 0:
                max_bar_y = y_start + (len(sig_comparisons) * bar_spacing) + bar_height * 0.6
                ax.set_ylim(y_min, max_bar_y)

        # Add all comparisons as text - one per line, no separators
        if len(comparison_texts) > 0:
            # Each comparison on its own line
            comparison_str = "\n".join(comparison_texts)

            # Calculate how much space we need for text (estimate)
            n_lines = len(comparison_texts)
            # Adjust bottom margin based on number of comparisons
            bottom_margin = max(0.15, 0.05 + (n_lines * 0.025))

            # Add text box at bottom - each comparison on separate line
            fig.text(
                0.5,
                0.01,
                comparison_str,
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

            # Adjust layout to accommodate text
            plt.tight_layout(rect=[0, bottom_margin, 1, 1])
        else:
            plt.tight_layout()

    ax.set_ylabel("Program Score (z-scored)", fontsize=11)
    ax.set_xlabel("Tumor Type", fontsize=11)
    ax.set_title(program_col.replace("sig:", "").replace("_z", ""), fontsize=12, fontweight="bold")
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    sns.despine(ax=ax)

    if output_path:
        plt.savefig(output_path, bbox_inches="tight", dpi=300)
        plt.close()
    else:
        plt.show()


# Create output directory
os.makedirs("figures/manuscript/tumor_differences", exist_ok=True)

# Plot tumor programs
print("\nGenerating boxplots for tumor programs...")
for program in tumor_programs:
    if program in available_programs:
        program_name = program.replace("sig:", "").replace("_z", "").replace("/", "_")
        output_pdf = f"figures/manuscript/tumor_differences/{program_name}_boxplot.pdf"
        output_png = f"figures/manuscript/tumor_differences/{program_name}_boxplot.png"

        plot_program_boxplot(
            adata, program, group_col="tumor_type", stats_df=stats_combined, output_path=output_pdf
        )
        plot_program_boxplot(
            adata, program, group_col="tumor_type", stats_df=stats_combined, output_path=output_png
        )
        print(f"  Saved: {program_name}")

# Plot immune programs
print("\nGenerating boxplots for immune programs...")
for program in immune_programs:
    if program in available_programs:
        program_name = program.replace("sig:", "").replace("_z", "").replace("/", "_")
        output_pdf = f"figures/manuscript/tumor_differences/{program_name}_boxplot.pdf"
        output_png = f"figures/manuscript/tumor_differences/{program_name}_boxplot.png"

        plot_program_boxplot(
            adata, program, group_col="tumor_type", stats_df=stats_combined, output_path=output_pdf
        )
        plot_program_boxplot(
            adata, program, group_col="tumor_type", stats_df=stats_combined, output_path=output_png
        )
        print(f"  Saved: {program_name}")

print("\n✅ Differential program analysis complete!")
print("   - Statistics saved to: figures/stats/tumor_differences_1vsall_stats.csv")
print("   - Plots saved to: figures/manuscript/tumor_differences/")
