from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
from matplotlib.patches import Wedge

# -----------------------------
# Config supplied by you
# -----------------------------
_ann_new = Path("results/adata.annotated.h5ad")
_ann_old = Path("results/adata.annotated.p2.h5ad")
ADATA_PATH = _ann_new if _ann_new.exists() else _ann_old
cluster_key: str = "cluster"
library_id_key: str = "library_id"


def prepare_table(adata, solidity_order, arch_order, drop_arch=None):
    df = (
        adata.obs[["solidity_type", "architecture_type"]]
        .value_counts()
        .to_frame("count")
        .reset_index()
    )

    if drop_arch:
        df = df[~df["architecture_type"].isin(drop_arch)]

    table = (
        df.pivot(index="solidity_type", columns="architecture_type", values="count")
        .reindex(index=solidity_order, columns=arch_order)
        .fillna(0)
    )

    return table


def nested_donut(table, title=None, cmap_outer="tab20c", cmap_inner="tab20b"):
    solidity = table.index
    arch = table.columns

    # Outer ring totals
    outer_counts = table.sum(axis=1).values
    total_outer = outer_counts.sum()

    # Colors
    outer_colors = plt.get_cmap(cmap_outer)(np.linspace(0.1, 0.9, len(solidity)))
    inner_colors = plt.get_cmap(cmap_inner)(np.linspace(0.1, 0.9, len(arch)))

    fig, ax = plt.subplots(figsize=(4, 3), dpi=200)

    # ---- Draw outer ring via pie (donut) ----
    outer_wedges, _ = ax.pie(
        outer_counts,
        radius=1.0,
        startangle=90,
        labels=None,
        colors=outer_colors,
        wedgeprops=dict(width=0.3, edgecolor="white"),
    )

    # ---- Draw inner wedges manually using Wedge patches ----
    for s_idx, (sol_name, wedge) in enumerate(zip(solidity, outer_wedges)):
        theta1, theta2 = wedge.theta1, wedge.theta2
        span = theta2 - theta1

        row = table.loc[sol_name].values
        row_total = row.sum()
        if row_total == 0:
            continue

        proportions = row / row_total

        # architecture slice angles inside the solidity wedge
        arch_angles = proportions * span
        cum_angles = np.cumsum(np.concatenate([[0], arch_angles[:-1]]))

        for a_idx, (arch_name, ang, w) in enumerate(zip(arch, cum_angles, arch_angles)):
            if row[a_idx] == 0:
                continue

            a1 = theta1 + ang
            a2 = a1 + w

            # Draw architecture slice using Wedge
            patch = Wedge(
                center=(0, 0),
                r=0.7,
                theta1=a1,
                theta2=a2,
                width=0.3,
                facecolor=inner_colors[a_idx],
                edgecolor="black",  # border visible
                linewidth=0.7,
            )
            ax.add_patch(patch)

            # Label inner slice
            ang_mid = (a1 + a2) / 2
            x = np.cos(np.deg2rad(ang_mid)) * 0.55
            y = np.sin(np.deg2rad(ang_mid)) * 0.55
            pct = row[a_idx] / row_total * 100

            ax.text(
                x,
                y,
                f"{arch_name}\n{row[a_idx]} ({pct:.2f}%)",
                ha="center",
                va="center",
                fontsize=7.5,
            )

        # Label outer slice
        outer_pct = outer_counts[s_idx] / total_outer * 100
        ang_mid = (theta1 + theta2) / 2
        x = np.cos(np.deg2rad(ang_mid)) * 1.18
        y = np.sin(np.deg2rad(ang_mid)) * 1.18

        ax.text(
            x,
            y,
            f"{sol_name}\n{outer_counts[s_idx]} ({outer_pct:.2f}%)",
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
        )

    ax.set(aspect="equal")
    if title:
        ax.set_title(title, fontsize=16)

    return fig, ax




def stacked_bar_with_labels(table, colors=None, title=None):
    """
    table:
        pd.DataFrame indexed by solidity_type, columns = architecture_type
        values = counts

    Produces one stacked barplot:
        - bar height = counts
        - each stack annotated with count and proportion within that solidity
    """

    solidity = table.index.tolist()
    arch = table.columns.tolist()

    counts = table.values
    proportions = counts / counts.sum(axis=1, keepdims=True)

    if colors is None:
        cmap = plt.get_cmap("tab20b")
        colors = cmap(np.linspace(0.1, 0.9, len(arch)))

    fig, ax = plt.subplots(figsize=(10, 6))

    bottom = np.zeros(len(solidity))

    # ---- Draw stacks ----
    for i, a in enumerate(arch):
        ax.bar(
            solidity,
            counts[:, i],
            bottom=bottom,
            color=colors[i],
            edgecolor="black",
            linewidth=0.6,
            label=a,
        )

        # Annotate each segment
        for j, sol in enumerate(solidity):
            c = counts[j, i]
            if c == 0:
                continue

            pct = proportions[j, i] * 100

            y_pos = bottom[j] + c / 2

            ax.text(
                j, y_pos, f"{c} ({pct:.2f}%)", ha="center", va="center", fontsize=8, color="black"
            )

        bottom += counts[:, i]

    # ---- Axis formatting ----
    ax.set_ylabel("Count")
    ax.set_title(title, fontsize=14)
    ax.legend(title="Architecture", bbox_to_anchor=(1.05, 1), loc="upper left")

    ax.set_xticks(np.arange(len(solidity)))
    ax.set_xticklabels(solidity, rotation=45, ha="right")

    return fig, ax


# -----------------------------
# Run
# -----------------------------
adata = sc.read(ADATA_PATH)

# # Ensure cluster is categorical and colors exist
# if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
#     adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

# if f"{cluster_key}_colors" not in adata.uns:
#     adata.uns[f"{cluster_key}_colors"] = ["#cccccc"] * len(adata.obs[cluster_key].cat.categories)

solidity_order = ["Normal", "Non-Solid", "Solid"]

arch_order_with_core = ["Core", "Scar Tissue", "Blood Vessel", "Bronchus", "TLS"]
arch_order_no_core = ["Scar Tissue", "Blood Vessel", "Bronchus", "TLS"]

table_with_core = prepare_table(
    adata, solidity_order=solidity_order, arch_order=arch_order_with_core, drop_arch=None
)

nested_donut(table_with_core, title="Structure Proportion (with core)")
plt.tight_layout()
sns.despine()
plt.savefig("figures/nested_solidity_structure_count_with_core_pie.png")
plt.savefig("figures/nested_solidity_structure_count_with_core_pie.pdf")
plt.close()

stacked_bar_with_labels(table_with_core, title="Structure Proportion (with core)")
plt.tight_layout()
sns.despine()
plt.savefig("figures/nested_solidity_structure_count_with_core_bar.png")
plt.savefig("figures/nested_solidity_structure_count_with_core_bar.pdf")
plt.close()

table_no_core = prepare_table(
    adata, solidity_order=solidity_order, arch_order=arch_order_no_core, drop_arch=["Core"]
)
plt.tight_layout()
sns.despine()
plt.savefig("figures/nested_solidity_structure_count_without_core_pie.png")
plt.savefig("figures/nested_solidity_structure_count_without_core_pie.pdf")
plt.close()

nested_donut(table_no_core, title="Structure Proportion")
plt.tight_layout()
sns.despine()
plt.savefig("figures/nested_solidity_structure_count_without_core_pie.png")
plt.savefig("figures/nested_solidity_structure_count_without_core_pie.pdf")
plt.close()

stacked_bar_with_labels(table_no_core, title="Structure Proportion")
plt.tight_layout()
sns.despine()
plt.savefig("figures/nested_solidity_structure_count_without_core_bar.png")
plt.savefig("figures/nested_solidity_structure_count_without_core_bar.pdf")
plt.close()
