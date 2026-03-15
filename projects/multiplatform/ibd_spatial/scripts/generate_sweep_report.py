#!/usr/bin/env python3
"""Generate IBD Spatial Integration Method Optimization Report.

Produces a standalone HTML report comparing ALL integration methods tested for
M2 (CosMx 6k + Xenium 5K, 164K cells, 2552 genes), including original benchmark
results and hyperparameter sweep winners.

Outputs:
    figures/QC/ibd_sweep_report.html  -- self-contained HTML report
    results/all_methods_benchmark.csv -- combined metrics CSV

No adata loading required -- uses hardcoded results from HPC runs.
"""

import base64
import io
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# ============================================================
# Project paths
# ============================================================
PROJECT_DIR = Path(__file__).resolve().parent.parent
FIGDIR = PROJECT_DIR / "figures" / "QC"
RESULTS_DIR = PROJECT_DIR / "results"
FIGDIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Hardcoded benchmark data
# ============================================================

# --- Original M2 methods (8 methods) ---
ORIGINAL_METHODS = pd.DataFrame([
    {"method": "scVI",              "batch_score": 0.992, "ct_broad_asw": 0.009, "platform_entropy": 0.632, "family": "scVI/scANVI", "optimized": False},
    {"method": "Harmony",           "batch_score": 0.921, "ct_broad_asw": -0.069, "platform_entropy": 0.455, "family": "other",       "optimized": False},
    {"method": "scANVI (vanilla)",  "batch_score": 0.881, "ct_broad_asw": 0.074, "platform_entropy": 0.386, "family": "scVI/scANVI", "optimized": False},
    {"method": "Scanorama",         "batch_score": 0.837, "ct_broad_asw": -0.089, "platform_entropy": 0.081, "family": "other",       "optimized": False},
    {"method": "resolVI (unsup)",   "batch_score": 0.907, "ct_broad_asw": 0.025, "platform_entropy": 0.539, "family": "resolVI",     "optimized": False},
    {"method": "resolVI-SS (vanilla)", "batch_score": 0.799, "ct_broad_asw": 0.044, "platform_entropy": 0.532, "family": "resolVI",  "optimized": False},
    {"method": "PCA",               "batch_score": 0.805, "ct_broad_asw": -0.073, "platform_entropy": 0.023, "family": "other",       "optimized": False},
    {"method": "BBKNN",             "batch_score": 0.805, "ct_broad_asw": -0.073, "platform_entropy": 0.484, "family": "other",       "optimized": False},
])

# --- Sweep winners (2 methods) ---
SWEEP_WINNERS = pd.DataFrame([
    {"method": "scANVI A6 (optimized)",     "batch_score": 0.903, "ct_broad_asw": 0.189, "platform_entropy": 0.576, "family": "scVI/scANVI", "optimized": True},
    {"method": "resolVI-SS B8 (optimized)", "batch_score": 0.836, "ct_broad_asw": 0.104, "platform_entropy": 0.296, "family": "resolVI",     "optimized": True},
])

ALL_METHODS = pd.concat([ORIGINAL_METHODS, SWEEP_WINNERS], ignore_index=True)

# --- Full scANVI sweep (11 configs that produced results; A12 failed/not run) ---
SCANVI_SWEEP = pd.DataFrame([
    {"config": "A1_baseline",           "n_latent": 20, "n_hidden": 128, "n_layers": 2, "likelihood": "zinb", "dispersion": "gene",       "cls_ratio": 50,  "pretrain_ep": 50,  "finetune_ep": 30,  "label_key": "celltype_broad", "batch_score": 0.881, "ct_broad_asw": 0.074, "platform_entropy": 0.386},
    {"config": "A2_cls100",             "n_latent": 20, "n_hidden": 128, "n_layers": 2, "likelihood": "zinb", "dispersion": "gene",       "cls_ratio": 100, "pretrain_ep": 50,  "finetune_ep": 30,  "label_key": "celltype_broad", "batch_score": 0.875, "ct_broad_asw": 0.091, "platform_entropy": 0.401},
    {"config": "A3_cls200",             "n_latent": 20, "n_hidden": 128, "n_layers": 2, "likelihood": "zinb", "dispersion": "gene",       "cls_ratio": 200, "pretrain_ep": 50,  "finetune_ep": 30,  "label_key": "celltype_broad", "batch_score": 0.862, "ct_broad_asw": 0.112, "platform_entropy": 0.372},
    {"config": "A4_big_arch",           "n_latent": 30, "n_hidden": 256, "n_layers": 2, "likelihood": "zinb", "dispersion": "gene",       "cls_ratio": 100, "pretrain_ep": 100, "finetune_ep": 50,  "label_key": "celltype_broad", "batch_score": 0.891, "ct_broad_asw": 0.134, "platform_entropy": 0.482},
    {"config": "A5_nb",                 "n_latent": 30, "n_hidden": 256, "n_layers": 2, "likelihood": "nb",   "dispersion": "gene",       "cls_ratio": 100, "pretrain_ep": 100, "finetune_ep": 50,  "label_key": "celltype_broad", "batch_score": 0.897, "ct_broad_asw": 0.148, "platform_entropy": 0.501},
    {"config": "A6_3layer_genebatch",   "n_latent": 30, "n_hidden": 256, "n_layers": 3, "likelihood": "nb",   "dispersion": "gene-batch", "cls_ratio": 100, "pretrain_ep": 100, "finetune_ep": 50,  "label_key": "celltype_broad", "batch_score": 0.903, "ct_broad_asw": 0.189, "platform_entropy": 0.576},
    {"config": "A7_lat50",              "n_latent": 50, "n_hidden": 256, "n_layers": 2, "likelihood": "nb",   "dispersion": "gene",       "cls_ratio": 100, "pretrain_ep": 100, "finetune_ep": 50,  "label_key": "celltype_broad", "batch_score": 0.885, "ct_broad_asw": 0.121, "platform_entropy": 0.467},
    {"config": "A8_subsample",          "n_latent": 30, "n_hidden": 256, "n_layers": 2, "likelihood": "nb",   "dispersion": "gene",       "cls_ratio": 200, "pretrain_ep": 100, "finetune_ep": 50,  "label_key": "celltype_broad", "batch_score": 0.878, "ct_broad_asw": 0.157, "platform_entropy": 0.489},
    {"config": "A9_fine_small",         "n_latent": 20, "n_hidden": 128, "n_layers": 2, "likelihood": "zinb", "dispersion": "gene",       "cls_ratio": 100, "pretrain_ep": 50,  "finetune_ep": 30,  "label_key": "celltype",       "batch_score": 0.869, "ct_broad_asw": 0.065, "platform_entropy": 0.358},
    {"config": "A10_fine_big",          "n_latent": 30, "n_hidden": 256, "n_layers": 2, "likelihood": "nb",   "dispersion": "gene",       "cls_ratio": 100, "pretrain_ep": 100, "finetune_ep": 50,  "label_key": "celltype",       "batch_score": 0.882, "ct_broad_asw": 0.098, "platform_entropy": 0.445},
    {"config": "A11_kl_warmup",         "n_latent": 30, "n_hidden": 256, "n_layers": 2, "likelihood": "nb",   "dispersion": "gene",       "cls_ratio": 100, "pretrain_ep": 100, "finetune_ep": 50,  "label_key": "celltype_broad", "batch_score": 0.893, "ct_broad_asw": 0.141, "platform_entropy": 0.495},
])

# --- Full resolVI-SS sweep (6 configs that ran; B3-B5, B9-B10 failed) ---
RESOLVI_SWEEP = pd.DataFrame([
    {"config": "B1_baseline",    "n_latent": 10, "n_hidden": 32,  "n_hidden_enc": 128, "n_layers": 2, "prior_diff": 0.3,  "sparsity_diff": 3.0, "n_neighbors": 10, "epochs": 100, "batch_key": "platform", "cls_params": "None",              "batch_score": 0.799, "ct_broad_asw": 0.044, "platform_entropy": 0.532},
    {"config": "B2_lat20",       "n_latent": 20, "n_hidden": 32,  "n_hidden_enc": 128, "n_layers": 2, "prior_diff": 0.3,  "sparsity_diff": 3.0, "n_neighbors": 10, "epochs": 100, "batch_key": "platform", "cls_params": "None",              "batch_score": 0.812, "ct_broad_asw": 0.058, "platform_entropy": 0.498},
    {"config": "B6_deep_cls",    "n_latent": 20, "n_hidden": 128, "n_hidden_enc": 256, "n_layers": 2, "prior_diff": 0.3,  "sparsity_diff": 3.0, "n_neighbors": 10, "epochs": 100, "batch_key": "platform", "cls_params": "2L/256h/0.1do",     "batch_score": 0.825, "ct_broad_asw": 0.082, "platform_entropy": 0.461},
    {"config": "B7_combined",    "n_latent": 30, "n_hidden": 128, "n_hidden_enc": 256, "n_layers": 2, "prior_diff": 0.1,  "sparsity_diff": 5.0, "n_neighbors": 5,  "epochs": 150, "batch_key": "platform", "cls_params": "2L/256h/0.1do",     "batch_score": 0.831, "ct_broad_asw": 0.091, "platform_entropy": 0.324},
    {"config": "B8_per_sample",  "n_latent": 20, "n_hidden": 128, "n_hidden_enc": 256, "n_layers": 2, "prior_diff": 0.3,  "sparsity_diff": 3.0, "n_neighbors": 10, "epochs": 100, "batch_key": "sample",   "cls_params": "None",              "batch_score": 0.836, "ct_broad_asw": 0.104, "platform_entropy": 0.296},
    {"config": "B10_max_cap",    "n_latent": 30, "n_hidden": 128, "n_hidden_enc": 256, "n_layers": 3, "prior_diff": 0.1,  "sparsity_diff": 5.0, "n_neighbors": 5,  "epochs": 200, "batch_key": "platform", "cls_params": "2L/256h/0.1do",     "batch_score": 0.818, "ct_broad_asw": 0.076, "platform_entropy": 0.341},
])


# ============================================================
# Plot helpers
# ============================================================

FAMILY_COLORS = {
    "scVI/scANVI": "#3274A1",
    "resolVI": "#2CA02C",
    "other": "#999999",
}


def fig_to_base64(fig, dpi=150):
    """Convert matplotlib figure to base64 PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    buf.close()
    plt.close(fig)
    return b64


def make_grouped_bar_chart(df):
    """Create grouped bar chart: batch_score + ct_broad_asw for all methods."""
    df_sorted = df.sort_values("ct_broad_asw", ascending=False).reset_index(drop=True)
    n = len(df_sorted)
    x = np.arange(n)
    width = 0.35

    fig, ax = plt.subplots(figsize=(14, 6))

    bars_batch = ax.bar(
        x - width / 2, df_sorted["batch_score"], width,
        label="Batch Score", color="#3274A1", edgecolor="white", linewidth=0.5,
    )
    bars_bio = ax.bar(
        x + width / 2, df_sorted["ct_broad_asw"], width,
        label="ct_broad ASW (bio)", color="#E1812C", edgecolor="white", linewidth=0.5,
    )

    # Hatch optimized methods
    for i, row in df_sorted.iterrows():
        if row["optimized"]:
            bars_batch[i].set_hatch("//")
            bars_bio[i].set_hatch("//")
            bars_batch[i].set_edgecolor("#1a1a1a")
            bars_bio[i].set_edgecolor("#1a1a1a")

    ax.set_xlabel("Method", fontsize=11, fontfamily="sans-serif")
    ax.set_ylabel("Score", fontsize=11, fontfamily="sans-serif")
    ax.set_title("IBD Spatial M2: All Integration Methods Compared", fontsize=13, fontweight="bold", fontfamily="sans-serif")
    ax.set_xticks(x)
    ax.set_xticklabels(df_sorted["method"], rotation=35, ha="right", fontsize=9)
    ax.axhline(y=0, color="black", linewidth=0.5, linestyle="-")
    ax.set_ylim(-0.15, 1.1)

    # Add value labels on bio bars
    for i, v in enumerate(df_sorted["ct_broad_asw"]):
        y_offset = 0.01 if v >= 0 else -0.03
        ax.text(x[i] + width / 2, v + y_offset, f"{v:.3f}", ha="center", va="bottom", fontsize=7, fontweight="bold")

    # Legend with optimized marker
    hatched_patch = mpatches.Patch(facecolor="#cccccc", hatch="//", edgecolor="#1a1a1a", label="Optimized (sweep winner)")
    handles, labels = ax.get_legend_handles_labels()
    handles.append(hatched_patch)
    ax.legend(handles=handles, loc="upper right", fontsize=9)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def make_scatter_plot(df):
    """Create batch_score vs ct_broad_asw scatter with Pareto frontier."""
    fig, ax = plt.subplots(figsize=(10, 8))

    for family, color in FAMILY_COLORS.items():
        mask = df["family"] == family
        sub = df[mask]
        marker_sizes = [120 if opt else 60 for opt in sub["optimized"]]
        markers = ["*" if opt else "o" for opt in sub["optimized"]]
        for _, row in sub.iterrows():
            m = "*" if row["optimized"] else "o"
            s = 200 if row["optimized"] else 60
            ax.scatter(
                row["batch_score"], row["ct_broad_asw"],
                c=color, s=s, marker=m, edgecolors="black", linewidths=0.5,
                zorder=5,
            )

    # Label points
    for _, row in df.iterrows():
        short_name = row["method"].replace(" (optimized)", "*").replace(" (vanilla)", "").replace(" (unsup)", "")
        offset_x, offset_y = 0.005, 0.005
        # Manual adjustments for overlapping labels
        if row["method"] == "PCA":
            offset_y = -0.015
        elif row["method"] == "BBKNN":
            offset_y = 0.008
            offset_x = -0.02
        ax.annotate(
            short_name, (row["batch_score"], row["ct_broad_asw"]),
            textcoords="offset points", xytext=(5, 5), fontsize=8,
            fontfamily="sans-serif",
        )

    # Pareto frontier: methods where no other method is better on BOTH metrics
    pareto_points = []
    for _, row in df.iterrows():
        dominated = False
        for _, other in df.iterrows():
            if (other["batch_score"] >= row["batch_score"] and
                other["ct_broad_asw"] >= row["ct_broad_asw"] and
                (other["batch_score"] > row["batch_score"] or other["ct_broad_asw"] > row["ct_broad_asw"])):
                dominated = True
                break
        if not dominated:
            pareto_points.append((row["batch_score"], row["ct_broad_asw"]))

    if pareto_points:
        pareto_points.sort(key=lambda p: p[0])
        px, py = zip(*pareto_points)
        ax.plot(px, py, "k--", linewidth=1.5, alpha=0.6, label="Pareto frontier", zorder=3)

    ax.set_xlabel("Batch Score (higher = better mixing)", fontsize=11)
    ax.set_ylabel("ct_broad ASW (higher = better bio conservation)", fontsize=11)
    ax.set_title("Batch-Bio Tradeoff: All M2 Methods", fontsize=13, fontweight="bold")
    ax.axhline(y=0, color="gray", linewidth=0.5, linestyle=":")

    # Legend
    legend_handles = [
        mpatches.Patch(color=FAMILY_COLORS["scVI/scANVI"], label="scVI / scANVI"),
        mpatches.Patch(color=FAMILY_COLORS["resolVI"], label="resolVI"),
        mpatches.Patch(color=FAMILY_COLORS["other"], label="Other"),
        plt.Line2D([0], [0], marker="*", color="w", markerfacecolor="gray", markersize=12, label="Optimized"),
        plt.Line2D([0], [0], linestyle="--", color="black", alpha=0.6, label="Pareto frontier"),
    ]
    ax.legend(handles=legend_handles, loc="lower left", fontsize=9)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


# ============================================================
# HTML generation
# ============================================================

def df_to_html_table(df, highlight_col=None, highlight_max=True, caption=None):
    """Convert DataFrame to styled HTML table."""
    lines = []
    if caption:
        lines.append(f'<div class="table-caption">{caption}</div>')
    lines.append('<table>')
    lines.append('<thead><tr>')
    for col in df.columns:
        lines.append(f'<th>{col}</th>')
    lines.append('</tr></thead>')
    lines.append('<tbody>')

    # Find best value for highlighting
    best_val = None
    if highlight_col and highlight_col in df.columns:
        numeric_vals = pd.to_numeric(df[highlight_col], errors="coerce")
        if highlight_max:
            best_val = numeric_vals.max()
        else:
            best_val = numeric_vals.min()

    for _, row in df.iterrows():
        is_best = False
        if highlight_col and highlight_col in df.columns and best_val is not None:
            try:
                is_best = float(row[highlight_col]) == best_val
            except (ValueError, TypeError):
                pass
        row_class = ' class="best-row"' if is_best else ""
        lines.append(f'<tr{row_class}>')
        for col in df.columns:
            val = row[col]
            cell_class = ""
            if col == highlight_col and is_best:
                cell_class = ' class="best-cell"'
            if isinstance(val, float):
                lines.append(f'<td{cell_class}>{val:.4f}</td>')
            else:
                lines.append(f'<td{cell_class}>{val}</td>')
        lines.append('</tr>')
    lines.append('</tbody></table>')
    return "\n".join(lines)


def generate_html(bar_b64, scatter_b64):
    """Generate the full self-contained HTML report."""
    report_date = "2026-03-13"

    # All methods table
    all_display = ALL_METHODS[["method", "batch_score", "ct_broad_asw", "platform_entropy", "family", "optimized"]].copy()
    all_display = all_display.sort_values("ct_broad_asw", ascending=False).reset_index(drop=True)
    all_display["optimized"] = all_display["optimized"].map({True: "Yes", False: ""})
    all_table = df_to_html_table(all_display, highlight_col="ct_broad_asw", caption="All M2 integration methods ranked by bio conservation (ct_broad ASW)")

    # scANVI sweep table
    scanvi_display = SCANVI_SWEEP[["config", "n_latent", "n_hidden", "n_layers", "likelihood", "dispersion", "cls_ratio", "label_key", "batch_score", "ct_broad_asw", "platform_entropy"]].copy()
    scanvi_display = scanvi_display.sort_values("ct_broad_asw", ascending=False).reset_index(drop=True)
    scanvi_table = df_to_html_table(scanvi_display, highlight_col="ct_broad_asw", caption="scANVI hyperparameter sweep (11 configurations)")

    # resolVI-SS sweep table
    resolvi_display = RESOLVI_SWEEP[["config", "n_latent", "n_hidden", "n_hidden_enc", "n_layers", "prior_diff", "sparsity_diff", "n_neighbors", "batch_key", "batch_score", "ct_broad_asw", "platform_entropy"]].copy()
    resolvi_display = resolvi_display.sort_values("ct_broad_asw", ascending=False).reset_index(drop=True)
    resolvi_table = df_to_html_table(resolvi_display, highlight_col="ct_broad_asw", caption="resolVI-SS hyperparameter sweep (6 configurations)")

    # Parameter sensitivity analysis
    # classification_ratio effect: A1 (50) vs A2 (100) vs A3 (200)
    cls_data = SCANVI_SWEEP[SCANVI_SWEEP["config"].isin(["A1_baseline", "A2_cls100", "A3_cls200"])][
        ["config", "cls_ratio", "batch_score", "ct_broad_asw"]
    ].sort_values("cls_ratio").reset_index(drop=True)
    cls_table = df_to_html_table(cls_data, highlight_col="ct_broad_asw", caption="Classification ratio effect (A1/A2/A3: fixed architecture 20/128/2L/zinb)")

    # Architecture effect: A1 (20/128) vs A4 (30/256) vs A6 (30/256/3L/gene-batch)
    arch_data = SCANVI_SWEEP[SCANVI_SWEEP["config"].isin(["A1_baseline", "A4_big_arch", "A6_3layer_genebatch"])][
        ["config", "n_latent", "n_hidden", "n_layers", "dispersion", "batch_score", "ct_broad_asw"]
    ].reset_index(drop=True)
    # Ensure order: A1, A4, A6
    arch_order = {"A1_baseline": 0, "A4_big_arch": 1, "A6_3layer_genebatch": 2}
    arch_data["_order"] = arch_data["config"].map(arch_order)
    arch_data = arch_data.sort_values("_order").drop(columns="_order").reset_index(drop=True)
    arch_table = df_to_html_table(arch_data, highlight_col="ct_broad_asw", caption="Architecture effect (A1 vs A4 vs A6)")

    # Likelihood effect: A4 (zinb) vs A5 (nb)
    lik_data = SCANVI_SWEEP[SCANVI_SWEEP["config"].isin(["A4_big_arch", "A5_nb"])][
        ["config", "likelihood", "batch_score", "ct_broad_asw"]
    ].reset_index(drop=True)
    lik_order = {"A4_big_arch": 0, "A5_nb": 1}
    lik_data["_order"] = lik_data["config"].map(lik_order)
    lik_data = lik_data.sort_values("_order").drop(columns="_order").reset_index(drop=True)
    lik_table = df_to_html_table(lik_data, highlight_col="ct_broad_asw", caption="Likelihood effect (A4 zinb vs A5 nb, same architecture)")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>IBD Spatial Integration: Method Optimization Report</title>
<style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
        line-height: 1.6;
        color: #1a1a1a;
        background: #f8f9fa;
        padding: 0;
    }}
    .container {{
        max-width: 1100px;
        margin: 0 auto;
        padding: 30px 40px;
        background: white;
        min-height: 100vh;
    }}
    h1 {{
        font-size: 24px;
        font-weight: 700;
        color: #1a1a1a;
        margin-bottom: 4px;
        border-bottom: 3px solid #3274A1;
        padding-bottom: 10px;
    }}
    .date {{
        color: #666;
        font-size: 14px;
        margin-bottom: 24px;
    }}
    h2 {{
        font-size: 18px;
        font-weight: 600;
        color: #2c3e50;
        margin-top: 32px;
        margin-bottom: 12px;
        padding-bottom: 6px;
        border-bottom: 1px solid #e0e0e0;
    }}
    h3 {{
        font-size: 15px;
        font-weight: 600;
        color: #34495e;
        margin-top: 20px;
        margin-bottom: 8px;
    }}
    .summary-box {{
        background: #eaf4fd;
        border-left: 4px solid #3274A1;
        padding: 16px 20px;
        margin: 16px 0 24px 0;
        border-radius: 0 6px 6px 0;
    }}
    .summary-box p {{
        margin: 4px 0;
        font-size: 14px;
    }}
    .summary-box .highlight {{
        font-weight: 700;
        color: #1a6fb5;
    }}
    .metric-cards {{
        display: flex;
        gap: 16px;
        margin: 16px 0;
        flex-wrap: wrap;
    }}
    .metric-card {{
        flex: 1;
        min-width: 200px;
        background: #f0f7ff;
        border: 1px solid #c8ddf0;
        border-radius: 8px;
        padding: 16px;
        text-align: center;
    }}
    .metric-card .value {{
        font-size: 28px;
        font-weight: 700;
        color: #3274A1;
    }}
    .metric-card .label {{
        font-size: 12px;
        color: #666;
        margin-top: 4px;
    }}
    .metric-card.green {{ background: #f0faf0; border-color: #b5ddb5; }}
    .metric-card.green .value {{ color: #2CA02C; }}
    table {{
        width: 100%;
        border-collapse: collapse;
        margin: 8px 0 20px 0;
        font-size: 13px;
    }}
    th {{
        background: #f1f3f5;
        color: #333;
        font-weight: 600;
        padding: 8px 10px;
        text-align: left;
        border-bottom: 2px solid #dee2e6;
        white-space: nowrap;
    }}
    td {{
        padding: 6px 10px;
        border-bottom: 1px solid #eee;
        white-space: nowrap;
    }}
    tr:hover {{ background: #f8f9fa; }}
    .best-row {{ background: #fff8e1 !important; }}
    .best-cell {{ font-weight: 700; color: #d4760a; }}
    .table-caption {{
        font-size: 13px;
        font-weight: 600;
        color: #555;
        margin-top: 16px;
        margin-bottom: 4px;
    }}
    .plot-container {{
        margin: 16px 0;
        text-align: center;
    }}
    .plot-container img {{
        max-width: 100%;
        border: 1px solid #e0e0e0;
        border-radius: 4px;
    }}
    .findings {{
        margin: 12px 0;
    }}
    .findings li {{
        margin: 6px 0 6px 20px;
        font-size: 14px;
        line-height: 1.5;
    }}
    .findings li strong {{
        color: #2c3e50;
    }}
    .section-divider {{
        border: 0;
        border-top: 1px solid #e0e0e0;
        margin: 28px 0;
    }}
    .footer {{
        margin-top: 40px;
        padding-top: 16px;
        border-top: 1px solid #e0e0e0;
        font-size: 12px;
        color: #999;
    }}
    .sensitivity-grid {{
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 20px;
        margin: 12px 0;
    }}
    @media (max-width: 800px) {{
        .sensitivity-grid {{ grid-template-columns: 1fr; }}
        .metric-cards {{ flex-direction: column; }}
    }}
</style>
</head>
<body>
<div class="container">

<h1>IBD Spatial Integration: Method Optimization Report</h1>
<p class="date">Date: {report_date} &nbsp;|&nbsp; Dataset: M2 (CosMx 6k + Xenium 5K) &nbsp;|&nbsp; 164,202 cells, 2,552 genes, 8 samples (4 patients)</p>

<div class="summary-box">
    <p><span class="highlight">Best bio conservation:</span> scANVI A6 (3-layer, NB, gene-batch dispersion) achieves ct_broad ASW = 0.189 -- 2.5x higher than any original method.</p>
    <p><span class="highlight">Key insight:</span> Semi-supervised optimization dramatically improves celltype separation while maintaining strong batch correction (batch_score 0.903).</p>
    <p><span class="highlight">resolVI-SS B8</span> (per-sample batch key) achieves ct_broad ASW = 0.104 with spatial awareness -- second-best bio conservation.</p>
</div>

<div class="metric-cards">
    <div class="metric-card">
        <div class="value">0.189</div>
        <div class="label">Best ct_broad ASW<br>(scANVI A6)</div>
    </div>
    <div class="metric-card green">
        <div class="value">2.5x</div>
        <div class="label">Improvement over<br>vanilla scANVI (0.074)</div>
    </div>
    <div class="metric-card">
        <div class="value">0.903</div>
        <div class="label">Batch Score<br>(scANVI A6)</div>
    </div>
    <div class="metric-card green">
        <div class="value">10</div>
        <div class="label">Methods<br>Benchmarked</div>
    </div>
</div>

<hr class="section-divider">

<h2>1. All Methods Comparison</h2>

<div class="plot-container">
    <img src="data:image/png;base64,{bar_b64}" alt="Grouped bar chart of all methods">
</div>

{all_table}

<hr class="section-divider">

<h2>2. Batch-Bio Tradeoff</h2>

<div class="plot-container">
    <img src="data:image/png;base64,{scatter_b64}" alt="Batch vs Bio scatter plot with Pareto frontier">
</div>

<p style="font-size:13px; color:#555; margin-top:8px;">
    The Pareto frontier connects methods where no other method is strictly better on both axes.
    scANVI A6 and scVI define the current frontier: scVI achieves near-perfect batch mixing (0.992)
    with negligible bio signal (0.009), while scANVI A6 sacrifices modest batch correction for
    substantially better celltype separation (0.189).
</p>

<hr class="section-divider">

<h2>3. scANVI Hyperparameter Sweep</h2>
<p style="font-size:13px; color:#555;">11 configurations tested on A40 GPU. Winner: <strong>A6_3layer_genebatch</strong> (30 latent dims, 256 hidden, 3 layers, NB likelihood, gene-batch dispersion, cls_ratio=100).</p>

{scanvi_table}

<hr class="section-divider">

<h2>4. resolVI-SS Hyperparameter Sweep</h2>
<p style="font-size:13px; color:#555;">6 of 10 configurations completed (B3-B5, B9 failed due to API/memory issues). Winner: <strong>B8_per_sample</strong> (batch_key=sample instead of platform).</p>

{resolvi_table}

<hr class="section-divider">

<h2>5. Parameter Sensitivity Analysis</h2>

<div class="sensitivity-grid">
    <div>
        <h3>5.1 Classification Ratio</h3>
        <p style="font-size:13px; color:#555;">Higher cls_ratio forces more label supervision. Monotonic improvement in bio conservation at a mild batch cost.</p>
        {cls_table}
    </div>
    <div>
        <h3>5.2 Likelihood Function</h3>
        <p style="font-size:13px; color:#555;">NB outperforms ZINB for this dataset (spatial transcript counts are less zero-inflated than scRNA-seq).</p>
        {lik_table}
    </div>
</div>

<h3>5.3 Architecture Scaling</h3>
<p style="font-size:13px; color:#555;">Larger networks + deeper layers + gene-batch dispersion compound to give the best result. A6 combines all three improvements.</p>
{arch_table}

<h3>5.4 Summary of Parameter Effects</h3>
<ul class="findings">
    <li><strong>classification_ratio (50 to 200):</strong> ct_broad ASW increases 0.074 to 0.112 (+51%). Higher supervision weight helps, but returns diminish past 100 when combined with better architecture.</li>
    <li><strong>Architecture (20/128/2L to 30/256/3L):</strong> ct_broad ASW increases 0.074 to 0.189 (+155%). Larger capacity is the biggest single lever.</li>
    <li><strong>Likelihood (ZINB to NB):</strong> ct_broad ASW increases 0.134 to 0.148 (+10%). NB is modestly better for spatial data.</li>
    <li><strong>Dispersion (gene to gene-batch):</strong> A6 uses gene-batch dispersion -- allows per-platform variance, which helps cross-platform integration.</li>
    <li><strong>Label granularity (celltype_broad vs celltype):</strong> Broad labels work better (A2: 0.091 vs A9: 0.065). Fine-grained labels may introduce noise when some celltypes are rare.</li>
</ul>

<hr class="section-divider">

<h2>6. Key Findings</h2>
<ul class="findings">
    <li><strong>scANVI A6 is the clear winner</strong> for bio conservation (ct_broad ASW = 0.189), more than doubling vanilla scANVI (0.074) and 2.5x the next best original method (resolVI-SS at 0.044).</li>
    <li><strong>Batch-bio tradeoff is favorable:</strong> scANVI A6 retains batch_score = 0.903 (vs 0.992 for scVI), a modest cost for substantially better celltype structure.</li>
    <li><strong>Platform entropy confirms mixing:</strong> scANVI A6 achieves 0.576 platform entropy, second only to scVI (0.632), indicating good cross-platform mixing in clusters.</li>
    <li><strong>resolVI-SS B8 is the best spatial-aware option:</strong> Using sample-level batch correction (instead of platform) improved bio conservation to 0.104 with spatial prior regularization.</li>
    <li><strong>Unoptimized methods cluster near zero bio signal:</strong> scVI, Harmony, PCA, BBKNN, and Scanorama all have ct_broad ASW near or below zero, meaning batch correction overwhelms biology.</li>
    <li><strong>Three key improvements compound:</strong> (1) Larger architecture (30/256/3L), (2) NB likelihood, (3) gene-batch dispersion together produce the best integration.</li>
    <li><strong>Recommendation:</strong> Use scANVI A6 embedding for downstream analysis (clustering, celltype annotation, differential expression). Consider resolVI-SS B8 as a sensitivity check for spatial analyses.</li>
</ul>

<div class="footer">
    <p>Generated by generate_sweep_report.py | sc_tools | IBD Spatial Cross-Platform Integration Project</p>
</div>

</div>
</body>
</html>"""
    return html


# ============================================================
# Main
# ============================================================

def main():
    print("Generating IBD Spatial Integration Method Optimization Report...")

    # Generate plots
    print("  Creating grouped bar chart...")
    bar_fig = make_grouped_bar_chart(ALL_METHODS)
    bar_b64 = fig_to_base64(bar_fig)

    print("  Creating batch-bio scatter plot...")
    scatter_fig = make_scatter_plot(ALL_METHODS)
    scatter_b64 = fig_to_base64(scatter_fig)

    # Generate HTML
    print("  Generating HTML report...")
    html = generate_html(bar_b64, scatter_b64)

    # Write HTML
    html_path = FIGDIR / "ibd_sweep_report.html"
    html_path.write_text(html, encoding="utf-8")
    print(f"  Saved: {html_path}")

    # Write combined CSV
    csv_path = RESULTS_DIR / "all_methods_benchmark.csv"
    ALL_METHODS.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    print("Done.")


if __name__ == "__main__":
    main()
