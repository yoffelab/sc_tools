#!/usr/bin/env python3
"""
Phase 3.1: Validate cell types by checking marker expression patterns.

For each panel, compute mean marker expression per cell type and generate
a heatmap. Check expected patterns (CD20+ B cells, CD3+CD4+ T cells, etc.).

Usage:
    python scripts/validate_celltypes.py [--panel immune|stromal|both]

Input:
    results/adata.immune.raw.p1.h5ad (or .annotated.p2.h5ad)
    results/adata.stromal.raw.p1.h5ad (or .annotated.p2.h5ad)

Output:
    outputs/celltype_validation_report.md
    figures/QC/celltype_marker_heatmap_immune.pdf
    figures/QC/celltype_marker_heatmap_stromal.pdf
"""

import argparse
import logging
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
RESULTS_DIR = PROJECT_DIR / "results"
FIGURES_DIR = PROJECT_DIR / "figures" / "QC"
OUTPUT_DIR = PROJECT_DIR / "outputs"

# Expected marker-celltype associations (IMC markers)
EXPECTED_PATTERNS = {
    "immune": {
        "B cell": ["CD20", "PAX5", "CD79a"],
        "T cell CD4+": ["CD3", "CD4"],
        "T cell CD8+": ["CD3", "CD8"],
        "Macrophage": ["CD68", "CD163"],
        "Myeloid": ["CD11c", "HLA-DR", "CD11b"],
    },
    "stromal": {
        "B cell": ["CD20"],
        "Stroma": ["PDPN", "VIM", "FN"],
        "Endothelial": ["CD31", "vWF"],
        "Macrophage": ["CD68", "CD163", "CD206"],
    },
}


def compute_mean_expression(adata: ad.AnnData, groupby: str) -> pd.DataFrame:
    """Compute mean expression per group."""
    if groupby not in adata.obs.columns:
        logger.warning(f"Column '{groupby}' not in obs")
        return pd.DataFrame()

    groups = adata.obs[groupby].unique()
    if len(groups) > 50:
        logger.warning(f"Too many groups ({len(groups)}); skipping")
        return pd.DataFrame()

    # Use layers['raw'] if X is empty (p4 may have zeroed X)
    data_matrix = adata.layers.get("raw", adata.X) if adata.layers else adata.X

    mean_expr = pd.DataFrame(index=groups, columns=adata.var_names, dtype=float)
    for group in groups:
        mask = adata.obs[groupby] == group
        if hasattr(data_matrix, "toarray"):
            mean_expr.loc[group] = np.array(data_matrix[mask].toarray().mean(axis=0)).flatten()
        else:
            mean_expr.loc[group] = np.array(data_matrix[mask].mean(axis=0)).flatten()

    return mean_expr


def plot_heatmap(mean_expr: pd.DataFrame, panel: str, output_path: Path):
    """Plot mean marker expression heatmap."""
    if mean_expr.empty:
        return

    # Z-score across cell types for visualization
    z_scored = mean_expr.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=0)

    fig_height = max(6, len(mean_expr) * 0.4)
    fig_width = max(10, len(mean_expr.columns) * 0.3)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        z_scored,
        cmap="RdBu_r",
        center=0,
        vmin=-2,
        vmax=2,
        xticklabels=True,
        yticklabels=True,
        ax=ax,
        cbar_kws={"label": "z-score"},
    )
    ax.set_title(f"Mean Marker Expression per Cell Type ({panel} panel)")
    ax.set_xlabel("Markers")
    ax.set_ylabel("Cell Type")
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(fontsize=8)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info(f"  Saved heatmap: {output_path}")


def validate_expected_patterns(
    mean_expr: pd.DataFrame, panel: str
) -> list[dict]:
    """Check if expected marker-celltype patterns hold."""
    results = []
    patterns = EXPECTED_PATTERNS.get(panel, {})

    for celltype, markers in patterns.items():
        for marker in markers:
            # Find matching marker in var names (case-insensitive partial match)
            matched_markers = [
                m for m in mean_expr.columns
                if marker.lower() in m.lower()
            ]
            if not matched_markers:
                results.append({
                    "celltype": celltype,
                    "marker": marker,
                    "status": "MARKER_NOT_FOUND",
                    "details": f"Marker {marker} not in var_names",
                })
                continue

            # Find matching celltype
            matched_celltypes = [
                ct for ct in mean_expr.index
                if celltype.lower() in str(ct).lower()
            ]
            if not matched_celltypes:
                results.append({
                    "celltype": celltype,
                    "marker": marker,
                    "status": "CELLTYPE_NOT_FOUND",
                    "details": f"Cell type {celltype} not in groups",
                })
                continue

            for m in matched_markers:
                for ct in matched_celltypes:
                    expr = mean_expr.loc[ct, m]
                    # Check if this celltype has above-median expression
                    median_expr = mean_expr[m].median()
                    if expr > median_expr:
                        results.append({
                            "celltype": ct,
                            "marker": m,
                            "status": "PASS",
                            "details": f"expr={expr:.3f} > median={median_expr:.3f}",
                        })
                    else:
                        results.append({
                            "celltype": ct,
                            "marker": m,
                            "status": "WARN",
                            "details": f"expr={expr:.3f} <= median={median_expr:.3f}",
                        })

    return results


def generate_report(all_results: dict) -> str:
    """Generate validation report."""
    lines = [
        "# Cell Type Validation Report\n",
        "## Summary\n",
    ]

    for panel, data in all_results.items():
        lines.append(f"### {panel.title()} Panel\n")

        if "mean_expr" in data and not data["mean_expr"].empty:
            me = data["mean_expr"]
            lines.append(f"- Cell types: {len(me)}")
            lines.append(f"- Markers: {len(me.columns)}")
            lines.append(f"- Cell type names: {sorted(me.index.tolist())}")

        if "pattern_checks" in data:
            checks = data["pattern_checks"]
            n_pass = sum(1 for c in checks if c["status"] == "PASS")
            n_warn = sum(1 for c in checks if c["status"] == "WARN")
            n_missing = sum(1 for c in checks if "NOT_FOUND" in c["status"])
            lines.append(f"\n**Pattern checks:** {n_pass} pass, {n_warn} warnings, {n_missing} not found\n")

            lines.append("| Cell Type | Marker | Status | Details |")
            lines.append("|-----------|--------|--------|---------|")
            for c in checks:
                lines.append(f"| {c['celltype']} | {c['marker']} | {c['status']} | {c['details']} |")

        lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Validate cell types")
    parser.add_argument("--panel", choices=["immune", "stromal", "both"], default="both")
    args = parser.parse_args()

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    panels = ["immune", "stromal"] if args.panel == "both" else [args.panel]
    all_results = {}

    for panel in panels:
        # Try annotated first, then raw
        input_path = RESULTS_DIR / f"adata.{panel}.annotated.p2.h5ad"
        if not input_path.exists():
            input_path = RESULTS_DIR / f"adata.{panel}.raw.p1.h5ad"
        if not input_path.exists():
            logger.warning(f"No checkpoint found for {panel} panel")
            continue

        logger.info(f"=== Validating {panel} panel: {input_path} ===")
        adata = ad.read_h5ad(input_path)

        # Find best celltype column
        celltype_col = None
        for col in ["celltype", "celltype_broad", "labels", "meta"]:
            if col in adata.obs.columns:
                n_unique = adata.obs[col].nunique()
                if 2 <= n_unique <= 50:
                    celltype_col = col
                    break

        if celltype_col is None:
            logger.warning(f"  No suitable cell type column found")
            all_results[panel] = {"error": "No cell type column"}
            continue

        logger.info(f"  Using cell type column: '{celltype_col}' ({adata.obs[celltype_col].nunique()} types)")

        # Compute mean expression
        mean_expr = compute_mean_expression(adata, celltype_col)

        # Plot heatmap
        heatmap_path = FIGURES_DIR / f"celltype_marker_heatmap_{panel}.pdf"
        plot_heatmap(mean_expr, panel, heatmap_path)

        # Validate patterns
        pattern_checks = validate_expected_patterns(mean_expr, panel)

        all_results[panel] = {
            "mean_expr": mean_expr,
            "pattern_checks": pattern_checks,
        }

    # Generate report
    report = generate_report(all_results)
    report_path = OUTPUT_DIR / "celltype_validation_report.md"
    report_path.write_text(report)
    logger.info(f"Report: {report_path}")


if __name__ == "__main__":
    main()
