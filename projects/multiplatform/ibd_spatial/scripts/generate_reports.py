#!/usr/bin/env python3
"""Generate integration benchmark HTML reports for M0 and M1.

Uses sc_tools.bm.integration.compare_integrations() for proper scib-style
metrics, then sc_tools.bm.report.generate_integration_report() for HTML output.

Usage:
    python generate_reports.py [--milestone m0|m1|all]
"""

import os
import sys
import argparse

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout.reconfigure(line_buffering=True)

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")


def generate_m0_report():
    """Generate M0 integration report (noseg vs withseg)."""
    print("=== M0 Report ===")
    adata_path = WORKDIR / "results" / "m0_benchmark" / "adata.m0.h5ad"
    if not adata_path.exists():
        print(f"  SKIP: {adata_path} not found")
        return

    adata = ad.read_h5ad(adata_path)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # Map method names to embedding keys
    embeddings = {}
    if "X_pca" in adata.obsm:
        embeddings["PCA (unintegrated)"] = "X_pca"
    if "X_harmony" in adata.obsm:
        embeddings["Harmony"] = "X_harmony"
    if "X_scvi" in adata.obsm:
        embeddings["scVI"] = "X_scvi"

    batch_key = "panel_variant"
    celltype_key = None  # noseg/withseg have no celltype annotations

    from sc_tools.bm.integration import compare_integrations

    print("  Computing metrics...")
    comparison_df = compare_integrations(
        adata,
        embeddings,
        batch_key=batch_key,
        celltype_key=celltype_key,
        include_unintegrated=False,
    )
    print(comparison_df.to_string(index=False))

    # Generate UMAP embeddings dict for the report
    umap_embeddings = {}
    for name, key in embeddings.items():
        umap_key = f"X_umap_{key.split('_', 1)[1] if '_' in key else key}"
        if umap_key in adata.obsm:
            umap_embeddings[name] = umap_key

    from sc_tools.bm.report import generate_integration_report

    output_path = WORKDIR / "figures" / "QC" / f"m0_integration_report.html"
    generate_integration_report(
        comparison_df,
        adata=adata,
        embeddings=umap_embeddings if umap_embeddings else None,
        batch_key=batch_key,
        celltype_key="patient_id",
        color_by="patient_id",
        output_path=output_path,
        title="M0: Xenium noseg vs withseg Integration Benchmark",
    )
    print(f"  Report saved: {output_path}")
    return comparison_df


def generate_m1_report():
    """Generate M1 integration report (CosMx 1k vs Xenium MT)."""
    print("\n=== M1 Report ===")
    adata_path = WORKDIR / "results" / "m1_benchmark" / "adata.m1.h5ad"
    if not adata_path.exists():
        print(f"  SKIP: {adata_path} not found")
        return

    adata = ad.read_h5ad(adata_path)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # Subsample for faster metric computation
    MAX_CELLS = 50000
    if adata.n_obs > MAX_CELLS:
        np.random.seed(42)
        idx = []
        for platform in adata.obs["platform"].unique():
            pmask = adata.obs["platform"] == platform
            pidx = np.where(pmask)[0]
            n_take = min(len(pidx), MAX_CELLS // 2)
            idx.extend(np.random.choice(pidx, n_take, replace=False))
        idx = sorted(idx)
        adata_sub = adata[idx].copy()
        print(f"  Subsampled to {adata_sub.n_obs} cells for metrics")
    else:
        adata_sub = adata

    embeddings = {}
    if "X_pca" in adata_sub.obsm:
        embeddings["PCA (unintegrated)"] = "X_pca"
    if "X_harmony" in adata_sub.obsm:
        embeddings["Harmony"] = "X_harmony"
    if "X_scvi" in adata_sub.obsm:
        embeddings["scVI"] = "X_scvi"

    batch_key = "platform"
    celltype_key = "celltype_broad"  # use broad celltypes (more reliable with 119 genes)

    from sc_tools.bm.integration import compare_integrations

    print("  Computing metrics...")
    comparison_df = compare_integrations(
        adata_sub,
        embeddings,
        batch_key=batch_key,
        celltype_key=celltype_key,
        include_unintegrated=False,
    )
    print(comparison_df.to_string(index=False))

    # For UMAP visualization, use full adata but subsample for the plot
    umap_embeddings = {}
    for name, key in embeddings.items():
        umap_key = f"X_umap_{key.split('_', 1)[1] if '_' in key else key}"
        if umap_key in adata_sub.obsm:
            umap_embeddings[name] = umap_key

    from sc_tools.bm.report import generate_integration_report

    output_path = WORKDIR / "figures" / "QC" / f"m1_integration_report.html"
    generate_integration_report(
        comparison_df,
        adata=adata_sub,
        embeddings=umap_embeddings if umap_embeddings else None,
        batch_key=batch_key,
        celltype_key=celltype_key,
        color_by="platform",
        output_path=output_path,
        title="M1: CosMx 1k vs Xenium MT Cross-Platform Integration",
    )
    print(f"  Report saved: {output_path}")

    # Also generate a celltype-colored version
    output_path2 = WORKDIR / "figures" / "QC" / f"m1_integration_report_celltype.html"
    generate_integration_report(
        comparison_df,
        adata=adata_sub,
        embeddings=umap_embeddings if umap_embeddings else None,
        batch_key=batch_key,
        celltype_key=celltype_key,
        color_by=celltype_key,
        output_path=output_path2,
        title="M1: CosMx 1k vs Xenium MT (colored by cell type)",
    )
    print(f"  Report saved: {output_path2}")

    return comparison_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--milestone", default="all", choices=["m0", "m1", "all"])
    args = parser.parse_args()

    if args.milestone in ("m0", "all"):
        generate_m0_report()
    if args.milestone in ("m1", "all"):
        generate_m1_report()

    print("\nDone!")
