#!/usr/bin/env python3
"""
Phase 1.4: Build spatial community AnnData object.

Load SO_k30_community_cluster.h5ad, map 20 communities to descriptive names
from manuscript, join with stromal expression + cell types.

Usage:
    python scripts/build_spatial_adata.py

Input:
    results/seurat_converted/stroma_spatial/SO_k30_community_cluster.h5ad
    results/adata.stromal.annotated.p2.h5ad (optional, for cell type enrichment)

Output:
    results/adata.stromal.spatial.communities.h5ad
"""

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
SEURAT_DIR = RESULTS_DIR / "seurat_converted"


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def main():
    config = load_config()

    # Load spatial community object
    spatial_path = PROJECT_DIR / config["objects"]["spatial_communities"]
    logger.info(f"Loading spatial community object: {spatial_path}")
    adata = ad.read_h5ad(spatial_path)
    logger.info(f"  Shape: {adata.shape}")
    logger.info(f"  obs columns: {sorted(adata.obs.columns.tolist())}")

    # Report on community and cluster columns
    for col in ["community", "cluster", "TME", "major_group", "recluster"]:
        if col in adata.obs.columns:
            vals = adata.obs[col].unique()
            logger.info(f"  {col}: {len(vals)} unique values")
            if len(vals) <= 30:
                logger.info(f"    Values: {sorted(vals.tolist())}")

    # Store raw
    adata.layers["raw"] = adata.X.copy()

    # Map community IDs to descriptive names (placeholder - refine after manuscript review)
    # The manuscript defines 20 spatial communities; names come from marker enrichment
    if "community" in adata.obs.columns:
        adata.obs["community_id"] = adata.obs["community"].astype(str)
        logger.info(f"  Communities: {adata.obs['community_id'].nunique()} unique")

    # Derive sample from orig.ident
    if "orig.ident" in adata.obs.columns:
        adata.obs["sample"] = adata.obs["orig.ident"].astype(str)

    # Try to join cell types from stromal P2 if available
    p2_path = RESULTS_DIR / "adata.stromal.annotated.p2.h5ad"
    if p2_path.exists():
        logger.info(f"Joining cell types from {p2_path}")
        p2 = ad.read_h5ad(p2_path, backed="r")
        for col in ["celltype", "celltype_broad"]:
            if col in p2.obs.columns:
                col_map = p2.obs[col].to_dict()
                adata.obs[col] = adata.obs.index.map(
                    lambda x, cm=col_map: cm.get(x, "Unknown")
                )
                n_mapped = sum(adata.obs[col] != "Unknown")
                logger.info(f"  Mapped {col}: {n_mapped:,}/{adata.n_obs:,}")
        p2.file.close()

    # Clean dtypes
    for col in adata.obs.columns:
        if adata.obs[col].dtype == "category":
            adata.obs[col] = adata.obs[col].astype(str)

    # Save
    out_path = RESULTS_DIR / "adata.stromal.spatial.communities.h5ad"
    logger.info(f"Saving: {out_path}")
    adata.write_h5ad(out_path)
    logger.info(f"  {adata.n_obs:,} cells, {adata.n_vars} markers")
    logger.info("Done.")


if __name__ == "__main__":
    main()
