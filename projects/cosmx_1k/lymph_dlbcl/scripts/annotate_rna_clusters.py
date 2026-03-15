#!/usr/bin/env python3
"""Annotate CosMx RNA clusters (a-j) with cell type labels.

Scores each cluster against lineage marker gene sets from the 1000-plex panel,
then assigns cell types based on highest-scoring signature per cluster.

Input:  results/cosmx_rna.h5ad  (464K cells, 1000 genes, clusters a-j)
Output: results/cosmx_rna_annotated.h5ad  (+ obs['celltype'])
"""

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"

# Lineage marker gene sets for cell type annotation
# Curated from manuscript + known CosMx 1000-plex panel markers
MARKER_GENES = {
    "B cell": ["MS4A1", "CD19", "PAX5", "BACH2", "BCL6", "CD79A", "CD79B", "BANK1"],
    "CD4 T cell": ["CD4", "IL7R", "FOXP3", "IL2RA", "CXCR5", "TCF7", "LEF1"],
    "CD8 T cell": ["CD8A", "CD8B", "GZMB", "PRF1", "GZMA", "GZMK", "LAG3", "NKG7"],
    "NK": ["NCAM1", "KLRD1", "NKG7", "GNLY", "KLRB1", "KLRC1"],
    "Macrophage": ["CD68", "CD163", "MRC1", "CSF1R", "CD14", "MSR1", "MARCO"],
    "DC": ["ITGAX", "CLEC9A", "CLEC4C", "FLT3", "IRF8", "BATF3"],
    "Fibroblast": ["PDPN", "FAP", "COL1A1", "COL1A2", "ACTA2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "ERG", "FLT1"],
    "Neutrophil": ["FCGR3B", "CSF3R", "CXCR2", "S100A8", "S100A9"],
    "Plasma cell": ["SDC1", "XBP1", "PRDM1", "IRF4", "TNFRSF17"],
}


def score_clusters(adata, cluster_col="cluster"):
    """Score each cluster against marker gene sets. Returns DataFrame of z-scores."""
    clusters = adata.obs[cluster_col].unique()
    genes_in_data = set(adata.var_names)

    scores = {}
    for celltype, markers in MARKER_GENES.items():
        present = [g for g in markers if g in genes_in_data]
        if not present:
            logger.warning("No markers found for %s", celltype)
            scores[celltype] = dict.fromkeys(clusters, 0.0)
            continue

        logger.info(
            "%s: %d/%d markers present (%s)",
            celltype,
            len(present),
            len(markers),
            ", ".join(present),
        )

        # Mean expression per cluster (use raw counts if available)
        cluster_means = {}
        for c in clusters:
            mask = adata.obs[cluster_col] == c
            subset = adata[mask, present].X
            if hasattr(subset, "toarray"):
                subset = subset.toarray()
            cluster_means[c] = np.mean(subset, axis=0).mean()

        scores[celltype] = cluster_means

    score_df = pd.DataFrame(scores, index=clusters)
    # Z-score per cell type (across clusters) for comparability
    score_z = (score_df - score_df.mean()) / score_df.std().replace(0, 1)
    return score_df, score_z


def assign_celltypes(score_z, cluster_col="cluster"):
    """Assign cell type to each cluster based on highest z-score."""
    assignments = {}
    for cluster in score_z.index:
        best = score_z.loc[cluster].idxmax()
        best_score = score_z.loc[cluster].max()
        assignments[cluster] = best
        logger.info(
            "Cluster %s -> %s (z=%.2f)",
            cluster,
            best,
            best_score,
        )
    return assignments


def main():
    rna_path = RESULTS_DIR / "cosmx_rna.h5ad"
    if not rna_path.exists():
        logger.error("cosmx_rna.h5ad not found at %s", rna_path)
        sys.exit(1)

    logger.info("Loading %s", rna_path)
    adata = ad.read_h5ad(rna_path)
    logger.info("Shape: %s", adata.shape)
    logger.info("Clusters: %s", sorted(adata.obs["cluster"].unique()))

    # Normalize for scoring (log1p of counts)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Score clusters
    score_df, score_z = score_clusters(adata, cluster_col="cluster")
    logger.info("\nCluster scores (z-scored):\n%s", score_z.round(2).to_string())

    # Save score table
    score_z.to_csv(RESULTS_DIR / "cosmx_rna_cluster_scores.csv")

    # Assign cell types
    assignments = assign_celltypes(score_z)

    # Map to adata
    adata.obs["celltype"] = adata.obs["cluster"].map(assignments).astype("category")
    logger.info("\nCell type distribution:")
    for ct, n in adata.obs["celltype"].value_counts().items():
        logger.info("  %s: %d (%.1f%%)", ct, n, 100 * n / len(adata))

    # Save
    out_path = RESULTS_DIR / "cosmx_rna_annotated.h5ad"
    logger.info("Saving %s", out_path)
    adata.write_h5ad(out_path)

    # Also save assignment map
    assign_df = pd.DataFrame(
        [(k, v) for k, v in assignments.items()],
        columns=["cluster", "celltype"],
    )
    assign_df.to_csv(RESULTS_DIR / "cosmx_rna_cluster_assignments.csv", index=False)
    logger.info("Done.")


if __name__ == "__main__":
    main()
