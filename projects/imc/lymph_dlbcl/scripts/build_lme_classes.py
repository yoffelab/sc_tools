#!/usr/bin/env python3
"""
Phase 3.2: Assign LME (Lymphoma Microenvironment) classes.

Load TME cluster assignments from downloaded CSVs and map to 5 LME classes:
- Cold (35.1%)
- Stromal (21.3%)
- Cytotoxic (20.7%)
- CD206 Enriched (8.2%)
- T cell Regulated (14.6%)

Add obs['LME_class'] to both panel AnnData objects.

Usage:
    python scripts/build_lme_classes.py [--panel immune|stromal|both]

Input:
    data/downloaded/metadata/2.17.22.v2stroma_C8_clusters.csv (or TME_alt6.csv)
    data/downloaded/metadata/1.13.21.merged.abundance.csv
    results/adata.immune.annotated.p2.h5ad
    results/adata.stromal.annotated.p2.h5ad

Output:
    results/adata.immune.celltyped.p4.h5ad (with LME_class)
    results/adata.stromal.celltyped.p4.h5ad (with LME_class)
    metadata/lme_class_assignments.csv
"""

import argparse
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
METADATA_DIR = PROJECT_DIR / "metadata"

# LME class mapping from k=10 k-means to 5 classes
# From DLBCL_case_clustering.ipynb analysis
LME_CLASS_MAP = {
    0: "Cold",
    1: "CD206_Enriched",
    2: "Cytotoxic",
    3: "Stroma",
    4: "Stroma",
    5: "T_cell_Regulated",
    6: "T_cell_Regulated",
    7: "Cold",
    8: "CD206_Enriched",
    9: "Cold",
    10: "Cytotoxic",
}

# Display names
LME_DISPLAY_NAMES = {
    "Cold": "Cold",
    "CD206_Enriched": "CD206 Enriched",
    "Cytotoxic": "Cytotoxic",
    "Stroma": "Stromal",
    "T_cell_Regulated": "T cell Regulated",
}

# Expected prevalences from manuscript
LME_EXPECTED_PCT = {
    "Cold": 35.1,
    "Stromal": 21.3,
    "Cytotoxic": 20.7,
    "CD206 Enriched": 8.2,
    "T cell Regulated": 14.6,
}


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def load_tme_assignments(config: dict) -> pd.DataFrame | None:
    """Load TME cluster assignments from downloaded CSVs."""

    # Try primary: v2stroma_C8_clusters.csv
    tme_path = PROJECT_DIR / config["metadata"].get("tme_clusters", "")
    if tme_path.exists():
        logger.info(f"Loading TME clusters: {tme_path}")
        df = pd.read_csv(tme_path)
        logger.info(f"  Shape: {df.shape}, columns: {df.columns.tolist()}")
        return df

    # Try alternative: TME_alt6.csv
    alt_path = PROJECT_DIR / config["metadata"].get("tme_alt6", "")
    if alt_path.exists():
        logger.info(f"Loading TME alt6: {alt_path}")
        df = pd.read_csv(alt_path)
        logger.info(f"  Shape: {df.shape}, columns: {df.columns.tolist()}")
        return df

    # Try patient TME
    patient_path = PROJECT_DIR / config["metadata"].get("patient_tme", "")
    if patient_path.exists():
        logger.info(f"Loading patient TME: {patient_path}")
        df = pd.read_csv(patient_path)
        logger.info(f"  Shape: {df.shape}, columns: {df.columns.tolist()}")
        return df

    logger.warning("No TME assignment CSV found")
    return None


def assign_lme_from_csv(tme_df: pd.DataFrame) -> pd.DataFrame:
    """Map TME cluster assignments to LME classes."""

    # Find the cluster column
    cluster_col = None
    for col in tme_df.columns:
        if "cluster" in col.lower() or "tme" in col.lower() or "class" in col.lower():
            cluster_col = col
            break

    if cluster_col is None:
        # Use last numeric column
        for col in reversed(tme_df.columns.tolist()):
            if tme_df[col].dtype in [int, float, np.int64, np.float64]:
                cluster_col = col
                break

    if cluster_col is None:
        logger.error(f"No cluster column found in TME CSV. Columns: {tme_df.columns.tolist()}")
        return pd.DataFrame()

    logger.info(f"  Cluster column: '{cluster_col}'")
    logger.info(f"  Unique clusters: {sorted(tme_df[cluster_col].unique().tolist())}")

    # Find sample ID column
    sample_col = None
    for col in tme_df.columns:
        if col.lower() in ["dlc_code", "sample", "patient", "case_id", "id"]:
            sample_col = col
            break
    if sample_col is None:
        sample_col = tme_df.columns[0]

    logger.info(f"  Sample column: '{sample_col}'")

    # Map clusters to LME classes
    result = tme_df[[sample_col, cluster_col]].copy()
    result.columns = ["sample", "tme_cluster"]
    result["tme_cluster_int"] = pd.to_numeric(result["tme_cluster"], errors="coerce")

    # Apply mapping
    result["LME_class"] = result["tme_cluster_int"].map(LME_CLASS_MAP)
    result["LME_display"] = result["LME_class"].map(LME_DISPLAY_NAMES)

    # Report
    if "LME_display" in result.columns:
        dist = result["LME_display"].value_counts(normalize=True) * 100
        logger.info("  LME class distribution:")
        for cls, pct in dist.items():
            expected = LME_EXPECTED_PCT.get(cls, "?")
            logger.info(f"    {cls}: {pct:.1f}% (expected: {expected}%)")

    return result


def assign_lme_from_kmeans(config: dict) -> pd.DataFrame | None:
    """Fallback: compute LME classes from abundance data using k-means."""
    abundance_path = PROJECT_DIR / config["metadata"].get("abundance", "")
    if not abundance_path.exists():
        logger.warning(f"No abundance file: {abundance_path}")
        return None

    logger.info(f"Computing LME from abundance: {abundance_path}")
    abundance = pd.read_csv(abundance_path, index_col=0)
    logger.info(f"  Shape: {abundance.shape}")

    # Normalize to proportions
    row_sums = abundance.sum(axis=1)
    abundance_norm = abundance.div(row_sums, axis=0)

    # Scale
    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import KMeans

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(abundance_norm)

    # K-means with k=10
    kmeans = KMeans(n_clusters=10, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(X_scaled)

    result = pd.DataFrame({
        "sample": abundance.index,
        "tme_cluster": clusters,
    })

    # Map to 5 classes
    result["LME_class"] = result["tme_cluster"].map(LME_CLASS_MAP)
    result["LME_display"] = result["LME_class"].map(LME_DISPLAY_NAMES)

    return result


def main():
    parser = argparse.ArgumentParser(description="Assign LME classes to panel AnnData")
    parser.add_argument("--panel", choices=["immune", "stromal", "both"], default="both")
    parser.add_argument("--method", choices=["csv", "kmeans", "auto"], default="auto",
                        help="How to assign LME classes")
    args = parser.parse_args()

    config = load_config()
    METADATA_DIR.mkdir(parents=True, exist_ok=True)

    # Get LME assignments
    lme_assignments = None

    if args.method in ("csv", "auto"):
        tme_df = load_tme_assignments(config)
        if tme_df is not None:
            lme_assignments = assign_lme_from_csv(tme_df)

    if lme_assignments is None or lme_assignments.empty:
        if args.method in ("kmeans", "auto"):
            lme_assignments = assign_lme_from_kmeans(config)

    if lme_assignments is None or lme_assignments.empty:
        logger.error("Could not determine LME assignments. Need TME CSVs or abundance data.")
        return

    # Save LME assignments
    lme_path = METADATA_DIR / "lme_class_assignments.csv"
    lme_assignments.to_csv(lme_path, index=False)
    logger.info(f"LME assignments saved: {lme_path} ({len(lme_assignments)} samples)")

    # Attach to panel AnnData
    panels = ["immune", "stromal"] if args.panel == "both" else [args.panel]

    for panel in panels:
        # Try p2 first, then p1
        input_path = RESULTS_DIR / f"adata.{panel}.annotated.p2.h5ad"
        if not input_path.exists():
            input_path = RESULTS_DIR / f"adata.{panel}.raw.p1.h5ad"
        if not input_path.exists():
            logger.warning(f"No checkpoint for {panel}")
            continue

        output_path = RESULTS_DIR / f"adata.{panel}.celltyped.p4.h5ad"
        logger.info(f"Attaching LME to {panel}: {input_path}")
        adata = ad.read_h5ad(input_path)

        # Join LME_class via sample ID
        if "sample" in adata.obs.columns:
            sample_to_lme = dict(
                zip(
                    lme_assignments["sample"].astype(str),
                    lme_assignments["LME_display"],
                )
            )
            adata.obs["LME_class"] = adata.obs["sample"].astype(str).map(sample_to_lme)
            n_mapped = adata.obs["LME_class"].notna().sum()
            logger.info(f"  LME mapped: {n_mapped:,}/{adata.n_obs:,} cells")

            if n_mapped > 0:
                dist = adata.obs["LME_class"].value_counts(normalize=True) * 100
                for cls, pct in dist.items():
                    logger.info(f"    {cls}: {pct:.1f}%")
        else:
            logger.warning(f"  No 'sample' column in {panel} obs")

        # Clean dtypes
        for col in adata.obs.columns:
            if adata.obs[col].dtype == "category":
                adata.obs[col] = adata.obs[col].astype(str)

        logger.info(f"Saving: {output_path}")
        adata.write_h5ad(output_path)
        logger.info(f"  {adata.n_obs:,} cells, {len(adata.obs.columns)} obs columns")

    logger.info("Done.")


if __name__ == "__main__":
    main()
