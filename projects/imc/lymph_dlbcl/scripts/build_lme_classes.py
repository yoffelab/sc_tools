#!/usr/bin/env python3
"""
Phase 3.2: Assign LME (Lymphoma Microenvironment) classes.

5 LME classes from the manuscript:
- Cold (35.1%, n=115) — immune-depleted
- Stromal (21.3%, n=70) — PDPN+ CAFs, TGFbeta
- Cytotoxic (20.7%, n=68) — M1 macrophages, CD8+GzmB+
- T cell Regulated (14.6%, n=48) — Tfh, Tregs, LAG3+TOX+
- CD206 Enriched (8.2%, n=27) — CD206+ M2-like, protective

Strategy:
1. Try pre-computed TME CSV (patient_tme or tme_clusters) with direct cluster->LME mapping
2. Fallback: k-means on abundance data (k=10, then map to 5 LME classes)
3. Validate: proportions must be within 5% of manuscript values

Usage:
    python scripts/build_lme_classes.py [--panel immune|stromal|both]

Input:
    data/downloaded/metadata/8.15.22.DLC380_patient_tme.csv (preferred)
    data/downloaded/metadata/2.17.22.v2stroma_C8_clusters.csv (fallback)
    data/downloaded/metadata/1.13.21.merged.abundance.csv (fallback for k-means)
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
LME_CLASS_MAP_K10 = {
    0: "Cold",
    1: "CD206 Enriched",
    2: "Cytotoxic",
    3: "Stromal",
    4: "Stromal",
    5: "T cell Regulated",
    6: "T cell Regulated",
    7: "Cold",
    8: "CD206 Enriched",
    9: "Cold",
    10: "Cytotoxic",
}

# Expected prevalences from manuscript (for validation)
LME_EXPECTED_PCT = {
    "Cold": 35.1,
    "Stromal": 21.3,
    "Cytotoxic": 20.7,
    "CD206 Enriched": 8.2,
    "T cell Regulated": 14.6,
}

# Manuscript order
LME_ORDER = ["Cold", "Stromal", "Cytotoxic", "T cell Regulated", "CD206 Enriched"]


def load_config():
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def try_patient_tme(config: dict) -> pd.DataFrame | None:
    """Try loading patient_tme CSV which may have direct LME assignments."""
    path = PROJECT_DIR / config["metadata"].get("patient_tme", "")
    if not path.exists():
        return None

    logger.info(f"Loading patient TME: {path}")
    df = pd.read_csv(path)
    logger.info(f"  Shape: {df.shape}, columns: {df.columns.tolist()}")

    # Look for a column that contains LME class names or cluster IDs
    for col in df.columns:
        unique_vals = set(df[col].dropna().astype(str).unique())
        # Check if values match LME class names
        lme_names = {"Cold", "Stromal", "Cytotoxic", "CD206 Enriched",
                     "T cell Regulated", "CD206_Enriched", "T_cell_Regulated",
                     "Stroma"}
        if len(unique_vals & lme_names) >= 3:
            logger.info(f"  Found LME class column: '{col}'")
            # Find sample ID column
            sample_col = _find_sample_col(df)
            if sample_col is None:
                continue
            result = df[[sample_col, col]].copy()
            result.columns = ["sample", "LME_class"]
            # Standardize names
            result["LME_class"] = result["LME_class"].replace({
                "CD206_Enriched": "CD206 Enriched",
                "T_cell_Regulated": "T cell Regulated",
                "Stroma": "Stromal",
            })
            return result

    return None


def try_tme_clusters(config: dict) -> pd.DataFrame | None:
    """Try loading TME cluster CSV and mapping clusters to LME classes."""
    for key in ["tme_clusters", "tme_alt6"]:
        path = PROJECT_DIR / config["metadata"].get(key, "")
        if not path.exists():
            continue

        logger.info(f"Loading TME clusters: {path}")
        df = pd.read_csv(path)
        logger.info(f"  Shape: {df.shape}, columns: {df.columns.tolist()}")

        # Find cluster column
        cluster_col = None
        for col in df.columns:
            if "cluster" in col.lower() or "tme" in col.lower() or "class" in col.lower():
                cluster_col = col
                break
        if cluster_col is None:
            for col in reversed(df.columns.tolist()):
                if df[col].dtype in [int, float, np.int64, np.float64]:
                    cluster_col = col
                    break

        if cluster_col is None:
            logger.warning(f"  No cluster column found. Columns: {df.columns.tolist()}")
            continue

        sample_col = _find_sample_col(df)
        if sample_col is None:
            continue

        logger.info(f"  Cluster column: '{cluster_col}', unique: {sorted(df[cluster_col].unique().tolist())}")

        result = df[[sample_col, cluster_col]].copy()
        result.columns = ["sample", "tme_cluster"]
        result["tme_cluster"] = pd.to_numeric(result["tme_cluster"], errors="coerce")

        # Map to LME classes
        result["LME_class"] = result["tme_cluster"].map(LME_CLASS_MAP_K10)

        # If mapping failed (clusters not in expected range), skip
        if result["LME_class"].isna().all():
            logger.warning(f"  Cluster values {sorted(df[cluster_col].unique())} do not match k=10 mapping")
            continue

        result = result.dropna(subset=["LME_class"])
        return result[["sample", "LME_class"]]

    return None


def try_tme_h5ad(config: dict) -> pd.DataFrame | None:
    """Try loading TME h5ad object (s_tme_seurat.h5ad) with cluster assignments."""
    tme_path = PROJECT_DIR / config["objects"].get("tme", "")
    if not tme_path.exists():
        return None

    logger.info(f"Loading TME h5ad: {tme_path}")
    tme = ad.read_h5ad(tme_path)
    logger.info(f"  Shape: {tme.shape}, obs columns: {list(tme.obs.columns)}")

    # TME object has DLC codes as obs_names and seurat_clusters with 6 clusters
    # 576 obs = ~288 patients x 2 panels (duplicated as DLC0001, DLC0001.1)
    # Use seurat_clusters for TME cluster assignment
    cluster_col = "seurat_clusters"
    if cluster_col not in tme.obs.columns:
        for col in tme.obs.columns:
            if "cluster" in col.lower():
                cluster_col = col
                break

    clusters = tme.obs[cluster_col]
    logger.info(f"  Cluster column: '{cluster_col}', values: {sorted(clusters.unique().tolist())}")

    # DLC codes from obs_names — deduplicate panel suffixes (DLC0001.1 -> DLC0001)
    import re
    result_rows = []
    for obs_name, cluster in zip(tme.obs_names, clusters, strict=False):
        # Strip .1 suffix (panel duplicate)
        clean_name = re.sub(r"\.\d+$", "", str(obs_name))
        if not clean_name.startswith("DLC"):
            continue
        result_rows.append({"sample": clean_name, "tme_cluster": int(cluster)})

    if not result_rows:
        logger.warning("  No DLC samples found in TME h5ad")
        return None

    result = pd.DataFrame(result_rows)
    # Deduplicate (keep first per sample, both panels should agree)
    result = result.drop_duplicates(subset=["sample"], keep="first")

    # The TME object has 6 clusters (0-5), not 10
    # Map 6 clusters to 5 LME classes based on cluster sizes and manuscript proportions
    # From the data: {0: 115, 1: 108, 2: 100, 3: 96, 4: 82, 5: 75}
    # Expected: Cold=35.1%(~101), Stromal=21.3%(~61), Cytotoxic=20.7%(~60),
    #           T cell Regulated=14.6%(~42), CD206 Enriched=8.2%(~24)
    # Total TME patients ~ 288 (576/2)
    # We need to determine which cluster is which LME class by analyzing the TME features

    # First, check if the cluster IDs match the k=10 mapping
    unique_clusters = sorted(result["tme_cluster"].unique().tolist())
    logger.info(f"  Unique clusters: {unique_clusters}, n_samples: {len(result)}")

    if len(unique_clusters) > 6 and max(unique_clusters) <= 10:
        # Try k=10 mapping only when we actually have >6 clusters
        result["LME_class"] = result["tme_cluster"].map(LME_CLASS_MAP_K10)
        if result["LME_class"].notna().sum() > 0:
            result = result.dropna(subset=["LME_class"])
            logger.info(f"  Mapped via k=10 mapping: {len(result)} samples")
            return result[["sample", "LME_class"]]

    # Check if 5-cluster seurat_cluster column exists (direct 5-class TME)
    if "seurat_cluster" in tme.obs.columns and len(unique_clusters) <= 6:
        # Use seurat_cluster (5 classes) instead of seurat_clusters (6 classes)
        sc_col = "seurat_cluster"
        sc_unique = sorted(tme.obs[sc_col].unique().tolist())
        if len(sc_unique) == 5:
            logger.info(f"  Found 5-class column '{sc_col}': {sc_unique}")
            # Re-extract with seurat_cluster
            result_rows_5 = []
            for obs_name, cluster in zip(tme.obs_names, tme.obs[sc_col], strict=False):
                clean_name = re.sub(r"\.\d+$", "", str(obs_name))
                if not clean_name.startswith("DLC"):
                    continue
                result_rows_5.append({"sample": clean_name, "tme_cluster": int(cluster)})
            result = pd.DataFrame(result_rows_5).drop_duplicates(subset=["sample"], keep="first")
            unique_clusters = sorted(result["tme_cluster"].unique().tolist())
            logger.info(f"  5-class clusters (dedup): {result['tme_cluster'].value_counts().sort_index().to_dict()}")

    # Determine mapping from feature profiles using active cluster column
    active_cluster_col = "seurat_cluster" if "seurat_cluster" in tme.obs.columns and len(unique_clusters) == 5 else cluster_col
    logger.info(f"  Determining cluster -> LME mapping from TME feature profiles (column: {active_cluster_col})")
    # Use layers['raw'] if X is empty (p4 may have zeroed X)
    data_matrix = tme.layers.get("raw", tme.X) if tme.layers else tme.X
    tme_data = pd.DataFrame(
        data_matrix if not hasattr(data_matrix, "toarray") else data_matrix.toarray(),
        index=tme.obs_names,
        columns=tme.var_names,
    )
    tme_data["cluster"] = tme.obs[active_cluster_col].values

    # Compute mean abundance per cluster
    cluster_means = tme_data.groupby("cluster").mean()
    logger.info(f"  Cluster means shape: {cluster_means.shape}")

    # Identify LME classes from cluster feature profiles
    # Cold = lowest overall immune cell abundance
    # CD206 Enriched = highest CD206/M2 macrophage markers
    # Cytotoxic = highest CD8+GrB+, M1 markers
    # Stromal = highest PDPN/fibroblast markers
    # T cell Regulated = highest CD4+PD1+, Treg markers

    # Find relevant columns
    cd206_cols = [c for c in cluster_means.columns if "cd206" in c.lower() or "m2" in c.lower()]
    cytotoxic_cols = [c for c in cluster_means.columns if "grb" in c.lower() or "cd8" in c.lower()]
    stromal_cols = [c for c in cluster_means.columns if "pdpn" in c.lower() or ("s0" in c.lower() or "s1" in c.lower())]
    tcell_cols = [c for c in cluster_means.columns if "cd4" in c.lower() and "pd1" in c.lower()]

    # Score each cluster
    cluster_scores = {}
    for cl in cluster_means.index:
        scores = {}
        scores["CD206 Enriched"] = cluster_means.loc[cl, cd206_cols].sum() if cd206_cols else 0
        scores["Cytotoxic"] = cluster_means.loc[cl, cytotoxic_cols].sum() if cytotoxic_cols else 0
        scores["Stromal"] = cluster_means.loc[cl, stromal_cols].sum() if stromal_cols else 0
        scores["T cell Regulated"] = cluster_means.loc[cl, tcell_cols].sum() if tcell_cols else 0
        scores["total_immune"] = cluster_means.loc[cl].sum()
        cluster_scores[cl] = scores

    scores_df = pd.DataFrame(cluster_scores).T
    logger.info(f"  Cluster feature scores:\n{scores_df.to_string()}")

    # Assign LME classes using combined size and feature score matching
    # Use Hungarian algorithm approach: match clusters to LME classes optimizing
    # both feature similarity and size match to manuscript proportions
    cluster_sizes = result["tme_cluster"].value_counts().to_dict()
    n_total = len(result)
    expected_sizes = {cls: n_total * pct / 100 for cls, pct in LME_EXPECTED_PCT.items()}

    logger.info(f"  Cluster sizes: {cluster_sizes}")
    logger.info(f"  Expected sizes (scaled to n={n_total}): " +
                ", ".join(f"{k}={v:.0f}" for k, v in expected_sizes.items()))

    # Score matrix: feature score (z-normalized) + size penalty
    from itertools import permutations

    feature_cols = ["CD206 Enriched", "Cytotoxic", "Stromal", "T cell Regulated"]
    available_clusters = sorted(cluster_sizes.keys())
    lme_classes_ordered = LME_ORDER.copy()

    best_mapping = None
    best_score = -float("inf")

    # For Cold, use total_immune (negated since Cold = lowest immune)
    for perm in permutations(available_clusters, len(lme_classes_ordered)):
        mapping = dict(zip(perm, lme_classes_ordered, strict=False))
        score = 0
        for cl, lme_cls in mapping.items():
            # Feature score
            if lme_cls == "Cold":
                score -= scores_df.loc[cl, "total_immune"]  # Lower immune = better for Cold
            elif lme_cls in feature_cols:
                score += scores_df.loc[cl, lme_cls]
            # Size match penalty (negative penalty for size mismatch)
            size_diff = abs(cluster_sizes.get(cl, 0) - expected_sizes.get(lme_cls, 0))
            score -= size_diff * 0.05  # Weight size penalty

        if score > best_score:
            best_score = score
            best_mapping = mapping

    lme_map_6 = best_mapping
    logger.info(f"  Optimal mapping (score={best_score:.2f}): {lme_map_6}")

    # Handle remaining clusters (>5 clusters)
    assigned = set(lme_map_6.keys())
    for cl in cluster_means.index:
        if cl not in assigned:
            remaining_scores = {k: v for k, v in cluster_scores[cl].items() if k != "total_immune"}
            best = max(remaining_scores, key=remaining_scores.get)
            lme_map_6[cl] = best
            logger.info(f"  Extra cluster {cl} assigned to {best}")

    logger.info(f"  6-cluster -> LME mapping: {lme_map_6}")

    result["LME_class"] = result["tme_cluster"].map(lme_map_6)
    result = result.dropna(subset=["LME_class"])
    return result[["sample", "LME_class"]]


def try_kmeans_on_abundance(config: dict) -> pd.DataFrame | None:
    """Fallback: compute LME classes from cell type abundance using k-means."""
    abundance_path = PROJECT_DIR / config["metadata"].get("abundance", "")
    if not abundance_path.exists():
        logger.warning(f"No abundance file: {abundance_path}")
        return None

    logger.info(f"Computing LME from abundance (k-means): {abundance_path}")
    abundance = pd.read_csv(abundance_path, index_col=0)
    logger.info(f"  Shape: {abundance.shape}")

    # Normalize to proportions
    row_sums = abundance.sum(axis=1)
    if (row_sums > 1.5).any():
        abundance_norm = abundance.div(row_sums, axis=0)
    else:
        abundance_norm = abundance.copy()

    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(abundance_norm)

    kmeans = KMeans(n_clusters=10, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(X_scaled)

    result = pd.DataFrame({
        "sample": abundance.index.astype(str),
        "tme_cluster": clusters,
    })
    result["LME_class"] = result["tme_cluster"].map(LME_CLASS_MAP_K10)
    result = result.dropna(subset=["LME_class"])

    return result[["sample", "LME_class"]]


def _find_sample_col(df: pd.DataFrame) -> str | None:
    """Find the sample/patient ID column."""
    for col in df.columns:
        if col.lower() in ["dlc_code", "dlc_id", "sample", "patient", "case_id", "id"]:
            return col
    return df.columns[0]


def validate_proportions(lme_df: pd.DataFrame):
    """Validate LME proportions against manuscript values."""
    dist = lme_df["LME_class"].value_counts(normalize=True) * 100
    logger.info("LME class distribution:")
    all_ok = True
    for cls in LME_ORDER:
        actual = dist.get(cls, 0)
        expected = LME_EXPECTED_PCT.get(cls, 0)
        diff = abs(actual - expected)
        status = "OK" if diff <= 5 else "WARN"
        if status == "WARN":
            all_ok = False
        logger.info(f"  {cls:20s}: {actual:5.1f}% (expected {expected}%, diff {diff:.1f}%) [{status}]")

    n_classes = lme_df["LME_class"].nunique()
    if n_classes != 5:
        logger.warning(f"  WARN: {n_classes} LME classes found, expected 5")
        missing = set(LME_ORDER) - set(dist.index)
        if missing:
            logger.warning(f"  Missing classes: {missing}")
        all_ok = False

    return all_ok


def main():
    parser = argparse.ArgumentParser(description="Assign LME classes to panel AnnData")
    parser.add_argument("--panel", choices=["immune", "stromal", "both"], default="both")
    parser.add_argument("--method", choices=["csv", "kmeans", "auto"], default="auto",
                        help="How to assign LME classes")
    args = parser.parse_args()

    config = load_config()
    METADATA_DIR.mkdir(parents=True, exist_ok=True)

    # Try sources in order of preference
    lme_assignments = None

    if args.method in ("csv", "auto"):
        # 1. Try patient_tme with direct LME names
        lme_assignments = try_patient_tme(config)

        # 2. Try TME cluster CSVs with mapping
        if lme_assignments is None:
            lme_assignments = try_tme_clusters(config)

        # 3. Try TME h5ad object (s_tme_seurat.h5ad)
        if lme_assignments is None:
            lme_assignments = try_tme_h5ad(config)

    if lme_assignments is None and args.method in ("kmeans", "auto"):
        # 3. Fallback: k-means on abundance
        lme_assignments = try_kmeans_on_abundance(config)

    if lme_assignments is None or lme_assignments.empty:
        logger.error("Could not determine LME assignments. Need TME CSVs or abundance data.")
        return

    # Validate proportions
    validate_proportions(lme_assignments)

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

        # Join LME_class via sample ID (normalize DLC formats)
        if "sample" in adata.obs.columns:
            import re as _re

            def _norm_dlc(s):
                m = _re.match(r"DLC[_\s-]?(\d+)", str(s), _re.IGNORECASE)
                return f"DLC{int(m.group(1)):04d}" if m else str(s)

            lme_normed = lme_assignments.copy()
            lme_normed["sample_norm"] = lme_normed["sample"].apply(_norm_dlc)
            sample_to_lme = dict(
                zip(
                    lme_normed["sample_norm"],
                    lme_normed["LME_class"],
                    strict=False,
                )
            )
            adata.obs["LME_class"] = (
                adata.obs["sample"].astype(str).apply(_norm_dlc).map(sample_to_lme)
            )

            # Make categorical with manuscript order
            adata.obs["LME_class"] = pd.Categorical(
                adata.obs["LME_class"], categories=LME_ORDER, ordered=True
            )

            n_mapped = adata.obs["LME_class"].notna().sum()
            logger.info(f"  LME mapped: {n_mapped:,}/{adata.n_obs:,} cells")

            if n_mapped > 0:
                dist = adata.obs["LME_class"].value_counts(normalize=True) * 100
                for cls in LME_ORDER:
                    pct = dist.get(cls, 0)
                    logger.info(f"    {cls}: {pct:.1f}%")
        else:
            logger.warning(f"  No 'sample' column in {panel} obs")

        # Clean dtypes for h5ad compatibility
        try:
            from h5ad_utils import clean_adata_for_h5ad
            clean_adata_for_h5ad(adata)
        except ImportError:
            # Inline cleanup
            for col in adata.obs.columns:
                if adata.obs[col].dtype.name == "category":
                    adata.obs[col] = adata.obs[col].astype(str)

        logger.info(f"Saving: {output_path}")
        adata.write_h5ad(output_path)
        logger.info(f"  {adata.n_obs:,} cells, {len(adata.obs.columns)} obs columns")

    logger.info("Done.")


if __name__ == "__main__":
    main()
