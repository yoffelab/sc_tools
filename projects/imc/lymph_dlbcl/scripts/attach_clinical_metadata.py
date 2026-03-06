#!/usr/bin/env python3
"""
Phase 1.2: Attach clinical metadata to panel AnnData objects.

Joins clinical data from downloaded CSVs:
- DLBCL_clinical_full.csv -> OS_time, OS_event, PFS_time, PFS_event, COO, LymphGen
- CTMA_121_punch_notes.csv -> sample-to-slide mapping
- CTMA121_mut_table.csv -> mutation status per gene

Usage:
    python scripts/attach_clinical_metadata.py [--panel immune|stromal|both]

Input:
    results/adata.immune.raw.p1.h5ad
    results/adata.stromal.raw.p1.h5ad

Output:
    results/adata.immune.annotated.p2.h5ad
    results/adata.stromal.annotated.p2.h5ad
"""

import argparse
import logging
from pathlib import Path

import anndata as ad
import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"


def load_config():
    """Load project config."""
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def load_clinical_data(config: dict) -> dict:
    """Load all clinical/metadata CSVs."""
    data = {}

    # Clinical full
    clinical_path = PROJECT_DIR / config["clinical"]["clinical_full"]
    if clinical_path.exists():
        df = pd.read_csv(clinical_path)
        logger.info(f"Clinical full: {df.shape}, columns: {df.columns.tolist()[:10]}")
        data["clinical"] = df
    else:
        logger.warning(f"Clinical full not found: {clinical_path}")

    # Punch notes (sample-to-slide mapping)
    punch_path = PROJECT_DIR / config["clinical"]["punch_notes"]
    if punch_path.exists():
        df = pd.read_csv(punch_path)
        logger.info(f"Punch notes: {df.shape}, columns: {df.columns.tolist()[:10]}")
        data["punch_notes"] = df
    else:
        logger.warning(f"Punch notes not found: {punch_path}")

    # Mutation table
    mut_path = PROJECT_DIR / config["clinical"]["mutation_table"]
    if mut_path.exists():
        df = pd.read_csv(mut_path)
        logger.info(f"Mutation table: {df.shape}, columns: {df.columns.tolist()[:10]}")
        data["mutations"] = df
    else:
        logger.warning(f"Mutation table not found: {mut_path}")

    # TME cluster assignments
    tme_path = PROJECT_DIR / config["metadata"].get("tme_clusters", "")
    if tme_path and Path(tme_path).exists():
        df = pd.read_csv(tme_path)
        logger.info(f"TME clusters: {df.shape}")
        data["tme_clusters"] = df

    # TME z-scores
    zscore_path = PROJECT_DIR / config["metadata"].get("tme_zscore", "")
    if zscore_path and Path(zscore_path).exists():
        df = pd.read_csv(zscore_path)
        logger.info(f"TME z-scores: {df.shape}")
        data["tme_zscore"] = df

    return data


def find_join_key(clinical_df: pd.DataFrame, adata_samples: set) -> str | None:
    """Find the column in clinical_df that maps to adata sample IDs."""
    for col in clinical_df.columns:
        vals = set(clinical_df[col].astype(str).unique())
        overlap = vals & adata_samples
        if len(overlap) > 10:  # Reasonable threshold
            logger.info(f"  Join key candidate: '{col}' ({len(overlap)} overlapping samples)")
            return col
    return None


def attach_clinical(adata: ad.AnnData, clinical_data: dict, panel: str) -> ad.AnnData:
    """Attach clinical metadata to AnnData."""
    logger.info(f"Attaching clinical metadata to {panel} panel ({adata.n_obs:,} cells)")

    sample_col = "sample"
    if sample_col not in adata.obs.columns:
        logger.error(f"  No '{sample_col}' column in adata.obs")
        return adata

    samples = set(adata.obs[sample_col].astype(str).unique())
    logger.info(f"  {len(samples)} unique samples in panel")

    # Join clinical data
    if "clinical" in clinical_data:
        clin = clinical_data["clinical"].copy()
        join_key = find_join_key(clin, samples)

        if join_key:
            clin[join_key] = clin[join_key].astype(str)
            clin = clin.set_index(join_key)

            # Target columns for clinical data
            target_cols = []
            for col_pattern in ["OS", "PFS", "COO", "LymphGen", "Age", "Sex",
                                "Stage", "IPI", "ECOG", "LDH", "Response",
                                "os_time", "os_event", "pfs_time", "pfs_event"]:
                for col in clin.columns:
                    if col_pattern.lower() in col.lower() and col not in target_cols:
                        target_cols.append(col)

            if not target_cols:
                target_cols = [c for c in clin.columns if c not in [join_key]][:20]

            logger.info(f"  Joining {len(target_cols)} clinical columns via '{join_key}'")

            sample_to_clin = {}
            for sample_id in samples:
                if sample_id in clin.index:
                    sample_to_clin[sample_id] = clin.loc[sample_id, target_cols]

            for col in target_cols:
                adata.obs[col] = adata.obs[sample_col].astype(str).map(
                    {s: row[col] for s, row in sample_to_clin.items() if col in row.index}
                )

            n_mapped = sum(adata.obs[sample_col].astype(str).isin(sample_to_clin.keys()))
            logger.info(f"  Mapped {n_mapped:,}/{adata.n_obs:,} cells to clinical data")
        else:
            logger.warning("  Could not find join key between clinical and sample IDs")
            logger.info(f"  Clinical columns: {clinical_data['clinical'].columns.tolist()[:10]}")
            logger.info(f"  Sample IDs (first 5): {sorted(list(samples))[:5]}")

    # Join TME clusters
    if "tme_clusters" in clinical_data:
        tme = clinical_data["tme_clusters"].copy()
        join_key = find_join_key(tme, samples)
        if join_key:
            tme[join_key] = tme[join_key].astype(str)
            tme_map = tme.set_index(join_key)
            for col in tme_map.columns:
                adata.obs[f"TME_{col}"] = adata.obs[sample_col].astype(str).map(
                    tme_map[col].to_dict()
                )
            logger.info(f"  Attached {len(tme_map.columns)} TME cluster columns")

    # Join mutations
    if "mutations" in clinical_data:
        mut = clinical_data["mutations"].copy()
        join_key = find_join_key(mut, samples)
        if join_key:
            mut[join_key] = mut[join_key].astype(str)
            mut_map = mut.set_index(join_key)
            for col in mut_map.columns:
                adata.obs[f"mut_{col}"] = adata.obs[sample_col].astype(str).map(
                    mut_map[col].to_dict()
                )
            logger.info(f"  Attached {len(mut_map.columns)} mutation columns")

    return adata


def main():
    parser = argparse.ArgumentParser(description="Attach clinical metadata to panel AnnData")
    parser.add_argument("--panel", choices=["immune", "stromal", "both"], default="both")
    args = parser.parse_args()

    config = load_config()
    clinical_data = load_clinical_data(config)

    panels = ["immune", "stromal"] if args.panel == "both" else [args.panel]

    for panel in panels:
        input_path = RESULTS_DIR / f"adata.{panel}.raw.p1.h5ad"
        output_path = RESULTS_DIR / f"adata.{panel}.annotated.p2.h5ad"

        if not input_path.exists():
            logger.error(f"Input not found: {input_path}")
            continue

        logger.info(f"Loading: {input_path}")
        adata = ad.read_h5ad(input_path)

        adata = attach_clinical(adata, clinical_data, panel)

        # Clean dtypes for h5ad compatibility — force all non-numeric to string
        for col in adata.obs.columns:
            s = adata.obs[col]
            numeric = pd.to_numeric(s, errors="coerce")
            if numeric.notna().sum() == s.notna().sum() and s.notna().any():
                adata.obs[col] = numeric
            else:
                adata.obs[col] = s.astype(str).replace(
                    {"nan": "", "None": "", "<NA>": ""}
                )
        for col in adata.var.columns:
            s = adata.var[col]
            numeric = pd.to_numeric(s, errors="coerce")
            if numeric.notna().sum() == s.notna().sum() and s.notna().any():
                adata.var[col] = numeric
            else:
                adata.var[col] = s.astype(str).replace(
                    {"nan": "", "None": "", "<NA>": ""}
                )
        # Fix _index reserved column
        if "_index" in adata.var.columns:
            adata.var = adata.var.drop(columns=["_index"])
        if adata.raw is not None and "_index" in adata.raw.var.columns:
            raw_var = adata.raw.var.copy()
            raw_var = raw_var.drop(columns=["_index"])
            from anndata import Raw
            adata._raw = Raw(
                adata, X=adata.raw.X, var=raw_var, varm=adata.raw.varm
            )

        logger.info(f"Saving: {output_path}")
        adata.write_h5ad(output_path)
        logger.info(f"  {adata.n_obs:,} cells, {len(adata.obs.columns)} obs columns")

    logger.info("Done.")


if __name__ == "__main__":
    main()
