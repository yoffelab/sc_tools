#!/usr/bin/env python3
"""
Phase 1.2: Attach clinical metadata to panel AnnData objects.

Uses metadata/DLC380_clinical.tsv directly (348 cases, proper survival/COO/IPI).
Joins via DLC_ID -> obs['sample'] (DLC_code format, e.g. DLC_0001).

Usage:
    python scripts/attach_clinical_metadata.py [--panel immune|stromal|both]

Input:
    results/adata.immune.raw.p1.h5ad
    results/adata.stromal.raw.p1.h5ad
    metadata/DLC380_clinical.tsv

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
METADATA_DIR = PROJECT_DIR / "metadata"

# Column renaming from DLC380_clinical.tsv to standard names
CLINICAL_RENAME = {
    "Overall survival (y)": "OS_time",
    "CODE_OS": "OS_event",
    "Progression free survival (y)": "PFS_time",
    "CODE_PFS": "PFS_event",
    "Disease specific survival (y)": "DSS_time",
    "CODE_DSS": "DSS_event",
    "Time to progression (y)": "TTP_time",
    "CODE_TTP": "TTP_event",
    "LYMPH2CX_COO": "COO",
    "AGE": "age",
    "SEX": "sex",
    "IPI": "IPI",
    "STAGE": "stage",
    "PS": "ECOG",
    "LDH": "LDH",
    "LDHRATIO": "LDH_ratio",
    "MYC_BA": "MYC_break",
    "BCL2_BA": "BCL2_break",
    "BCL6_BA": "BCL6_break",
}


def load_config():
    """Load project config."""
    with open(PROJECT_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def load_clinical() -> pd.DataFrame:
    """Load and standardize DLC380_clinical.tsv."""
    clinical_path = METADATA_DIR / "DLC380_clinical.tsv"
    if not clinical_path.exists():
        raise FileNotFoundError(f"Clinical metadata not found: {clinical_path}")

    df = pd.read_csv(clinical_path, sep="\t")
    logger.info(f"Loaded clinical: {df.shape[0]} cases, {df.shape[1]} columns")

    # Filter to final cohort
    if "FINAL_COHORT" in df.columns:
        before = len(df)
        df = df[df["FINAL_COHORT"] == "YES"].copy()
        logger.info(f"  Filtered FINAL_COHORT=YES: {before} -> {len(df)} cases")

    # Rename columns to standard names
    rename_map = {k: v for k, v in CLINICAL_RENAME.items() if k in df.columns}
    df = df.rename(columns=rename_map)
    logger.info(f"  Renamed {len(rename_map)} columns to standard names")

    # Set DLC_ID as index for joining
    df["DLC_ID"] = df["DLC_ID"].astype(str)
    df = df.set_index("DLC_ID")

    return df


def load_mutations(config: dict) -> pd.DataFrame | None:
    """Load mutation table if available."""
    mut_path = PROJECT_DIR / config["clinical"]["mutation_table"]
    if not mut_path.exists():
        logger.info(f"Mutation table not found: {mut_path}")
        return None

    mut = pd.read_csv(mut_path)
    logger.info(f"Mutation table: {mut.shape}")
    return mut


def normalize_dlc_id(dlc_id: str) -> str:
    """Normalize DLC ID formats: DLC0002, DLC_0002, DLC_002 -> DLC_0002."""
    s = str(dlc_id).strip()
    # Extract numeric part from DLC prefix
    import re
    m = re.match(r"DLC[_\s-]?(\d+)", s, re.IGNORECASE)
    if m:
        return f"DLC_{int(m.group(1)):04d}"
    return s


def attach_clinical(adata: ad.AnnData, clinical: pd.DataFrame, panel: str) -> ad.AnnData:
    """Attach clinical metadata to AnnData via sample column."""
    logger.info(f"Attaching clinical metadata to {panel} panel ({adata.n_obs:,} cells)")

    if "sample" not in adata.obs.columns:
        logger.error("  No 'sample' column in adata.obs")
        return adata

    samples = set(adata.obs["sample"].astype(str).unique())
    logger.info(f"  {len(samples)} unique samples in panel")

    # Normalize clinical index
    clinical_norm = clinical.copy()
    clinical_norm.index = clinical_norm.index.map(normalize_dlc_id)

    # Normalize sample IDs and check overlap
    sample_norm_map = {s: normalize_dlc_id(s) for s in samples}
    samples_norm = set(sample_norm_map.values())
    overlap = samples_norm & set(clinical_norm.index)
    logger.info(f"  {len(overlap)}/{len(samples)} samples match clinical data (after normalization)")

    if len(overlap) == 0:
        logger.warning("  No overlap even after normalization")
        sample_examples = sorted(samples)[:5]
        clinical_examples = sorted(clinical.index.tolist())[:5]
        logger.info(f"  Sample IDs: {sample_examples}")
        logger.info(f"  Clinical IDs: {clinical_examples}")
        return adata

    # Target columns to attach
    target_cols = [c for c in clinical_norm.columns if c != "FINAL_COHORT"]
    logger.info(f"  Attaching {len(target_cols)} clinical columns")

    # Map sample -> normalized -> clinical data
    sample_str = adata.obs["sample"].astype(str)
    sample_normalized = sample_str.map(sample_norm_map)
    for col in target_cols:
        col_map = clinical_norm[col].to_dict()
        adata.obs[col] = sample_normalized.map(col_map)

    n_mapped = sample_normalized.isin(clinical_norm.index).sum()
    logger.info(f"  Mapped {n_mapped:,}/{adata.n_obs:,} cells to clinical data")

    # Report coverage
    for col in ["OS_time", "PFS_time", "COO", "age", "IPI"]:
        if col in adata.obs.columns:
            n_valid = adata.obs[col].notna().sum()
            logger.info(f"    {col}: {n_valid:,} non-null")

    return adata


def attach_mutations(adata: ad.AnnData, mut_df: pd.DataFrame) -> ad.AnnData:
    """Attach mutation status to AnnData."""
    if mut_df is None:
        return adata

    samples = set(adata.obs["sample"].astype(str).unique())

    # Find the sample ID column in mutation table
    join_col = None
    for col in mut_df.columns:
        overlap = set(mut_df[col].astype(str)) & samples
        if len(overlap) > 10:
            join_col = col
            break

    if join_col is None:
        logger.warning("  Could not find join key for mutations")
        return adata

    mut_df[join_col] = mut_df[join_col].astype(str)
    mut_map = mut_df.set_index(join_col)
    gene_cols = [c for c in mut_map.columns if c != join_col]

    sample_str = adata.obs["sample"].astype(str)
    for col in gene_cols:
        adata.obs[f"mut_{col}"] = sample_str.map(mut_map[col].to_dict())

    logger.info(f"  Attached {len(gene_cols)} mutation columns")
    return adata


def clean_adata_dtypes(adata: ad.AnnData):
    """Clean obs/var dtypes for h5ad compatibility."""
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


def main():
    parser = argparse.ArgumentParser(description="Attach clinical metadata to panel AnnData")
    parser.add_argument("--panel", choices=["immune", "stromal", "both"], default="both")
    args = parser.parse_args()

    config = load_config()
    clinical = load_clinical()
    mut_df = load_mutations(config)

    panels = ["immune", "stromal"] if args.panel == "both" else [args.panel]

    for panel in panels:
        input_path = RESULTS_DIR / f"adata.{panel}.raw.p1.h5ad"
        output_path = RESULTS_DIR / f"adata.{panel}.annotated.p2.h5ad"

        if not input_path.exists():
            logger.error(f"Input not found: {input_path}")
            continue

        logger.info(f"Loading: {input_path}")
        adata = ad.read_h5ad(input_path)

        adata = attach_clinical(adata, clinical, panel)
        adata = attach_mutations(adata, mut_df)
        clean_adata_dtypes(adata)

        logger.info(f"Saving: {output_path}")
        adata.write_h5ad(output_path)
        logger.info(f"  {adata.n_obs:,} cells, {len(adata.obs.columns)} obs columns")

    logger.info("Done.")


if __name__ == "__main__":
    main()
