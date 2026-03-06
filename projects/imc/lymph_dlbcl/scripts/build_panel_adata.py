#!/usr/bin/env python3
"""
Phase 1.1: Build panel-specific AnnData checkpoints from Seurat-converted h5ad.

For each panel (immune T2, stromal S2):
1. Load full preprocessing object as base (expression matrix + markers)
2. Map sample IDs (DLC_code from T2; orig.ident->DLC_code mapping for S2)
3. Transfer cell type labels from merged/subset objects
4. Create obs['celltype'] (30 subpops) + obs['celltype_broad'] (12 major)
5. Store raw intensities in layers['raw']

Usage:
    python scripts/build_panel_adata.py [--panel immune|stromal|both]
    python scripts/build_panel_adata.py --panel immune
    python scripts/build_panel_adata.py --panel stromal

Output:
    results/adata.immune.raw.p1.h5ad
    results/adata.stromal.raw.p1.h5ad
"""

import argparse
import json
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
METADATA_DIR = PROJECT_DIR / "metadata"


def load_config():
    """Load project config."""
    config_path = PROJECT_DIR / "config.yaml"
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_celltype_map(panel: str) -> dict:
    """Load cell type mapping JSON for a panel."""
    map_path = METADATA_DIR / f"celltype_map_{panel}.json"
    if map_path.exists():
        with open(map_path) as f:
            return json.load(f)
    logger.warning(f"No celltype map found at {map_path}; will use labels/meta directly")
    return {}


def build_immune_panel(config: dict) -> ad.AnnData:
    """Build immune panel (T2) AnnData."""
    logger.info("=== Building immune panel (T2) ===")

    # 1. Load T2 full object as base
    t2_path = PROJECT_DIR / config["objects"]["immune_full"]
    logger.info(f"Loading T2 full object: {t2_path}")
    adata = ad.read_h5ad(t2_path)
    logger.info(f"  Shape: {adata.shape}")

    # Store raw intensities
    adata.layers["raw"] = adata.X.copy()

    # 2. DLC_code is present in T2
    if "DLC_code" in adata.obs.columns:
        adata.obs["sample"] = adata.obs["DLC_code"].astype(str)
        n_samples = adata.obs["sample"].nunique()
        logger.info(f"  DLC_code present: {n_samples} unique samples")
    else:
        logger.warning("  DLC_code not found in T2 full object; using orig.ident")
        adata.obs["sample"] = adata.obs["orig.ident"].astype(str)

    # 3. Labels: T2 full has 'labels' column
    if "labels" in adata.obs.columns:
        labels = adata.obs["labels"]
        n_types = labels.nunique()
        logger.info(f"  labels column present: {n_types} unique types")
        logger.info(f"  Label values: {sorted(labels.unique().tolist())}")
    else:
        logger.warning("  No 'labels' column in T2 full; will transfer from subsets")
        # Transfer from merged object
        merged_path = PROJECT_DIR / config["objects"]["immune_merged"]
        if merged_path.exists():
            logger.info(f"  Loading immune merged for label transfer: {merged_path}")
            merged = ad.read_h5ad(merged_path)
            if "labels" in merged.obs.columns:
                # Join on index
                label_map = merged.obs["labels"].to_dict()
                adata.obs["labels"] = adata.obs.index.map(
                    lambda x: label_map.get(x, "Unknown")
                )
                logger.info(f"  Transferred {sum(adata.obs['labels'] != 'Unknown')} labels")

    # Also check 'meta' column for broader categories
    if "meta" in adata.obs.columns:
        logger.info(f"  meta column: {sorted(adata.obs['meta'].unique().tolist())}")

    # 4. Map to manuscript nomenclature
    celltype_map = load_celltype_map("immune")
    if celltype_map:
        # Map labels -> celltype (fine) and celltype_broad
        if "label_to_celltype" in celltype_map:
            adata.obs["celltype"] = adata.obs["labels"].map(
                celltype_map["label_to_celltype"]
            ).fillna("Unknown")
        if "label_to_celltype_broad" in celltype_map:
            adata.obs["celltype_broad"] = adata.obs["labels"].map(
                celltype_map["label_to_celltype_broad"]
            ).fillna("Unknown")
    else:
        # Use labels directly as celltype; derive broad from meta or label prefix
        if "labels" in adata.obs.columns:
            adata.obs["celltype"] = adata.obs["labels"].astype(str)
        if "meta" in adata.obs.columns:
            adata.obs["celltype_broad"] = adata.obs["meta"].astype(str)
        elif "labels" in adata.obs.columns:
            # Derive broad category from label prefix (before underscore or number)
            adata.obs["celltype_broad"] = adata.obs["labels"].astype(str)

    # 5. Ensure standard obs columns
    adata.obs["panel"] = "immune"
    adata.obs["panel_version"] = "T2"

    # Store orig.ident for traceability
    if "orig.ident" in adata.obs.columns:
        adata.obs["orig_ident"] = adata.obs["orig.ident"].astype(str)

    logger.info(f"  Final shape: {adata.shape}")
    logger.info(f"  obs columns: {sorted(adata.obs.columns.tolist())}")

    return adata


def build_stromal_panel(config: dict) -> ad.AnnData:
    """Build stromal panel (S2) AnnData."""
    logger.info("=== Building stromal panel (S2) ===")

    # 1. Load S2 full object as base
    s2_path = PROJECT_DIR / config["objects"]["stromal_full"]
    logger.info(f"Loading S2 full object: {s2_path}")
    adata = ad.read_h5ad(s2_path)
    logger.info(f"  Shape: {adata.shape}")

    # Store raw intensities
    adata.layers["raw"] = adata.X.copy()

    # 2. DLC_code mapping for S2 (not present in S2 preprocessing objects)
    # Try cell ID CSVs first
    s2_cellid_path = PROJECT_DIR / config["metadata"].get("s2_cellid", "")
    s1_cellid_path = PROJECT_DIR / config["metadata"].get("s1_cellid", "")

    if s2_cellid_path.exists():
        logger.info(f"  Loading S2 cell ID mapping: {s2_cellid_path}")
        cellid_df = pd.read_csv(s2_cellid_path)
        logger.info(f"  Cell ID columns: {cellid_df.columns.tolist()}")
        # Map cell barcodes to DLC_code
        if "DLC_code" in cellid_df.columns:
            barcode_col = cellid_df.columns[0]  # First col is usually barcode/index
            id_map = dict(zip(cellid_df[barcode_col], cellid_df["DLC_code"]))
            adata.obs["sample"] = adata.obs.index.map(
                lambda x: id_map.get(x, "Unknown")
            )
        else:
            logger.warning("  No DLC_code in cell ID CSV; using orig.ident")
            adata.obs["sample"] = adata.obs["orig.ident"].astype(str)
    else:
        logger.info("  No S2 cell ID CSV found; using orig.ident as sample")
        adata.obs["sample"] = adata.obs["orig.ident"].astype(str)

    n_samples = adata.obs["sample"].nunique()
    logger.info(f"  Samples: {n_samples} unique")

    # 3. Transfer labels from merged/subset objects
    # S2 subsets have 'meta' and 'recluster' but not always 'labels'
    # The merged stroma object has 'meta' column with cell type categories
    merged_path = PROJECT_DIR / config["objects"]["stromal_merged"]
    if merged_path.exists():
        logger.info(f"  Loading stromal merged for label transfer: {merged_path}")
        merged = ad.read_h5ad(merged_path)
        for col in ["meta", "recluster", "labels"]:
            if col in merged.obs.columns:
                col_map = merged.obs[col].to_dict()
                adata.obs[f"merged_{col}"] = adata.obs.index.map(
                    lambda x, cm=col_map: cm.get(x, "Unknown")
                )
                n_mapped = sum(adata.obs[f"merged_{col}"] != "Unknown")
                logger.info(f"  Transferred merged.{col}: {n_mapped}/{adata.n_obs} mapped")

    # Also transfer from S2 subsets
    subset_dirs = {
        "bcell": "stroma_2_preprocessing/S1_seurat_bcell.h5ad",
        "other": "stroma_2_preprocessing/S1_seurat_other.h5ad",
        "stroma": "stroma_2_preprocessing/S1_seurat_stroma.h5ad",
        "tcell": "stroma_2_preprocessing/S1_seurat_tcell.h5ad",
    }

    # Build celltype from subset membership
    celltype_from_subset = pd.Series("Unknown", index=adata.obs.index)
    for subset_name, subset_rel in subset_dirs.items():
        subset_path = SEURAT_DIR / subset_rel
        if subset_path.exists():
            subset_adata = ad.read_h5ad(subset_path, backed="r")
            subset_idx = set(subset_adata.obs.index)
            mask = adata.obs.index.isin(subset_idx)
            celltype_from_subset[mask] = subset_name

            # If subset has 'meta' or 'recluster', use for finer annotation
            if "meta" in subset_adata.obs.columns:
                meta_map = subset_adata.obs["meta"].to_dict()
                for idx in adata.obs.index[mask]:
                    if idx in meta_map:
                        celltype_from_subset[idx] = f"{subset_name}:{meta_map[idx]}"

            subset_adata.file.close()
            logger.info(f"  Subset {subset_name}: {mask.sum()} cells")

    adata.obs["celltype_from_subset"] = celltype_from_subset

    # 4. Map to manuscript nomenclature
    celltype_map = load_celltype_map("stromal")
    if celltype_map:
        if "label_to_celltype" in celltype_map:
            src_col = "merged_meta" if "merged_meta" in adata.obs.columns else "celltype_from_subset"
            adata.obs["celltype"] = adata.obs[src_col].map(
                celltype_map["label_to_celltype"]
            ).fillna("Unknown")
        if "label_to_celltype_broad" in celltype_map:
            src_col = "merged_meta" if "merged_meta" in adata.obs.columns else "celltype_from_subset"
            adata.obs["celltype_broad"] = adata.obs[src_col].map(
                celltype_map["label_to_celltype_broad"]
            ).fillna("Unknown")
    else:
        # Use best available annotation
        if "merged_meta" in adata.obs.columns:
            adata.obs["celltype"] = adata.obs["merged_meta"].astype(str)
            adata.obs["celltype_broad"] = adata.obs["celltype_from_subset"].str.split(":").str[0]
        else:
            adata.obs["celltype"] = adata.obs["celltype_from_subset"].astype(str)
            adata.obs["celltype_broad"] = adata.obs["celltype_from_subset"].str.split(":").str[0]

    # 5. Standard obs columns
    adata.obs["panel"] = "stromal"
    adata.obs["panel_version"] = "S2"

    if "orig.ident" in adata.obs.columns:
        adata.obs["orig_ident"] = adata.obs["orig.ident"].astype(str)

    logger.info(f"  Final shape: {adata.shape}")
    logger.info(f"  obs columns: {sorted(adata.obs.columns.tolist())}")

    return adata


def save_checkpoint(adata: ad.AnnData, panel: str, phase: str = "p1"):
    """Save AnnData checkpoint."""
    out_path = RESULTS_DIR / f"adata.{panel}.raw.{phase}.h5ad"
    logger.info(f"Saving: {out_path}")

    # Clean up for h5ad compatibility — force all non-numeric obs/var to string
    for col in adata.obs.columns:
        s = adata.obs[col]
        # Try numeric first
        numeric = pd.to_numeric(s, errors="coerce")
        if numeric.notna().sum() == s.notna().sum() and s.notna().any():
            adata.obs[col] = numeric
        else:
            # Force to string — handles categories, mixed types, objects
            adata.obs[col] = s.astype(str).replace({"nan": "", "None": "", "<NA>": ""})
    for col in adata.var.columns:
        s = adata.var[col]
        numeric = pd.to_numeric(s, errors="coerce")
        if numeric.notna().sum() == s.notna().sum() and s.notna().any():
            adata.var[col] = numeric
        else:
            adata.var[col] = s.astype(str).replace({"nan": "", "None": "", "<NA>": ""})

    # Fix _index column in var/raw.var (reserved name in anndata)
    if "_index" in adata.var.columns:
        adata.var = adata.var.drop(columns=["_index"])
    if adata.raw is not None and "_index" in adata.raw.var.columns:
        raw_var = adata.raw.var.copy()
        raw_var = raw_var.drop(columns=["_index"])
        # Rebuild raw without _index
        from anndata import Raw
        adata._raw = Raw(adata, X=adata.raw.X, var=raw_var, varm=adata.raw.varm)

    adata.write_h5ad(out_path)
    logger.info(f"  Saved: {adata.shape[0]:,} cells x {adata.shape[1]} markers")
    return out_path


def main():
    parser = argparse.ArgumentParser(description="Build panel-specific AnnData checkpoints")
    parser.add_argument("--panel", choices=["immune", "stromal", "both"], default="both",
                        help="Which panel to build (default: both)")
    args = parser.parse_args()

    config = load_config()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if args.panel in ("immune", "both"):
        adata_immune = build_immune_panel(config)
        save_checkpoint(adata_immune, "immune")

    if args.panel in ("stromal", "both"):
        adata_stromal = build_stromal_panel(config)
        save_checkpoint(adata_stromal, "stromal")

    logger.info("Done.")


if __name__ == "__main__":
    main()
