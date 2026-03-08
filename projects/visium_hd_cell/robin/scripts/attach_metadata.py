"""Phase 2 (metadata_attach): Attach clinical metadata to robin cell data.

Reads results/adata.filtered.h5ad and metadata/sample_metadata.xlsx,
maps sample IDs to clinical annotations (patient, tumor status, organ, timepoint),
saves results/adata.annotated.h5ad.

Run from project root or via SLURM.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = PROJECT_ROOT.parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


# Sample metadata mapping (from sample_metadata.xlsx)
# Batch2 sample IDs have suffixes (_S1, _S2, etc.) stripped to match the ROBIN Sample info column
SAMPLE_METADATA = {
    "PT01-1_NAT": {"patient": 1, "batch": "B1", "tissue_type": "Non-Tumor (biopsy)", "tumor": "Non_Tumor", "organ": "Colon", "timepoint": "BL", "fixation": "FF"},
    "PT01-2_TUM": {"patient": 1, "batch": "B1", "tissue_type": "Tumor (biopsy)", "tumor": "Tumor", "organ": "Colon", "timepoint": "BL", "fixation": "FF"},
    "PT01-3_NAT": {"patient": 1, "batch": "B1", "tissue_type": "Non-Tumor (biopsy)", "tumor": "Non_Tumor", "organ": "Colon", "timepoint": "PRT", "fixation": "FF"},
    "PT01-4_TUM": {"patient": 1, "batch": "B1", "tissue_type": "Tumor (biopsy)", "tumor": "Tumor", "organ": "Colon", "timepoint": "PRT", "fixation": "FF"},
    "PT01-5_NAT": {"patient": 1, "batch": "B1", "tissue_type": "Non-Tumor", "tumor": "Non_Tumor", "organ": "Colon", "timepoint": "S", "fixation": "FF"},
    "PT01-6_TUM": {"patient": 1, "batch": "B1", "tissue_type": "Tumor", "tumor": "Tumor", "organ": "Colon", "timepoint": "S", "fixation": "FF"},
    "PT05-2_TUM": {"patient": 5, "batch": "B1", "tissue_type": "Tumor (biopsy)", "tumor": "Tumor", "organ": "Colon", "timepoint": "BL", "fixation": "FFPE"},
    "PT05-4_TUM": {"patient": 5, "batch": "B1", "tissue_type": "Tumor (biopsy)", "tumor": "Tumor", "organ": "Colon", "timepoint": "PRT", "fixation": "FF"},
    "P3_S1_BL_NT_1_iF10_S7": {"patient": 3, "batch": "B2", "tissue_type": "Baseline-NonTumor", "tumor": "Non_Tumor", "organ": "Colon", "timepoint": "BL", "fixation": "FF"},
    "P3_S2_BL_T_2_iG10_S8": {"patient": 3, "batch": "B2", "tissue_type": "Baseline-Tumor", "tumor": "Tumor", "organ": "Colon", "timepoint": "BL", "fixation": "FF"},
    "Pat_1_LN_PRT_Samp8_D1_iC9_S3": {"patient": 1, "batch": "B2", "tissue_type": "LN (RT)", "tumor": None, "organ": "LN", "timepoint": "PRT", "fixation": "FF"},
    "Pat_3_LN_PRT_Samp_8_D1_iB9_S2": {"patient": 3, "batch": "B2", "tissue_type": "LN (RT)", "tumor": None, "organ": "LN", "timepoint": "PRT", "fixation": "FF"},
    "Pat_3_LN_Samp_7_D1_iA9_S1": {"patient": 3, "batch": "B2", "tissue_type": "LN (non-RT)", "tumor": None, "organ": "LN", "timepoint": "BL", "fixation": "FF"},
    "Pat_5_BL_Tu_Samp2_A1_iD9_S4": {"patient": 5, "batch": "B2", "tissue_type": "Tumor (biopsy)", "tumor": "Tumor", "organ": "Colon", "timepoint": "BL", "fixation": "FF"},
    "Pat5_Samp1_BL_2_iD10_S5": {"patient": 5, "batch": "B2", "tissue_type": "Baseline-NonTumor", "tumor": "Non_Tumor", "organ": "Colon", "timepoint": "BL", "fixation": "FF"},
}


def main():
    import pandas as pd
    import scanpy as sc

    raw_path = PROJECT_ROOT / "results" / "adata.filtered.h5ad"
    out_path = PROJECT_ROOT / "results" / "adata.annotated.h5ad"

    if not raw_path.exists():
        logger.error("adata.filtered.h5ad not found at %s", raw_path)
        sys.exit(1)

    logger.info("Reading %s", raw_path)
    adata = sc.read_h5ad(raw_path)
    logger.info("Loaded: %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Ensure sample column exists
    if "sample" not in adata.obs.columns:
        logger.error("obs['sample'] not found")
        sys.exit(1)

    samples = adata.obs["sample"].unique().tolist()
    logger.info("Samples found: %s", samples)

    # Map metadata
    meta_df = pd.DataFrame.from_dict(SAMPLE_METADATA, orient="index")
    meta_df.index.name = "sample"

    for col in meta_df.columns:
        mapping = meta_df[col].to_dict()
        adata.obs[col] = adata.obs["sample"].map(mapping)
        n_mapped = adata.obs[col].notna().sum()
        logger.info("Mapped %s: %d/%d cells", col, n_mapped, adata.n_obs)

    # Ensure library_id matches sample if not already set
    if "library_id" not in adata.obs.columns:
        adata.obs["library_id"] = adata.obs["sample"].copy()

    # Convert categoricals
    for col in ["patient", "batch", "tissue_type", "tumor", "organ", "timepoint", "fixation"]:
        if col in adata.obs.columns:
            adata.obs[col] = adata.obs[col].astype("category")

    adata.write_h5ad(out_path)
    logger.info("Saved annotated adata: %s (%d cells x %d genes)", out_path, adata.n_obs, adata.n_vars)


if __name__ == "__main__":
    main()
