"""Phase 0b ingest: Load per-sample Visium HD cell segmentation data into AnnData.

Reads SpaceRanger 4 segmented outputs for all robin samples, extracts cell
centroids from geojson, saves per-sample AnnData to data/{sample_id}/adata.h5ad,
then concatenates into results/adata.raw.h5ad.

Run from project root (projects/visium_hd_cell/robin/) or via SLURM on brb.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# Batch1 and Batch2 SR4 output base directories on brb
BATCH1_BASE = Path("/athena/elementolab/scratch/juk4007/robin/Robin_batch1/spaceranger_outputs")
BATCH2_BASE = Path("/athena/elementolab/scratch/juk4007/robin/Robin_batch2/spaceranger_outputs")

# Sample IDs per batch (from all_samples.tsv)
BATCH1_SAMPLES = [
    "PT01-1_NAT",
    "PT01-2_TUM",
    "PT01-3_NAT",
    "PT01-4_TUM",
    "PT01-5_NAT",
    "PT01-6_TUM",
    "PT05-2_TUM",
    "PT05-4_TUM",
]

BATCH2_SAMPLES = [
    "P3_S1_BL_NT_1_iF10_S7",
    "P3_S2_BL_T_2_iG10_S8",
    "Pat_1_LN_PRT_Samp8_D1_iC9_S3",
    "Pat_3_LN_PRT_Samp_8_D1_iB9_S2",
    "Pat_3_LN_Samp_7_D1_iA9_S1",
    "Pat_5_BL_Tu_Samp2_A1_iD9_S4",
    "Pat5_Samp1_BL_2_iD10_S5",
    # Pat5_Samp3_PRT_1_iE10_S6 excluded (failed QC)
]

# Project root (this script is at scripts/ingest_phase0b.py)
PROJECT_ROOT = Path(__file__).resolve().parent.parent


def main():
    from sc_tools.ingest.loaders import concat_samples, load_visium_hd_cell_sample

    adatas = []
    failed = []

    all_samples = [(s, BATCH1_BASE / s) for s in BATCH1_SAMPLES] + [
        (s, BATCH2_BASE / s) for s in BATCH2_SAMPLES
    ]

    for sample_id, sr_dir in all_samples:
        logger.info("Loading cell seg %s from %s", sample_id, sr_dir)
        try:
            adata = load_visium_hd_cell_sample(
                sr_dir,
                sample_id,
                load_images=False,
            )

            # Save per-sample checkpoint
            out_dir = PROJECT_ROOT / "data" / sample_id
            out_dir.mkdir(parents=True, exist_ok=True)
            out_path = out_dir / "adata.h5ad"
            adata.write_h5ad(out_path)
            logger.info(
                "Saved %s: %d cells x %d genes -> %s",
                sample_id,
                adata.n_obs,
                adata.n_vars,
                out_path,
            )
            adatas.append(adata)
        except Exception as exc:
            logger.error("FAILED %s: %s", sample_id, exc)
            failed.append(sample_id)

    if failed:
        logger.warning("Failed samples: %s", failed)

    if not adatas:
        logger.error("No samples loaded successfully")
        sys.exit(1)

    # Concatenate
    logger.info("Concatenating %d samples...", len(adatas))
    adata_concat = concat_samples(adatas)

    # Save concatenated
    results_dir = PROJECT_ROOT / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / "adata.raw.h5ad"
    adata_concat.write_h5ad(out_path)
    logger.info(
        "Saved concatenated: %d cells x %d genes -> %s",
        adata_concat.n_obs,
        adata_concat.n_vars,
        out_path,
    )

    if failed:
        logger.warning("Completed with %d failures: %s", len(failed), failed)
        sys.exit(1)

    logger.info(
        "Phase 0b ingest complete: %d samples, %d total cells", len(adatas), adata_concat.n_obs
    )


if __name__ == "__main__":
    main()
