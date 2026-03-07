#!/usr/bin/env python3
"""Generic Phase 0b ingestion script.

Reads a batch manifest (all_samples.tsv), loads per-sample AnnData objects
using modality-specific loaders, optionally saves per-sample checkpoints,
and concatenates into a single AnnData for Phase 1.

Usage:
  # From project directory:
  python scripts/ingest.py --modality imc --manifest metadata/phase0/all_samples.tsv

  # Auto-collect manifest from batch TSVs:
  python scripts/ingest.py --modality visium --phase0-dir metadata/phase0/

  # Save per-sample checkpoints:
  python scripts/ingest.py --modality imc --save-per-sample
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def _get_loader(modality: str):
    """Return (loader_func, dir_column) for the given modality."""
    from sc_tools.ingest import (
        load_imc_sample,
        load_visium_hd_cell_sample,
        load_visium_hd_sample,
        load_visium_sample,
        load_xenium_sample,
    )

    loaders = {
        "visium": (load_visium_sample, "fastq_dir"),
        "visium_hd": (load_visium_hd_sample, "fastq_dir"),
        "visium_hd_cell": (load_visium_hd_cell_sample, "fastq_dir"),
        "xenium": (load_xenium_sample, "xenium_dir"),
        "imc": (load_imc_sample, "processed_dir"),
    }
    if modality == "cosmx":
        raise NotImplementedError("CosMx loader not yet implemented")
    if modality not in loaders:
        raise ValueError(f"Unknown modality: {modality}. Valid: {list(loaders)}")
    return loaders[modality]


def load_single_sample(row, modality: str, load_images: bool = False):
    """Load one sample from a manifest row."""
    loader_func, dir_col = _get_loader(modality)
    sample_id = row["sample_id"]
    sample_dir = row[dir_col]

    kwargs = {"load_images": load_images} if modality == "imc" else {}
    if modality in ("visium", "visium_hd", "visium_hd_cell"):
        kwargs["load_images"] = load_images

    return loader_func(sample_dir, sample_id, **kwargs)


def run_ingestion(
    manifest,
    modality: str,
    load_images: bool = False,
    save_per_sample: bool = False,
    per_sample_dir: Path | None = None,
):
    """Load all samples from manifest, return (adatas, failed_ids)."""
    adatas = []
    failed = []

    for _, row in manifest.iterrows():
        sample_id = row["sample_id"]
        try:
            adata = load_single_sample(row, modality, load_images=load_images)
            adatas.append(adata)

            if save_per_sample and per_sample_dir is not None:
                out_dir = per_sample_dir / sample_id
                out_dir.mkdir(parents=True, exist_ok=True)
                out_path = out_dir / "adata.h5ad"  # ingest_load canonical name
                adata.write_h5ad(out_path)
                logger.info("Saved per-sample checkpoint: %s", out_path)

        except Exception:
            logger.warning("Failed to load sample %s", sample_id, exc_info=True)
            failed.append(sample_id)

    return adatas, failed


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Phase 0b ingestion: manifest -> per-sample AnnData -> concatenated h5ad"
    )
    parser.add_argument(
        "--modality",
        required=True,
        choices=["visium", "visium_hd", "visium_hd_cell", "xenium", "imc", "cosmx"],
        help="Data modality.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Path to all_samples.tsv (default: auto-collect from --phase0-dir).",
    )
    parser.add_argument(
        "--phase0-dir",
        type=Path,
        default=Path("metadata/phase0"),
        help="Directory with batch TSVs (default: metadata/phase0/).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/adata.raw.h5ad"),
        help="Output path for concatenated h5ad (default: results/adata.raw.h5ad).",
    )
    parser.add_argument(
        "--save-per-sample",
        action="store_true",
        help="Save per-sample adata.h5ad checkpoints.",
    )
    parser.add_argument(
        "--per-sample-dir",
        type=Path,
        default=Path("data"),
        help="Directory for per-sample checkpoints (default: data/).",
    )
    parser.add_argument(
        "--load-images",
        action="store_true",
        help="Load images (IMC TIFFs, Visium H&E, etc.).",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    from sc_tools.ingest import (
        collect_all_batches,
        concat_samples,
        load_batch_manifest,
        validate_manifest,
    )

    # Load or collect manifest
    if args.manifest is not None:
        manifest = load_batch_manifest(args.manifest)
    else:
        all_samples = args.phase0_dir / "all_samples.tsv"
        if all_samples.exists():
            manifest = load_batch_manifest(all_samples)
        else:
            manifest = collect_all_batches(args.phase0_dir)

    if manifest.empty:
        logger.error("Empty manifest. Nothing to ingest.")
        sys.exit(1)

    # Validate
    issues = validate_manifest(manifest, args.modality)
    if issues:
        for issue in issues:
            logger.error("Manifest validation: %s", issue)
        sys.exit(1)

    logger.info("Manifest: %d samples, modality=%s", len(manifest), args.modality)

    # Load samples
    adatas, failed = run_ingestion(
        manifest,
        args.modality,
        load_images=args.load_images,
        save_per_sample=args.save_per_sample,
        per_sample_dir=args.per_sample_dir,
    )

    if not adatas:
        logger.error("All %d samples failed to load.", len(failed))
        sys.exit(1)

    # Concatenate
    adata = concat_samples(adatas, calculate_qc=True)

    # Save
    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)

    # Summary
    logger.info(
        "Ingestion complete: %d/%d samples loaded, %d failed -> %d obs x %d vars -> %s",
        len(adatas),
        len(manifest),
        len(failed),
        adata.n_obs,
        adata.n_vars,
        args.output,
    )
    if failed:
        logger.warning("Failed samples: %s", failed)


if __name__ == "__main__":
    main()
