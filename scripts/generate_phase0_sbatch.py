#!/usr/bin/env python3
"""Generate SLURM sbatch scripts for Phase 0 HPC execution.

Reads a batch manifest TSV and optional config.yaml, then generates
one sbatch script per sample for SpaceRanger, Xenium Ranger, or IMC.

Usage:
    python scripts/generate_phase0_sbatch.py \\
        --manifest projects/visium_hd/robin/metadata/phase0/batch2_samples.tsv \\
        --modality visium_hd \\
        --config projects/visium_hd/robin/config.yaml \\
        --sbatch-dir projects/visium_hd/robin/scripts/sbatch

    python scripts/generate_phase0_sbatch.py \\
        --manifest projects/imc/lymph_dlbcl/metadata/phase0/batch1_samples.tsv \\
        --modality imc \\
        --sbatch-dir projects/imc/lymph_dlbcl/scripts/sbatch
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate SLURM sbatch scripts from batch manifest"
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to batch manifest TSV (e.g. metadata/phase0/batch2_samples.tsv)",
    )
    parser.add_argument(
        "--modality",
        required=True,
        choices=["visium", "visium_hd", "visium_hd_cell", "xenium", "imc"],
        help="Data modality",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Path to config.yaml with transcriptome, probe_set, slurm settings",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="SpaceRanger/tool output directory (overrides config)",
    )
    parser.add_argument(
        "--sbatch-dir",
        default=None,
        help="Directory to write generated sbatch scripts",
    )
    parser.add_argument(
        "--inventory",
        action="store_true",
        help="Generate phase0_inventory.md alongside or instead of sbatch scripts",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing sbatch scripts",
    )
    parser.add_argument(
        "--transcriptome",
        default=None,
        help="Path to reference transcriptome (overrides config)",
    )
    parser.add_argument(
        "--probe-set",
        default=None,
        help="Path to probe set CSV (overrides config)",
    )
    parser.add_argument(
        "--spaceranger-path",
        default=None,
        help="Path to spaceranger binary (overrides config)",
    )

    args = parser.parse_args()

    # Late imports to keep CLI fast on --help
    import pandas as pd

    try:
        import yaml
    except ImportError:
        yaml = None

    from sc_tools.ingest.slurm import (
        build_batch_sbatch,
        generate_phase0_inventory,
        load_phase0_status,
        write_sbatch_script,
    )

    # Load manifest
    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        logger.error("Manifest not found: %s", manifest_path)
        sys.exit(1)

    manifest = pd.read_csv(manifest_path, sep="\t")
    logger.info("Loaded manifest: %d samples from %s", len(manifest), manifest_path)

    # Load config if provided
    config = {}
    if args.config:
        config_path = Path(args.config)
        if config_path.exists():
            if yaml is None:
                logger.error("PyYAML required to read config.yaml: pip install pyyaml")
                sys.exit(1)
            with open(config_path) as f:
                config = yaml.safe_load(f) or {}
            logger.info("Loaded config: %s", config_path)
        else:
            logger.warning("Config file not found: %s (using defaults)", config_path)

    # Validate: need at least one output mode
    if not args.inventory and not args.sbatch_dir:
        logger.error("Must specify --sbatch-dir, --inventory, or both")
        sys.exit(1)

    # Resolve parameters (CLI > config > defaults)
    transcriptome = args.transcriptome or config.get("transcriptome")
    probe_set = args.probe_set or config.get("probe_set")
    spaceranger_path = args.spaceranger_path or config.get("spaceranger_path", "spaceranger")
    output_dir = args.output_dir or config.get("output_dir", "spaceranger_outputs")
    slurm = config.get("slurm", {})

    # Generate inventory if requested
    if args.inventory:
        manifest_dir = manifest_path.parent
        status_path = manifest_dir / "phase0_status.tsv"
        status = load_phase0_status(status_path)

        inventory_path = manifest_dir / "phase0_inventory.md"
        generate_phase0_inventory(
            manifest=manifest,
            modality=args.modality,
            config=config,
            status=status,
            output_path=inventory_path,
        )
        logger.info("Generated inventory: %s", inventory_path)

    # Generate sbatch scripts if requested
    if args.sbatch_dir:
        scripts = build_batch_sbatch(
            manifest=manifest,
            modality=args.modality,
            output_dir=output_dir,
            transcriptome=transcriptome,
            probe_set=probe_set,
            spaceranger_path=spaceranger_path,
            slurm=slurm if slurm else None,
        )

        if not scripts:
            logger.warning("No scripts generated. Check manifest and modality.")
            sys.exit(1)

        sbatch_dir = Path(args.sbatch_dir)
        for sample_id, script_text in scripts:
            script_path = sbatch_dir / f"sr4_{sample_id}.sbatch"
            if args.modality == "xenium":
                script_path = sbatch_dir / f"xr_{sample_id}.sbatch"
            elif args.modality == "imc":
                script_path = sbatch_dir / f"imc_{sample_id}.sbatch"

            write_sbatch_script(script_text, script_path, overwrite=args.overwrite)
            logger.info("  %s", script_path)

        logger.info("Generated %d sbatch scripts in %s", len(scripts), sbatch_dir)


if __name__ == "__main__":
    main()
