"""Generate panel CSV files for Plan B IMC reprocessing pipeline.

Extracts marker names from the Seurat-converted h5ad var_names and writes
imctools-compatible panel CSV files for the immune (T2) and stromal (S2) panels.

Usage
-----
    python scripts/generate_panel_csv.py [--output-dir metadata/]

The generated CSVs have format required by ``sc_tools.ingest.imc.build_imc_pipeline_cmd``:

    channel,Target,Metal_Tag,full,ilastik
    DNA1,DNA1,,1,1
    CD3,CD3,,1,0
    ...

NOTE: The ``Metal_Tag`` column is left empty because isotope tag information is
not stored in the h5ad files.  Fill it in using the MCD file headers or the
original Hyperion acquisition panel document once raw files are accessible on
cayuga.  After filling Metal_Tag, change the ``channel`` column to the
``MarkerName(IsotopeTag)`` format (e.g. ``CD3(Er170)``) required by imctools.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
from pathlib import Path

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Morphological / housekeeping var names to exclude (CellProfiler output)
# ---------------------------------------------------------------------------
_MORPH_PREFIXES = (
    "filled-area",
    "solidity",
    "eccentricity",
    "perimeter",
    "full-index",
    "n-s",   # normalization standard
    "p-a",   # placeholder / Palladium-A channel
)

# Compartment suffixes added by CellProfiler
_COMPARTMENT_SUFFIXES = ("-mem", "-nuc")

# Markers flagged ilastik=1 (used for cell/nucleus segmentation)
_ILASTIK_MARKERS = {
    "dna1",
    "dna2",
    "total h3",
    "h3",
    "histone h3",
    "ir191",
    "ir193",
}


def _strip_compartment(var_name: str) -> str | None:
    """Return protein name with compartment suffix removed, or None if morph."""
    vl = var_name.lower()
    for prefix in _MORPH_PREFIXES:
        if vl.startswith(prefix):
            return None
    for suffix in _COMPARTMENT_SUFFIXES:
        if vl.endswith(suffix):
            return var_name[: -len(suffix)]
    # No known suffix -- return as-is (e.g. "N-S", "full-index" already excluded)
    return None


def _is_ilastik(protein_name: str) -> int:
    return 1 if protein_name.lower() in _ILASTIK_MARKERS else 0


def extract_panel_from_varnames(var_names: list[str]) -> list[dict]:
    """Convert h5ad var_names to panel records (deduped by protein name).

    Parameters
    ----------
    var_names:
        List of var_names from an AnnData object (e.g. ``['CD38-mem', 'DNA1-nuc', ...]``).

    Returns
    -------
    List of dicts with keys: channel, Target, Metal_Tag, full, ilastik.
    """
    seen: dict[str, dict] = {}
    for var in var_names:
        protein = _strip_compartment(var)
        if protein is None:
            continue
        key = protein.lower()
        if key in seen:
            continue
        seen[key] = {
            "channel": protein,      # placeholder; update to Protein(IsotopeTag) after MCD review
            "Target": protein,
            "Metal_Tag": "",         # FILL IN from MCD file headers
            "full": 1,
            "ilastik": _is_ilastik(protein),
        }
    return list(seen.values())


def write_panel_csv(records: list[dict], path: Path) -> None:
    """Write panel records to CSV."""
    import csv
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["channel", "Target", "Metal_Tag", "full", "ilastik"]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    logger.info("Wrote %d channels to %s", len(records), path)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--immune-h5ad",
        default="results/seurat_converted/tcell_2_preprocessing/t2_SO_seurat.h5ad",
        help="Path to immune panel h5ad (default: %(default)s)",
    )
    parser.add_argument(
        "--stromal-h5ad",
        default="results/seurat_converted/stroma_2_preprocessing/S1_seurat_SO.h5ad",
        help="Path to stromal panel h5ad (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        default="metadata",
        help="Directory to write panel CSVs (default: %(default)s)",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    try:
        import anndata
    except ImportError:
        logger.error("anndata is required: pip install anndata")
        sys.exit(1)

    out_dir = Path(args.output_dir)

    # --- Immune panel (T2) ---
    immune_path = Path(args.immune_h5ad)
    if not immune_path.exists():
        logger.error("Immune h5ad not found: %s", immune_path)
        sys.exit(1)
    logger.info("Loading immune h5ad from %s", immune_path)
    adata_immune = anndata.read_h5ad(immune_path, backed="r")
    immune_records = extract_panel_from_varnames(list(adata_immune.var_names))
    immune_out = out_dir / "panel_immune_t2.csv"
    write_panel_csv(immune_records, immune_out)

    # --- Stromal panel (S2) ---
    stromal_path = Path(args.stromal_h5ad)
    if not stromal_path.exists():
        logger.error("Stromal h5ad not found: %s", stromal_path)
        sys.exit(1)
    logger.info("Loading stromal h5ad from %s", stromal_path)
    adata_stromal = anndata.read_h5ad(stromal_path, backed="r")
    stromal_records = extract_panel_from_varnames(list(adata_stromal.var_names))
    stromal_out = out_dir / "panel_stromal_s2.csv"
    write_panel_csv(stromal_records, stromal_out)

    print(f"Immune panel: {len(immune_records)} channels -> {immune_out}")
    print(f"Stromal panel: {len(stromal_records)} channels -> {stromal_out}")
    print()
    print("Next step: fill in 'Metal_Tag' column in both CSVs using MCD file headers.")
    print("Then update 'channel' to 'MarkerName(IsotopeTag)' format, e.g. 'CD3(Er170)'.")


if __name__ == "__main__":
    main()
