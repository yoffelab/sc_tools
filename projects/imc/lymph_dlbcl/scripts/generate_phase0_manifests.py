"""Generate Phase 0 batch manifests for Plan B IMC reprocessing.

Reads the DLC clinical TSV and the h5ad var_names to produce per-panel batch
manifest TSVs under ``metadata/phase0/``.  The ``mcd_file`` column contains
placeholder paths derived from the expected layout on cayuga:

    /home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2/{sample_id}.mcd

**Step 1 (before running this script):** SSH to cayuga and verify that the
actual MCD/txt file paths match the expected pattern.  The file_path_sheet CSVs
in the DLBCLv2 directory (``6.22.22.file_path_sheet_s1.csv``,
``3.28.23.DLBCL_image_sheet.csv``) list the authoritative paths -- copy these
to ``metadata/phase0/`` and pass them via ``--path-sheet-immune`` /
``--path-sheet-stromal`` to use real paths instead of placeholders.

Usage
-----
    # With placeholder paths (default, cayuga not required)
    python scripts/generate_phase0_manifests.py

    # With real path sheets from cayuga (preferred)
    python scripts/generate_phase0_manifests.py \\
        --path-sheet-immune metadata/phase0/file_path_sheet_s2.csv \\
        --path-sheet-stromal metadata/phase0/file_path_sheet_s1.csv

Output
------
    metadata/phase0/batch1_immune.tsv
    metadata/phase0/batch1_stromal.tsv
    metadata/phase0/all_samples.tsv  (combined; created by collect_all_batches)
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# Raw data backup root on cayuga
CAYUGA_RAW_DIR = "/home/fs01/juk4007/elementolab/backup/dylan/hyperion/DLBCLv2"

# Expected per-panel panel CSV paths (relative to project root, used at pipeline runtime)
PANEL_IMMUNE = "metadata/panel_immune_t2.csv"
PANEL_STROMAL = "metadata/panel_stromal_s2.csv"


# ---------------------------------------------------------------------------
# DLC code normalisation helpers
# ---------------------------------------------------------------------------

def _dlc_h5ad_to_clinical(code: str) -> str | None:
    """Convert h5ad DLC code (``DLC0002``) to clinical TSV format (``DLC_0002``).

    Returns None for non-DLC codes (e.g. ``CTMA121``, ``0``).
    """
    m = re.match(r"^DLC(\d+)$", code, re.IGNORECASE)
    if m:
        return f"DLC_{int(m.group(1)):04d}"
    return None


def _clinical_to_h5ad(code: str) -> str | None:
    """Convert clinical DLC code (``DLC_0002``) to h5ad format (``DLC0002``)."""
    m = re.match(r"^DLC_(\d+)$", code, re.IGNORECASE)
    if m:
        return f"DLC{int(m.group(1)):04d}"
    return None


# ---------------------------------------------------------------------------
# Load sample IDs from h5ad
# ---------------------------------------------------------------------------

def _get_dlc_codes_from_h5ad(h5ad_path: Path) -> list[str]:
    """Return sorted list of DLC codes from h5ad obs.

    Tries ``DLC_code`` first, then falls back to ``orig.ident`` and index.
    """
    try:
        import anndata
    except ImportError:
        logger.error("anndata is required: pip install anndata")
        sys.exit(1)

    adata = anndata.read_h5ad(h5ad_path, backed="r")

    # Try known DLC code columns in priority order
    for col in ("DLC_code", "orig.ident", "sample", "library_id"):
        if col in adata.obs.columns:
            codes = adata.obs[col].dropna().unique().tolist()
            dlc_codes = sorted(
                str(c) for c in codes if re.match(r"^DLC\d+$", str(c), re.IGNORECASE)
            )
            if dlc_codes:
                logger.info("Using column '%s' from %s", col, h5ad_path.name)
                return dlc_codes

    # Last resort: try obs_names (cell barcodes) if they encode DLC
    index_codes = sorted(
        set(
            re.search(r"DLC\d+", str(idx)).group()
            for idx in adata.obs_names
            if re.search(r"DLC\d+", str(idx), re.IGNORECASE)
        )
    )
    if index_codes:
        logger.info("Extracted DLC codes from obs_names of %s", h5ad_path.name)
        return index_codes

    logger.warning("No DLC codes found in %s -- returning empty list", h5ad_path)
    return []


# ---------------------------------------------------------------------------
# Build manifest
# ---------------------------------------------------------------------------

def build_manifest(
    sample_ids: list[str],
    panel: str,
    panel_csv: str,
    raw_dir: str = CAYUGA_RAW_DIR,
    path_sheet: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build a Phase 0 batch manifest DataFrame.

    Parameters
    ----------
    sample_ids:
        List of sample IDs in h5ad format (e.g. ``['DLC0002', 'DLC0006']``).
    panel:
        Panel label: ``'immune'`` or ``'stromal'``.
    panel_csv:
        Relative path to the panel CSV file.
    raw_dir:
        Root directory containing MCD files on cayuga.
    path_sheet:
        Optional DataFrame loaded from a path_sheet CSV.  If provided, MCD
        paths are looked up from this table instead of inferred from raw_dir.

    Returns
    -------
    DataFrame with columns: sample_id, mcd_file, panel_csv, panel, batch.
    """
    rows = []
    for sid in sample_ids:
        mcd_file = _resolve_mcd_path(sid, raw_dir, path_sheet)
        rows.append({
            "sample_id": sid,
            "mcd_file": mcd_file,
            "panel_csv": panel_csv,
            "panel": panel,
            "batch": "batch1",
        })
    return pd.DataFrame(rows, columns=["sample_id", "mcd_file", "panel_csv", "panel", "batch"])


def _resolve_mcd_path(
    sample_id: str,
    raw_dir: str,
    path_sheet: pd.DataFrame | None,
) -> str:
    """Return MCD path for a sample_id, using path_sheet when available."""
    if path_sheet is not None:
        # Look for sample_id in path_sheet (try various column name patterns)
        for col in ("sample_id", "SampleID", "sample", "DLC_code", "dlc_code"):
            if col in path_sheet.columns:
                matches = path_sheet[path_sheet[col].astype(str) == sample_id]
                if not matches.empty:
                    for path_col in ("mcd_file", "mcd_path", "path", "file_path"):
                        if path_col in path_sheet.columns:
                            return str(matches.iloc[0][path_col])
    # Fallback: infer from expected directory layout
    return f"{raw_dir}/{sample_id}.mcd"


# ---------------------------------------------------------------------------
# Collect all batches (thin wrapper around sc_tools.ingest.config)
# ---------------------------------------------------------------------------

def collect_all_batches(phase0_dir: Path) -> Path:
    """Merge all batch*.tsv files in phase0_dir into all_samples.tsv."""
    batch_files = sorted(phase0_dir.glob("batch*.tsv"))
    if not batch_files:
        logger.warning("No batch*.tsv files found in %s", phase0_dir)
        out = phase0_dir / "all_samples.tsv"
        pd.DataFrame(columns=["sample_id", "mcd_file", "panel_csv", "panel", "batch"]).to_csv(
            out, sep="\t", index=False
        )
        return out
    dfs = [pd.read_csv(f, sep="\t") for f in batch_files]
    combined = pd.concat(dfs, ignore_index=True)
    out = phase0_dir / "all_samples.tsv"
    combined.to_csv(out, sep="\t", index=False)
    logger.info("Wrote %d samples to %s", len(combined), out)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--immune-h5ad",
        default="results/seurat_converted/tcell_2_preprocessing/t2_SO_seurat.h5ad",
        help="Immune panel h5ad (to extract sample IDs)",
    )
    parser.add_argument(
        "--stromal-h5ad",
        default="results/seurat_converted/stroma_2_preprocessing/S1_seurat_SO.h5ad",
        help="Stromal panel h5ad (to extract sample IDs)",
    )
    parser.add_argument(
        "--path-sheet-immune",
        default=None,
        help="Optional path sheet CSV for immune panel (from cayuga backup)",
    )
    parser.add_argument(
        "--path-sheet-stromal",
        default=None,
        help="Optional path sheet CSV for stromal panel (from cayuga backup)",
    )
    parser.add_argument(
        "--raw-dir",
        default=CAYUGA_RAW_DIR,
        help="Root directory of raw MCD files on cayuga (default: %(default)s)",
    )
    parser.add_argument(
        "--clinical-tsv",
        default="metadata/DLC380_clinical.tsv",
        help="Clinical TSV with DLC_ID column (used as stromal sample fallback)",
    )
    parser.add_argument(
        "--use-clinical-for-stromal",
        action="store_true",
        default=False,
        help="Use clinical TSV DLC IDs as stromal sample IDs when h5ad has none",
    )
    parser.add_argument(
        "--output-dir",
        default="metadata/phase0",
        help="Output directory for manifest TSVs (default: %(default)s)",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load optional path sheets
    immune_sheet = None
    if args.path_sheet_immune:
        immune_sheet = pd.read_csv(args.path_sheet_immune)
        logger.info("Loaded immune path sheet: %d rows", len(immune_sheet))

    stromal_sheet = None
    if args.path_sheet_stromal:
        stromal_sheet = pd.read_csv(args.path_sheet_stromal)
        logger.info("Loaded stromal path sheet: %d rows", len(stromal_sheet))

    # Extract sample IDs from h5ad
    immune_path = Path(args.immune_h5ad)
    stromal_path = Path(args.stromal_h5ad)

    immune_samples = []
    if immune_path.exists():
        immune_samples = _get_dlc_codes_from_h5ad(immune_path)
        logger.info("Immune panel: %d DLC samples from h5ad", len(immune_samples))
    else:
        logger.warning("Immune h5ad not found: %s -- manifest will be empty", immune_path)

    stromal_samples = []
    if stromal_path.exists():
        stromal_samples = _get_dlc_codes_from_h5ad(stromal_path)
        logger.info("Stromal panel: %d DLC samples from h5ad", len(stromal_samples))
    else:
        logger.warning("Stromal h5ad not found: %s -- manifest will be empty", stromal_path)

    # Fallback: use clinical DLC IDs for stromal panel when h5ad has no patient info
    if not stromal_samples and args.use_clinical_for_stromal:
        clinical_path = Path(args.clinical_tsv)
        if clinical_path.exists():
            clin = pd.read_csv(clinical_path, sep="\t")
            # Convert DLC_XXXX -> DLCXXXX to match h5ad format
            all_ids = clin["DLC_ID"].dropna().tolist()
            stromal_samples = sorted(
                _clinical_to_h5ad(c) for c in all_ids if _clinical_to_h5ad(c)
            )
            logger.info(
                "Stromal panel: %d DLC samples from clinical TSV (fallback)", len(stromal_samples)
            )
        else:
            logger.warning("Clinical TSV not found: %s", clinical_path)

    # Build manifests
    immune_df = build_manifest(
        immune_samples, "immune", PANEL_IMMUNE, args.raw_dir, immune_sheet
    )
    stromal_df = build_manifest(
        stromal_samples, "stromal", PANEL_STROMAL, args.raw_dir, stromal_sheet
    )

    # Write per-panel batch TSVs
    immune_out = out_dir / "batch1_immune.tsv"
    stromal_out = out_dir / "batch1_stromal.tsv"
    immune_df.to_csv(immune_out, sep="\t", index=False)
    stromal_df.to_csv(stromal_out, sep="\t", index=False)
    logger.info("Wrote immune manifest: %s (%d samples)", immune_out, len(immune_df))
    logger.info("Wrote stromal manifest: %s (%d samples)", stromal_out, len(stromal_df))

    # Merge into all_samples.tsv
    all_out = collect_all_batches(out_dir)
    logger.info("Merged manifest: %s", all_out)

    print(f"\nImmune manifest:  {immune_out} ({len(immune_df)} samples)")
    print(f"Stromal manifest: {stromal_out} ({len(stromal_df)} samples)")
    print(f"All samples:      {all_out}")
    if immune_sheet is None and stromal_sheet is None:
        print(
            "\nWARNING: mcd_file paths are PLACEHOLDERS."
            "\nSSH to cayuga and verify paths, or re-run with --path-sheet-* flags."
        )


if __name__ == "__main__":
    main()
