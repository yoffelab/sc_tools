"""Batch manifest parsing and collection for Phase 0.

Supports per-batch TSV files under metadata/phase0/ with modality-specific
column schemas. Concatenates all batch files into a collected manifest.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# Required columns per modality
REQUIRED_COLUMNS = {
    "visium": {"sample_id", "fastq_dir", "image", "slide", "area"},
    "visium_hd": {"sample_id", "fastq_dir", "cytaimage", "slide", "area"},
    "visium_hd_cell": {"sample_id", "fastq_dir", "cytaimage", "slide", "area"},
    "xenium": {"sample_id", "xenium_dir"},
    "imc": {"sample_id", "mcd_file", "panel_csv"},
    "cosmx": set(),  # CosMx: no Phase 0 (data assumed processed)
}


def load_batch_manifest(path: str | Path) -> pd.DataFrame:
    """Load a single batch TSV manifest.

    Parameters
    ----------
    path
        Path to a tab-separated manifest file.

    Returns
    -------
    DataFrame with manifest rows.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")
    df = pd.read_csv(path, sep="\t")
    logger.info("Loaded %d samples from %s", len(df), path.name)
    return df


def collect_all_batches(
    phase0_dir: str | Path,
    output: str | Path | None = None,
) -> pd.DataFrame:
    """Glob metadata/phase0/*_samples.tsv, concatenate, write all_samples.tsv.

    Parameters
    ----------
    phase0_dir
        Directory containing batch TSV files (e.g., metadata/phase0/).
    output
        Path to write the collected manifest. Defaults to
        ``phase0_dir/all_samples.tsv``.

    Returns
    -------
    Concatenated DataFrame of all batch manifests.
    """
    phase0_dir = Path(phase0_dir)
    batch_files = sorted(phase0_dir.glob("*_samples.tsv"))

    if not batch_files:
        logger.warning("No *_samples.tsv files found in %s", phase0_dir)
        return pd.DataFrame()

    dfs = []
    for f in batch_files:
        if f.name == "all_samples.tsv":
            continue
        df = load_batch_manifest(f)
        if "batch" not in df.columns:
            # Infer batch name from filename (e.g., batch1_samples.tsv -> batch1)
            batch_name = f.stem.replace("_samples", "")
            df["batch"] = batch_name
        dfs.append(df)

    if not dfs:
        return pd.DataFrame()

    combined = pd.concat(dfs, ignore_index=True)

    # Check for duplicate sample_ids
    if "sample_id" in combined.columns:
        dupes = combined["sample_id"][combined["sample_id"].duplicated()]
        if len(dupes) > 0:
            logger.warning(
                "Duplicate sample_ids across batches: %s",
                list(dupes.unique()),
            )

    if output is None:
        output = phase0_dir / "all_samples.tsv"
    output = Path(output)
    combined.to_csv(output, sep="\t", index=False)
    logger.info(
        "Collected %d samples from %d batch file(s) -> %s",
        len(combined),
        len(dfs),
        output,
    )

    return combined


def validate_manifest(
    df: pd.DataFrame,
    modality: str,
) -> list[str]:
    """Check that required columns exist for the modality.

    Parameters
    ----------
    df
        Manifest DataFrame.
    modality
        One of: visium, visium_hd, xenium, imc, cosmx.

    Returns
    -------
    List of validation issue messages. Empty means valid.
    """
    if modality not in REQUIRED_COLUMNS:
        return [f"Unknown modality '{modality}'. Must be one of: {list(REQUIRED_COLUMNS.keys())}"]

    required = REQUIRED_COLUMNS[modality]
    if not required:
        return []

    missing = required - set(df.columns)
    issues = []
    if missing:
        issues.append(f"Missing required columns for {modality}: {sorted(missing)}")

    if "sample_id" in df.columns and df["sample_id"].isna().any():
        issues.append("sample_id contains null values")

    return issues
