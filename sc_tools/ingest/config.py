"""Batch manifest parsing and collection for Phase 0.

Supports per-batch TSV files under metadata/phase0/ with modality-specific
column schemas. Concatenates all batch files into a collected manifest.

IMC manifest note: ``REQUIRED_COLUMNS["imc"]`` lists the minimum columns for
Phase 0b loading (``processed_dir``). Running the pipeline via Phase 0a also
requires ``mcd_file`` and ``panel_csv``; include them in the TSV when the
pipeline has not yet run.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def _read_tsv(path: str | os.PathLike) -> pd.DataFrame:
    """Read a TSV file from a local path or remote URI."""
    uri = str(path)
    if "://" in uri:
        try:
            from sc_tools.storage import smart_read_csv

            return smart_read_csv(uri, sep="\t")
        except ImportError:
            pass
    return pd.read_csv(uri, sep="\t")


# Required columns per modality
REQUIRED_COLUMNS = {
    "visium": {"sample_id", "fastq_dir", "image", "slide", "area"},
    "visium_hd": {"sample_id", "fastq_dir", "cytaimage", "slide", "area"},
    "visium_hd_cell": {"sample_id", "fastq_dir", "cytaimage", "slide", "area"},
    "xenium": {"sample_id", "xenium_dir"},
    # IMC Phase 0b minimum: processed_dir points to processed/{sample}/.
    # Also include mcd_file + panel_csv when running the pipeline (Phase 0a).
    "imc": {"sample_id", "processed_dir"},
    "cosmx": {"sample_id", "cosmx_dir"},  # flat CSV/Parquet or RDS output dir
}


def load_batch_manifest(path: str | os.PathLike) -> pd.DataFrame:
    """Load a single batch TSV manifest from a local path or remote URI.

    Parameters
    ----------
    path
        Path or URI to a tab-separated manifest file.

    Returns
    -------
    DataFrame with manifest rows.
    """
    uri = str(path)
    if "://" not in uri:
        # Local path — check existence before reading
        local = Path(uri)
        if not local.exists():
            raise FileNotFoundError(f"Manifest not found: {local}")
    df = _read_tsv(uri)
    logger.info("Loaded %d samples from %s", len(df), Path(uri).name)
    return df


def collect_all_batches(
    phase0_dir: str | os.PathLike,
    output: str | os.PathLike | None = None,
) -> pd.DataFrame:
    """Glob ``metadata/phase0/*_samples.tsv``, concatenate, write ``all_samples.tsv``.

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
    phase0_dir_str = str(phase0_dir)
    if "://" in phase0_dir_str:
        # Remote directory: use fsspec glob
        try:
            from sc_tools.storage import resolve_fs

            fs, remote_path = resolve_fs(phase0_dir_str)
            raw_matches = fs.glob(remote_path.rstrip("/") + "/*_samples.tsv")
            batch_files_strs = sorted(raw_matches)
            # Rebuild full URIs for remote batch files
            protocol = phase0_dir_str.split("://")[0]
            batch_files_uris = [f"{protocol}://{p}" for p in batch_files_strs]
        except ImportError:
            batch_files_uris = []
        batch_files = batch_files_uris  # type: ignore[assignment]
        phase0_dir = Path(phase0_dir_str)  # for output path fallback
    else:
        phase0_dir = Path(phase0_dir_str)
        batch_files = sorted(phase0_dir.glob("*_samples.tsv"))  # type: ignore[assignment]

    if not batch_files:
        logger.warning("No *_samples.tsv files found in %s", phase0_dir)
        return pd.DataFrame()

    dfs = []
    for f in batch_files:
        fname = f if isinstance(f, str) else f.name
        if str(fname).endswith("all_samples.tsv"):
            continue
        df = load_batch_manifest(f)
        if "batch" not in df.columns:
            # Infer batch name from filename (e.g., batch1_samples.tsv -> batch1)
            stem = Path(str(f)).stem
            batch_name = stem.replace("_samples", "")
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
