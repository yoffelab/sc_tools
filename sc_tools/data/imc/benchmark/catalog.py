"""Dataset discovery and ROI cataloging for IMC segmentation benchmark.

Scans HPC directories for internal datasets, builds a flat DataFrame
of all ROIs with metadata (dataset, tissue, paths to TIFF/mask/prob files).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

__all__ = [
    "IMCDatasetEntry",
    "ROIRecord",
    "discover_internal_datasets",
    "discover_rois",
    "build_benchmark_catalog",
]

logger = logging.getLogger(__name__)

# Known internal IMC datasets with expected locations relative to base_dir
INTERNAL_DATASETS: dict[str, dict] = {
    "ggo-imc": {
        "subdir": "ggo-imc/processed",
        "tissue": "Lung (GGO)",
        "panel_pattern": "PANEL_*",
    },
    "aapc": {
        "subdir": "aapc/processed",
        "tissue": "Colon (adenoma)",
        "panel_pattern": None,
    },
    "clonevo": {
        "subdir": "clonevo/processed",
        "tissue": "Bladder cancer",
        "panel_pattern": None,
    },
    "colon_liver_meta": {
        "subdir": "colon_liver_meta/processed",
        "tissue": "Colon/Liver",
        "panel_pattern": None,
    },
    "imc_normal_lung": {
        "subdir": "imc_normal_lung/processed",
        "tissue": "Normal lung",
        "panel_pattern": None,
    },
    "ggo_human": {
        "subdir": "ggo_human/processed",
        "tissue": "Lung (GGO)",
        "panel_pattern": None,
    },
    "lung_covid": {
        "subdir": "lung_covid/processed",
        "tissue": "COVID lung",
        "panel_pattern": None,
    },
    "bcg_bladder": {
        "subdir": "bcg_bladder/processed",
        "tissue": "Bladder (BCG)",
        "panel_pattern": None,
    },
}


@dataclass
class ROIRecord:
    """Metadata for a single ROI in the benchmark catalog."""

    dataset: str
    tissue: str
    sample_id: str
    roi_id: str
    tiff_path: str
    channel_csv_path: str | None = None
    mask_path: str | None = None
    prob_path: str | None = None
    panel_csv_path: str | None = None
    source: str = "internal"  # "internal" or "public"


@dataclass
class IMCDatasetEntry:
    """A discovered IMC dataset with its ROIs."""

    name: str
    base_dir: str
    tissue: str
    source: str = "internal"
    rois: list[ROIRecord] = field(default_factory=list)

    @property
    def n_rois(self) -> int:
        return len(self.rois)


def discover_rois(
    processed_dir: Path,
    dataset_name: str,
    tissue: str,
    source: str = "internal",
) -> list[ROIRecord]:
    """Discover ROIs in a processed IMC directory.

    Looks for ``tiffs/`` subdirectories containing ``*_full.tiff`` files
    (steinbock/ElementoLab convention).

    Parameters
    ----------
    processed_dir
        Path to the processed/ directory (or a sample subdirectory).
    dataset_name
        Name of the dataset.
    tissue
        Tissue type label.
    source
        "internal" or "public".

    Returns
    -------
    List of ROIRecord objects.
    """
    processed_dir = Path(processed_dir)
    records = []

    # Pattern 1: processed/{sample}/tiffs/{roi}_full.tiff
    # Pattern 2: processed/tiffs/{roi}_full.tiff (flat structure)
    tiff_dirs = list(processed_dir.glob("*/tiffs")) + (
        [processed_dir / "tiffs"] if (processed_dir / "tiffs").is_dir() else []
    )

    for tiff_dir in tiff_dirs:
        sample_id = tiff_dir.parent.name if tiff_dir.parent != processed_dir else "default"

        for tiff_path in sorted(tiff_dir.glob("*_full.tiff")):
            roi_id = tiff_path.stem.replace("_full", "")

            # Look for companion files
            channel_csv = tiff_path.with_name(f"{roi_id}_full.csv")
            mask_file = tiff_path.with_name(f"{roi_id}_full_mask.tiff")
            prob_file = tiff_path.with_name(f"{roi_id}_Probabilities.tiff")

            # Look for panel CSV at sample or dataset level
            panel_csv = None
            for candidate in [
                tiff_dir.parent / "channel_labels.csv",
                processed_dir / "channel_labels.csv",
                tiff_dir.parent / f"{sample_id}.channel_labels.csv",
            ]:
                if candidate.is_file():
                    panel_csv = str(candidate)
                    break

            records.append(
                ROIRecord(
                    dataset=dataset_name,
                    tissue=tissue,
                    sample_id=sample_id,
                    roi_id=roi_id,
                    tiff_path=str(tiff_path),
                    channel_csv_path=str(channel_csv) if channel_csv.is_file() else None,
                    mask_path=str(mask_file) if mask_file.is_file() else None,
                    prob_path=str(prob_file) if prob_file.is_file() else None,
                    panel_csv_path=panel_csv,
                    source=source,
                )
            )

    return records


def discover_internal_datasets(
    base_dir: str | Path,
    datasets: list[str] | None = None,
) -> list[IMCDatasetEntry]:
    """Scan HPC base directory for known internal IMC datasets.

    Parameters
    ----------
    base_dir
        Root directory containing dataset subdirectories
        (e.g., ``/athena/elementolab/scratch/juk4007/``).
    datasets
        Specific dataset names to look for. None = all known.

    Returns
    -------
    List of IMCDatasetEntry objects (only those found on disk).
    """
    base_dir = Path(base_dir)
    entries = []

    targets = datasets if datasets else list(INTERNAL_DATASETS.keys())

    for name in targets:
        if name not in INTERNAL_DATASETS:
            logger.warning("Unknown dataset: %s, skipping", name)
            continue

        info = INTERNAL_DATASETS[name]
        dataset_dir = base_dir / info["subdir"]

        if not dataset_dir.is_dir():
            logger.debug("Dataset %s not found at %s", name, dataset_dir)
            continue

        # If panel_pattern is set, scan each panel subdirectory
        if info.get("panel_pattern"):
            panel_dirs = sorted(dataset_dir.glob(info["panel_pattern"]))
            all_rois = []
            for pd_ in panel_dirs:
                all_rois.extend(discover_rois(pd_, name, info["tissue"], source="internal"))
            rois = all_rois
        else:
            rois = discover_rois(dataset_dir, name, info["tissue"], source="internal")

        if rois:
            entry = IMCDatasetEntry(
                name=name,
                base_dir=str(dataset_dir),
                tissue=info["tissue"],
                source="internal",
                rois=rois,
            )
            entries.append(entry)
            logger.info("Found %d ROIs in %s", len(rois), name)
        else:
            logger.debug("No ROIs found in %s", name)

    return entries


def build_benchmark_catalog(
    base_dir: str | Path | None = None,
    datasets: list[str] | None = None,
    public_dir: str | Path | None = None,
    include_public: bool = False,
) -> pd.DataFrame:
    """Build a flat DataFrame catalog of all benchmark ROIs.

    Parameters
    ----------
    base_dir
        Root directory for internal datasets.
    datasets
        Specific internal dataset names. None = all.
    public_dir
        Directory containing downloaded public datasets.
    include_public
        Whether to include public datasets.

    Returns
    -------
    DataFrame with columns: dataset, tissue, sample_id, roi_id,
    tiff_path, channel_csv_path, mask_path, prob_path, panel_csv_path,
    source, has_gt (bool).
    """
    all_rois: list[ROIRecord] = []

    # Internal datasets
    if base_dir is not None:
        entries = discover_internal_datasets(base_dir, datasets)
        for entry in entries:
            all_rois.extend(entry.rois)

    # Public datasets
    if include_public and public_dir is not None:
        public_dir = Path(public_dir)
        if public_dir.is_dir():
            for ds_dir in sorted(public_dir.iterdir()):
                if ds_dir.is_dir():
                    rois = discover_rois(ds_dir, ds_dir.name, tissue="public", source="public")
                    all_rois.extend(rois)

    if not all_rois:
        logger.warning("No ROIs found in catalog")
        return pd.DataFrame(
            columns=[
                "dataset",
                "tissue",
                "sample_id",
                "roi_id",
                "tiff_path",
                "channel_csv_path",
                "mask_path",
                "prob_path",
                "panel_csv_path",
                "source",
                "has_gt",
            ]
        )

    rows = []
    for r in all_rois:
        rows.append(
            {
                "dataset": r.dataset,
                "tissue": r.tissue,
                "sample_id": r.sample_id,
                "roi_id": r.roi_id,
                "tiff_path": r.tiff_path,
                "channel_csv_path": r.channel_csv_path,
                "mask_path": r.mask_path,
                "prob_path": r.prob_path,
                "panel_csv_path": r.panel_csv_path,
                "source": r.source,
                "has_gt": r.mask_path is not None,
            }
        )

    df = pd.DataFrame(rows)
    logger.info(
        "Catalog: %d ROIs from %d datasets (%d with GT)",
        len(df),
        df["dataset"].nunique(),
        df["has_gt"].sum(),
    )
    return df
