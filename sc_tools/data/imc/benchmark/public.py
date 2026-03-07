"""Public IMC dataset download and standardization.

Downloads datasets from Zenodo/Mendeley and converts to steinbock convention
(``tiffs/{roi}_full.tiff``, ``tiffs/{roi}_full.csv``, ``tiffs/{roi}_full_mask.tiff``).
"""

from __future__ import annotations

import logging
import shutil
import zipfile
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "PUBLIC_DATASETS",
    "download_public_dataset",
    "standardize_public_dataset",
]

logger = logging.getLogger(__name__)


@dataclass
class PublicDatasetInfo:
    """Registry entry for a public IMC dataset."""

    name: str
    url: str
    n_images_approx: int
    n_channels: int
    tissue: str
    has_gt: bool
    gt_type: str
    description: str


PUBLIC_DATASETS: dict[str, PublicDatasetInfo] = {
    "steinbock_example": PublicDatasetInfo(
        name="steinbock_example",
        url="https://zenodo.org/records/7624451/files/steinbock_example.zip",
        n_images_approx=10,
        n_channels=40,
        tissue="Mixed",
        has_gt=True,
        gt_type="masks + prob maps",
        description="Steinbock example dataset with masks and probability maps",
    ),
    "jackson_2020_breast": PublicDatasetInfo(
        name="jackson_2020_breast",
        url="https://zenodo.org/records/3518284",
        n_images_approx=720,
        n_channels=35,
        tissue="Breast TMA",
        has_gt=True,
        gt_type="CellProfiler masks",
        description="Jackson et al. 2020 breast cancer TMA dataset",
    ),
    "immuncan_2024": PublicDatasetInfo(
        name="immuncan_2024",
        url="https://zenodo.org/records/12912567",
        n_images_approx=179,
        n_channels=40,
        tissue="5 cancer types",
        has_gt=True,
        gt_type="340K annotated cells",
        description="ImmuCan 2024 pan-cancer IMC dataset",
    ),
    "damond_2019_pancreas": PublicDatasetInfo(
        name="damond_2019_pancreas",
        url="https://data.mendeley.com/datasets/cydmwsfztj/2",
        n_images_approx=845,
        n_channels=38,
        tissue="Pancreas T1D",
        has_gt=True,
        gt_type="CellProfiler masks",
        description="Damond et al. 2019 pancreas T1D dataset",
    ),
    "rendeiro_2021_covid": PublicDatasetInfo(
        name="rendeiro_2021_covid",
        url="https://zenodo.org/records/4110560",
        n_images_approx=240,
        n_channels=36,
        tissue="COVID lung",
        has_gt=True,
        gt_type="Steinbock masks",
        description="Rendeiro et al. 2021 COVID lung IMC dataset",
    ),
}


def download_public_dataset(
    name: str,
    output_dir: str | Path,
    force: bool = False,
) -> Path:
    """Download a public IMC dataset.

    Parameters
    ----------
    name
        Dataset name (must be in PUBLIC_DATASETS).
    output_dir
        Directory to download into.
    force
        Re-download even if already present.

    Returns
    -------
    Path to the downloaded/extracted dataset directory.
    """
    if name not in PUBLIC_DATASETS:
        raise ValueError(f"Unknown dataset: {name!r}. Available: {list(PUBLIC_DATASETS.keys())}")

    info = PUBLIC_DATASETS[name]
    output_dir = Path(output_dir)
    ds_dir = output_dir / name
    done_sentinel = ds_dir / ".download_complete"

    if done_sentinel.is_file() and not force:
        logger.info("Dataset %s already downloaded at %s", name, ds_dir)
        return ds_dir

    ds_dir.mkdir(parents=True, exist_ok=True)

    try:
        import urllib.request
    except ImportError as e:
        raise ImportError("urllib is required for downloading datasets") from e

    url = info.url
    if url.endswith(".zip"):
        zip_path = ds_dir / f"{name}.zip"
        logger.info("Downloading %s from %s ...", name, url)
        urllib.request.urlretrieve(url, str(zip_path))

        logger.info("Extracting %s ...", zip_path.name)
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(ds_dir)
        zip_path.unlink()
    else:
        logger.warning(
            "Dataset %s URL is not a direct zip download. Manual download may be required from: %s",
            name,
            url,
        )
        # Create a README with download instructions
        readme = ds_dir / "README_DOWNLOAD.txt"
        readme.write_text(
            f"Download dataset manually from:\n{url}\n\n"
            f"Extract into this directory: {ds_dir}\n"
            f"Then run: python -m sc_tools.bm.cli prepare --dataset {name}\n"
        )
        return ds_dir

    done_sentinel.touch()
    logger.info("Download complete: %s (%s)", name, ds_dir)
    return ds_dir


def standardize_public_dataset(
    ds_dir: str | Path,
    name: str,
    output_dir: str | Path | None = None,
) -> Path:
    """Convert a downloaded public dataset to steinbock convention.

    Produces a ``processed/`` directory with:
    - ``tiffs/{roi}_full.tiff`` — multi-channel intensity stack
    - ``tiffs/{roi}_full.csv`` — channel index CSV
    - ``tiffs/{roi}_full_mask.tiff`` — labeled segmentation mask (if available)

    Parameters
    ----------
    ds_dir
        Path to the downloaded dataset.
    name
        Dataset name (for format-specific conversion).
    output_dir
        Where to write standardized output. Default: ``ds_dir/processed/``.

    Returns
    -------
    Path to the processed/ directory.
    """
    ds_dir = Path(ds_dir)
    if output_dir is None:
        output_dir = ds_dir / "processed"
    output_dir = Path(output_dir)

    if name == "steinbock_example":
        return _standardize_steinbock_example(ds_dir, output_dir)
    elif name == "jackson_2020_breast":
        return _standardize_jackson(ds_dir, output_dir)
    else:
        logger.warning(
            "No standardization implemented for %s. Attempting generic steinbock format detection.",
            name,
        )
        return _standardize_generic(ds_dir, output_dir)


def _standardize_steinbock_example(ds_dir: Path, output_dir: Path) -> Path:
    """Steinbock example: already in correct format, just symlink or copy."""
    # Look for the steinbock output structure
    candidates = [
        ds_dir / "steinbock_example",
        ds_dir,
    ]
    for candidate in candidates:
        tiff_dir = candidate / "img"
        mask_dir = candidate / "masks"
        if tiff_dir.is_dir():
            out_tiff = output_dir / "tiffs"
            out_tiff.mkdir(parents=True, exist_ok=True)

            try:
                import tifffile  # noqa: F401
            except ImportError as e:
                raise ImportError(
                    "tifffile required for standardization: pip install tifffile"
                ) from e

            for img_file in sorted(tiff_dir.glob("*.tiff")):
                roi_id = img_file.stem
                # Copy TIFF as _full.tiff
                dst = out_tiff / f"{roi_id}_full.tiff"
                if not dst.exists():
                    shutil.copy2(img_file, dst)

                # Create channel CSV if panel.csv exists
                panel = candidate / "panel.csv"
                if panel.is_file():
                    _create_channel_csv_from_panel(panel, out_tiff / f"{roi_id}_full.csv")

                # Copy mask if available
                mask_file = mask_dir / f"{roi_id}.tiff" if mask_dir.is_dir() else None
                if mask_file and mask_file.is_file():
                    shutil.copy2(mask_file, out_tiff / f"{roi_id}_full_mask.tiff")

            logger.info(
                "Standardized steinbock_example: %d ROIs", len(list(out_tiff.glob("*_full.tiff")))
            )
            return output_dir

    logger.warning("Could not find steinbock example data in %s", ds_dir)
    return output_dir


def _standardize_jackson(ds_dir: Path, output_dir: Path) -> Path:
    """Jackson 2020 breast: OME-TIFF per ROI with CellProfiler masks."""
    logger.info(
        "Jackson 2020 standardization: dataset is large (~720 images). Conversion may take time."
    )
    # This is a placeholder — full implementation depends on exact download structure
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "README.txt").write_text(
        "Jackson 2020 conversion requires manual setup.\nSee: https://zenodo.org/records/3518284\n"
    )
    return output_dir


def _standardize_generic(ds_dir: Path, output_dir: Path) -> Path:
    """Attempt generic steinbock format detection."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if already in steinbock format
    for subdir in [ds_dir, ds_dir / "processed"]:
        tiff_dir = subdir / "tiffs"
        if tiff_dir.is_dir() and list(tiff_dir.glob("*_full.tiff")):
            logger.info("Data already in steinbock format at %s", tiff_dir)
            if subdir != output_dir:
                out_tiff = output_dir / "tiffs"
                out_tiff.mkdir(parents=True, exist_ok=True)
                for f in tiff_dir.glob("*"):
                    dst = out_tiff / f.name
                    if not dst.exists():
                        shutil.copy2(f, dst)
            return output_dir

    logger.warning("Could not auto-detect format in %s", ds_dir)
    return output_dir


def _create_channel_csv_from_panel(panel_csv: Path, output_csv: Path) -> None:
    """Create a steinbock channel CSV from a panel.csv file."""
    import pandas as pd

    panel = pd.read_csv(panel_csv)

    # Try common column names for channel names
    name_col = None
    for candidate in ["name", "Name", "target", "Target", "channel_name", "marker"]:
        if candidate in panel.columns:
            name_col = candidate
            break

    if name_col is None:
        logger.warning("Cannot determine channel name column in %s", panel_csv)
        return

    tag_col = None
    for candidate in ["channel", "Channel", "Metal_Tag", "metal_tag"]:
        if candidate in panel.columns:
            tag_col = candidate
            break

    rows = []
    for idx, row in panel.iterrows():
        name = str(row[name_col])
        tag = str(row[tag_col]) if tag_col else ""
        if tag and not name.endswith(f"({tag})"):
            channel_str = f"{name}({tag})"
        else:
            channel_str = name
        rows.append({"index": idx, "channel": channel_str})

    out_df = pd.DataFrame(rows)
    out_df.to_csv(output_csv, index=False, header=False)
