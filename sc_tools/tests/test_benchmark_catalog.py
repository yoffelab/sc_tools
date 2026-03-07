"""Tests for sc_tools.data.imc.benchmark.catalog — dataset discovery and ROI cataloging."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def _create_mock_roi(base_dir: Path, sample: str, roi_id: str) -> None:
    """Create a minimal mock ROI directory structure."""
    tiff_dir = base_dir / sample / "tiffs"
    tiff_dir.mkdir(parents=True, exist_ok=True)

    # Create mock TIFF (3 channels, 10x10)
    import tifffile

    img = np.random.rand(3, 10, 10).astype(np.float32)
    tifffile.imwrite(str(tiff_dir / f"{roi_id}_full.tiff"), img)

    # Create channel CSV
    csv_path = tiff_dir / f"{roi_id}_full.csv"
    csv_path.write_text("0,DNA1(Ir191)\n1,CD3(Er170)\n2,PanCK(Pt195)\n")

    # Create mask
    mask = np.random.randint(0, 5, (10, 10), dtype=np.int32)
    tifffile.imwrite(str(tiff_dir / f"{roi_id}_full_mask.tiff"), mask)


class TestROIRecord:
    def test_dataclass_fields(self):
        from sc_tools.data.imc.benchmark.catalog import ROIRecord

        r = ROIRecord(
            dataset="test",
            tissue="lung",
            sample_id="s1",
            roi_id="roi1",
            tiff_path="/tmp/test.tiff",
        )
        assert r.dataset == "test"
        assert r.source == "internal"
        assert r.mask_path is None


class TestDiscoverROIs:
    def test_discovers_rois(self, tmp_path):
        from sc_tools.data.imc.benchmark.catalog import discover_rois

        _create_mock_roi(tmp_path, "sample1", "roi1")
        _create_mock_roi(tmp_path, "sample1", "roi2")

        rois = discover_rois(tmp_path, "test_ds", "lung")
        assert len(rois) == 2
        assert rois[0].dataset == "test_ds"
        assert rois[0].tissue == "lung"
        assert rois[0].channel_csv_path is not None
        assert rois[0].mask_path is not None

    def test_empty_directory(self, tmp_path):
        from sc_tools.data.imc.benchmark.catalog import discover_rois

        rois = discover_rois(tmp_path, "empty", "none")
        assert len(rois) == 0


class TestBuildCatalog:
    def test_builds_catalog_from_dir(self, tmp_path):
        from sc_tools.data.imc.benchmark.catalog import build_benchmark_catalog

        # Create mock dataset
        ds_dir = tmp_path / "ggo-imc" / "processed" / "PANEL_1"
        _create_mock_roi(ds_dir, "sample1", "roi1")
        _create_mock_roi(ds_dir, "sample1", "roi2")

        df = build_benchmark_catalog(base_dir=tmp_path, datasets=["ggo-imc"])
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2
        assert "tiff_path" in df.columns
        assert "has_gt" in df.columns
        assert all(df["has_gt"])

    def test_empty_catalog(self, tmp_path):
        from sc_tools.data.imc.benchmark.catalog import build_benchmark_catalog

        df = build_benchmark_catalog(base_dir=tmp_path, datasets=["nonexistent"])
        assert len(df) == 0


class TestIMCDatasetEntry:
    def test_n_rois(self):
        from sc_tools.data.imc.benchmark.catalog import IMCDatasetEntry, ROIRecord

        entry = IMCDatasetEntry(
            name="test",
            base_dir="/tmp",
            tissue="lung",
            rois=[
                ROIRecord("test", "lung", "s1", "r1", "/tmp/r1.tiff"),
                ROIRecord("test", "lung", "s1", "r2", "/tmp/r2.tiff"),
            ],
        )
        assert entry.n_rois == 2
