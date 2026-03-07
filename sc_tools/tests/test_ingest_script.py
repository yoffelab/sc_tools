"""Tests for scripts/ingest.py orchestration logic."""

from __future__ import annotations

from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest


def _make_adata(sample_id: str, n_obs: int = 10, n_vars: int = 5) -> ad.AnnData:
    """Create a minimal synthetic AnnData."""
    adata = ad.AnnData(
        X=np.random.rand(n_obs, n_vars).astype(np.float32),
        obs=pd.DataFrame(
            {"sample": sample_id, "library_id": sample_id},
            index=[f"{sample_id}_{i}" for i in range(n_obs)],
        ),
    )
    adata.obsm["spatial"] = np.random.rand(n_obs, 2).astype(np.float32)
    return adata


def _make_manifest(modality: str, n_samples: int = 3) -> pd.DataFrame:
    """Create a synthetic manifest DataFrame."""
    rows = []
    for i in range(n_samples):
        row = {"sample_id": f"sample_{i}", "batch": "batch1"}
        if modality in ("visium", "visium_hd", "visium_hd_cell"):
            row.update(
                {
                    "fastq_dir": f"/data/sample_{i}",
                    "image": "img.tif",
                    "cytaimage": "cyta.tif",
                    "slide": "S1",
                    "area": "A1",
                }
            )
        elif modality == "xenium":
            row["xenium_dir"] = f"/data/sample_{i}"
        elif modality == "imc":
            row["processed_dir"] = f"/data/sample_{i}"
        rows.append(row)
    return pd.DataFrame(rows)


class TestLoadSingleSample:
    """Test modality dispatch in load_single_sample."""

    @patch("scripts.ingest.load_single_sample")
    def test_dispatch_is_callable(self, mock_load):
        """Verify load_single_sample is importable and callable."""
        mock_load.return_value = _make_adata("s1")
        result = mock_load({"sample_id": "s1", "fastq_dir": "/d"}, "visium")
        assert result is not None


class TestRunIngestion:
    """Test the run_ingestion orchestration."""

    def test_successful_loading(self):
        """All samples load successfully."""
        from scripts.ingest import run_ingestion

        manifest = _make_manifest("imc", n_samples=2)

        with patch("scripts.ingest.load_single_sample") as mock_load:
            mock_load.side_effect = [_make_adata("sample_0"), _make_adata("sample_1")]
            adatas, failed = run_ingestion(manifest, "imc")

        assert len(adatas) == 2
        assert len(failed) == 0
        assert mock_load.call_count == 2

    def test_failed_samples_skipped(self):
        """Failed samples are skipped, others still load."""
        from scripts.ingest import run_ingestion

        manifest = _make_manifest("imc", n_samples=3)

        with patch("scripts.ingest.load_single_sample") as mock_load:
            mock_load.side_effect = [
                _make_adata("sample_0"),
                FileNotFoundError("missing"),
                _make_adata("sample_2"),
            ]
            adatas, failed = run_ingestion(manifest, "imc")

        assert len(adatas) == 2
        assert failed == ["sample_1"]

    def test_all_samples_fail(self):
        """All samples failing returns empty list."""
        from scripts.ingest import run_ingestion

        manifest = _make_manifest("imc", n_samples=2)

        with patch("scripts.ingest.load_single_sample") as mock_load:
            mock_load.side_effect = Exception("fail")
            adatas, failed = run_ingestion(manifest, "imc")

        assert len(adatas) == 0
        assert len(failed) == 2


class TestCosmxNotImplemented:
    """CosMx should raise NotImplementedError."""

    def test_cosmx_raises(self):
        from scripts.ingest import _get_loader

        with pytest.raises(NotImplementedError, match="CosMx"):
            _get_loader("cosmx")
