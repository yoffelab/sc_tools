"""Unit tests for sc_tools.ingest module."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from sc_tools.ingest.config import (
    collect_all_batches,
    load_batch_manifest,
    validate_manifest,
)
from sc_tools.ingest.imc import build_imc_pipeline_cmd
from sc_tools.ingest.loaders import (
    concat_samples,
    load_imc_sample,
    load_visium_hd_sample,
    load_visium_sample,
    load_xenium_sample,
)
from sc_tools.ingest.spaceranger import (
    build_batch_commands,
    build_spaceranger_count_cmd,
)
from sc_tools.ingest.xenium import build_xenium_ranger_cmd

# ---------------------------------------------------------------------------
# config.py tests
# ---------------------------------------------------------------------------


class TestLoadBatchManifest:
    def test_load_tsv(self, tmp_path):
        tsv = tmp_path / "batch1_samples.tsv"
        tsv.write_text(
            "sample_id\tfastq_dir\tcytaimage\tslide\tarea\n"
            "S1\t/data/S1\t/img/S1.tif\tH1-ABC\tD1\n"
            "S2\t/data/S2\t/img/S2.tif\tH1-ABC\tA1\n"
        )
        df = load_batch_manifest(tsv)
        assert len(df) == 2
        assert list(df.columns) == ["sample_id", "fastq_dir", "cytaimage", "slide", "area"]

    def test_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_batch_manifest(tmp_path / "nonexistent.tsv")


class TestCollectAllBatches:
    def test_collect_two_batches(self, tmp_path):
        phase0 = tmp_path / "phase0"
        phase0.mkdir()

        (phase0 / "batch1_samples.tsv").write_text(
            "sample_id\tfastq_dir\nS1\t/data/S1\nS2\t/data/S2\n"
        )
        (phase0 / "batch2_samples.tsv").write_text("sample_id\tfastq_dir\nS3\t/data/S3\n")

        df = collect_all_batches(phase0)
        assert len(df) == 3
        assert "batch" in df.columns
        assert set(df["batch"]) == {"batch1", "batch2"}
        assert (phase0 / "all_samples.tsv").exists()

    def test_empty_dir(self, tmp_path):
        phase0 = tmp_path / "phase0"
        phase0.mkdir()
        df = collect_all_batches(phase0)
        assert len(df) == 0

    def test_skips_all_samples_tsv(self, tmp_path):
        phase0 = tmp_path / "phase0"
        phase0.mkdir()
        (phase0 / "batch1_samples.tsv").write_text("sample_id\nS1\n")
        (phase0 / "all_samples.tsv").write_text("sample_id\nOLD\n")
        df = collect_all_batches(phase0)
        assert len(df) == 1
        assert df.iloc[0]["sample_id"] == "S1"

    def test_preserves_existing_batch_column(self, tmp_path):
        phase0 = tmp_path / "phase0"
        phase0.mkdir()
        (phase0 / "custom_samples.tsv").write_text("sample_id\tbatch\nS1\tRobin_batch2\n")
        df = collect_all_batches(phase0)
        assert df.iloc[0]["batch"] == "Robin_batch2"

    def test_warns_duplicate_sample_ids(self, tmp_path, caplog):
        phase0 = tmp_path / "phase0"
        phase0.mkdir()
        (phase0 / "batch1_samples.tsv").write_text("sample_id\nS1\n")
        (phase0 / "batch2_samples.tsv").write_text("sample_id\nS1\n")
        import logging

        with caplog.at_level(logging.WARNING):
            df = collect_all_batches(phase0)
        assert len(df) == 2
        assert "Duplicate" in caplog.text


class TestValidateManifest:
    def test_valid_visium_hd(self):
        df = pd.DataFrame(
            {
                "sample_id": ["S1"],
                "fastq_dir": ["/data"],
                "cytaimage": ["/img.tif"],
                "slide": ["H1"],
                "area": ["D1"],
            }
        )
        assert validate_manifest(df, "visium_hd") == []

    def test_missing_columns(self):
        df = pd.DataFrame({"sample_id": ["S1"]})
        issues = validate_manifest(df, "visium_hd")
        assert len(issues) == 1
        assert "Missing required columns" in issues[0]

    def test_unknown_modality(self):
        df = pd.DataFrame()
        issues = validate_manifest(df, "merfish")
        assert "Unknown modality" in issues[0]

    def test_cosmx_no_requirements(self):
        df = pd.DataFrame({"anything": [1]})
        assert validate_manifest(df, "cosmx") == []

    def test_null_sample_id(self):
        df = pd.DataFrame(
            {
                "sample_id": ["S1", None],
                "xenium_dir": ["/a", "/b"],
            }
        )
        issues = validate_manifest(df, "xenium")
        assert any("null" in i for i in issues)


# ---------------------------------------------------------------------------
# spaceranger.py tests
# ---------------------------------------------------------------------------


class TestBuildSpacerangerCmd:
    def test_visium_hd_cmd(self):
        cmd = build_spaceranger_count_cmd(
            sample_id="S1",
            fastqs="/data/fastqs",
            transcriptome="/ref/GRCh38",
            cytaimage="/img/cyto.tif",
            slide="H1-ABC",
            area="D1",
        )
        assert "spaceranger count" in cmd
        assert "--id=S1" in cmd
        assert "--cytaimage=/img/cyto.tif" in cmd
        assert "--slide=H1-ABC" in cmd

    def test_visium_cmd(self):
        cmd = build_spaceranger_count_cmd(
            sample_id="S1",
            fastqs="/data/fastqs",
            transcriptome="/ref/GRCh38",
            image="/img/he.tif",
        )
        assert "--image=/img/he.tif" in cmd

    def test_no_image_raises(self):
        with pytest.raises(ValueError, match="image"):
            build_spaceranger_count_cmd(
                sample_id="S1",
                fastqs="/data",
                transcriptome="/ref",
            )

    def test_custom_resources(self):
        cmd = build_spaceranger_count_cmd(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            cytaimage="/img.tif",
            localcores=32,
            localmem=128,
        )
        assert "--localcores=32" in cmd
        assert "--localmem=128" in cmd


class TestBuildBatchCommands:
    def test_batch_commands(self):
        manifest = pd.DataFrame(
            {
                "sample_id": ["S1", "S2"],
                "fastq_dir": ["/data/S1", "/data/S2"],
                "cytaimage": ["/img/S1.tif", "/img/S2.tif"],
                "slide": ["H1", "H1"],
                "area": ["D1", "A1"],
            }
        )
        cmds = build_batch_commands(manifest, transcriptome="/ref/GRCh38")
        assert len(cmds) == 2
        assert cmds[0][0] == "S1"
        assert "--id=S1" in cmds[0][1]

    def test_skips_invalid_rows(self):
        manifest = pd.DataFrame(
            {
                "sample_id": ["S1"],
                "fastq_dir": ["/data/S1"],
                # No image or cytaimage -> should skip
            }
        )
        cmds = build_batch_commands(manifest, transcriptome="/ref")
        assert len(cmds) == 0


# ---------------------------------------------------------------------------
# xenium.py tests
# ---------------------------------------------------------------------------


class TestXeniumCmd:
    def test_basic_cmd(self):
        cmd = build_xenium_ranger_cmd("X1", "/raw/xenium")
        assert "xeniumranger resegment" in cmd
        assert "--id=X1" in cmd


# ---------------------------------------------------------------------------
# imc.py tests
# ---------------------------------------------------------------------------


class TestImcCmd:
    def test_basic_cmd(self):
        cmd = build_imc_pipeline_cmd("/data/sample.mcd", "/panel.csv", "/output")
        assert "--mcd /data/sample.mcd" in cmd
        assert "--panel /panel.csv" in cmd

    def test_custom_pipeline_dir(self):
        cmd = build_imc_pipeline_cmd(
            "/data/sample.mcd", "/panel.csv", "/output", pipeline_dir="/custom/imc"
        )
        assert "/custom/imc/run_pipeline.py" in cmd


# ---------------------------------------------------------------------------
# loaders.py tests
# ---------------------------------------------------------------------------


class TestConcatSamples:
    def _make_adata(self, sample_id, n_obs=20, n_vars=10):
        X = np.random.poisson(5, (n_obs, n_vars)).astype(np.float32)
        obs = pd.DataFrame(
            {"sample": sample_id},
            index=[f"{sample_id}_cell_{i}" for i in range(n_obs)],
        )
        var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.obsm["spatial"] = np.random.rand(n_obs, 2)
        return adata

    def test_concat_two(self):
        a1 = self._make_adata("S1")
        a2 = self._make_adata("S2")
        combined = concat_samples([a1, a2], calculate_qc=False)
        assert combined.n_obs == 40
        assert "spatial" in combined.obsm

    def test_single_sample(self):
        a1 = self._make_adata("S1")
        combined = concat_samples([a1], calculate_qc=False)
        assert combined.n_obs == 20

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="No AnnData"):
            concat_samples([])

    def test_qc_metrics_calculated(self):
        a1 = self._make_adata("S1", n_vars=600)
        combined = concat_samples([a1], calculate_qc=True)
        assert "total_counts" in combined.obs.columns


class TestLoadImcSample:
    def test_load_from_h5ad(self, tmp_path):
        n_obs, n_vars = 30, 20
        X = np.random.rand(n_obs, n_vars).astype(np.float32)
        obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
        var = pd.DataFrame(index=[f"marker_{i}" for i in range(n_vars)])
        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.obsm["spatial"] = np.random.rand(n_obs, 2)

        path = tmp_path / "roi1.h5ad"
        adata.write_h5ad(path)

        loaded = load_imc_sample(path, "ROI_1")
        assert loaded.obs["sample"].iloc[0] == "ROI_1"
        assert loaded.obs["raw_data_dir"].iloc[0] == str(path)
        assert "spatial" in loaded.obsm

    def test_x_spatial_renamed(self, tmp_path):
        """obs['X_spatial'] should be copied to obsm['spatial'] if spatial missing."""
        n_obs, n_vars = 10, 5
        adata = ad.AnnData(
            X=np.random.rand(n_obs, n_vars).astype(np.float32),
            obs=pd.DataFrame(index=[f"c_{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=[f"m_{i}" for i in range(n_vars)]),
        )
        adata.obsm["X_spatial"] = np.random.rand(n_obs, 2)

        path = tmp_path / "roi2.h5ad"
        adata.write_h5ad(path)

        loaded = load_imc_sample(path, "ROI_2")
        assert "spatial" in loaded.obsm
        np.testing.assert_array_equal(loaded.obsm["spatial"], loaded.obsm["X_spatial"])


# ---------------------------------------------------------------------------
# Mock-based loader tests: Visium, Visium HD, Xenium
# ---------------------------------------------------------------------------


def _make_mock_adata(n_obs=30, n_vars=50):
    """Create a synthetic AnnData that mimics what scanpy/squidpy would return."""
    X = np.random.poisson(5, (n_obs, n_vars)).astype(np.float32)
    obs = pd.DataFrame(index=[f"BARCODE_{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=[f"Gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["spatial"] = np.random.rand(n_obs, 2)
    return adata


class TestLoadVisiumSample:
    """Mock-based tests for load_visium_sample (Visium)."""

    def test_sets_required_obs_keys(self, tmp_path):
        sr_dir = tmp_path / "sample1"
        (sr_dir / "outs").mkdir(parents=True)

        mock_adata = _make_mock_adata()
        with patch("scanpy.read_visium", return_value=mock_adata) as mock_read:
            result = load_visium_sample(sr_dir, "VIS_01")
            mock_read.assert_called_once_with(sr_dir / "outs", load_images=True)

        assert (result.obs["sample"] == "VIS_01").all()
        assert (result.obs["library_id"] == "VIS_01").all()
        assert (result.obs["raw_data_dir"] == str(sr_dir)).all()

    def test_spatial_preserved(self, tmp_path):
        sr_dir = tmp_path / "sample2"
        (sr_dir / "outs").mkdir(parents=True)

        mock_adata = _make_mock_adata()
        with patch("scanpy.read_visium", return_value=mock_adata):
            result = load_visium_sample(sr_dir, "VIS_02")

        assert "spatial" in result.obsm
        assert result.obsm["spatial"].shape == (30, 2)

    def test_no_outs_dir_uses_root(self, tmp_path):
        """If outs/ does not exist, use spaceranger_dir directly."""
        sr_dir = tmp_path / "sample3"
        sr_dir.mkdir(parents=True)
        # No outs/ subdirectory

        mock_adata = _make_mock_adata()
        with patch("scanpy.read_visium", return_value=mock_adata) as mock_read:
            load_visium_sample(sr_dir, "VIS_03")
            mock_read.assert_called_once_with(sr_dir, load_images=True)

    def test_load_images_false(self, tmp_path):
        sr_dir = tmp_path / "sample4"
        (sr_dir / "outs").mkdir(parents=True)

        mock_adata = _make_mock_adata()
        with patch("scanpy.read_visium", return_value=mock_adata) as mock_read:
            load_visium_sample(sr_dir, "VIS_04", load_images=False)
            mock_read.assert_called_once_with(sr_dir / "outs", load_images=False)

    def test_var_names_made_unique(self, tmp_path):
        sr_dir = tmp_path / "sample5"
        (sr_dir / "outs").mkdir(parents=True)

        mock_adata = _make_mock_adata(n_vars=3)
        mock_adata.var_names = pd.Index(["GeneA", "GeneA", "GeneB"])
        with patch("scanpy.read_visium", return_value=mock_adata):
            result = load_visium_sample(sr_dir, "VIS_05")

        assert result.var_names.is_unique


class TestLoadVisiumHdSample:
    """Mock-based tests for load_visium_hd_sample (Visium HD)."""

    def test_sets_required_obs_keys(self, tmp_path):
        sr_dir = tmp_path / "hd_sample1"
        bin_dir = sr_dir / "square_008um"
        bin_dir.mkdir(parents=True)

        mock_adata = _make_mock_adata()
        with patch("squidpy.read.visium", return_value=mock_adata):
            result = load_visium_hd_sample(sr_dir, "HD_01")

        assert (result.obs["sample"] == "HD_01").all()
        assert (result.obs["library_id"] == "HD_01").all()
        assert (result.obs["raw_data_dir"] == str(sr_dir)).all()

    def test_parquet_positions_loaded(self, tmp_path):
        sr_dir = tmp_path / "hd_sample2"
        bin_dir = sr_dir / "square_008um"
        spatial_dir = bin_dir / "spatial"
        spatial_dir.mkdir(parents=True)

        n_obs = 20
        barcodes = [f"BARCODE_{i}" for i in range(n_obs)]
        mock_adata = _make_mock_adata(n_obs=n_obs)
        mock_adata.obs_names = pd.Index(barcodes)

        # Write parquet with tissue positions
        pos_df = pd.DataFrame(
            {
                "barcode": barcodes,
                "pxl_col_in_fullres": np.random.randint(0, 5000, n_obs),
                "pxl_row_in_fullres": np.random.randint(0, 5000, n_obs),
                "in_tissue": [1] * n_obs,
            }
        )
        pos_df.to_parquet(spatial_dir / "tissue_positions.parquet")

        with patch("squidpy.read.visium", return_value=mock_adata):
            result = load_visium_hd_sample(sr_dir, "HD_02")

        assert "spatial" in result.obsm
        assert result.obsm["spatial"].shape == (n_obs, 2)
        assert "pxl_col_in_fullres" in result.obs.columns

    def test_falls_back_to_outs_bin_dir(self, tmp_path):
        sr_dir = tmp_path / "hd_sample3"
        outs_bin = sr_dir / "outs" / "square_008um"
        outs_bin.mkdir(parents=True)

        mock_adata = _make_mock_adata()
        with patch("squidpy.read.visium", return_value=mock_adata) as mock_read:
            load_visium_hd_sample(sr_dir, "HD_03")
            mock_read.assert_called_once_with(
                outs_bin,
                load_images=False,
                counts_file="filtered_feature_bc_matrix.h5",
            )

    def test_custom_bin_size(self, tmp_path):
        sr_dir = tmp_path / "hd_sample4"
        bin_dir = sr_dir / "square_016um"
        bin_dir.mkdir(parents=True)

        mock_adata = _make_mock_adata()
        with patch("squidpy.read.visium", return_value=mock_adata) as mock_read:
            load_visium_hd_sample(sr_dir, "HD_04", bin_size="square_016um")
            mock_read.assert_called_once_with(
                bin_dir,
                load_images=False,
                counts_file="filtered_feature_bc_matrix.h5",
            )

    def test_bin_dir_not_found_raises(self, tmp_path):
        sr_dir = tmp_path / "hd_sample5"
        sr_dir.mkdir(parents=True)
        # No bin directory at all

        with pytest.raises(FileNotFoundError, match="Bin directory not found"):
            load_visium_hd_sample(sr_dir, "HD_05")

    def test_no_parquet_still_works(self, tmp_path):
        """If parquet does not exist, spatial from squidpy is kept."""
        sr_dir = tmp_path / "hd_sample6"
        bin_dir = sr_dir / "square_008um"
        bin_dir.mkdir(parents=True)

        mock_adata = _make_mock_adata()
        with patch("squidpy.read.visium", return_value=mock_adata):
            result = load_visium_hd_sample(sr_dir, "HD_06")

        assert (result.obs["sample"] == "HD_06").all()


class TestLoadXeniumSample:
    """Mock-based tests for load_xenium_sample (Xenium)."""

    def test_spatialdata_io_path(self, tmp_path):
        """Test the spatialdata-io loading path."""
        xenium_dir = tmp_path / "xenium1"
        xenium_dir.mkdir()

        mock_adata = _make_mock_adata()
        mock_sdata = MagicMock()
        mock_sdata.tables = {"table": mock_adata}

        with patch.dict("sys.modules", {"spatialdata_io": MagicMock()}) as _:
            import sys

            mock_sdio = sys.modules["spatialdata_io"]
            mock_sdio.xenium.return_value = mock_sdata

            result = load_xenium_sample(xenium_dir, "XEN_01")

        assert (result.obs["sample"] == "XEN_01").all()
        assert (result.obs["library_id"] == "XEN_01").all()
        assert (result.obs["raw_data_dir"] == str(xenium_dir)).all()

    def test_scanpy_fallback_with_h5(self, tmp_path):
        """When spatialdata-io not available, falls back to scanpy read_10x_h5."""
        xenium_dir = tmp_path / "xenium2"
        xenium_dir.mkdir()
        # Create the h5 file path so the exists() check passes
        (xenium_dir / "cell_feature_matrix.h5").touch()

        mock_adata = _make_mock_adata()

        with (
            patch.dict("sys.modules", {"spatialdata_io": None}),
            patch("scanpy.read_10x_h5", return_value=mock_adata) as mock_read,
        ):
            result = load_xenium_sample(xenium_dir, "XEN_02")
            mock_read.assert_called_once_with(str(xenium_dir / "cell_feature_matrix.h5"))

        assert (result.obs["sample"] == "XEN_02").all()

    def test_scanpy_fallback_with_coords_csv(self, tmp_path):
        """When cells.csv exists, spatial coordinates are loaded from it."""
        xenium_dir = tmp_path / "xenium3"
        xenium_dir.mkdir()
        (xenium_dir / "cell_feature_matrix.h5").touch()

        n_obs = 25
        mock_adata = _make_mock_adata(n_obs=n_obs)
        # Remove spatial so we can test CSV loading
        del mock_adata.obsm["spatial"]

        coords_df = pd.DataFrame(
            {
                "cell_id": range(n_obs),
                "x_centroid": np.random.rand(n_obs) * 1000,
                "y_centroid": np.random.rand(n_obs) * 1000,
            }
        )
        coords_df.to_csv(xenium_dir / "cells.csv", index=False)

        with (
            patch.dict("sys.modules", {"spatialdata_io": None}),
            patch("scanpy.read_10x_h5", return_value=mock_adata),
        ):
            result = load_xenium_sample(xenium_dir, "XEN_03")

        assert "spatial" in result.obsm
        assert result.obsm["spatial"].shape == (n_obs, 2)

    def test_scanpy_fallback_mtx(self, tmp_path):
        """When no h5, falls back to read_10x_mtx."""
        xenium_dir = tmp_path / "xenium4"
        xenium_dir.mkdir()
        # No h5 file exists

        mock_adata = _make_mock_adata()

        with (
            patch.dict("sys.modules", {"spatialdata_io": None}),
            patch("scanpy.read_10x_mtx", return_value=mock_adata) as mock_mtx,
        ):
            result = load_xenium_sample(xenium_dir, "XEN_04")
            mock_mtx.assert_called_once_with(str(xenium_dir / "cell_feature_matrix"))

        assert (result.obs["sample"] == "XEN_04").all()

    def test_var_names_made_unique(self, tmp_path):
        xenium_dir = tmp_path / "xenium5"
        xenium_dir.mkdir()
        (xenium_dir / "cell_feature_matrix.h5").touch()

        mock_adata = _make_mock_adata(n_vars=3)
        mock_adata.var_names = pd.Index(["BRCA1", "BRCA1", "TP53"])

        with (
            patch.dict("sys.modules", {"spatialdata_io": None}),
            patch("scanpy.read_10x_h5", return_value=mock_adata),
        ):
            result = load_xenium_sample(xenium_dir, "XEN_05")

        assert result.var_names.is_unique
