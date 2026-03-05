"""Unit tests for sc_tools.ingest module."""

from __future__ import annotations

from pathlib import Path
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
from sc_tools.ingest.imc import IMCPanelMapper, build_imc_composite, build_imc_pipeline_cmd
from sc_tools.ingest.loaders import (
    concat_samples,
    load_he_image,
    load_imc_sample,
    load_visium_hd_sample,
    load_visium_sample,
    load_xenium_sample,
)
from sc_tools.ingest.slurm import (
    build_batch_sbatch,
    build_imc_sbatch,
    build_sbatch_header,
    build_spaceranger_sbatch,
    build_xenium_sbatch,
    write_sbatch_script,
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

    def test_cosmx_requires_cosmx_dir(self):
        df = pd.DataFrame({"sample_id": ["S1"], "cosmx_dir": ["/path/to/cosmx"]})
        assert validate_manifest(df, "cosmx") == []

    def test_cosmx_missing_cosmx_dir(self):
        df = pd.DataFrame({"sample_id": ["S1"]})
        issues = validate_manifest(df, "cosmx")
        assert len(issues) == 1
        assert "cosmx_dir" in issues[0]

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

    def test_probe_set(self):
        cmd = build_spaceranger_count_cmd(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            cytaimage="/img.tif",
            probe_set="/probes/v2.csv",
        )
        assert "--probe-set=/probes/v2.csv" in cmd

    def test_sample_filter(self):
        cmd = build_spaceranger_count_cmd(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            cytaimage="/img.tif",
            sample_filter="PT01-1_NAT",
        )
        assert "--sample=PT01-1_NAT" in cmd

    def test_create_bam_true(self):
        cmd = build_spaceranger_count_cmd(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            cytaimage="/img.tif",
            create_bam=True,
        )
        assert "--create-bam=true" in cmd

    def test_create_bam_false_default(self):
        cmd = build_spaceranger_count_cmd(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            cytaimage="/img.tif",
        )
        assert "--create-bam=false" in cmd

    def test_custom_spaceranger_path(self):
        cmd = build_spaceranger_count_cmd(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            cytaimage="/img.tif",
            spaceranger_path="/opt/sr4/spaceranger",
        )
        assert cmd.startswith("/opt/sr4/spaceranger count")
        assert "spaceranger count" not in cmd.replace("/opt/sr4/spaceranger count", "")


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
    def test_load_from_processed_dir(self, tmp_path):
        n_obs, n_vars = 30, 20
        X = np.random.rand(n_obs, n_vars).astype(np.float32)
        obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
        var = pd.DataFrame(index=[f"marker_{i}" for i in range(n_vars)])
        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.obsm["spatial"] = np.random.rand(n_obs, 2)

        processed_dir = tmp_path / "processed" / "sample1"
        processed_dir.mkdir(parents=True)
        adata.write_h5ad(processed_dir / "cells.h5ad")

        loaded = load_imc_sample(processed_dir, "ROI_1")
        assert loaded.obs["sample"].iloc[0] == "ROI_1"
        assert loaded.obs["raw_data_dir"].iloc[0] == str(processed_dir)
        assert "spatial" in loaded.obsm

    def test_x_spatial_renamed(self, tmp_path):
        """obsm['X_spatial'] should be copied to obsm['spatial'] if spatial missing."""
        n_obs, n_vars = 10, 5
        adata = ad.AnnData(
            X=np.random.rand(n_obs, n_vars).astype(np.float32),
            obs=pd.DataFrame(index=[f"c_{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=[f"m_{i}" for i in range(n_vars)]),
        )
        adata.obsm["X_spatial"] = np.random.rand(n_obs, 2)

        processed_dir = tmp_path / "processed" / "sample2"
        processed_dir.mkdir(parents=True)
        adata.write_h5ad(processed_dir / "cells.h5ad")

        loaded = load_imc_sample(processed_dir, "ROI_2")
        assert "spatial" in loaded.obsm
        np.testing.assert_array_equal(loaded.obsm["spatial"], loaded.obsm["X_spatial"])

    def test_missing_cells_h5ad_raises(self, tmp_path):
        """FileNotFoundError when cells.h5ad is absent."""
        processed_dir = tmp_path / "processed" / "empty"
        processed_dir.mkdir(parents=True)
        with pytest.raises(FileNotFoundError, match="cells.h5ad"):
            load_imc_sample(processed_dir, "ROI_3")


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


# ---------------------------------------------------------------------------
# IMCPanelMapper tests
# ---------------------------------------------------------------------------


def _write_full_csv(tmp_path, rows: list[tuple[int, str]]) -> Path:
    """Write a *_full.csv style channel index file."""
    csv_path = tmp_path / "roi_full.csv"
    lines = [",channel"] + [f"{idx},{name}" for idx, name in rows]
    csv_path.write_text("\n".join(lines))
    return csv_path


class TestIMCPanelMapper:
    CHANNELS = [
        (0, "HH3(In113)"),
        (1, "CD45RA(Nd143)"),
        (2, "CD163(Sm147)"),
        (3, "CD14(Nd148)"),
        (4, "CD3(Er170)"),
        (5, "DNA1(Ir191)"),
        (6, "PanCK(Pt195)"),
        (7, "CD68(Pt196)"),
        (8, "<EMPTY>(In115)"),
    ]

    def test_resolve_exact_protein(self, tmp_path):
        csv = _write_full_csv(tmp_path, self.CHANNELS)
        mapper = IMCPanelMapper(full_csv=csv)
        assert mapper.resolve("CD3") == 4
        assert mapper.resolve("PanCK") == 6
        assert mapper.resolve("DNA1") == 5

    def test_resolve_case_insensitive(self, tmp_path):
        csv = _write_full_csv(tmp_path, self.CHANNELS)
        mapper = IMCPanelMapper(full_csv=csv)
        assert mapper.resolve("cd3") == 4
        assert mapper.resolve("PANCK") == 6

    def test_resolve_isotope_tag(self, tmp_path):
        csv = _write_full_csv(tmp_path, self.CHANNELS)
        mapper = IMCPanelMapper(full_csv=csv)
        assert mapper.resolve("Er170") == 4
        assert mapper.resolve("Ir191") == 5

    def test_resolve_full_string(self, tmp_path):
        csv = _write_full_csv(tmp_path, self.CHANNELS)
        mapper = IMCPanelMapper(full_csv=csv)
        assert mapper.resolve("CD3(Er170)") == 4

    def test_resolve_partial_match(self, tmp_path):
        csv = _write_full_csv(tmp_path, self.CHANNELS)
        mapper = IMCPanelMapper(full_csv=csv)
        # CD45RA matches partial 'CD45'
        assert mapper.resolve("CD45") == 1

    def test_resolve_missing_returns_none(self, tmp_path):
        csv = _write_full_csv(tmp_path, self.CHANNELS)
        mapper = IMCPanelMapper(full_csv=csv)
        assert mapper.resolve("NONEXISTENT") is None

    def test_resolve_empty_channel_returns_none(self, tmp_path):
        csv = _write_full_csv(tmp_path, self.CHANNELS)
        mapper = IMCPanelMapper(full_csv=csv)
        assert mapper.resolve("<EMPTY>") is None

    def test_build_rgb_indices(self, tmp_path):
        csv = _write_full_csv(tmp_path, self.CHANNELS)
        mapper = IMCPanelMapper(full_csv=csv)
        rgb = mapper.build_rgb_indices(r="PanCK", g="CD3", b="DNA1")
        assert rgb == {"R": 6, "G": 4, "B": 5}

    def test_build_rgb_indices_missing_channel(self, tmp_path):
        csv = _write_full_csv(tmp_path, self.CHANNELS)
        mapper = IMCPanelMapper(full_csv=csv)
        rgb = mapper.build_rgb_indices(r="NONEXISTENT", g="CD3", b="DNA1")
        assert rgb["R"] is None
        assert rgb["G"] == 4

    def test_channel_names(self, tmp_path):
        csv = _write_full_csv(tmp_path, self.CHANNELS)
        mapper = IMCPanelMapper(full_csv=csv)
        names = mapper.channel_names()
        assert names[4] == "CD3"
        assert names[5] == "DNA1"

    def test_n_channels(self, tmp_path):
        csv = _write_full_csv(tmp_path, self.CHANNELS)
        mapper = IMCPanelMapper(full_csv=csv)
        assert mapper.n_channels() == len(self.CHANNELS)

    def test_set_from_var_names(self):
        mapper = IMCPanelMapper()
        mapper.set_from_var_names(["CD3(Er170)", "PanCK(Pt195)", "DNA1(Ir191)"])
        assert mapper.resolve("CD3") == 0
        assert mapper.resolve("PanCK") == 1
        assert mapper.resolve("DNA1") == 2

    def test_panel_csv_ilastik_flags(self, tmp_path):
        panel_path = tmp_path / "channel_labels.csv"
        panel_path.write_text(
            "channel,Target,Metal_Tag,Atom,full,ilastik\n"
            "CD3(Er170),CD3,Er170,170,1,1\n"
            "PanCK(Pt195),PanCK,Pt195,195,1,0\n"
            "DNA1(Ir191),DNA1,Ir191,191,1,1\n"
        )
        mapper = IMCPanelMapper(panel_csv=panel_path)
        ilastik_idx = mapper.get_ilastik_indices()
        assert len(ilastik_idx) == 2


# ---------------------------------------------------------------------------
# build_imc_composite tests
# ---------------------------------------------------------------------------


def _write_synthetic_tiff(tmp_path: Path, n_ch: int, h: int, w: int) -> Path:
    """Write a synthetic (C, H, W) uint16 TIFF for testing."""
    try:
        import tifffile
    except ImportError:
        pytest.skip("tifffile not installed")
    arr = (np.random.rand(n_ch, h, w) * 1000).astype(np.uint16)
    tiff_path = tmp_path / "roi_full.tiff"
    tifffile.imwrite(str(tiff_path), arr)
    return tiff_path


class TestBuildImcComposite:
    CHANNELS = [
        (0, "HH3(In113)"),
        (1, "CD45RA(Nd143)"),
        (2, "CD3(Er170)"),
        (3, "DNA1(Ir191)"),
        (4, "PanCK(Pt195)"),
        (5, "CD68(Pt196)"),
    ]

    def test_basic_composite(self, tmp_path):
        tiff_path = _write_synthetic_tiff(tmp_path, n_ch=6, h=64, w=64)
        csv_path = _write_full_csv(tmp_path, self.CHANNELS)

        result = build_imc_composite(
            tiff_path=tiff_path,
            channel_csv=csv_path,
            r="PanCK",
            g="CD3",
            b="DNA1",
        )

        assert "images" in result
        assert "hires" in result["images"]
        assert "full" in result["images"]
        assert result["images"]["hires"].shape == (64, 64, 3)
        assert result["images"]["hires"].dtype == np.uint8
        assert result["images"]["full"].shape == (6, 64, 64)
        assert result["images"]["full"].dtype == np.float32

    def test_rgb_indices_stored(self, tmp_path):
        tiff_path = _write_synthetic_tiff(tmp_path, n_ch=6, h=32, w=32)
        csv_path = _write_full_csv(tmp_path, self.CHANNELS)
        result = build_imc_composite(tiff_path, csv_path, r="PanCK", g="CD3", b="DNA1")

        rgb_ch = result["metadata"]["rgb_channels"]
        assert rgb_ch["R"] == "PanCK"
        assert rgb_ch["G"] == "CD3"
        assert rgb_ch["B"] == "DNA1"

        rgb_idx = result["metadata"]["rgb_indices"]
        assert rgb_idx["R"] == 4
        assert rgb_idx["G"] == 2
        assert rgb_idx["B"] == 3

    def test_scalefactors_no_downsample(self, tmp_path):
        tiff_path = _write_synthetic_tiff(tmp_path, n_ch=6, h=32, w=32)
        csv_path = _write_full_csv(tmp_path, self.CHANNELS)
        result = build_imc_composite(tiff_path, csv_path)
        assert result["scalefactors"]["tissue_hires_scalef"] == 1.0

    def test_downsample(self, tmp_path):
        tiff_path = _write_synthetic_tiff(tmp_path, n_ch=6, h=64, w=64)
        csv_path = _write_full_csv(tmp_path, self.CHANNELS)
        result = build_imc_composite(tiff_path, csv_path, downsample=2)
        assert result["images"]["hires"].shape == (32, 32, 3)
        assert result["scalefactors"]["tissue_hires_scalef"] == 0.5

    def test_missing_channel_in_rgb(self, tmp_path):
        """Missing RGB channel produces a zero plane, not an error."""
        tiff_path = _write_synthetic_tiff(tmp_path, n_ch=6, h=32, w=32)
        csv_path = _write_full_csv(tmp_path, self.CHANNELS)
        result = build_imc_composite(tiff_path, csv_path, r="NONEXISTENT", g="CD3", b="DNA1")
        assert result["images"]["hires"][..., 0].max() == 0

    def test_arcsinh_normalization(self, tmp_path):
        tiff_path = _write_synthetic_tiff(tmp_path, n_ch=6, h=32, w=32)
        csv_path = _write_full_csv(tmp_path, self.CHANNELS)
        result = build_imc_composite(tiff_path, csv_path)
        assert result["images"]["full"].min() >= 0

    def test_tiff_not_found_raises(self, tmp_path):
        csv_path = _write_full_csv(tmp_path, self.CHANNELS)
        with pytest.raises(FileNotFoundError, match="TIFF stack"):
            build_imc_composite(tmp_path / "missing.tiff", csv_path)

    def test_csv_not_found_raises(self, tmp_path):
        tiff_path = _write_synthetic_tiff(tmp_path, n_ch=6, h=32, w=32)
        with pytest.raises(FileNotFoundError, match="Channel CSV"):
            build_imc_composite(tiff_path, tmp_path / "missing.csv")


# ---------------------------------------------------------------------------
# load_imc_sample with load_images tests
# ---------------------------------------------------------------------------


class TestLoadImcSampleImages:
    """Tests for load_imc_sample with load_images=True."""

    def _make_imc_processed_dir(self, tmp_path, sample_id: str, n_markers: int = 5):
        """Create a minimal processed IMC directory with cells.h5ad and tiffs/."""
        processed_dir = tmp_path / "processed" / sample_id.split("-")[0]
        tiff_dir = processed_dir / "tiffs"
        tiff_dir.mkdir(parents=True)

        # Write cells.h5ad
        n_obs = 20
        X = np.random.rand(n_obs, n_markers).astype(np.float32)
        markers = ["CD3", "PanCK", "DNA1", "CD68", "CD45RA"][:n_markers]
        adata_inner = ad.AnnData(
            X=X,
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=markers),
        )
        adata_inner.obsm["spatial"] = np.random.rand(n_obs, 2) * 100
        adata_inner.write_h5ad(processed_dir / "cells.h5ad")

        # Write channel CSV
        channel_strings = [
            "CD3(Er170)",
            "PanCK(Pt195)",
            "DNA1(Ir191)",
            "CD68(Pt196)",
            "CD45RA(Nd143)",
        ][:n_markers]
        csv_path = tiff_dir / f"{sample_id}_full.csv"
        csv_path.write_text(
            ",channel\n" + "\n".join(f"{i},{n}" for i, n in enumerate(channel_strings))
        )

        # Write synthetic TIFF
        try:
            import tifffile

            arr = (np.random.rand(n_markers, 32, 32) * 1000).astype(np.uint16)
            tifffile.imwrite(str(tiff_dir / f"{sample_id}_full.tiff"), arr)
        except ImportError:
            (tiff_dir / f"{sample_id}_full.tiff").write_bytes(b"dummy")

        return processed_dir

    def test_load_images_populates_uns_spatial(self, tmp_path):
        try:
            import tifffile  # noqa: F401
        except ImportError:
            pytest.skip("tifffile not installed")

        sample_id = "SAMPLE-01"
        processed_dir = self._make_imc_processed_dir(tmp_path, sample_id)
        result = load_imc_sample(processed_dir, sample_id, load_images=True)

        assert "spatial" in result.uns
        assert sample_id in result.uns["spatial"]
        spatial = result.uns["spatial"][sample_id]
        assert "images" in spatial
        assert "hires" in spatial["images"]
        assert "full" in spatial["images"]
        assert spatial["images"]["hires"].dtype == np.uint8
        assert spatial["images"]["hires"].shape[-1] == 3

    def test_load_images_false_no_spatial(self, tmp_path):
        sample_id = "SAMPLE-02"
        processed_dir = self._make_imc_processed_dir(tmp_path, sample_id)
        result = load_imc_sample(processed_dir, sample_id, load_images=False)
        assert not result.uns.get("spatial", {})

    def test_load_images_missing_tiffs_dir_warns(self, tmp_path, caplog):
        """Missing tiffs/ directory logs a warning and does not raise."""
        import logging

        n_obs, n_vars = 10, 3
        processed_dir = tmp_path / "processed" / "sample_no_tiffs"
        processed_dir.mkdir(parents=True)
        adata_inner = ad.AnnData(
            X=np.random.rand(n_obs, n_vars).astype(np.float32),
            obs=pd.DataFrame(index=[f"c_{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=[f"m_{i}" for i in range(n_vars)]),
        )
        adata_inner.write_h5ad(processed_dir / "cells.h5ad")

        with caplog.at_level(logging.WARNING):
            result = load_imc_sample(processed_dir, "NO_TIFFS", load_images=True)

        assert "tiffs/" in caplog.text
        assert result.n_obs == n_obs


# ---------------------------------------------------------------------------
# load_he_image tests
# ---------------------------------------------------------------------------


class TestLoadHeImage:
    def test_inject_into_existing_adata(self, tmp_path):
        try:
            import tifffile
        except ImportError:
            pytest.skip("tifffile not installed")

        img = (np.random.rand(64, 64, 3) * 255).astype(np.uint8)
        he_path = tmp_path / "sample.tiff"
        tifffile.imwrite(str(he_path), img)

        adata = ad.AnnData(
            X=np.zeros((5, 3)),
            obs=pd.DataFrame(index=[f"c_{i}" for i in range(5)]),
            var=pd.DataFrame(index=["A", "B", "C"]),
        )

        load_he_image(he_path, "LIB_01", adata)

        assert "spatial" in adata.uns
        assert "LIB_01" in adata.uns["spatial"]
        spatial = adata.uns["spatial"]["LIB_01"]
        assert "images" in spatial
        assert "hires" in spatial["images"]
        assert spatial["images"]["hires"].dtype == np.uint8
        assert spatial["images"]["hires"].shape == (64, 64, 3)

    def test_grayscale_converted_to_rgb(self, tmp_path):
        try:
            import tifffile
        except ImportError:
            pytest.skip("tifffile not installed")

        img = (np.random.rand(32, 32) * 255).astype(np.uint8)
        he_path = tmp_path / "gray.tiff"
        tifffile.imwrite(str(he_path), img)

        adata = ad.AnnData(X=np.zeros((2, 2)))
        load_he_image(he_path, "LIB_GRAY", adata)

        out = adata.uns["spatial"]["LIB_GRAY"]["images"]["hires"]
        assert out.shape == (32, 32, 3)

    def test_downsample(self, tmp_path):
        try:
            import tifffile
        except ImportError:
            pytest.skip("tifffile not installed")

        img = (np.random.rand(64, 64, 3) * 255).astype(np.uint8)
        he_path = tmp_path / "he.tiff"
        tifffile.imwrite(str(he_path), img)

        adata = ad.AnnData(X=np.zeros((2, 2)))
        load_he_image(he_path, "LIB_DS", adata, downsample=2)

        out = adata.uns["spatial"]["LIB_DS"]["images"]["hires"]
        assert out.shape == (32, 32, 3)

    def test_file_not_found_raises(self, tmp_path):
        adata = ad.AnnData(X=np.zeros((2, 2)))
        with pytest.raises(FileNotFoundError):
            load_he_image(tmp_path / "nonexistent.tiff", "LIB_X", adata)

    def test_custom_image_key(self, tmp_path):
        try:
            import tifffile
        except ImportError:
            pytest.skip("tifffile not installed")

        img = (np.random.rand(16, 16, 3) * 255).astype(np.uint8)
        he_path = tmp_path / "he_low.tiff"
        tifffile.imwrite(str(he_path), img)

        adata = ad.AnnData(X=np.zeros((2, 2)))
        load_he_image(he_path, "LIB_KEY", adata, image_key="lowres")

        assert "lowres" in adata.uns["spatial"]["LIB_KEY"]["images"]


# ---------------------------------------------------------------------------
# SLURM sbatch generation tests
# ---------------------------------------------------------------------------


class TestBuildSbatchHeader:
    def test_basic_header(self):
        header = build_sbatch_header("test_job")
        assert "#SBATCH --job-name=test_job" in header
        assert "#SBATCH --partition=scu-cpu" in header
        assert "#SBATCH --cpus-per-task=32" in header
        assert "#SBATCH --mem=240G" in header
        assert "#SBATCH --time=2-00:00:00" in header
        assert "#SBATCH --output=logs/test_job_%j.out" in header
        assert "#SBATCH --error=logs/test_job_%j.err" in header

    def test_custom_resources(self):
        header = build_sbatch_header(
            "myjob",
            partition="gpu",
            cpus_per_task=8,
            mem="64G",
            time="4:00:00",
        )
        assert "#SBATCH --partition=gpu" in header
        assert "#SBATCH --cpus-per-task=8" in header
        assert "#SBATCH --mem=64G" in header
        assert "#SBATCH --time=4:00:00" in header

    def test_custom_log_dir(self):
        header = build_sbatch_header("j1", log_dir="output/logs")
        assert "#SBATCH --output=output/logs/j1_%j.out" in header

    def test_extra_directives(self):
        header = build_sbatch_header("j1", extra_directives={"account": "mylab"})
        assert "#SBATCH --account=mylab" in header


class TestBuildSpacerangerSbatch:
    def test_basic_visium_hd(self):
        script = build_spaceranger_sbatch(
            sample_id="PT01-2_TUM",
            fastqs="/data/batch1",
            transcriptome="/ref/GRCh38",
            output_dir="/scratch/outputs",
            cytaimage="/data/batch1/PT01-2_TUM/cyto.tif",
            slide="H1-ABC",
            area="D1",
        )
        assert "#!/bin/bash" in script
        assert "#SBATCH --job-name=sr4_PT01-2_TUM" in script
        assert "set -euo pipefail" in script
        assert 'SR="spaceranger"' in script
        assert 'TRANSCRIPTOME="/ref/GRCh38"' in script
        assert 'CYTAIMAGE="/data/batch1/PT01-2_TUM/cyto.tif"' in script
        assert "--create-bam=false" in script
        assert "test -x" in script
        assert "test -d" in script
        assert "test -f" in script
        assert "[VERIFY]" in script

    def test_probe_set_included(self):
        script = build_spaceranger_sbatch(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            output_dir="/out",
            cytaimage="/img.tif",
            probe_set="/probes/v2.csv",
        )
        assert 'PROBE_SET="/probes/v2.csv"' in script
        assert '--probe-set="${PROBE_SET}"' in script
        assert "Probe set not found" in script

    def test_sample_filter_included(self):
        script = build_spaceranger_sbatch(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            output_dir="/out",
            cytaimage="/img.tif",
            sample_filter="PT01-1",
        )
        assert 'SAMPLE_FILTER="PT01-1"' in script
        assert '--sample="${SAMPLE_FILTER}"' in script

    def test_custom_spaceranger_path(self):
        script = build_spaceranger_sbatch(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            output_dir="/out",
            cytaimage="/img.tif",
            spaceranger_path="/opt/sr4/spaceranger",
        )
        assert 'SR="/opt/sr4/spaceranger"' in script

    def test_create_bam_true(self):
        script = build_spaceranger_sbatch(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            output_dir="/out",
            cytaimage="/img.tif",
            create_bam=True,
        )
        assert "--create-bam=true" in script

    def test_slurm_overrides(self):
        script = build_spaceranger_sbatch(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            output_dir="/out",
            cytaimage="/img.tif",
            slurm={"partition": "gpu", "mem": "128G"},
        )
        assert "#SBATCH --partition=gpu" in script
        assert "#SBATCH --mem=128G" in script

    def test_no_image_raises(self):
        with pytest.raises(ValueError, match="image"):
            build_spaceranger_sbatch(
                sample_id="S1",
                fastqs="/data",
                transcriptome="/ref",
                output_dir="/out",
            )

    def test_both_images(self):
        script = build_spaceranger_sbatch(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            output_dir="/out",
            cytaimage="/cyto.tif",
            image="/he.tif",
        )
        assert 'CYTAIMAGE="/cyto.tif"' in script
        assert 'IMAGE="/he.tif"' in script
        assert "CytAssist image not found" in script
        assert "H&E image not found" in script

    def test_cleanup_block(self):
        script = build_spaceranger_sbatch(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            output_dir="/out",
            cytaimage="/img.tif",
        )
        assert "[CLEANUP]" in script
        assert "rm -rf" in script

    def test_post_verification_block(self):
        script = build_spaceranger_sbatch(
            sample_id="S1",
            fastqs="/data",
            transcriptome="/ref",
            output_dir="/out",
            cytaimage="/img.tif",
        )
        assert "square_008um" in script
        assert "cell_segmentation" in script
        assert "cloupe" in script
        assert "[FAIL]" in script


class TestBuildXeniumSbatch:
    def test_basic(self):
        script = build_xenium_sbatch(
            sample_id="XEN_01",
            xenium_bundle="/data/xenium/run1",
            output_dir="/out",
        )
        assert "#!/bin/bash" in script
        assert "#SBATCH --job-name=xr_XEN_01" in script
        assert 'XENIUM_BUNDLE="/data/xenium/run1"' in script
        assert "xeniumranger" in script
        assert "--xenium-bundle" in script
        assert "test -d" in script

    def test_custom_binary(self):
        script = build_xenium_sbatch(
            sample_id="XEN_01",
            xenium_bundle="/data/xenium/run1",
            output_dir="/out",
            xenium_ranger_path="/opt/xr/xeniumranger",
        )
        assert 'XR="/opt/xr/xeniumranger"' in script


class TestBuildImcSbatch:
    def test_basic(self):
        script = build_imc_sbatch(
            sample_id="IMC_01",
            mcd_file="/data/sample.mcd",
            panel_csv="/panel.csv",
            output_dir="/out",
        )
        assert "#!/bin/bash" in script
        assert "#SBATCH --job-name=imc_IMC_01" in script
        assert 'MCD_FILE="/data/sample.mcd"' in script
        assert 'PANEL_CSV="/panel.csv"' in script
        assert "test -f" in script
        assert "tiffs/" in script
        assert "cells.h5ad" in script


class TestBuildBatchSbatch:
    def test_spaceranger_batch(self):
        manifest = pd.DataFrame(
            {
                "sample_id": ["S1", "S2"],
                "fastq_dir": ["/data/S1", "/data/S2"],
                "cytaimage": ["/img/S1.tif", "/img/S2.tif"],
                "slide": ["H1", "H1"],
                "area": ["D1", "A1"],
            }
        )
        scripts = build_batch_sbatch(
            manifest,
            modality="visium_hd",
            output_dir="/out",
            transcriptome="/ref",
        )
        assert len(scripts) == 2
        assert scripts[0][0] == "S1"
        assert scripts[1][0] == "S2"
        assert "sr4_S1" in scripts[0][1]
        assert "sr4_S2" in scripts[1][1]

    def test_xenium_batch(self):
        manifest = pd.DataFrame(
            {
                "sample_id": ["X1"],
                "xenium_dir": ["/data/xenium/X1"],
            }
        )
        scripts = build_batch_sbatch(
            manifest,
            modality="xenium",
            output_dir="/out",
        )
        assert len(scripts) == 1
        assert "xr_X1" in scripts[0][1]

    def test_imc_batch(self):
        manifest = pd.DataFrame(
            {
                "sample_id": ["I1"],
                "mcd_file": ["/data/I1.mcd"],
                "panel_csv": ["/panel.csv"],
            }
        )
        scripts = build_batch_sbatch(
            manifest,
            modality="imc",
            output_dir="/out",
        )
        assert len(scripts) == 1
        assert "imc_I1" in scripts[0][1]

    def test_skips_missing_transcriptome(self):
        manifest = pd.DataFrame(
            {
                "sample_id": ["S1"],
                "fastq_dir": ["/data/S1"],
                "cytaimage": ["/img/S1.tif"],
            }
        )
        scripts = build_batch_sbatch(
            manifest,
            modality="visium_hd",
            output_dir="/out",
            # No transcriptome
        )
        assert len(scripts) == 0

    def test_unsupported_modality(self):
        manifest = pd.DataFrame({"sample_id": ["S1"]})
        scripts = build_batch_sbatch(
            manifest,
            modality="cosmx",
            output_dir="/out",
        )
        assert len(scripts) == 0

    def test_probe_set_from_manifest(self):
        manifest = pd.DataFrame(
            {
                "sample_id": ["S1"],
                "fastq_dir": ["/data/S1"],
                "cytaimage": ["/img.tif"],
                "probe_set": ["/probes/custom.csv"],
            }
        )
        scripts = build_batch_sbatch(
            manifest,
            modality="visium_hd",
            output_dir="/out",
            transcriptome="/ref",
            probe_set="/probes/default.csv",
        )
        assert len(scripts) == 1
        # Manifest probe_set should override the function arg
        assert "/probes/custom.csv" in scripts[0][1]


class TestWriteSbatchScript:
    def test_creates_file(self, tmp_path):
        script = "#!/bin/bash\necho hello\n"
        path = write_sbatch_script(script, tmp_path / "test.sbatch")
        assert path.exists()
        assert path.read_text() == script
        # Check executable
        assert path.stat().st_mode & 0o111

    def test_creates_parent_dirs(self, tmp_path):
        path = write_sbatch_script("#!/bin/bash\n", tmp_path / "sub" / "dir" / "test.sbatch")
        assert path.exists()

    def test_no_overwrite_by_default(self, tmp_path):
        script_path = tmp_path / "test.sbatch"
        script_path.write_text("old")
        with pytest.raises(FileExistsError):
            write_sbatch_script("new", script_path)

    def test_overwrite_flag(self, tmp_path):
        script_path = tmp_path / "test.sbatch"
        script_path.write_text("old")
        write_sbatch_script("new", script_path, overwrite=True)
        assert script_path.read_text() == "new"
