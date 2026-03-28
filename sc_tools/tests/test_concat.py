"""Unit tests for sct concat command (CONCAT-01 through CONCAT-04).

Tests cover:
- CLI command invocation and validation
- Spatial metadata preservation during concatenation
- Provenance sidecar creation with SHA256 checksums
- Pipeline DAG registration of concat phase
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def spatial_adata_pair(tmp_path):
    """Two minimal AnnDatas with spatial metadata for concat testing."""
    adatas = []
    for _i, name in enumerate(["sample_A", "sample_B"]):
        n_obs, n_var = 50, 20
        rng = np.random.default_rng(42 + _i)
        X = rng.random((n_obs, n_var)).astype(np.float32)
        obs = pd.DataFrame(
            {"sample": name},
            index=[f"{name}_cell_{j}" for j in range(n_obs)],
        )
        var = pd.DataFrame(index=[f"gene_{j}" for j in range(n_var)])
        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.obsm["spatial"] = rng.random((n_obs, 2)).astype(np.float64)
        adata.uns["spatial"] = {name: {"scalefactors": {"spot_diameter": 1.0}}}
        path = tmp_path / f"{name}.h5ad"
        adata.write_h5ad(path)
        adatas.append((str(path), name))
    return adatas, tmp_path


# ---------------------------------------------------------------------------
# TestConcatCommand (CONCAT-01)
# ---------------------------------------------------------------------------


class TestConcatCommand:
    """Tests for the sct concat CLI command."""

    def test_concat_two_samples(self, spatial_adata_pair):
        """Two minimal AnnDatas merge into output with correct n_obs."""
        from sc_tools.cli.concat import register_concat

        import typer
        from typer.testing import CliRunner

        pairs, tmp = spatial_adata_pair
        output = str(tmp / "merged.h5ad")

        test_app = typer.Typer()
        register_concat(test_app)
        runner = CliRunner()
        result = runner.invoke(
            test_app,
            ["concat", "--input", pairs[0][0], "--input", pairs[1][0], "--output", output],
        )
        assert result.exit_code == 0, f"CLI failed: {result.output}"
        assert Path(output).exists()
        merged = ad.read_h5ad(output)
        assert merged.n_obs == 100  # 50 + 50

    def test_concat_dry_run(self, spatial_adata_pair):
        """With --dry-run, no output file is created."""
        from sc_tools.cli.concat import register_concat

        import typer
        from typer.testing import CliRunner

        pairs, tmp = spatial_adata_pair
        output = str(tmp / "dry_output.h5ad")

        test_app = typer.Typer()
        register_concat(test_app)
        runner = CliRunner()
        result = runner.invoke(
            test_app,
            [
                "concat",
                "--input", pairs[0][0],
                "--input", pairs[1][0],
                "--output", output,
                "--dry-run",
            ],
        )
        assert result.exit_code == 0
        assert not Path(output).exists()

    def test_concat_single_sample_error(self, spatial_adata_pair):
        """Single input returns an error exit code."""
        from sc_tools.cli.concat import register_concat

        import typer
        from typer.testing import CliRunner

        pairs, tmp = spatial_adata_pair
        output = str(tmp / "single.h5ad")

        test_app = typer.Typer()
        register_concat(test_app)
        runner = CliRunner()
        result = runner.invoke(
            test_app,
            ["concat", "--input", pairs[0][0], "--output", output],
        )
        # Should fail with exit code 1 (user error)
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# TestSpatialPreservation (CONCAT-02)
# ---------------------------------------------------------------------------


class TestSpatialPreservation:
    """Tests for spatial metadata preservation during concatenation."""

    def test_spatial_keys_preserved(self, spatial_adata_pair):
        """merged.uns['spatial'] contains all input library_ids."""
        from sc_tools.cli.concat import register_concat

        import typer
        from typer.testing import CliRunner

        pairs, tmp = spatial_adata_pair
        output = str(tmp / "merged_spatial.h5ad")

        test_app = typer.Typer()
        register_concat(test_app)
        runner = CliRunner()
        result = runner.invoke(
            test_app,
            ["concat", "--input", pairs[0][0], "--input", pairs[1][0], "--output", output],
        )
        assert result.exit_code == 0, f"CLI failed: {result.output}"
        merged = ad.read_h5ad(output)
        spatial_keys = set(merged.uns.get("spatial", {}).keys())
        assert "sample_A" in spatial_keys
        assert "sample_B" in spatial_keys

    def test_spatial_coords_preserved(self, spatial_adata_pair):
        """merged.obsm['spatial'] has correct shape after concat."""
        from sc_tools.cli.concat import register_concat

        import typer
        from typer.testing import CliRunner

        pairs, tmp = spatial_adata_pair
        output = str(tmp / "merged_coords.h5ad")

        test_app = typer.Typer()
        register_concat(test_app)
        runner = CliRunner()
        result = runner.invoke(
            test_app,
            ["concat", "--input", pairs[0][0], "--input", pairs[1][0], "--output", output],
        )
        assert result.exit_code == 0, f"CLI failed: {result.output}"
        merged = ad.read_h5ad(output)
        assert "spatial" in merged.obsm
        assert merged.obsm["spatial"].shape == (100, 2)


# ---------------------------------------------------------------------------
# TestConcatProvenance (CONCAT-03)
# ---------------------------------------------------------------------------


class TestConcatProvenance:
    """Tests for provenance sidecar creation."""

    def test_provenance_sidecar_created(self, spatial_adata_pair):
        """.provenance.json exists alongside output."""
        from sc_tools.cli.concat import register_concat

        import typer
        from typer.testing import CliRunner

        pairs, tmp = spatial_adata_pair
        output = str(tmp / "merged_prov.h5ad")

        test_app = typer.Typer()
        register_concat(test_app)
        runner = CliRunner()
        result = runner.invoke(
            test_app,
            ["concat", "--input", pairs[0][0], "--input", pairs[1][0], "--output", output],
        )
        assert result.exit_code == 0, f"CLI failed: {result.output}"
        prov_path = Path(output).with_suffix(".provenance.json")
        assert prov_path.exists(), f"Provenance sidecar not found at {prov_path}"

    def test_provenance_contains_sha256(self, spatial_adata_pair):
        """Sidecar has sha256 for each input."""
        from sc_tools.cli.concat import register_concat

        import typer
        from typer.testing import CliRunner

        pairs, tmp = spatial_adata_pair
        output = str(tmp / "merged_sha.h5ad")

        test_app = typer.Typer()
        register_concat(test_app)
        runner = CliRunner()
        result = runner.invoke(
            test_app,
            ["concat", "--input", pairs[0][0], "--input", pairs[1][0], "--output", output],
        )
        assert result.exit_code == 0, f"CLI failed: {result.output}"
        prov_path = Path(output).with_suffix(".provenance.json")
        prov_data = json.loads(prov_path.read_text())
        # Should have input files with sha256
        inputs = prov_data.get("inputs", [])
        assert len(inputs) >= 2
        for inp in inputs:
            assert "sha256" in inp, f"Input record missing sha256: {inp}"


# ---------------------------------------------------------------------------
# TestConcatPipeline (CONCAT-04)
# ---------------------------------------------------------------------------


class TestConcatPipeline:
    """Tests for concat pipeline phase registration."""

    def test_concat_in_standard_phases(self):
        """'concat' key exists in STANDARD_PHASES."""
        from sc_tools.pipeline import STANDARD_PHASES

        assert "concat" in STANDARD_PHASES

    def test_concat_depends_on_ingest_load(self):
        """concat.depends_on includes ingest_load."""
        from sc_tools.pipeline import STANDARD_PHASES, _dp

        spec = STANDARD_PHASES["concat"]
        assert _dp("ingest_load") in spec.depends_on

    def test_concat_visible_in_status(self):
        """concat appears between ingest_load and qc_filter in DAG."""
        from sc_tools.pipeline import get_available_next

        available = get_available_next(["ingest_raw", "ingest_load"])
        slugs = [k[1] for k in available]
        assert "concat" in slugs
        assert "qc_filter" in slugs
