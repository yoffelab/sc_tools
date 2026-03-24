"""Tests for the sct assemble CLI command group."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

mudata = pytest.importorskip("mudata")

import anndata as ad
from typer.testing import CliRunner

from sc_tools.cli import app

runner = CliRunner()


def _make_h5ad(path: Path, n_obs: int, n_var: int, subject_ids: list[str], prefix: str) -> None:
    """Write a minimal h5ad file for testing."""
    rng = np.random.default_rng(42)
    X = rng.random((n_obs, n_var)).astype(np.float32)
    cells_per_patient = n_obs // len(subject_ids)
    sids = []
    celltypes = []
    for pat in subject_ids:
        sids.extend([pat] * cells_per_patient)
        celltypes.extend([f"type_{j % 3}" for j in range(cells_per_patient)])
    obs = pd.DataFrame(
        {
            "subject_id": pd.Categorical(sids),
            "celltype": pd.Categorical(celltypes),
        },
        index=[f"{prefix}_cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"{prefix}_gene_{i}" for i in range(n_var)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.write_h5ad(str(path))


class TestAssembleBuild:
    """Test sct assemble build command."""

    def test_build_command(self, tmp_path):
        """Build h5mu from 2 h5ad files succeeds."""
        rna_path = tmp_path / "rna.h5ad"
        imc_path = tmp_path / "imc.h5ad"
        out_path = tmp_path / "atlas.h5mu"

        _make_h5ad(rna_path, 30, 100, ["PAT1", "PAT2", "PAT3"], "rna")
        _make_h5ad(imc_path, 20, 40, ["PAT1", "PAT2"], "imc")

        result = runner.invoke(
            app,
            [
                "assemble", "build",
                str(rna_path), str(imc_path),
                "-m", "rna", "-m", "imc",
                "-o", str(out_path),
            ],
        )

        assert result.exit_code == 0, f"CLI failed: {result.output}"
        assert out_path.exists()

        output = json.loads(result.output)
        assert output["status"] == "success"
        assert "rna" in output["data"]["modalities"]
        assert "imc" in output["data"]["modalities"]

    def test_build_mismatched_args(self, tmp_path):
        """Mismatched inputs/modalities count gives exit code 1."""
        rna_path = tmp_path / "rna.h5ad"
        _make_h5ad(rna_path, 30, 100, ["PAT1", "PAT2", "PAT3"], "rna")

        result = runner.invoke(
            app,
            [
                "assemble", "build",
                str(rna_path),
                "-m", "rna", "-m", "imc",  # 2 modalities for 1 input
                "-o", str(tmp_path / "out.h5mu"),
            ],
        )

        assert result.exit_code == 1
        output = json.loads(result.output)
        assert output["status"] == "error"

    def test_build_missing_file(self, tmp_path):
        """Non-existent input file gives exit code 1."""
        result = runner.invoke(
            app,
            [
                "assemble", "build",
                str(tmp_path / "nonexistent.h5ad"),
                "-m", "rna",
                "-o", str(tmp_path / "out.h5mu"),
            ],
        )

        assert result.exit_code == 1
        output = json.loads(result.output)
        assert output["status"] == "error"


class TestAssembleQuery:
    """Test sct assemble query command."""

    def test_query_command(self, tmp_path):
        """Query celltype_proportions from built h5mu returns valid JSON."""
        rna_path = tmp_path / "rna.h5ad"
        imc_path = tmp_path / "imc.h5ad"
        out_path = tmp_path / "atlas.h5mu"

        _make_h5ad(rna_path, 30, 100, ["PAT1", "PAT2", "PAT3"], "rna")
        _make_h5ad(imc_path, 20, 40, ["PAT1", "PAT2"], "imc")

        # Build first
        runner.invoke(
            app,
            [
                "assemble", "build",
                str(rna_path), str(imc_path),
                "-m", "rna", "-m", "imc",
                "-o", str(out_path),
            ],
        )

        # Now query
        result = runner.invoke(
            app,
            ["assemble", "query", str(out_path)],
        )

        assert result.exit_code == 0, f"CLI failed: {result.output}"
        output = json.loads(result.output)
        assert output["status"] == "success"
        assert "proportions" in output["data"]
        assert len(output["data"]["proportions"]) > 0


class TestAssembleEmbed:
    """Test sct assemble embed command."""

    def test_embed_command(self, tmp_path):
        """Embed with mofa on built h5mu succeeds."""
        pytest.importorskip("muon")
        pytest.importorskip("mofapy2")

        rna_path = tmp_path / "rna.h5ad"
        imc_path = tmp_path / "imc.h5ad"
        out_path = tmp_path / "atlas.h5mu"

        _make_h5ad(rna_path, 30, 100, ["PAT1", "PAT2", "PAT3"], "rna")
        _make_h5ad(imc_path, 20, 40, ["PAT1", "PAT2"], "imc")

        # Build first
        runner.invoke(
            app,
            [
                "assemble", "build",
                str(rna_path), str(imc_path),
                "-m", "rna", "-m", "imc",
                "-o", str(out_path),
            ],
        )

        # Now embed
        result = runner.invoke(
            app,
            [
                "assemble", "embed", str(out_path),
                "--method", "mofa",
                "--n-factors", "5",
            ],
        )

        assert result.exit_code == 0, f"CLI failed: {result.output}"
        output = json.loads(result.output)
        assert output["status"] == "success"
        assert output["data"]["method"] == "mofa"
        assert output["data"]["n_factors"] == 5


class TestAssembleHelp:
    """Test sct assemble --help."""

    def test_assemble_help(self):
        """sct assemble --help lists build, embed, query subcommands."""
        result = runner.invoke(app, ["assemble", "--help"])

        assert result.exit_code == 0
        assert "build" in result.output
        assert "embed" in result.output
        assert "query" in result.output
