"""Tests for sct estimate CLI command (07-02 Task 1)."""

from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner


@pytest.fixture()
def tmp_h5ad(tmp_path):
    """Create a small sparse h5ad for estimate testing."""
    n_obs, n_vars = 50, 20
    rng = np.random.default_rng(42)
    dense = rng.standard_normal((n_obs, n_vars)).astype(np.float32)
    dense[dense < 0.5] = 0
    X = sp.csr_matrix(dense)
    obs_df = pd.DataFrame(
        {"batch": [f"b{i % 3}" for i in range(n_obs)]},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var_df = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=X, obs=obs_df, var=var_df)
    path = tmp_path / "test.h5ad"
    adata.write_h5ad(path)
    return path


@pytest.fixture()
def runner():
    return CliRunner()


def test_estimate_cli_returns_json(runner, tmp_h5ad):
    """Estimate command returns valid JSON with estimated_peak_mb key."""
    from sc_tools.cli import app

    result = runner.invoke(app, ["estimate", "preprocess", str(tmp_h5ad)])
    assert result.exit_code == 0, f"Exit code {result.exit_code}: {result.output}"
    data = json.loads(result.output)
    assert "estimated_peak_mb" in data.get("data", {})


def test_estimate_cli_has_n_obs(runner, tmp_h5ad):
    """Estimate output JSON contains n_obs matching test file cell count."""
    from sc_tools.cli import app

    result = runner.invoke(app, ["estimate", "preprocess", str(tmp_h5ad)])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["data"]["n_obs"] == 50


def test_estimate_cli_method_multiplier(runner, tmp_h5ad):
    """Preprocess estimate has higher peak than qc estimate (2.5 vs 1.5 multiplier)."""
    from sc_tools.cli import app

    res_preprocess = runner.invoke(app, ["estimate", "preprocess", str(tmp_h5ad)])
    res_qc = runner.invoke(app, ["estimate", "qc", str(tmp_h5ad)])
    assert res_preprocess.exit_code == 0
    assert res_qc.exit_code == 0

    peak_preprocess = json.loads(res_preprocess.output)["data"]["estimated_peak_mb"]
    peak_qc = json.loads(res_qc.output)["data"]["estimated_peak_mb"]
    assert peak_preprocess > peak_qc


def test_estimate_cli_missing_file(runner):
    """Estimate on nonexistent file exits with code 1 and error JSON."""
    from sc_tools.cli import app

    result = runner.invoke(app, ["estimate", "preprocess", "/nonexistent/file.h5ad"])
    assert result.exit_code == 1
    data = json.loads(result.output)
    assert data["status"] == "error"


def test_estimate_cli_unknown_command(runner, tmp_h5ad):
    """Estimate for unknown command name uses default multiplier, still returns valid JSON."""
    from sc_tools.cli import app

    result = runner.invoke(app, ["estimate", "unknown_method", str(tmp_h5ad)])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert "estimated_peak_mb" in data.get("data", {})
    # Should use default multiplier (2.0)
    assert data["data"]["method_multiplier"] == 2.0
