"""Tests for --dry-run and --force flags on data-touching CLI commands (07-02 Task 2)."""

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
    """Create a small sparse h5ad for dry-run testing."""
    n_obs, n_vars = 50, 20
    rng = np.random.default_rng(42)
    dense = rng.standard_normal((n_obs, n_vars)).astype(np.float32)
    dense[dense < 0.5] = 0
    X = sp.csr_matrix(dense)
    obs_df = pd.DataFrame(
        {
            "batch": [f"b{i % 3}" for i in range(n_obs)],
            "celltype": [f"ct{i % 5}" for i in range(n_obs)],
            "library_id": [f"lib{i % 2}" for i in range(n_obs)],
        },
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


def test_dryrun_qc_run(runner, tmp_h5ad):
    """--dry-run on qc run exits 0 with status=skipped."""
    from sc_tools.cli import app

    result = runner.invoke(app, ["qc", "run", str(tmp_h5ad), "--dry-run"])
    assert result.exit_code == 0, f"Exit code {result.exit_code}: {result.output}"
    data = json.loads(result.output)
    assert data["status"] == "skipped"


def test_dryrun_preprocess_run(runner, tmp_h5ad):
    """--dry-run on preprocess run exits 0 with status=skipped and planned_command."""
    from sc_tools.cli import app

    result = runner.invoke(app, ["preprocess", "run", str(tmp_h5ad), "--dry-run"])
    assert result.exit_code == 0, f"Exit code {result.exit_code}: {result.output}"
    data = json.loads(result.output)
    assert data["status"] == "skipped"
    assert "planned_command" in data.get("data", {})


def test_dryrun_no_side_effects(runner, tmp_h5ad, tmp_path):
    """--dry-run on preprocess run does NOT create an output file."""
    from sc_tools.cli import app

    out_path = tmp_path / "output.h5ad"
    result = runner.invoke(
        app,
        ["preprocess", "run", str(tmp_h5ad), "--output", str(out_path), "--dry-run"],
    )
    assert result.exit_code == 0
    assert not out_path.exists(), "Output file should not exist after --dry-run"


def test_dryrun_missing_file(runner):
    """--dry-run on qc run with missing file exits 1 with error JSON."""
    from sc_tools.cli import app

    result = runner.invoke(app, ["qc", "run", "/nonexistent.h5ad", "--dry-run"])
    assert result.exit_code == 1
    data = json.loads(result.output)
    assert data["status"] == "error"


def test_dryrun_reports_estimate(runner, tmp_h5ad):
    """--dry-run on qc run includes estimate data with n_obs."""
    from sc_tools.cli import app

    result = runner.invoke(app, ["qc", "run", str(tmp_h5ad), "--dry-run"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    estimate = data["data"].get("estimate", {})
    assert estimate is not None
    assert estimate["n_obs"] == 50


def test_force_flag_accepted(runner, tmp_h5ad):
    """--force flag is accepted by qc run (does not cause flag parsing error)."""
    from sc_tools.cli import app

    # This will likely fail on actual QC execution, but the flag itself should be accepted
    # The --force flag is popped by cli_handler before the function runs
    result = runner.invoke(app, ["qc", "run", str(tmp_h5ad), "--force"])
    # Exit code may be non-zero due to QC deps, but should NOT be 2 (usage error)
    assert result.exit_code != 2, f"Flag parsing error: {result.output}"
