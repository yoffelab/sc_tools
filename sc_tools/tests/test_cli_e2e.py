"""End-to-end CLI tests with real data on HPC (TST-06).

These tests run actual CLI commands on real project data. They are skipped
when real data is not available (CI, local dev). Run on HPC with:

    pytest sc_tools/tests/test_cli_e2e.py -x -q

Data paths are configured via environment variables:
    SCT_TEST_DATA_DIR: Directory with real h5ad files
    SCT_TEST_EMBEDDING_DIR: Directory with per-method h5ad embeddings
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest
from typer.testing import CliRunner

from sc_tools.cli import app

runner = CliRunner()

# --- Skip guards ---

_DATA_DIR = os.environ.get("SCT_TEST_DATA_DIR", "")
_EMBEDDING_DIR = os.environ.get("SCT_TEST_EMBEDDING_DIR", "")

_HAS_REAL_DATA = bool(_DATA_DIR) and Path(_DATA_DIR).is_dir()
_HAS_EMBEDDINGS = bool(_EMBEDDING_DIR) and Path(_EMBEDDING_DIR).is_dir()

skip_no_data = pytest.mark.skipif(not _HAS_REAL_DATA, reason="SCT_TEST_DATA_DIR not set or not a directory")
skip_no_embeddings = pytest.mark.skipif(not _HAS_EMBEDDINGS, reason="SCT_TEST_EMBEDDING_DIR not set or not a directory")


@skip_no_data
class TestQCRunE2E:
    """sct qc run on real data."""

    def test_qc_run_real_data(self):
        """Run QC on a real h5ad file."""
        h5ad_files = list(Path(_DATA_DIR).glob("*.h5ad"))
        assert h5ad_files, f"No h5ad files in {_DATA_DIR}"
        result = runner.invoke(app, ["qc", "run", str(h5ad_files[0])])
        assert result.exit_code == 0

    def test_validate_real_data(self):
        """Validate a real checkpoint."""
        h5ad_files = list(Path(_DATA_DIR).glob("*.h5ad"))
        assert h5ad_files, f"No h5ad files in {_DATA_DIR}"
        result = runner.invoke(app, ["validate", "run", "qc_filter", str(h5ad_files[0])])
        # May pass or fail (exit 0 or 2), but should not crash (exit 1 or 3)
        assert result.exit_code in (0, 2)


@skip_no_embeddings
class TestBenchmarkE2E:
    """sct benchmark integration on real embeddings."""

    def test_benchmark_real_embeddings(self):
        """Run benchmark on real per-method h5ad files."""
        result = runner.invoke(app, [
            "benchmark", "integration",
            "--from-dir", _EMBEDDING_DIR,
            "--subsample-n", "10000",
            "--seed", "42",
        ])
        assert result.exit_code == 0
