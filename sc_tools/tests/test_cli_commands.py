"""CLI command integration tests (TST-05).

Uses 100-cell AnnData fixtures from conftest.py to test actual command execution
through Typer's CliRunner. Tests verify CLIResult JSON output structure and
correctness without requiring real project data.
"""
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from sc_tools.cli import app

runner = CliRunner()


def _parse_result(output: str) -> dict:
    """Parse CLIResult JSON from command stdout."""
    return json.loads(output)


# ---------------------------------------------------------------------------
# sct validate (CMD-03)
# ---------------------------------------------------------------------------


class TestValidate:
    """sct validate run <phase> <file>."""

    def test_validate_qc_filter_passes(self, adata_100_h5ad):
        """adata_100 has sample, raw_data_dir, spatial -> qc_filter passes."""
        result = runner.invoke(app, ["validate", "run", "qc_filter", adata_100_h5ad])
        assert result.exit_code == 0
        data = _parse_result(result.output)
        assert data["status"] == "success"
        assert data["data"]["n_issues"] == 0

    def test_validate_preprocess_passes(self, adata_100_preprocess_checkpoint):
        """adata_100_preprocess_checkpoint has X_scvi, leiden, raw -> preprocess passes."""
        result = runner.invoke(app, ["validate", "run", "preprocess", adata_100_preprocess_checkpoint])
        assert result.exit_code == 0
        data = _parse_result(result.output)
        assert data["status"] == "success"

    def test_validate_missing_file(self, tmp_path):
        """Missing file raises user error (exit 1)."""
        result = runner.invoke(app, ["validate", "run", "qc_filter", str(tmp_path / "nope.h5ad")])
        assert result.exit_code == 1

    def test_validate_bad_phase(self, adata_100_h5ad):
        """Unknown phase slug raises user error (exit 1)."""
        result = runner.invoke(app, ["validate", "run", "invalid_phase", adata_100_h5ad])
        assert result.exit_code == 1

    def test_validate_json_structure(self, adata_100_h5ad):
        """Output contains all CLIResult fields."""
        result = runner.invoke(app, ["validate", "run", "qc_filter", adata_100_h5ad])
        data = _parse_result(result.output)
        for field in ("status", "command", "data", "provenance", "message"):
            assert field in data, f"Missing field: {field}"
        assert "phase" in data["data"]
        assert "issues" in data["data"]


# ---------------------------------------------------------------------------
# sct status (CMD-05)
# ---------------------------------------------------------------------------


class TestStatus:
    """sct status show."""

    def test_status_returns_dag(self):
        """Status shows pipeline phases even without registry (D-09)."""
        result = runner.invoke(app, ["status", "show"])
        assert result.exit_code == 0
        data = _parse_result(result.output)
        assert data["status"] == "success"
        assert "phases" in data["data"]
        assert len(data["data"]["phases"]) > 0

    def test_status_has_available_next(self):
        """Status lists available_next phases."""
        result = runner.invoke(app, ["status", "show"])
        data = _parse_result(result.output)
        assert "available_next" in data["data"]
        # With no completed phases, ingest_raw should be available
        assert len(data["data"]["available_next"]) > 0

    def test_status_registry_unavailable_message(self):
        """When registry unavailable, message indicates so (D-09)."""
        result = runner.invoke(app, ["status", "show"])
        data = _parse_result(result.output)
        if not data["data"]["registry_available"]:
            assert "registry unavailable" in data["message"].lower()


# ---------------------------------------------------------------------------
# sct qc run (CMD-01)
# ---------------------------------------------------------------------------


class TestQCRun:
    """sct qc run <file>."""

    def test_qc_run_produces_metrics(self, adata_100_h5ad):
        """QC run outputs JSON with per-sample metrics."""
        result = runner.invoke(app, ["qc", "run", adata_100_h5ad])
        assert result.exit_code == 0
        data = _parse_result(result.output)
        assert data["status"] == "success"
        assert data["data"]["n_cells"] == 100
        assert data["data"]["n_samples"] == 2  # L0, L1
        assert "metrics" in data["data"]

    def test_qc_run_missing_file(self, tmp_path):
        """Missing file gives exit 1."""
        result = runner.invoke(app, ["qc", "run", str(tmp_path / "nope.h5ad")])
        assert result.exit_code == 1

    def test_qc_run_json_structure(self, adata_100_h5ad):
        """Output has all CLIResult fields."""
        result = runner.invoke(app, ["qc", "run", adata_100_h5ad])
        data = _parse_result(result.output)
        for field in ("status", "command", "data", "provenance", "message"):
            assert field in data
        assert data["command"] == "qc run"


# ---------------------------------------------------------------------------
# sct preprocess run (CMD-02)
# ---------------------------------------------------------------------------


class TestPreprocessRun:
    """sct preprocess run <file> (CMD-02)."""

    def test_preprocess_help_shows_options(self):
        """Preprocess --help shows run subcommand with expected options."""
        result = runner.invoke(app, ["preprocess", "run", "--help"])
        assert result.exit_code == 0
        assert "--modality" in result.output
        assert "--integration" in result.output
        assert "--project-dir" in result.output

    def test_preprocess_missing_file(self, tmp_path):
        """Missing file gives exit 1."""
        result = runner.invoke(app, ["preprocess", "run", str(tmp_path / "nope.h5ad")])
        assert result.exit_code == 1

    def test_preprocess_bad_integration(self, adata_100_h5ad):
        """Unknown integration method gives exit 1."""
        result = runner.invoke(app, [
            "preprocess", "run", adata_100_h5ad,
            "--integration", "nonexistent_method",
        ])
        assert result.exit_code == 1

    def test_detect_modality_from_uns(self, adata_100):
        """_detect_modality reads adata.uns['modality'] (D-08)."""
        from sc_tools.cli.preprocess import _detect_modality
        assert _detect_modality(adata_100) == "visium"

    def test_detect_modality_raises_on_unknown(self, adata_100):
        """_detect_modality raises SCToolsDataError when cannot determine (D-08)."""
        from sc_tools.cli.preprocess import _detect_modality
        from sc_tools.errors import SCToolsDataError
        # Remove modality hint and use large panel (n_vars=200 > 1000 check won't trigger xenium)
        # But 200 < 1000 so it will detect as xenium. We need n_vars >= 1000.
        # Use a modified adata to test the error path.
        import numpy as np
        from anndata import AnnData

        # Create adata with >1000 genes and no modality hint
        rng = np.random.default_rng(99)
        big_adata = AnnData(rng.negative_binomial(5, 0.3, (10, 1500)).astype("float32"))
        with pytest.raises(SCToolsDataError, match="Cannot auto-detect"):
            _detect_modality(big_adata)


# ---------------------------------------------------------------------------
# sct benchmark integration (CMD-04)
# ---------------------------------------------------------------------------


class TestBenchmarkIntegration:
    """sct benchmark integration --from-dir <dir> (CMD-04)."""

    def test_benchmark_help_shows_options(self):
        """Benchmark --help shows integration subcommand with expected options."""
        result = runner.invoke(app, ["benchmark", "integration", "--help"])
        assert result.exit_code == 0
        assert "--from-dir" in result.output
        assert "--subsample-n" in result.output
        assert "--seed" in result.output
        assert "--report" in result.output
        assert "--batch-key" in result.output

    def test_benchmark_missing_dir(self, tmp_path):
        """Missing --from-dir gives exit 1."""
        result = runner.invoke(app, [
            "benchmark", "integration",
            "--from-dir", str(tmp_path / "nonexistent"),
        ])
        assert result.exit_code == 1

    def test_benchmark_empty_dir(self, tmp_path):
        """Directory with no h5ad files gives exit 1."""
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        result = runner.invoke(app, [
            "benchmark", "integration",
            "--from-dir", str(empty_dir),
        ])
        assert result.exit_code == 1


# ---------------------------------------------------------------------------
# CMD-08: fast-fail dependency check
# ---------------------------------------------------------------------------


class TestDepCheck:
    """_check_deps fast-fail (CMD-08)."""

    def test_check_deps_missing_package(self):
        """Missing dep raises SCToolsUserError with install instructions."""
        from sc_tools.cli import _check_deps
        from sc_tools.errors import SCToolsUserError

        with pytest.raises(SCToolsUserError, match="Missing required dependencies"):
            _check_deps(["nonexistent_package_xyz"])

    def test_check_deps_available_package(self):
        """Available dep does not raise."""
        from sc_tools.cli import _check_deps
        _check_deps(["json"])  # stdlib, always available


# ---------------------------------------------------------------------------
# CMD-07: Shared Result type (CLI + MCP)
# ---------------------------------------------------------------------------


class TestSharedResult:
    """CLI and MCP tools share CLIResult type (CMD-07)."""

    def test_mcp_validate_returns_cli_result(self, adata_100_h5ad):
        """MCP validate_checkpoint returns valid CLIResult JSON."""
        from sc_tools.mcp.tools_server import validate_checkpoint

        json_str = validate_checkpoint(adata_100_h5ad, "qc_filter")
        data = json.loads(json_str)
        for field in ("status", "command", "data", "provenance", "message"):
            assert field in data, f"MCP result missing field: {field}"

    def test_cli_and_mcp_same_result_type(self, adata_100_h5ad):
        """Both CLI and MCP produce CLIResult with same fields."""
        # CLI
        cli_result = runner.invoke(app, ["validate", "run", "qc_filter", adata_100_h5ad])
        cli_data = _parse_result(cli_result.output)

        # MCP
        from sc_tools.mcp.tools_server import validate_checkpoint
        mcp_data = json.loads(validate_checkpoint(adata_100_h5ad, "qc_filter"))

        # Both should have the same top-level fields
        cli_fields = set(cli_data.keys())
        mcp_fields = set(mcp_data.keys())
        assert cli_fields == mcp_fields, f"Field mismatch: CLI={cli_fields - mcp_fields}, MCP={mcp_fields - cli_fields}"
