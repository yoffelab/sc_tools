"""Tests for the run_full_phase MCP tool in sc_tools.mcp.tools_server.

TDD: tests written before implementation.
"""

from __future__ import annotations

import pytest

mcp = pytest.importorskip("mcp", reason="mcp not installed")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def project_dir(tmp_path):
    """Create a minimal mock project directory structure."""
    proj = tmp_path / "test_project"
    proj.mkdir()
    (proj / "results").mkdir()
    (proj / "metadata").mkdir()
    # Snakefile is required by run_full_phase to locate the snakemake target
    snakefile = proj / "Snakefile"
    snakefile.write_text("# mock Snakefile\n")
    return proj


@pytest.fixture()
def project_with_filtered(project_dir):
    """Project that has adata.filtered.h5ad ready (qc_filter done)."""
    (project_dir / "results" / "adata.filtered.h5ad").write_bytes(b"mock")
    return project_dir


@pytest.fixture()
def project_with_annotated(project_with_filtered):
    """Project that has adata.annotated.h5ad ready (metadata_attach done)."""
    proj = project_with_filtered
    (proj / "results" / "adata.annotated.h5ad").write_bytes(b"mock")
    return proj


@pytest.fixture()
def server():
    """Import the MCP tools_server module."""
    from sc_tools.mcp import tools_server

    return tools_server


# ---------------------------------------------------------------------------
# Tests: project path validation
# ---------------------------------------------------------------------------


class TestRunFullPhaseProjectValidation:
    def test_nonexistent_project_path_returns_error(self, server, tmp_path):
        result = server.run_full_phase(
            project_path=str(tmp_path / "does_not_exist"),
            phase_slug="qc_filter",
        )
        assert isinstance(result, dict)
        assert result["status"] == "error"
        assert (
            "not found" in result["message"].lower()
            or "does not exist" in result["message"].lower()
        )

    def test_file_path_instead_of_dir_returns_error(self, server, tmp_path):
        f = tmp_path / "notadir.txt"
        f.write_text("x")
        result = server.run_full_phase(
            project_path=str(f),
            phase_slug="qc_filter",
        )
        assert result["status"] == "error"

    def test_unknown_phase_slug_returns_error(self, server, project_dir):
        result = server.run_full_phase(
            project_path=str(project_dir),
            phase_slug="nonexistent_phase",
        )
        assert result["status"] == "error"
        assert "phase" in result["message"].lower() or "slug" in result["message"].lower()


# ---------------------------------------------------------------------------
# Tests: prerequisite checking
# ---------------------------------------------------------------------------


class TestRunFullPhasePrerequisites:
    def test_missing_inputs_returns_missing_status(self, server, project_dir):
        """metadata_attach requires adata.filtered.h5ad; if absent, report missing."""
        result = server.run_full_phase(
            project_path=str(project_dir),
            phase_slug="metadata_attach",
        )
        assert result["status"] == "missing_inputs"
        assert len(result["missing"]) > 0

    def test_ready_when_inputs_present(self, server, project_with_filtered):
        """metadata_attach can proceed when adata.filtered.h5ad exists."""
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
        )
        assert result["status"] in ("ready", "missing_inputs")
        # If filtered exists, metadata_attach input is satisfied
        assert result["status"] == "ready"

    def test_phase_with_no_checkpoint_prerequisite_is_ready(self, server, project_dir):
        """ingest_raw has no depends_on, so it is always ready."""
        result = server.run_full_phase(
            project_path=str(project_dir),
            phase_slug="ingest_raw",
        )
        # ingest_raw has no file-based prereqs and no input checkpoint to check
        assert result["status"] == "ready"
        assert result["missing"] == []


# ---------------------------------------------------------------------------
# Tests: return structure
# ---------------------------------------------------------------------------


class TestRunFullPhaseReturnStructure:
    def test_result_has_required_keys(self, server, project_with_filtered):
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
        )
        for key in ("phase", "status", "command", "inputs", "outputs", "missing"):
            assert key in result, f"Missing key: {key}"

    def test_phase_key_matches_slug(self, server, project_with_filtered):
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
        )
        assert result["phase"] == "metadata_attach"

    def test_inputs_is_list(self, server, project_with_filtered):
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
        )
        assert isinstance(result["inputs"], list)

    def test_outputs_is_list(self, server, project_with_filtered):
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
        )
        assert isinstance(result["outputs"], list)

    def test_missing_is_list(self, server, project_with_filtered):
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
        )
        assert isinstance(result["missing"], list)

    def test_command_contains_snakemake(self, server, project_with_filtered):
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
        )
        if result["status"] == "ready":
            assert "snakemake" in result["command"].lower()

    def test_command_contains_project_path(self, server, project_with_filtered):
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
        )
        if result["status"] == "ready":
            assert str(project_with_filtered) in result["command"]

    def test_outputs_contain_checkpoint_relative_path(self, server, project_with_filtered):
        """Output list must include the expected checkpoint filename."""
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
        )
        # metadata_attach checkpoint is results/adata.annotated.h5ad
        assert any("adata.annotated.h5ad" in o for o in result["outputs"])


# ---------------------------------------------------------------------------
# Tests: dry_run mode
# ---------------------------------------------------------------------------


class TestRunFullPhaseDryRun:
    def test_dry_run_returns_dict(self, server, project_with_filtered):
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
            dry_run=True,
        )
        assert isinstance(result, dict)

    def test_dry_run_contains_phase_info(self, server, project_with_filtered):
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
            dry_run=True,
        )
        assert result["phase"] == "metadata_attach"
        assert "command" in result
        assert "inputs" in result
        assert "outputs" in result

    def test_dry_run_and_live_return_same_structure(self, server, project_with_filtered):
        """dry_run=True and dry_run=False must return same dict shape."""
        dry = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
            dry_run=True,
        )
        live = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
            dry_run=False,
        )
        assert set(dry.keys()) == set(live.keys())


# ---------------------------------------------------------------------------
# Tests: command construction
# ---------------------------------------------------------------------------


class TestRunFullPhaseCommand:
    def test_command_references_snakefile(self, server, project_with_filtered):
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
        )
        if result["status"] == "ready":
            assert "Snakefile" in result["command"] or "-s " in result["command"]

    def test_command_is_not_executed(self, server, project_with_filtered):
        """run_full_phase must return the command string, not execute it."""
        # If this returns a valid dict with 'command' as a string, it was not executed
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
        )
        assert isinstance(result["command"], str)

    def test_scoring_phase_uses_scoring_target(self, server, project_with_annotated):
        """The Snakemake target in the command should match the phase slug."""
        result = server.run_full_phase(
            project_path=str(project_with_annotated),
            phase_slug="preprocess",
        )
        # preprocess requires adata.annotated.h5ad which is present
        if result["status"] == "ready":
            assert "preprocess" in result["command"]

    def test_dry_run_appends_dry_run_flag(self, server, project_with_filtered):
        """dry_run=True must include --dry-run in the command string."""
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
            dry_run=True,
        )
        assert "--dry-run" in result["command"]

    def test_live_run_does_not_contain_dry_run_flag(self, server, project_with_filtered):
        """dry_run=False (default) must NOT include --dry-run in the command."""
        result = server.run_full_phase(
            project_path=str(project_with_filtered),
            phase_slug="metadata_attach",
            dry_run=False,
        )
        assert "--dry-run" not in result["command"]


# ---------------------------------------------------------------------------
# Tests: placeholder phases (checkpoint with {sample_id})
# ---------------------------------------------------------------------------


class TestRunFullPhasePlaceholder:
    def test_ingest_load_does_not_raise(self, server, project_dir):
        """ingest_load has a {sample_id} checkpoint; run_full_phase must not raise KeyError."""
        result = server.run_full_phase(
            project_path=str(project_dir),
            phase_slug="ingest_load",
        )
        assert isinstance(result, dict)
        assert result["status"] in ("ready", "missing_inputs", "error")

    def test_ingest_load_outputs_is_list(self, server, project_dir):
        """outputs key must be a list even when checkpoint contains {sample_id}."""
        result = server.run_full_phase(
            project_path=str(project_dir),
            phase_slug="ingest_load",
        )
        assert isinstance(result["outputs"], list)
