"""CLI argument parsing tests (TST-04).

No data is loaded -- these test the CLI framework, not pipeline logic.
"""

from __future__ import annotations

import json

import pytest
from typer.testing import CliRunner

from sc_tools.cli import app

runner = CliRunner()


class TestHelp:
    """sct --help and subcommand help."""

    def test_sct_help_returns_zero(self):
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0

    def test_sct_help_shows_command_groups(self):
        result = runner.invoke(app, ["--help"])
        for group in ("qc", "preprocess", "integrate", "benchmark", "celltype"):
            assert group in result.output, f"{group} missing from help"

    def test_sct_help_shows_human_flag(self):
        result = runner.invoke(app, ["--help"])
        assert "--human" in result.output

    def test_subcommand_help(self):
        for group in ("qc", "preprocess", "integrate", "benchmark", "celltype"):
            result = runner.invoke(app, [group, "--help"])
            assert result.exit_code == 0, f"{group} --help failed"


class TestVersion:
    """sct --version."""

    def test_version_returns_zero(self):
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0

    def test_version_output_not_empty(self):
        result = runner.invoke(app, ["--version"])
        assert result.output.strip() != ""

    def test_version_command(self):
        result = runner.invoke(app, ["version"])
        assert result.exit_code == 0
        assert result.output.strip() != ""


class TestHumanFlag:
    """--human flag is accepted globally."""

    def test_human_flag_accepted(self):
        result = runner.invoke(app, ["--human", "--help"])
        assert result.exit_code == 0

    def test_human_flag_with_subcommand(self):
        result = runner.invoke(app, ["--human", "qc", "--help"])
        assert result.exit_code == 0


class TestUnknownCommand:
    """Unknown commands fail gracefully."""

    def test_unknown_command_exits_nonzero(self):
        result = runner.invoke(app, ["nonexistent_command_xyz"])
        assert result.exit_code != 0


class TestCLIResult:
    """CLIResult model tests (no CLI invocation, just model validation)."""

    def test_cli_result_json_valid(self):
        from sc_tools.models.result import CLIResult, Provenance, Status

        r = CLIResult(
            status=Status.success,
            command="test",
            provenance=Provenance(command="test"),
            message="ok",
        )
        d = json.loads(r.model_dump_json())
        assert d["status"] == "success"
        assert "command" in d
        assert "data" in d
        assert "artifacts" in d
        assert "provenance" in d
        assert "message" in d

    def test_cli_result_error_envelope(self):
        from sc_tools.models.result import CLIResult, ErrorInfo, Provenance, Status

        r = CLIResult(
            status=Status.error,
            command="fail",
            provenance=Provenance(command="fail"),
            error=ErrorInfo(category="fixable", suggestion="try again"),
        )
        d = json.loads(r.model_dump_json())
        assert d["error"]["category"] == "fixable"
        assert d["error"]["suggestion"] == "try again"

    def test_cli_result_partial_failures(self):
        from sc_tools.models.result import CLIResult, Provenance, Status

        r = CLIResult(
            status=Status.partial,
            command="qc",
            provenance=Provenance(command="qc"),
            partial_failures=[{"sample": "s1", "reason": "missing file"}],
        )
        d = json.loads(r.model_dump_json())
        assert d["partial_failures"][0]["sample"] == "s1"

    def test_status_enum_values(self):
        from sc_tools.models.result import Status

        assert set(s.value for s in Status) == {"success", "error", "partial", "skipped"}

    def test_dual_serialization(self):
        """model_dump_json for CLI, model_dump(mode='json') for MCP."""
        from sc_tools.models.result import CLIResult, Provenance, Status

        r = CLIResult(
            status=Status.success,
            command="test",
            provenance=Provenance(command="test"),
        )
        json_str = r.model_dump_json()
        json_dict = r.model_dump(mode="json")
        assert isinstance(json_str, str)
        assert isinstance(json_dict, dict)
        assert json_dict["status"] == "success"


class TestErrorTaxonomy:
    """Exception hierarchy maps to correct exit codes."""

    def test_user_error_exit_code(self):
        from sc_tools.errors import SCToolsUserError

        assert SCToolsUserError.exit_code == 1

    def test_data_error_exit_code(self):
        from sc_tools.errors import SCToolsDataError

        assert SCToolsDataError.exit_code == 2

    def test_runtime_error_exit_code(self):
        from sc_tools.errors import SCToolsRuntimeError

        assert SCToolsRuntimeError.exit_code == 3

    def test_fatal_error_exit_code(self):
        from sc_tools.errors import SCToolsFatalError

        assert SCToolsFatalError.exit_code == 3

    def test_error_categories(self):
        from sc_tools.errors import (
            SCToolsDataError,
            SCToolsFatalError,
            SCToolsRuntimeError,
            SCToolsUserError,
        )

        assert SCToolsUserError.category == "fixable"
        assert SCToolsDataError.category == "fixable"
        assert SCToolsRuntimeError.category == "retryable"
        assert SCToolsFatalError.category == "fatal"

    def test_error_with_suggestion(self):
        from sc_tools.errors import SCToolsUserError

        e = SCToolsUserError("bad input", suggestion="use --flag")
        assert str(e) == "bad input"
        assert e.suggestion == "use --flag"


class TestNonInteractive:
    """CLI-07: No interactive prompts."""

    def test_no_prompt_in_cli_module(self):
        """cli.py must not use typer.confirm or typer.prompt."""
        import inspect

        import sc_tools.cli as cli_mod

        source = inspect.getsource(cli_mod)
        assert "typer.confirm" not in source, "cli.py contains typer.confirm (violates CLI-07)"
        assert "typer.prompt" not in source, "cli.py contains typer.prompt (violates CLI-07)"


class TestLazyImports:
    """CLI-08: No heavy imports at startup."""

    def test_no_heavy_imports_in_source(self):
        """cli.py source must not have top-level heavy imports."""
        import inspect

        import sc_tools.cli as cli_mod

        source = inspect.getsource(cli_mod)
        # Check for top-level import statements (not inside functions)
        # We look for lines that start with "import scanpy" or "from scanpy"
        lines = source.split("\n")
        for i, line in enumerate(lines):
            stripped = line.strip()
            # Skip lines inside function bodies (indented)
            if line and not line[0].isspace():
                assert not stripped.startswith("import scanpy"), (
                    f"Line {i}: top-level 'import scanpy' in cli.py"
                )
                assert not stripped.startswith("import torch"), (
                    f"Line {i}: top-level 'import torch' in cli.py"
                )
                assert not stripped.startswith("import scvi"), (
                    f"Line {i}: top-level 'import scvi' in cli.py"
                )
                assert not stripped.startswith("from scanpy"), (
                    f"Line {i}: top-level 'from scanpy' in cli.py"
                )
                assert not stripped.startswith("from torch"), (
                    f"Line {i}: top-level 'from torch' in cli.py"
                )
                assert not stripped.startswith("from scvi"), (
                    f"Line {i}: top-level 'from scvi' in cli.py"
                )

    def test_startup_note(self):
        """Placeholder: authoritative check is `time sct --help` < 500ms."""
        # This test is a documentation anchor. The real validation is manual:
        # time sct --help  (should complete in < 500ms)
        pass
