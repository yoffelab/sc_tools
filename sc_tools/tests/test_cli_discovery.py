"""CLI discovery command tests (DSC-01, DSC-02, DSC-03).

Tests for `sct list-commands`, `sct describe`, and `sct schema` commands
that provide machine-readable CLI introspection for agents.
"""

from __future__ import annotations

import json
import sys

import pytest
from typer.testing import CliRunner

from sc_tools.cli import app

runner = CliRunner()


class TestListCommands:
    """DSC-01: sct list-commands --json returns machine-readable catalog."""

    def test_returns_success_json(self):
        result = runner.invoke(app, ["list-commands", "--json"])
        assert result.exit_code == 0, f"exit_code={result.exit_code}, output={result.output}"
        data = json.loads(result.output)
        assert data["status"] == "success"

    def test_contains_known_commands(self):
        result = runner.invoke(app, ["list-commands", "--json"])
        data = json.loads(result.output)
        commands = data["data"]["commands"]
        names = {c["name"] for c in commands}
        for expected in ("qc run", "benchmark integration", "preprocess run", "validate run", "status show", "report generate"):
            assert expected in names, f"Missing command: {expected}"

    def test_excludes_meta_commands(self):
        result = runner.invoke(app, ["list-commands", "--json"])
        data = json.loads(result.output)
        commands = data["data"]["commands"]
        names = {c["name"] for c in commands}
        for excluded in ("version", "list-commands", "describe", "schema"):
            assert excluded not in names, f"Meta command should be excluded: {excluded}"

    def test_param_types(self):
        result = runner.invoke(app, ["list-commands", "--json"])
        data = json.loads(result.output)
        commands = data["data"]["commands"]
        qc_run = next(c for c in commands if c["name"] == "qc run")
        params = qc_run["params"]
        assert params["type"] == "object"
        assert "file" in params["properties"]
        assert params["properties"]["file"]["type"] == "string"
        assert "file" in params["required"]

    def test_param_defaults(self):
        result = runner.invoke(app, ["list-commands", "--json"])
        data = json.loads(result.output)
        commands = data["data"]["commands"]
        qc_run = next(c for c in commands if c["name"] == "qc run")
        modality_prop = qc_run["params"]["properties"]["modality"]
        assert modality_prop["default"] == "visium"

    def test_count_matches(self):
        result = runner.invoke(app, ["list-commands", "--json"])
        data = json.loads(result.output)
        assert data["data"]["count"] == len(data["data"]["commands"])


class TestDescribe:
    """DSC-02: sct describe <cmd> returns JSON schema for command."""

    def test_known_command(self):
        result = runner.invoke(app, ["describe", "qc run"])
        assert result.exit_code == 0, f"exit_code={result.exit_code}, output={result.output}"
        data = json.loads(result.output)
        assert data["status"] == "success"
        assert data["data"]["name"] == "qc run"
        assert "params" in data["data"]
        assert "properties" in data["data"]["params"]
        assert "output_schema" in data["data"]

    def test_unknown_command(self):
        result = runner.invoke(app, ["describe", "nonexistent"])
        assert result.exit_code == 1, f"exit_code={result.exit_code}, output={result.output}"
        data = json.loads(result.output)
        assert data["status"] == "error"
        assert data["error"]["suggestion"]  # should list available commands

    def test_output_schema_has_defs(self):
        result = runner.invoke(app, ["describe", "qc run"])
        data = json.loads(result.output)
        output_schema = data["data"]["output_schema"]
        defs = output_schema.get("$defs", {})
        for expected_def in ("Status", "ErrorInfo", "Provenance"):
            assert expected_def in defs, f"Missing $def: {expected_def}"


class TestSchema:
    """DSC-03: sct schema returns full CLI contract."""

    def test_returns_success(self):
        result = runner.invoke(app, ["schema"])
        assert result.exit_code == 0, f"exit_code={result.exit_code}, output={result.output}"
        data = json.loads(result.output)
        assert data["status"] == "success"

    def test_commands_keyed_by_name(self):
        result = runner.invoke(app, ["schema"])
        data = json.loads(result.output)
        commands = data["data"]["commands"]
        assert isinstance(commands, dict)
        assert "qc run" in commands

    def test_defs_present(self):
        result = runner.invoke(app, ["schema"])
        data = json.loads(result.output)
        defs = data["data"]["$defs"]
        for expected in ("CLIResult", "Status", "ErrorInfo", "Provenance"):
            assert expected in defs, f"Missing $def: {expected}"

    def test_schema_version(self):
        result = runner.invoke(app, ["schema"])
        data = json.loads(result.output)
        assert isinstance(data["data"]["schema_version"], str)
        assert len(data["data"]["schema_version"]) > 0

    def test_sc_tools_version(self):
        result = runner.invoke(app, ["schema"])
        data = json.loads(result.output)
        assert isinstance(data["data"]["sc_tools_version"], str)
        assert len(data["data"]["sc_tools_version"]) > 0


class TestOutputFormat:
    """All discovery commands return CLIResult envelope."""

    def test_all_discovery_commands_return_cli_result(self):
        invocations = [
            ["list-commands", "--json"],
            ["describe", "qc run"],
            ["schema"],
        ]
        for args in invocations:
            result = runner.invoke(app, args)
            assert result.exit_code == 0, f"Failed for {args}: exit_code={result.exit_code}, output={result.output}"
            data = json.loads(result.output)
            for key in ("status", "command", "provenance"):
                assert key in data, f"Missing key '{key}' in response for {args}"


class TestLazyImports:
    """Discovery module must not import heavy dependencies."""

    def test_no_heavy_imports(self):
        # Remove discovery module if already imported so we can test cleanly
        mods_to_remove = [k for k in sys.modules if k.startswith("sc_tools.cli.discovery")]
        for m in mods_to_remove:
            del sys.modules[m]

        # Record heavy modules before import
        heavy = {"scanpy", "torch", "scvi"}
        before = {m for m in heavy if m in sys.modules}

        import sc_tools.cli.discovery  # noqa: F401

        after = {m for m in heavy if m in sys.modules}
        newly_imported = after - before
        assert not newly_imported, f"Heavy modules imported by discovery: {newly_imported}"
