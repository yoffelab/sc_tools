"""CLI discovery commands (DSC-01, DSC-02, DSC-03).

Provides machine-readable introspection for agents:

- ``sct list-commands --json`` -- catalog of all leaf commands with params
- ``sct describe <cmd>`` -- JSON schema for a specific command
- ``sct schema`` -- full CLI contract as a single JSON document

No heavy imports (scanpy, torch, scvi) -- only Typer/Click introspection.
"""

from __future__ import annotations

import click
import typer

from sc_tools.errors import SCToolsUserError
from sc_tools.models.result import CLIResult, Provenance, Status, _get_version

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_CLICK_TYPE_MAP: dict[str, str] = {
    "text": "string",
    "TEXT": "string",
    "integer": "integer",
    "INT": "integer",
    "float": "number",
    "FLOAT": "number",
    "boolean": "boolean",
    "BOOL": "boolean",
    "path": "string",
    "PATH": "string",
    "uuid": "string",
    "UUID": "string",
    "choice": "string",
    "Choice": "string",
    "STRING": "string",
    "FILENAME": "string",
}

# Commands to exclude from discovery output (meta/discovery commands, not pipeline operations)
_EXCLUDED_COMMANDS = {"version", "list-commands", "describe", "schema"}


# ---------------------------------------------------------------------------
# Introspection helpers
# ---------------------------------------------------------------------------


def _walk_commands(target_app: typer.Typer) -> dict[str, click.Command]:
    """Walk the Typer app tree and return {name: click.Command} for leaf commands.

    Skips stub groups with no subcommands and commands in _EXCLUDED_COMMANDS.
    """
    click_app = typer.main.get_command(target_app)
    ctx = click.Context(click_app)
    result: dict[str, click.Command] = {}

    for name in click_app.list_commands(ctx):
        if name in _EXCLUDED_COMMANDS:
            continue
        cmd = click_app.get_command(ctx, name)
        if cmd is None:
            continue
        if isinstance(cmd, click.Group):
            grp_ctx = click.Context(cmd, parent=ctx)
            sub_names = cmd.list_commands(grp_ctx)
            if not sub_names:
                # Stub group with no subcommands -- skip
                continue
            for sub_name in sub_names:
                sub_cmd = cmd.get_command(grp_ctx, sub_name)
                if sub_cmd is not None:
                    full_name = f"{name} {sub_name}"
                    if full_name not in _EXCLUDED_COMMANDS:
                        result[full_name] = sub_cmd
        else:
            result[name] = cmd

    return result


def _param_to_schema(param: click.Parameter) -> dict:
    """Convert a Click parameter to a JSON Schema property dict."""
    type_name = getattr(param.type, "name", "text") or "text"
    schema: dict = {"type": _CLICK_TYPE_MAP.get(type_name, "string")}

    if hasattr(param, "help") and param.help is not None:
        schema["description"] = param.help

    if param.default is not None:
        schema["default"] = param.default

    if hasattr(param.type, "choices") and param.type.choices:
        schema["enum"] = list(param.type.choices)

    return schema


def _command_to_entry(name: str, cmd: click.Command) -> dict:
    """Build a discovery catalog entry for a command."""
    properties: dict[str, dict] = {}
    required: list[str] = []

    for param in cmd.params:
        # Skip help and eager non-value params
        if param.name in ("help",):
            continue
        if getattr(param, "is_eager", False) and not getattr(param, "expose_value", True):
            continue

        param_name = param.name or ""
        # Convert underscores to hyphens for CLI-style naming, but keep
        # the Python-style name for JSON schema consistency
        properties[param_name] = _param_to_schema(param)

        if isinstance(param, click.Argument) and param.required:
            required.append(param_name)
        elif isinstance(param, click.Option) and param.required:
            required.append(param_name)

    params_schema: dict = {
        "type": "object",
        "properties": properties,
        "required": required,
    }

    group = name.split()[0] if " " in name else None
    description = cmd.help or cmd.short_help or ""

    return {
        "name": name,
        "group": group,
        "description": description,
        "params": params_schema,
    }


# ---------------------------------------------------------------------------
# Registration function
# ---------------------------------------------------------------------------


def register_discovery(target_app: typer.Typer) -> None:
    """Register discovery commands on the given Typer app."""
    from sc_tools.cli import cli_handler

    @target_app.command("list-commands")
    @cli_handler
    def list_commands(
        json_output: bool = typer.Option(False, "--json", help="JSON output (always JSON, flag for consistency)"),
    ) -> None:
        """List all available CLI commands with their parameters (DSC-01)."""
        from sc_tools.cli import app as main_app

        commands = _walk_commands(main_app)
        entries = [_command_to_entry(n, c) for n, c in sorted(commands.items())]

        return CLIResult(
            status=Status.success,
            command="list-commands",
            data={"commands": entries, "count": len(entries)},
            provenance=Provenance(command="list-commands"),
            message=f"{len(entries)} commands available",
        )

    @target_app.command("describe")
    @cli_handler
    def describe(
        command_name: str = typer.Argument(help="Command name (space-separated, e.g. 'qc run')"),
    ) -> None:
        """Describe a specific CLI command's parameters and output schema (DSC-02)."""
        from sc_tools.cli import app as main_app

        commands = _walk_commands(main_app)

        if command_name not in commands:
            available = sorted(commands.keys())
            raise SCToolsUserError(
                f"Unknown command: {command_name}",
                suggestion=f"Available commands: {', '.join(available)}",
            )

        entry = _command_to_entry(command_name, commands[command_name])
        entry["output_schema"] = CLIResult.model_json_schema()

        return CLIResult(
            status=Status.success,
            command="describe",
            data=entry,
            provenance=Provenance(command="describe"),
            message=f"Schema for '{command_name}'",
        )

    @target_app.command("schema")
    @cli_handler
    def schema() -> None:
        """Return the full CLI contract as a single JSON document (DSC-03)."""
        from sc_tools.cli import app as main_app

        commands = _walk_commands(main_app)
        commands_dict: dict = {}
        for name, cmd in sorted(commands.items()):
            entry = _command_to_entry(name, cmd)
            # Remove "name" key -- redundant since name is the dict key
            del entry["name"]
            commands_dict[name] = entry

        cli_schema = CLIResult.model_json_schema()
        defs = dict(cli_schema.get("$defs", {}))
        # Add CLIResult itself to $defs
        cli_result_def = {k: v for k, v in cli_schema.items() if k != "$defs"}
        defs["CLIResult"] = cli_result_def

        data = {
            "schema_version": "1.0",
            "sc_tools_version": _get_version(),
            "commands": commands_dict,
            "output": {"$ref": "#/$defs/CLIResult"},
            "$defs": defs,
        }

        return CLIResult(
            status=Status.success,
            command="schema",
            data=data,
            provenance=Provenance(command="schema"),
            message=f"Full CLI schema: {len(commands_dict)} commands, {len(defs)} shared types",
        )
