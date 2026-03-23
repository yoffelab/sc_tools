"""Pipeline status commands (CMD-05)."""

from __future__ import annotations

import typer

status_app = typer.Typer(help="Pipeline status commands")


@status_app.callback(invoke_without_command=True)
def status_callback(ctx: typer.Context) -> None:
    """Pipeline status commands. Run 'sct status --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
