"""Validation commands (CMD-03)."""

from __future__ import annotations

import typer

validate_app = typer.Typer(help="Validation commands")


@validate_app.callback(invoke_without_command=True)
def validate_callback(ctx: typer.Context) -> None:
    """Validation commands. Run 'sct validate --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
