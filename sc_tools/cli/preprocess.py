"""Preprocessing commands (CMD-02)."""

from __future__ import annotations

import typer

preprocess_app = typer.Typer(help="Preprocessing commands")


@preprocess_app.callback(invoke_without_command=True)
def preprocess_callback(ctx: typer.Context) -> None:
    """Preprocessing commands. Run 'sct preprocess --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
