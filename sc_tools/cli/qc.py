"""QC and report commands (CMD-01, CMD-06)."""

from __future__ import annotations

import typer

qc_app = typer.Typer(help="Quality control commands")


@qc_app.callback(invoke_without_command=True)
def qc_callback(ctx: typer.Context) -> None:
    """Quality control commands. Run 'sct qc --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
