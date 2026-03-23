"""Benchmarking commands (CMD-04)."""

from __future__ import annotations

import typer

benchmark_app = typer.Typer(help="Benchmarking commands")


@benchmark_app.callback(invoke_without_command=True)
def benchmark_callback(ctx: typer.Context) -> None:
    """Benchmarking commands. Run 'sct benchmark --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
