"""Typer CLI application for sc_tools.

cli/ package (D-01) providing the ``sct`` command with:

- Global ``--human`` flag for Rich output to stderr (D-04, CLI-04)
- Error handler mapping exceptions to exit codes 0/1/2/3 (D-12, D-15)
- Dependency check utility _check_deps (CMD-08, D-11)
- Seven command groups: qc, preprocess, validate, status, integrate, benchmark, celltype

No heavy imports (scanpy, torch, scvi) at module level (CLI-08).
"""

from __future__ import annotations

import functools
import sys

import typer

from sc_tools.errors import (
    SCToolsDataError,
    SCToolsFatalError,
    SCToolsRuntimeError,
    SCToolsUserError,
)
from sc_tools.models.result import CLIResult, ErrorInfo, Provenance, Status

# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------

app = typer.Typer(
    name="sct",
    help="sct -- sc_tools command-line interface.",
    pretty_exceptions_enable=False,  # We handle errors ourselves (Pitfall 3)
    add_completion=False,
)

# ---------------------------------------------------------------------------
# Global state & callback
# ---------------------------------------------------------------------------

_state: dict[str, bool] = {"human": False}


def _version_callback(value: bool) -> None:
    """Print version and exit when --version is passed."""
    if value:
        from importlib.metadata import version as get_version

        try:
            v = get_version("sci-sc-tools")
        except Exception:
            v = "unknown"
        typer.echo(v)
        raise typer.Exit()


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    human: bool = typer.Option(False, "--human", help="Rich-formatted output to stderr"),
    version: bool = typer.Option(  # noqa: ARG001
        False,
        "--version",
        help="Show version and exit",
        is_eager=True,
        callback=_version_callback,
    ),
) -> None:
    """sct -- sc_tools command-line interface."""
    _state["human"] = human
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------


def _emit(result: CLIResult) -> None:
    """Write CLIResult JSON to stdout. If --human, also render Rich to stderr."""
    sys.stdout.write(result.model_dump_json(indent=2) + "\n")
    sys.stdout.flush()
    if _state["human"]:
        _render_rich(result)


def _render_rich(result: CLIResult) -> None:
    """Render CLIResult as Rich formatted output to stderr."""
    from rich.console import Console
    from rich.panel import Panel

    console = Console(stderr=True)
    style = (
        "green"
        if result.status == Status.success
        else "red"
        if result.status == Status.error
        else "yellow"
    )
    console.print(
        Panel(
            result.message or str(result.status.value if isinstance(result.status, Status) else result.status),
            title=result.command,
            style=style,
        )
    )
    if result.error:
        console.print(f"[red]Error:[/red] {result.error.suggestion}")
    if result.data:
        from rich.table import Table

        table = Table(title="Data")
        table.add_column("Key")
        table.add_column("Value")
        for k, v in result.data.items():
            table.add_row(str(k), str(v))
        console.print(table)


# ---------------------------------------------------------------------------
# Dependency check utility (CMD-08, D-11)
# ---------------------------------------------------------------------------

_DEP_INSTALL: dict[str, str] = {
    "scanpy": "pip install scanpy",
    "scvi-tools": "pip install scvi-tools",
    "scib_metrics": "pip install scib-metrics",
    "h5py": "pip install h5py",
    "rapids_singlecell": "pip install rapids-singlecell (requires CUDA)",
}


def _check_deps(deps: list[str]) -> None:
    """Check optional dependencies before loading data (CMD-08). Raises SCToolsUserError."""
    missing = []
    for dep in deps:
        try:
            __import__(dep.replace("-", "_"))
        except ImportError:
            install = _DEP_INSTALL.get(dep, f"pip install {dep}")
            missing.append(f"  {dep}: {install}")
    if missing:
        raise SCToolsUserError(
            "Missing required dependencies:\n" + "\n".join(missing),
            suggestion="Install missing packages and retry",
        )


# ---------------------------------------------------------------------------
# Error handler decorator
# ---------------------------------------------------------------------------


def _make_error_result(exc: Exception, command: str) -> CLIResult:
    """Build a CLIResult from an sc_tools exception."""
    from sc_tools.errors import SCToolsError

    category = getattr(exc, "category", "fatal")
    suggestion = getattr(exc, "suggestion", "This is an unexpected error. Please report it.")
    details = getattr(exc, "details", None)

    return CLIResult(
        status=Status.error,
        command=command,
        provenance=Provenance(command=command),
        message=str(exc),
        error=ErrorInfo(
            category=category,
            suggestion=suggestion,
            details=details,
        ),
    )


def cli_handler(func):  # noqa: ANN001, ANN201
    """Wrap a CLI command with error handling and JSON output."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):  # noqa: ANN002, ANN003, ANN202
        try:
            result = func(*args, **kwargs)
            _emit(result)
            raise SystemExit(0)
        except SystemExit:
            raise
        except SCToolsUserError as e:
            _emit(_make_error_result(e, func.__name__))
            raise SystemExit(1)
        except SCToolsDataError as e:
            _emit(_make_error_result(e, func.__name__))
            raise SystemExit(2)
        except SCToolsRuntimeError as e:
            _emit(_make_error_result(e, func.__name__))
            raise SystemExit(3)
        except MemoryError:
            err = SCToolsRuntimeError(
                "Out of memory",
                suggestion="Retry with --subsample-n or reduce dataset size",
            )
            _emit(_make_error_result(err, func.__name__))
            raise SystemExit(3)
        except Exception as e:
            err = SCToolsFatalError(
                str(e),
                suggestion="This is an unexpected error. Please report it.",
            )
            _emit(_make_error_result(err, func.__name__))
            raise SystemExit(3)

    return wrapper


# ---------------------------------------------------------------------------
# Stub command groups kept in __init__.py (no Phase 3 commands for these)
# ---------------------------------------------------------------------------

integrate_app = typer.Typer(help="Integration commands")
celltype_app = typer.Typer(help="Cell typing commands")


@integrate_app.callback(invoke_without_command=True)
def integrate_callback(ctx: typer.Context) -> None:
    """Integration commands. Run 'sct integrate --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


@celltype_app.callback(invoke_without_command=True)
def celltype_callback(ctx: typer.Context) -> None:
    """Cell typing commands. Run 'sct celltype --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


# ---------------------------------------------------------------------------
# Version command
# ---------------------------------------------------------------------------


@app.command()
def version() -> None:
    """Show sc_tools version."""
    from importlib.metadata import version as get_version

    try:
        v = get_version("sci-sc-tools")
    except Exception:
        v = "unknown"
    typer.echo(v)


# ---------------------------------------------------------------------------
# Register subcommands (after all helpers are defined)
# ---------------------------------------------------------------------------

from sc_tools.cli.qc import qc_app, report_app  # noqa: E402
from sc_tools.cli.preprocess import preprocess_app  # noqa: E402
from sc_tools.cli.validate import validate_app  # noqa: E402
from sc_tools.cli.benchmark import benchmark_app  # noqa: E402
from sc_tools.cli.status import status_app  # noqa: E402

app.add_typer(qc_app, name="qc")
app.add_typer(report_app, name="report")
app.add_typer(preprocess_app, name="preprocess")
app.add_typer(validate_app, name="validate")
app.add_typer(status_app, name="status")
app.add_typer(integrate_app, name="integrate")
app.add_typer(benchmark_app, name="benchmark")
app.add_typer(celltype_app, name="celltype")

from sc_tools.cli.discovery import register_discovery  # noqa: E402
register_discovery(app)
