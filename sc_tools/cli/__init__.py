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
import logging
import sys
import time

import typer

from sc_tools.errors import (
    SCToolsDataError,
    SCToolsFatalError,
    SCToolsRuntimeError,
    SCToolsUserError,
)
from sc_tools.models.result import CLIResult, ErrorInfo, Provenance, Status

_cli_logger = logging.getLogger(__name__)

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
    "pydeseq2": "pip install 'pydeseq2>=0.5.0'",
    "mudata": "pip install 'mudata>=0.3.1'",
    "muon": "pip install 'muon>=0.1.6'",
    "mofapy2": "pip install 'mofapy2>=0.7.0'",
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


def _write_provenance_sidecars(
    result: CLIResult,
    kwargs: dict,
    input_files: list[str],
    runtime_s: float,
    peak_mb: float,
) -> None:
    """Write .provenance.json sidecars for each artifact (D-01).

    Additionally embeds provenance in adata.uns for h5ad artifacts (D-04).
    """
    from sc_tools.provenance.sidecar import (
        build_provenance_record,
        embed_provenance_in_adata,
        write_sidecar,
    )

    # Clean kwargs: remove typer internals
    params = {k: v for k, v in kwargs.items() if not k.startswith("_")}

    record = build_provenance_record(
        command=result.command,
        params=params,
        input_files=input_files,
        runtime_s=runtime_s,
        peak_memory_mb=peak_mb,
    )

    for artifact in result.artifacts:
        try:
            write_sidecar(artifact, record)
            # D-04: embed in adata.uns for h5ad artifacts
            if str(artifact).endswith(".h5ad"):
                embed_provenance_in_adata(artifact, record)
        except Exception:
            _cli_logger.warning(
                "Failed to write provenance sidecar for %s", artifact, exc_info=True
            )


def cli_handler(func=None, *, tier=None):  # noqa: ANN001, ANN201
    """Wrap a CLI command with error handling, JSON output, and dry-run support.

    Can be used as ``@cli_handler`` (no args) or ``@cli_handler(tier=DataTier.T3_FULL)``.
    When ``dry_run=True`` is passed as a kwarg, the wrapper validates inputs,
    reports planned operations and memory estimate, then exits 0 with
    ``status=skipped`` without executing the wrapped function.
    """

    def decorator(fn):  # noqa: ANN001, ANN202
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):  # noqa: ANN002, ANN003, ANN202
            dry_run = kwargs.pop("dry_run", False)
            _force = kwargs.pop("force", False)
            # Always pass dry_run and force explicitly so inner functions never
            # fall back to typer.Option(...) defaults, which are truthy OptionInfo objects.
            kwargs["dry_run"] = False
            kwargs["force"] = False
            start_time = time.monotonic()

            try:
                if dry_run:
                    # Build dry-run result without executing the command
                    from pathlib import Path

                    file_arg = (
                        kwargs.get("file")
                        or kwargs.get("input_file")
                        or kwargs.get("from_dir")
                    )
                    tier_label = tier.value if tier else "full"
                    dry_data: dict = {"planned_command": fn.__name__, "tier": tier_label}

                    if file_arg and Path(str(file_arg)).exists():
                        try:
                            from sc_tools.io.estimate import estimate_from_h5

                            est = estimate_from_h5(str(file_arg))
                            dry_data["estimate"] = est
                        except Exception:
                            dry_data["estimate"] = None
                    elif file_arg:
                        raise SCToolsUserError(
                            f"Input file not found: {file_arg}",
                            suggestion="Check the file path",
                        )

                    result = CLIResult(
                        status=Status.skipped,
                        command=fn.__name__,
                        data=dry_data,
                        provenance=Provenance(command=fn.__name__),
                        message="Dry run -- no data modified",
                    )
                    _emit(result)
                    raise SystemExit(0)

                result = fn(*args, **kwargs)

                runtime_s = time.monotonic() - start_time

                # Extract _input_files convention key before emitting (D-07)
                input_files = result.data.pop("_input_files", [])

                _emit(result)

                # Sidecar writing: only on success with artifacts (D-01, D-03)
                if result.status == Status.success and result.artifacts:
                    try:
                        from sc_tools.provenance.sidecar import get_peak_memory_mb

                        peak_mb = get_peak_memory_mb()
                        _write_provenance_sidecars(result, kwargs, input_files, runtime_s, peak_mb)
                    except Exception:
                        _cli_logger.warning("Provenance sidecar writing failed", exc_info=True)

                raise SystemExit(0)
            except SystemExit:
                raise
            except SCToolsUserError as e:
                _emit(_make_error_result(e, fn.__name__))
                raise SystemExit(1)  # noqa: B904
            except SCToolsDataError as e:
                _emit(_make_error_result(e, fn.__name__))
                raise SystemExit(2)  # noqa: B904
            except SCToolsRuntimeError as e:
                _emit(_make_error_result(e, fn.__name__))
                raise SystemExit(3)  # noqa: B904
            except MemoryError:
                err = SCToolsRuntimeError(
                    "Out of memory",
                    suggestion="Retry with --subsample-n or reduce dataset size",
                )
                _emit(_make_error_result(err, fn.__name__))
                raise SystemExit(3)  # noqa: B904
            except Exception as e:
                err = SCToolsFatalError(
                    str(e),
                    suggestion="This is an unexpected error. Please report it.",
                )
                _emit(_make_error_result(err, fn.__name__))
                raise SystemExit(3)  # noqa: B904

        return wrapper

    if func is not None:
        # Called as @cli_handler (no parentheses)
        return decorator(func)
    # Called as @cli_handler(tier=...) -- return the decorator
    return decorator


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

from sc_tools.cli.benchmark import benchmark_app  # noqa: E402
from sc_tools.cli.de import de_app  # noqa: E402
from sc_tools.cli.preprocess import preprocess_app  # noqa: E402
from sc_tools.cli.qc import qc_app, report_app  # noqa: E402
from sc_tools.cli.status import status_app  # noqa: E402
from sc_tools.cli.validate import validate_app  # noqa: E402

app.add_typer(qc_app, name="qc")
app.add_typer(report_app, name="report")
app.add_typer(preprocess_app, name="preprocess")
app.add_typer(validate_app, name="validate")
app.add_typer(status_app, name="status")
app.add_typer(integrate_app, name="integrate")
app.add_typer(benchmark_app, name="benchmark")
app.add_typer(celltype_app, name="celltype")
app.add_typer(de_app, name="de")

from sc_tools.cli.discovery import register_discovery  # noqa: E402

register_discovery(app)

from sc_tools.cli.provenance import register_provenance  # noqa: E402

register_provenance(app)

from sc_tools.cli.estimate import register_estimate  # noqa: E402

register_estimate(app)

from sc_tools.cli.assemble import register_assemble  # noqa: E402

register_assemble(app)

from sc_tools.cli.concat import register_concat  # noqa: E402

register_concat(app)
