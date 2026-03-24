"""Validation commands (CMD-03)."""

from __future__ import annotations

import typer

validate_app = typer.Typer(help="Validation commands")


@validate_app.callback(invoke_without_command=True)
def validate_callback(ctx: typer.Context) -> None:
    """Validation commands. Run 'sct validate --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


@validate_app.command("run")
def validate_run(
    phase: str = typer.Argument(..., help="Phase slug: qc_filter, metadata_attach, preprocess, scoring, celltype_manual"),
    file: str = typer.Argument(..., help="Path to .h5ad checkpoint file"),
    fix: bool = typer.Option(False, "--fix", help="Attempt auto-fixes"),
    project_dir: str = typer.Option(".", "--project-dir", help="Project directory (per D-03)"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate inputs and report plan without executing"),
) -> None:
    """Validate a checkpoint file against phase requirements."""
    from sc_tools.cli import _check_deps, _emit
    from sc_tools.errors import SCToolsUserError
    from sc_tools.models.result import CLIResult, Provenance, Status

    _check_deps(["anndata"])

    from pathlib import Path

    path = Path(file)
    if not path.exists():
        raise SCToolsUserError(f"File not found: {path}", suggestion="Check the file path and try again")

    if dry_run:
        # Dry-run: validate file exists, read T1 metadata, return early
        dry_data: dict = {"planned_command": "validate_run", "tier": "metadata", "phase": phase, "file": str(path)}
        try:
            from sc_tools.io.estimate import estimate_from_h5

            est = estimate_from_h5(str(path))
            dry_data["estimate"] = est
        except Exception:
            dry_data["estimate"] = None
        result = CLIResult(
            status=Status.skipped,
            command=f"validate {phase}",
            data=dry_data,
            provenance=Provenance(command=f"validate {phase}"),
            message="Dry run -- no data modified",
        )
        _emit(result)
        raise SystemExit(0)

    from sc_tools.validate import validate_file

    try:
        issues = validate_file(str(path), phase=phase, strict=False, fix=fix)
    except ValueError as e:
        raise SCToolsUserError(str(e), suggestion="Use a valid phase slug: qc_filter, metadata_attach, preprocess, scoring, celltype_manual")

    status = Status.success if not issues else Status.error
    result = CLIResult(
        status=status,
        command=f"validate {phase}",
        data={"phase": phase, "file": str(path), "issues": issues, "n_issues": len(issues), "fix_applied": fix},
        artifacts=[],
        provenance=Provenance(command=f"validate {phase}"),
        message=f"Validation {'passed' if not issues else 'failed'}: {len(issues)} issue(s) found",
    )

    _emit(result)
    raise SystemExit(0 if status == Status.success else 2)
