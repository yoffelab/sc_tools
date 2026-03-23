"""Pipeline status commands (CMD-05)."""

from __future__ import annotations

import typer

from sc_tools.cli import cli_handler

status_app = typer.Typer(help="Pipeline status commands")


@status_app.callback(invoke_without_command=True)
def status_callback(ctx: typer.Context) -> None:
    """Pipeline status commands. Run 'sct status --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        # Default: run the status command directly
        status_show()


@status_app.command("show")
@cli_handler
def status_show(
    project_dir: str = typer.Option(".", "--project-dir", help="Project directory (per D-03)"),
) -> None:
    """Show pipeline phase DAG state: completed phases, available next, checkpoint paths."""
    from sc_tools.models.result import CLIResult, Provenance, Status
    from sc_tools.pipeline import get_available_next, get_dag, tuple_to_display

    dag = get_dag()
    completed: list[tuple[str, str]] = []
    registry_available = False

    # D-09: graceful fallback when registry unavailable
    try:
        from sc_tools.db import get_engine  # noqa: F401

        # Attempt registry query -- if it works, populate completed
        # For now, registry query is best-effort
        registry_available = True
    except Exception:
        pass  # Registry unavailable -- show DAG without completion status

    available_next = get_available_next(completed)

    phases_data = {}
    for key, spec in dag.items():
        display = tuple_to_display(key)
        phases_data[display] = {
            "label": spec.label,
            "checkpoint": spec.checkpoint,
            "branch": spec.branch,
            "depends_on": [tuple_to_display(d) if isinstance(d, tuple) else d for d in spec.depends_on],
            "completed": key in {tuple(c) if isinstance(c, list) else c for c in completed},
            "optional": spec.optional,
        }

    data = {
        "phases": phases_data,
        "completed": [tuple_to_display(c) for c in completed],
        "available_next": [tuple_to_display(a) for a in available_next],
        "registry_available": registry_available,
        "total_phases": len(dag),
    }

    n_completed = len(completed)
    n_available = len(available_next)
    msg = f"{n_completed}/{len(dag)} phases complete, {n_available} available next"
    if not registry_available:
        msg += " (registry unavailable -- showing DAG structure only)"

    return CLIResult(
        status=Status.success,
        command="status",
        data=data,
        provenance=Provenance(command="status"),
        message=msg,
    )
