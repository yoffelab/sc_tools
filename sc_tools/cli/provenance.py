"""Provenance CLI commands (PRV-03, PRV-04).

Provides:

- ``sct provenance show <file>`` -- display provenance from sidecar or adata.uns
- ``sct provenance trace <file>`` -- walk lineage chain to origins

Follows the register_discovery(app) pattern from discovery.py (D-11).
No heavy imports at module level (CLI-08).
"""

from __future__ import annotations

import typer

from sc_tools.errors import SCToolsUserError
from sc_tools.models.result import CLIResult, Provenance, Status


def register_provenance(target_app: typer.Typer) -> None:
    """Register provenance commands on the given Typer app."""
    from sc_tools.cli import cli_handler

    provenance_app = typer.Typer(help="Provenance tracking commands")

    @provenance_app.command("show")
    @cli_handler
    def show(
        file: str = typer.Argument(..., help="Path to artifact file"),
    ) -> None:
        """Display provenance record for a file (PRV-03, D-13).

        Reads from .provenance.json sidecar first; falls back to
        adata.uns['sct_provenance'] for h5ad files.
        """
        from sc_tools.provenance.sidecar import read_sidecar

        prov = read_sidecar(file)

        # Fallback to adata.uns for h5ad files (D-05)
        if prov is None and file.endswith(".h5ad"):
            from sc_tools.provenance.trace import _read_uns_provenance
            from pathlib import Path

            prov = _read_uns_provenance(Path(file))

        if prov is None:
            raise SCToolsUserError(
                "No provenance found",
                suggestion="File may be a raw input or provenance sidecar was deleted",
            )

        return CLIResult(
            status=Status.success,
            command="provenance show",
            data=prov,
            artifacts=[],
            provenance=Provenance(command="provenance show"),
            message=f"Provenance for {file}",
        )

    @provenance_app.command("trace")
    @cli_handler
    def trace(
        file: str = typer.Argument(..., help="Path to artifact file"),
        project_dir: str = typer.Option(
            ".", "--project-dir", help="Project root for relative path resolution"
        ),
    ) -> None:
        """Walk lineage chain to origins (PRV-04, D-11, D-12).

        Follows provenance sidecar input references recursively.
        Output is chronological (oldest first).
        """
        from sc_tools.provenance.trace import trace_lineage

        steps = trace_lineage(file, project_dir=project_dir)

        return CLIResult(
            status=Status.success,
            command="provenance trace",
            data={"file": file, "lineage": steps, "depth": len(steps)},
            artifacts=[],
            provenance=Provenance(command="provenance trace"),
            message=f"Lineage trace: {len(steps)} step(s)",
        )

    target_app.add_typer(provenance_app, name="provenance")
