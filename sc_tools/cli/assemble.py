"""Multi-omic assembly commands (MOM-03, MOM-04).

Provides ``sct assemble build``, ``sct assemble embed``, and ``sct assemble query``
subcommands for building, embedding, and querying multi-omic atlases.
"""

from __future__ import annotations

import logging

import typer

logger = logging.getLogger(__name__)


def register_assemble(app: typer.Typer) -> None:
    """Register the assemble command group on the given Typer app."""
    from sc_tools.cli import _check_deps, cli_handler
    from sc_tools.models.result import CLIResult, Provenance, Status

    assemble_app = typer.Typer(help="Multi-omic assembly commands")

    @assemble_app.callback(invoke_without_command=True)
    def assemble_callback(ctx: typer.Context) -> None:
        """Multi-omic assembly commands. Run 'sct assemble --help' for subcommands."""
        if ctx.invoked_subcommand is None:
            typer.echo(ctx.get_help())

    # --- build subcommand ---

    @assemble_app.command("build")
    @cli_handler
    def assemble_build(
        inputs: list[str] = typer.Argument(..., help="Per-modality h5ad file paths"),
        modalities: list[str] = typer.Option(
            ..., "--modality", "-m", help="Modality label for each input (rna, imc, visium, xenium)"
        ),
        output: str = typer.Option("atlas.h5mu", "--output", "-o", help="Output h5mu file path"),
        subject_key: str = typer.Option("subject_id", "--subject-key", help="Column for subject ID join"),
        dry_run: bool = typer.Option(False, "--dry-run", help="Validate without building"),
        force: bool = typer.Option(False, "--force", help="Bypass memory guard"),
    ) -> None:
        """Build a multi-omic atlas from per-modality h5ad files."""
        from pathlib import Path

        from sc_tools.errors import SCToolsUserError

        _check_deps(["mudata"])

        # Validate inputs/modalities match
        if len(inputs) != len(modalities):
            raise SCToolsUserError(
                f"Number of inputs ({len(inputs)}) must match number of modalities ({len(modalities)})",
                suggestion="Provide one --modality/-m flag per input file",
            )

        # Validate files exist
        for inp in inputs:
            if not Path(inp).exists():
                raise SCToolsUserError(
                    f"Input file not found: {inp}",
                    suggestion="Check the file path",
                )

        # Lazy imports per CLI-08
        import anndata as ad

        from sc_tools.assembly._atlas import MultiOmicAtlas

        # Load modalities
        adatas = {}
        for inp, mod_name in zip(inputs, modalities):
            adatas[mod_name] = ad.read_h5ad(inp)

        atlas = MultiOmicAtlas.from_modalities(adatas, subject_key=subject_key)
        atlas.save(output)

        cell_counts = {mod: adata.n_obs for mod, adata in adatas.items()}

        return CLIResult(
            status=Status.success,
            command="assemble build",
            data={
                "modalities": list(adatas.keys()),
                "cell_counts": cell_counts,
                "total_cells": sum(cell_counts.values()),
                "_input_files": inputs,
            },
            artifacts=[output],
            provenance=Provenance(command="assemble build"),
            message=f"Built atlas with {len(adatas)} modalities, {sum(cell_counts.values())} total cells",
        )

    # --- embed subcommand ---

    @assemble_app.command("embed")
    @cli_handler
    def assemble_embed(
        input: str = typer.Argument(..., help="Path to h5mu atlas file"),
        method: str = typer.Option("mofa", "--method", help="Embedding method (mofa, multivi, totalvi)"),
        n_factors: int = typer.Option(15, "--n-factors", help="Number of latent factors"),
        output: str = typer.Option(None, "--output", "-o", help="Output path (default: overwrite input)"),
        dry_run: bool = typer.Option(False, "--dry-run", help="Validate without running"),
        force: bool = typer.Option(False, "--force", help="Bypass memory guard"),
    ) -> None:
        """Run joint embedding on an assembled multi-omic atlas."""
        from pathlib import Path

        from sc_tools.errors import SCToolsUserError

        deps = ["mudata", "muon"]
        if method in ("multivi", "totalvi"):
            deps.append("scvi-tools")
        _check_deps(deps)

        if not Path(input).exists():
            raise SCToolsUserError(
                f"Input file not found: {input}",
                suggestion="Check the file path",
            )

        from sc_tools.assembly._atlas import MultiOmicAtlas

        atlas = MultiOmicAtlas.load(input)
        embedding = atlas.embed(method=method, n_factors=n_factors)

        out_path = output or input
        atlas.save(out_path)

        return CLIResult(
            status=Status.success,
            command="assemble embed",
            data={
                "method": method,
                "n_factors": n_factors,
                "embedding_shape": list(embedding.shape),
                "_input_files": [input],
            },
            artifacts=[out_path],
            provenance=Provenance(command="assemble embed"),
            message=f"Embedded atlas with {method} ({embedding.shape[0]} cells, {embedding.shape[1]} factors)",
        )

    # --- query subcommand ---

    @assemble_app.command("query")
    @cli_handler
    def assemble_query(
        input: str = typer.Argument(..., help="Path to h5mu atlas file"),
        query_type: str = typer.Option(
            "celltype_proportions", "--query-type", help="Query type (celltype_proportions)"
        ),
        celltype_key: str = typer.Option("celltype", "--celltype-key", help="Column for cell type annotations"),
        group_by: str = typer.Option("subject_id", "--group-by", help="Column to group by"),
        patient: str = typer.Option(None, "--patient", help="Filter to specific patient"),
        dry_run: bool = typer.Option(False, "--dry-run", help="Validate without querying"),
        force: bool = typer.Option(False, "--force", help="Bypass memory guard"),
    ) -> None:
        """Query a multi-omic atlas for cross-modal statistics."""
        from pathlib import Path

        from sc_tools.errors import SCToolsUserError

        _check_deps(["mudata"])

        if not Path(input).exists():
            raise SCToolsUserError(
                f"Input file not found: {input}",
                suggestion="Check the file path",
            )

        from sc_tools.assembly._atlas import MultiOmicAtlas

        atlas = MultiOmicAtlas.load(input)

        if patient:
            # Filter to specific patient first
            view_mdata = atlas.patient_view(patient)
            from sc_tools.assembly._atlas import MultiOmicAtlas as MA

            atlas = MA(view_mdata)

        props = atlas.celltype_proportions(celltype_key=celltype_key, group_by=group_by)

        return CLIResult(
            status=Status.success,
            command="assemble query",
            data={
                "query_type": query_type,
                "proportions": props.to_dict(orient="records"),
                "n_groups": props[group_by].nunique() if group_by in props.columns else 0,
                "_input_files": [input],
            },
            provenance=Provenance(command="assemble query"),
            message=f"Query returned {len(props)} rows",
        )

    app.add_typer(assemble_app, name="assemble")
