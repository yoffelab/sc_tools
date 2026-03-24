"""Differential expression commands (SCI-01).

Provides ``sct de run`` for pseudobulk DE via PyDESeq2.
"""

from __future__ import annotations

import logging

import typer

from sc_tools.cli import _check_deps, cli_handler
from sc_tools.io.gateway import DataTier

logger = logging.getLogger(__name__)

de_app = typer.Typer(help="Differential expression commands")


@de_app.callback(invoke_without_command=True)
def de_callback(ctx: typer.Context) -> None:
    """Differential expression commands. Run 'sct de --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


@de_app.command("run")
@cli_handler(tier=DataTier.T3_FULL)
def de_run(
    file: str = typer.Argument(..., help="Path to cell-typed h5ad"),
    condition: str = typer.Option(..., "--condition", "-c", help="Condition column in obs"),
    project_dir: str = typer.Option(".", "--project-dir", help="Project directory for output"),
    subject_col: str = typer.Option("subject_id", "--subject-col", help="Subject ID column"),
    celltype_col: str = typer.Option("celltype", "--celltype-col", help="Cell type column"),
    reference: str = typer.Option(None, "--reference", help="Reference level for contrast"),
    formula: str = typer.Option(None, "--formula", help="Custom design formula (e.g. '~ condition + sex')"),
    min_subjects: int = typer.Option(3, "--min-subjects", help="Min subjects per condition group"),
    min_cells: int = typer.Option(10, "--min-cells", help="Min cells per subject+celltype"),
    layer: str = typer.Option(None, "--layer", help="Count layer name override"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate inputs and report plan without executing"),
    force: bool = typer.Option(False, "--force", help="Override memory guard for large files"),
) -> None:
    """Run pseudobulk differential expression per celltype (SCI-01).

    Aggregates raw counts by subject + celltype, runs PyDESeq2 with automatic
    design formula inference (including batch covariate detection with collinearity guard).
    Writes per-celltype CSV files to {project_dir}/results/de/.
    """
    from pathlib import Path

    from sc_tools.errors import SCToolsDataError
    from sc_tools.models.result import CLIResult, Provenance, Status

    _check_deps(["pydeseq2", "scanpy", "anndata"])

    import scanpy as sc

    from sc_tools.tl.de import run_pseudobulk_de

    path = Path(file)
    if not path.exists():
        from sc_tools.errors import SCToolsUserError

        raise SCToolsUserError(f"File not found: {path}", suggestion="Check the file path")

    logger.info("Loading %s", path.name)
    adata = sc.read_h5ad(path)
    logger.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Validate required columns exist
    if subject_col not in adata.obs.columns:
        raise SCToolsDataError(
            f"Subject column '{subject_col}' not found in obs. "
            f"Available columns: {list(adata.obs.columns)}",
            suggestion=f"Set --subject-col to a valid column name",
        )
    if condition not in adata.obs.columns:
        raise SCToolsDataError(
            f"Condition column '{condition}' not found in obs. "
            f"Available columns: {list(adata.obs.columns)}",
            suggestion=f"Set --condition to a valid column name",
        )

    # Run pseudobulk DE
    results = run_pseudobulk_de(
        adata,
        condition_key=condition,
        subject_key=subject_col,
        celltype_key=celltype_col,
        formula=formula,
        reference=reference,
        min_subjects_per_group=min_subjects,
        min_cells_per_combo=min_cells,
        layer=layer,
    )

    # Write per-celltype CSVs
    output_dir = Path(project_dir) / "results" / "de"
    output_dir.mkdir(parents=True, exist_ok=True)

    artifacts: list[str] = []
    total_de_genes = 0
    n_celltypes_tested = len(results)

    # Count how many celltypes were in the data but not tested
    all_celltypes = adata.obs[celltype_col].nunique()
    n_celltypes_skipped = all_celltypes - n_celltypes_tested

    for ct, df in results.items():
        # Sanitize celltype name for filename
        safe_ct = ct.replace("/", "_").replace(" ", "_")
        csv_path = output_dir / f"{safe_ct}.csv"
        df.to_csv(csv_path, index=False)
        artifacts.append(str(csv_path))

        if "padj" in df.columns:
            total_de_genes += int((df["padj"] < 0.05).sum())

    logger.info(
        "DE complete: %d celltypes tested, %d skipped, %d total DE genes (padj<0.05)",
        n_celltypes_tested, n_celltypes_skipped, total_de_genes,
    )

    return CLIResult(
        status=Status.success,
        command="de run",
        data={
            "n_celltypes_tested": n_celltypes_tested,
            "n_celltypes_skipped": n_celltypes_skipped,
            "total_de_genes": total_de_genes,
            "_input_files": [str(file)],
        },
        artifacts=artifacts,
        provenance=Provenance(command="de run"),
        message=(
            f"Pseudobulk DE complete: {n_celltypes_tested} celltypes tested, "
            f"{n_celltypes_skipped} skipped, {total_de_genes} DE genes (padj<0.05). "
            f"Results in {output_dir}"
        ),
    )
