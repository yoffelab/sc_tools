"""Preprocessing commands (CMD-02)."""
from __future__ import annotations

import logging

import typer

from sc_tools.cli import _check_deps, cli_handler
from sc_tools.io.gateway import DataTier

logger = logging.getLogger(__name__)

preprocess_app = typer.Typer(help="Preprocessing commands")


@preprocess_app.callback(invoke_without_command=True)
def preprocess_callback(ctx: typer.Context) -> None:
    """Preprocessing commands. Run 'sct preprocess --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


def _detect_modality(adata) -> str:
    """Auto-detect modality from adata.uns or panel size (D-08).

    Checks adata.uns for explicit modality keys, falls back to panel size heuristic.
    Raises SCToolsDataError if modality cannot be determined.
    """
    from sc_tools.pp.recipes import VALID_MODALITIES
    from sc_tools.errors import SCToolsDataError

    # Check adata.uns for explicit modality hint
    for key in ("modality", "platform", "data_type"):
        if key in adata.uns:
            val = str(adata.uns[key]).lower()
            if val in VALID_MODALITIES:
                logger.info("Auto-detected modality from adata.uns['%s']: %s", key, val)
                return val

    # Heuristic: small panel -> targeted (xenium/cosmx)
    if adata.n_vars < 1000:
        logger.info("Auto-detected targeted panel (n_vars=%d < 1000), using 'xenium'", adata.n_vars)
        return "xenium"

    # Cannot determine -- per D-08, raise SCToolsDataError
    raise SCToolsDataError(
        f"Cannot auto-detect modality (n_vars={adata.n_vars}, no modality key in adata.uns)",
        suggestion="Specify --modality visium|visium_hd|xenium|cosmx|imc",
    )


@preprocess_app.command("run")
@cli_handler(tier=DataTier.T3_FULL)
def preprocess_run(
    file: str = typer.Argument(..., help="Path to input h5ad file (raw counts)"),
    project_dir: str = typer.Option(".", "--project-dir", help="Project directory (per D-03)"),
    output: str = typer.Option(None, "--output", "-o", help="Output h5ad path (default: {project_dir}/results/adata.normalized.h5ad)"),
    modality: str = typer.Option("auto", "--modality", "-m", help="Data modality: auto, visium, visium_hd, xenium, cosmx, imc (D-08)"),
    batch_key: str = typer.Option("library_id", "--batch-key", help="Batch column in obs"),
    integration: str = typer.Option("scvi", "--integration", "-i", help="Integration method: scvi, harmony, cytovi, none"),
    n_top_genes: int = typer.Option(2000, "--n-top-genes", help="Number of HVGs"),
    resolution: float = typer.Option(0.8, "--resolution", help="Leiden clustering resolution"),
    use_gpu: str = typer.Option("auto", "--use-gpu", help="GPU usage: auto, true, false"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate inputs and report plan without executing"),
    force: bool = typer.Option(False, "--force", help="Override memory guard for large files"),
):  # returns CLIResult
    """Run modality-aware preprocessing on an h5ad file (CMD-02, D-06, D-08, D-13)."""
    from sc_tools.models.result import CLIResult, Provenance, Status
    from sc_tools.errors import SCToolsUserError
    from sc_tools.pp.recipes import VALID_MODALITIES, VALID_INTEGRATIONS

    _check_deps(["scanpy", "anndata"])

    from pathlib import Path
    import scanpy as sc

    path = Path(file)
    if not path.exists():
        raise SCToolsUserError(f"File not found: {path}", suggestion="Check the file path")

    if integration not in VALID_INTEGRATIONS:
        raise SCToolsUserError(
            f"Unknown integration method '{integration}'",
            suggestion=f"Use one of: {', '.join(sorted(VALID_INTEGRATIONS))}",
        )

    # Check integration-specific deps
    if integration == "scvi":
        _check_deps(["scvi"])
    elif integration == "harmony":
        _check_deps(["scanpy"])

    logger.info("Loading %s", path.name)
    adata = sc.read_h5ad(path)
    logger.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)

    # D-08: auto-detect modality
    resolved_modality = modality
    if modality == "auto":
        resolved_modality = _detect_modality(adata)
    elif modality not in VALID_MODALITIES:
        raise SCToolsUserError(
            f"Unknown modality '{modality}'",
            suggestion=f"Use one of: auto, {', '.join(sorted(VALID_MODALITIES))}",
        )

    # Convert use_gpu string to appropriate type
    gpu_val: str | bool = use_gpu
    if use_gpu.lower() == "true":
        gpu_val = True
    elif use_gpu.lower() == "false":
        gpu_val = False

    # D-13: communicate what will happen
    logger.info(
        "Preprocessing: modality=%s, integration=%s, batch_key=%s, n_top_genes=%d, resolution=%.2f",
        resolved_modality, integration, batch_key, n_top_genes, resolution,
    )

    from sc_tools.pp.recipes import preprocess

    adata = preprocess(
        adata,
        modality=resolved_modality,
        batch_key=batch_key,
        integration=integration,
        n_top_genes=n_top_genes,
        resolution=resolution,
        use_gpu=gpu_val,
    )

    # Determine output path
    proj = Path(project_dir)
    if output is None:
        output_path = proj / "results" / "adata.normalized.h5ad"
    else:
        output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Writing output to %s", output_path)
    adata.write_h5ad(output_path)

    return CLIResult(
        status=Status.success,
        command="preprocess run",
        data={
            "input": str(path),
            "output": str(output_path),
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "modality": resolved_modality,
            "modality_auto_detected": modality == "auto",
            "integration": integration,
            "batch_key": batch_key,
            "n_top_genes": n_top_genes,
            "resolution": resolution,
        },
        artifacts=[str(output_path)],
        provenance=Provenance(command="preprocess run"),
        message=f"Preprocessing complete: {adata.n_obs} cells, modality={resolved_modality}, integration={integration} -> {output_path}",
    )
