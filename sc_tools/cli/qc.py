"""QC and report commands (CMD-01, CMD-06)."""

from __future__ import annotations

import logging

import typer

from sc_tools.cli import _check_deps, cli_handler
from sc_tools.io.gateway import DataTier

logger = logging.getLogger(__name__)

qc_app = typer.Typer(help="Quality control commands")


@qc_app.callback(invoke_without_command=True)
def qc_callback(ctx: typer.Context) -> None:
    """Quality control commands. Run 'sct qc --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


@qc_app.command("run")
@cli_handler(tier=DataTier.T3_FULL)
def qc_run(
    file: str = typer.Argument(..., help="Path to h5ad checkpoint file"),
    project_dir: str = typer.Option(".", "--project-dir", help="Project directory"),
    modality: str = typer.Option("visium", "--modality", "-m", help="Data modality (visium, visium_hd, xenium, cosmx, imc)"),
    sample_col: str = typer.Option("library_id", "--sample-col", help="Sample column in obs"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate inputs and report plan without executing"),
    force: bool = typer.Option(False, "--force", help="Override memory guard for large files"),
) -> None:
    """Run QC metrics on a checkpoint file, output JSON summary (CMD-01)."""
    from sc_tools.models.result import CLIResult, Provenance, Status

    _check_deps(["scanpy", "anndata"])

    from pathlib import Path

    import scanpy as sc

    from sc_tools.errors import SCToolsUserError
    from sc_tools.qc.sample_qc import classify_samples, compute_sample_metrics

    path = Path(file)
    if not path.exists():
        raise SCToolsUserError(f"File not found: {path}", suggestion="Check the file path")

    logger.info("Loading %s (%s)", path.name, modality)
    adata = sc.read_h5ad(path)
    logger.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Ensure QC metrics exist (use sc_tools wrapper that creates mt/hb boolean columns first)
    if "n_genes_by_counts" not in adata.obs.columns:
        from sc_tools.qc.metrics import calculate_qc_metrics
        calculate_qc_metrics(adata, modality=modality, inplace=True, log1p=False)

    metrics = compute_sample_metrics(adata, sample_col=sample_col, modality=modality)
    classified = classify_samples(metrics, modality=modality)

    n_samples = len(metrics)
    n_pass = int(classified["qc_pass"].sum()) if "qc_pass" in classified.columns else n_samples
    n_fail = n_samples - n_pass

    # Build JSON-serializable metrics summary
    metrics_dict = metrics.to_dict(orient="index")
    classified_dict = classified.to_dict(orient="index")

    return CLIResult(
        status=Status.success,
        command="qc run",
        data={
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_samples": n_samples,
            "n_pass": n_pass,
            "n_fail": n_fail,
            "modality": modality,
            "sample_col": sample_col,
            "metrics": metrics_dict,
            "classified": classified_dict,
        },
        provenance=Provenance(command="qc run"),
        message=f"QC metrics computed: {adata.n_obs} cells, {n_samples} samples ({n_pass} pass, {n_fail} fail)",
    )


# ---------------------------------------------------------------------------
# sct report <type> (CMD-06)
# ---------------------------------------------------------------------------

_VALID_REPORT_TYPES = ("pre_filter", "post_filter", "post_integration", "post_celltyping")

report_app = typer.Typer(help="Report generation commands")


@report_app.callback(invoke_without_command=True)
def report_callback(ctx: typer.Context) -> None:
    """Report generation commands. Run 'sct report --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


@report_app.command("generate")
@cli_handler
def report_generate(
    report_type: str = typer.Argument(..., help="Report type: pre_filter, post_filter, post_integration, post_celltyping"),
    adata_path: str = typer.Option(..., "--adata", help="Path to primary h5ad file"),
    project_dir: str = typer.Option(".", "--project-dir", help="Project directory"),
    adata_post_path: str = typer.Option(None, "--adata-post", help="Path to post-filter h5ad (for post_filter report)"),
    modality: str = typer.Option("visium", "--modality", "-m", help="Data modality"),
    sample_col: str = typer.Option("library_id", "--sample-col", help="Sample column in obs"),
    batch_key: str = typer.Option("sample", "--batch-key", help="Batch key (for post_integration)"),
) -> None:
    """Generate an HTML QC report (CMD-06, D-04, D-05)."""
    from sc_tools.errors import SCToolsUserError
    from sc_tools.models.result import CLIResult, Provenance, Status

    if report_type not in _VALID_REPORT_TYPES:
        raise SCToolsUserError(
            f"Unknown report type '{report_type}'. Must be one of: {', '.join(_VALID_REPORT_TYPES)}",
            suggestion=f"Use one of: {', '.join(_VALID_REPORT_TYPES)}",
        )

    _check_deps(["scanpy", "anndata", "matplotlib"])

    from datetime import datetime
    from pathlib import Path

    import scanpy as sc

    from sc_tools.qc.report import (
        generate_post_celltyping_report,
        generate_post_filter_report,
        generate_post_integration_report,
        generate_pre_filter_report,
    )
    from sc_tools.qc.sample_qc import classify_samples, compute_sample_metrics

    adata_file = Path(adata_path)
    if not adata_file.exists():
        raise SCToolsUserError(f"File not found: {adata_file}", suggestion="Check --adata path")

    proj = Path(project_dir)
    # D-04: output path convention
    date_stamp = datetime.now().strftime("%y%m%d")
    output_dir = proj / "figures" / "reports"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Loading %s for %s report", adata_file.name, report_type)
    adata = sc.read_h5ad(adata_file)
    logger.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Ensure QC metrics for pre/post filter
    if report_type in ("pre_filter", "post_filter"):
        if "n_genes_by_counts" not in adata.obs.columns:
            sc.pp.calculate_qc_metrics(adata, qc_vars=["MT-", "mt-"], percent_top=None, log1p=False, inplace=True)
        metrics = compute_sample_metrics(adata, sample_col=sample_col, modality=modality)
        classified = classify_samples(metrics, modality=modality)

    report_path: Path
    if report_type == "pre_filter":
        report_path = generate_pre_filter_report(
            adata, metrics, classified, output_dir,
            sample_col=sample_col, modality=modality,
            date_stamp=date_stamp,
        )
    elif report_type == "post_filter":
        if adata_post_path is None:
            raise SCToolsUserError(
                "post_filter report requires --adata-post (post-filter h5ad)",
                suggestion="Provide --adata-post with the filtered AnnData path",
            )
        adata_post_file = Path(adata_post_path)
        if not adata_post_file.exists():
            raise SCToolsUserError(f"File not found: {adata_post_file}", suggestion="Check --adata-post path")
        adata_post = sc.read_h5ad(adata_post_file)
        report_path = generate_post_filter_report(
            adata, adata_post, metrics, classified, output_dir,
            sample_col=sample_col, modality=modality,
            date_stamp=date_stamp,
        )
    elif report_type == "post_integration":
        report_path = generate_post_integration_report(
            adata, output_dir,
            batch_key=batch_key,
            date_stamp=date_stamp,
        )
    elif report_type == "post_celltyping":
        report_path = generate_post_celltyping_report(
            adata, output_dir,
            date_stamp=date_stamp,
        )

    return CLIResult(
        status=Status.success,
        command=f"report {report_type}",
        data={"report_type": report_type, "report_path": str(report_path), "n_cells": int(adata.n_obs)},
        artifacts=[str(report_path)],
        provenance=Provenance(command=f"report {report_type}"),
        message=f"{report_type} report generated: {report_path}",
    )
