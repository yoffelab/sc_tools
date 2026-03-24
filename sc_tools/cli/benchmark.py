"""Benchmarking commands (CMD-04)."""
from __future__ import annotations

import logging

import typer

from sc_tools.cli import _check_deps, cli_handler
from sc_tools.io.gateway import DataTier

logger = logging.getLogger(__name__)

benchmark_app = typer.Typer(help="Benchmarking commands")


@benchmark_app.callback(invoke_without_command=True)
def benchmark_callback(ctx: typer.Context) -> None:
    """Benchmarking commands. Run 'sct benchmark --help' for available subcommands."""
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())


@benchmark_app.command("integration")
@cli_handler(tier=DataTier.T2_SUMMARY)
def benchmark_integration(
    from_dir: str = typer.Option(..., "--from-dir", help="Directory with per-method h5ad embedding files"),
    project_dir: str = typer.Option(".", "--project-dir", help="Project directory (per D-03)"),
    batch_key: str = typer.Option("sample", "--batch-key", help="Batch column in obs"),
    celltype_key: str = typer.Option(None, "--celltype-key", help="Cell type column in obs"),
    bio_key: str = typer.Option(None, "--bio-key", help="Bio conservation column (defaults to celltype_key)"),
    subsample_n: int = typer.Option(500_000, "--subsample-n", help="Subsample to N cells (D-07, default 500000)"),
    seed: int = typer.Option(0, "--seed", help="Random seed for subsampling (D-07, default 0)"),
    batch_weight: float = typer.Option(0.4, "--batch-weight", help="Weight for batch removal score"),
    bio_weight: float = typer.Option(0.6, "--bio-weight", help="Weight for bio conservation score"),
    resolution: float = typer.Option(1.0, "--resolution", help="Leiden resolution for ARI/NMI"),
    generate_report: bool = typer.Option(False, "--report", help="Generate HTML benchmark report (D-12)"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate inputs and report plan without executing"),
    force: bool = typer.Option(False, "--force", help="Override memory guard for large files"),
):  # returns CLIResult
    """Compare integration methods from pre-computed h5ad embeddings (CMD-04, D-07, D-12, D-13).

    Discovers h5ad files in --from-dir, loads embeddings via h5py, computes
    comparison metrics, and outputs ranked JSON. Subsampling is applied
    automatically when dataset exceeds --subsample-n cells.

    With --report, generates an HTML benchmark report that includes subsampling
    parameters per D-12 item 3.
    """
    from sc_tools.models.result import CLIResult, Provenance, Status
    from sc_tools.errors import SCToolsUserError

    _check_deps(["h5py", "scanpy"])

    from pathlib import Path

    from_path = Path(from_dir)
    if not from_path.is_dir():
        raise SCToolsUserError(f"Directory not found: {from_path}", suggestion="Check --from-dir path")

    # Discover h5ad files in the directory
    h5ad_files = sorted(from_path.glob("*.h5ad"))
    if not h5ad_files:
        raise SCToolsUserError(
            f"No .h5ad files found in {from_path}",
            suggestion="Ensure --from-dir contains per-method h5ad embedding files",
        )

    # Build embedding_files dict: method name = stem of filename
    embedding_files: dict[str, str] = {}
    for f in h5ad_files:
        method_name = f.stem
        embedding_files[method_name] = str(f)
    logger.info("Discovered %d methods: %s", len(embedding_files), list(embedding_files.keys()))

    # Load one file to check cell count for subsampling decision
    import h5py
    with h5py.File(h5ad_files[0], "r") as hf:
        n_cells = hf["obs"].attrs.get("_index_length", 0)
        if n_cells == 0 and "obs" in hf:
            # Fallback: try shape
            for key in hf["obs"]:
                n_cells = len(hf[f"obs/{key}"])
                break

    # D-07, D-12: subsampling decision and transparency
    actual_subsample_n = None
    subsampled = False
    if n_cells > subsample_n:
        actual_subsample_n = subsample_n
        subsampled = True
        logger.info("Subsampling from %d to %d cells (seed=%d)", n_cells, subsample_n, seed)
    else:
        logger.info("Dataset has %d cells (below subsample threshold %d), no subsampling", n_cells, subsample_n)

    from sc_tools.bm.integration import compare_integrations

    logger.info(
        "Running benchmark: batch_key=%s, celltype_key=%s, batch_weight=%.2f, bio_weight=%.2f, seed=%d",
        batch_key, celltype_key, batch_weight, bio_weight, seed,
    )

    results_df = compare_integrations(
        embedding_files=embedding_files,
        batch_key=batch_key,
        celltype_key=celltype_key,
        bio_key=bio_key,
        batch_weight=batch_weight,
        bio_weight=bio_weight,
        subsample_n=actual_subsample_n,
        seed=seed,
        resolution=resolution,
    )

    # Convert results to JSON-serializable dict
    results_dict = results_df.reset_index().to_dict(orient="records")

    artifacts = []

    # D-12: record subsampling in provenance data
    provenance_data = {
        "from_dir": str(from_path),
        "methods": list(embedding_files.keys()),
        "n_cells_original": int(n_cells),
        "subsampled": subsampled,
        "subsample_n": int(subsample_n),
        "seed": int(seed),
        "batch_key": batch_key,
        "celltype_key": celltype_key,
        "bio_key": bio_key,
        "batch_weight": float(batch_weight),
        "bio_weight": float(bio_weight),
        "resolution": float(resolution),
    }

    # D-12 item 3: generate HTML benchmark report with subsampling params
    if generate_report:
        _check_deps(["matplotlib"])
        from sc_tools.bm.report import generate_benchmark_report

        proj = Path(project_dir)
        report_output_dir = proj / "figures" / "reports"
        report_output_dir.mkdir(parents=True, exist_ok=True)

        # Pass subsampling params so they appear in the report header
        benchmark_params = {
            "n_cells_original": int(n_cells),
            "subsampled": subsampled,
            "subsample_n": int(subsample_n),
            "seed": int(seed),
            "batch_key": batch_key,
            "celltype_key": celltype_key,
            "batch_weight": float(batch_weight),
            "bio_weight": float(bio_weight),
            "resolution": float(resolution),
            "methods": list(embedding_files.keys()),
        }

        report_path = generate_benchmark_report(
            results_df,
            output_path=report_output_dir / "integration_benchmark_report.html",
            benchmark_params=benchmark_params,
        )
        artifacts.append(str(report_path))
        logger.info("Benchmark report generated: %s", report_path)

    # D-13: clear message about what was done
    top_method = results_dict[0]["method"] if results_dict and "method" in results_dict[0] else (results_dict[0].get("index", "?") if results_dict else "?")
    msg = f"Benchmark complete: {len(embedding_files)} methods compared"
    if subsampled:
        msg += f" (subsampled {n_cells} -> {subsample_n} cells, seed={seed})"
    msg += f". Top method: {top_method}"
    if generate_report and artifacts:
        msg += f". Report: {artifacts[0]}"

    return CLIResult(
        status=Status.success,
        command="benchmark integration",
        data={
            "results": results_dict,
            "n_methods": len(embedding_files),
            **provenance_data,
        },
        artifacts=artifacts,
        provenance=Provenance(command="benchmark integration"),
        message=msg,
    )
