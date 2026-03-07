"""CLI entry point for IMC segmentation benchmark.

Usage: ``python -m sc_tools.bm.cli <subcommand> [options]``

Subcommands:
- ``catalog``: Build ROI catalog from HPC directories
- ``download``: Download public datasets
- ``prepare``: Validate and prepare ROIs
- ``run``: Run full benchmark
- ``run-chunk``: Run a single chunk (for SLURM)
- ``train-vit``: Train SegFormer model (Strategy 4)
- ``evaluate``: Generate reports from results
- ``slurm``: Generate and submit SLURM jobs
"""

from __future__ import annotations

import argparse
import logging
import sys

logger = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> int:
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="sc_tools.bm.cli",
        description="IMC segmentation benchmark pipeline",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging")
    subparsers = parser.add_subparsers(dest="command", help="Subcommand")

    # catalog
    p_catalog = subparsers.add_parser("catalog", help="Build ROI catalog")
    p_catalog.add_argument("--base-dir", required=True, help="HPC base directory")
    p_catalog.add_argument("--datasets", nargs="*", help="Specific datasets")
    p_catalog.add_argument("--output", default="catalog.csv", help="Output CSV")
    p_catalog.add_argument("--include-public", action="store_true")
    p_catalog.add_argument("--public-dir", help="Public datasets directory")

    # download
    p_download = subparsers.add_parser("download", help="Download public datasets")
    p_download.add_argument("--dataset", required=True, help="Dataset name")
    p_download.add_argument("--output-dir", default="public_data", help="Output directory")
    p_download.add_argument("--force", action="store_true")

    # prepare
    p_prepare = subparsers.add_parser("prepare", help="Validate and prepare ROIs")
    p_prepare.add_argument("--catalog", required=True, help="Catalog CSV")
    p_prepare.add_argument("--output", default="validation_report.csv")

    # run
    p_run = subparsers.add_parser("run", help="Run full benchmark")
    p_run.add_argument("--catalog", required=True, help="Catalog CSV")
    p_run.add_argument("--config", required=True, help="Config YAML")
    p_run.add_argument("--output-dir", default="benchmark_results")

    # run-chunk
    p_chunk = subparsers.add_parser("run-chunk", help="Run a single SLURM chunk")
    p_chunk.add_argument("--catalog", required=True)
    p_chunk.add_argument("--config", required=True)
    p_chunk.add_argument("--chunk-id", type=int, required=True)
    p_chunk.add_argument("--chunk-size", type=int, default=20)
    p_chunk.add_argument("--output-dir", required=True)
    p_chunk.add_argument("--strategies", nargs="*", type=int)

    # train-vit
    p_vit = subparsers.add_parser("train-vit", help="Train SegFormer model")
    p_vit.add_argument("--catalog", required=True)
    p_vit.add_argument("--config", required=True)
    p_vit.add_argument("--output-dir", default="segformer_model")

    # evaluate
    p_eval = subparsers.add_parser("evaluate", help="Generate reports")
    p_eval.add_argument("--results", required=True, help="Results CSV")
    p_eval.add_argument("--output", default="benchmark_report.html")
    p_eval.add_argument("--title", default="IMC Segmentation Benchmark Report")

    # slurm
    p_slurm = subparsers.add_parser("slurm", help="Generate SLURM jobs")
    p_slurm.add_argument("--catalog", required=True)
    p_slurm.add_argument("--config", required=True)
    p_slurm.add_argument("--output-dir", default="benchmark_slurm")
    p_slurm.add_argument("--submit", action="store_true")
    p_slurm.add_argument("--dry-run", action="store_true")

    args = parser.parse_args(argv)

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    if args.command is None:
        parser.print_help()
        return 1

    return _dispatch(args)


def _dispatch(args: argparse.Namespace) -> int:
    """Dispatch to subcommand handler."""
    try:
        if args.command == "catalog":
            return _cmd_catalog(args)
        elif args.command == "download":
            return _cmd_download(args)
        elif args.command == "prepare":
            return _cmd_prepare(args)
        elif args.command == "run":
            return _cmd_run(args)
        elif args.command == "run-chunk":
            return _cmd_run_chunk(args)
        elif args.command == "train-vit":
            return _cmd_train_vit(args)
        elif args.command == "evaluate":
            return _cmd_evaluate(args)
        elif args.command == "slurm":
            return _cmd_slurm(args)
        else:
            logger.error("Unknown command: %s", args.command)
            return 1
    except Exception as e:
        logger.error("Error: %s", e, exc_info=True)
        return 1


def _cmd_catalog(args: argparse.Namespace) -> int:
    from sc_tools.data.imc.benchmark.catalog import build_benchmark_catalog

    df = build_benchmark_catalog(
        base_dir=args.base_dir,
        datasets=args.datasets,
        public_dir=args.public_dir,
        include_public=args.include_public,
    )
    df.to_csv(args.output, index=False)
    print(f"Catalog saved: {len(df)} ROIs -> {args.output}")
    return 0


def _cmd_download(args: argparse.Namespace) -> int:
    from sc_tools.data.imc.benchmark.public import download_public_dataset

    ds_dir = download_public_dataset(args.dataset, args.output_dir, force=args.force)
    print(f"Downloaded to: {ds_dir}")
    return 0


def _cmd_prepare(args: argparse.Namespace) -> int:
    import pandas as pd

    from sc_tools.data.imc.benchmark.prepare import validate_roi_files

    catalog = pd.read_csv(args.catalog)
    results = []
    for _, row in catalog.iterrows():
        val = validate_roi_files(row["tiff_path"], row.get("channel_csv_path"))
        val["roi_id"] = row["roi_id"]
        val["dataset"] = row["dataset"]
        results.append(val)

    report = pd.DataFrame(results)
    report.to_csv(args.output, index=False)
    n_valid = report["valid"].sum()
    print(f"Validation: {n_valid}/{len(report)} valid ROIs -> {args.output}")
    return 0


def _cmd_run(args: argparse.Namespace) -> int:
    import pandas as pd

    from sc_tools.bm.runner import run_benchmark
    from sc_tools.data.imc.benchmark.config import BenchmarkConfig

    catalog = pd.read_csv(args.catalog)
    config = BenchmarkConfig.from_yaml(args.config)
    config.output_dir = args.output_dir

    df = run_benchmark(catalog, config, output_dir=args.output_dir)
    print(f"Benchmark complete: {len(df)} results")
    return 0


def _cmd_run_chunk(args: argparse.Namespace) -> int:
    import pandas as pd

    from sc_tools.bm.runner import run_benchmark
    from sc_tools.data.imc.benchmark.config import BenchmarkConfig

    catalog = pd.read_csv(args.catalog)
    config = BenchmarkConfig.from_yaml(args.config)

    # Slice catalog for this chunk
    start = args.chunk_id * args.chunk_size
    end = min(start + args.chunk_size, len(catalog))
    chunk = catalog.iloc[start:end].reset_index(drop=True)

    if args.strategies:
        config.strategies = args.strategies

    chunk_output = f"{args.output_dir}/chunk_{args.chunk_id}"
    df = run_benchmark(chunk, config, output_dir=chunk_output)
    print(f"Chunk {args.chunk_id}: {len(df)} results")
    return 0


def _cmd_train_vit(args: argparse.Namespace) -> int:
    import pandas as pd

    from sc_tools.data.imc.benchmark.config import BenchmarkConfig

    catalog = pd.read_csv(args.catalog)
    config = BenchmarkConfig.from_yaml(args.config)

    from sc_tools.bm.strategy_vit import train_segformer

    train_segformer(catalog, config, output_dir=args.output_dir)
    print(f"SegFormer training complete -> {args.output_dir}")
    return 0


def _cmd_evaluate(args: argparse.Namespace) -> int:
    import pandas as pd

    from sc_tools.bm.report import generate_benchmark_report
    from sc_tools.bm.runner import aggregate_results

    results = pd.read_csv(args.results)
    aggregated = aggregate_results(results)

    generate_benchmark_report(
        results_df=results,
        aggregated=aggregated,
        output_path=args.output,
        title=args.title,
    )
    print(f"Report saved: {args.output}")
    return 0


def _cmd_slurm(args: argparse.Namespace) -> int:
    import pandas as pd

    from sc_tools.bm.slurm import generate_benchmark_sbatch, submit_benchmark_jobs
    from sc_tools.data.imc.benchmark.config import BenchmarkConfig

    catalog = pd.read_csv(args.catalog)
    config = BenchmarkConfig.from_yaml(args.config)

    scripts = generate_benchmark_sbatch(catalog, config, output_dir=args.output_dir)
    print(f"Generated {len(scripts)} sbatch scripts")

    if args.submit:
        job_ids = submit_benchmark_jobs(scripts, dry_run=args.dry_run)
        print(f"Submitted {len(job_ids)} jobs: {job_ids}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
