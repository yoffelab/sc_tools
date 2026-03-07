"""SLURM job generation for distributed IMC segmentation benchmark.

Splits the ROI catalog into chunks and generates sbatch scripts
(GPU for strategies 3-4, CPU for 1-2).
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

import pandas as pd

__all__ = [
    "generate_benchmark_sbatch",
    "submit_benchmark_jobs",
    "collect_benchmark_results",
]

logger = logging.getLogger(__name__)

_SBATCH_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --time={time_limit}
#SBATCH --mem={mem}
#SBATCH --cpus-per-task={cpus}
{gpu_line}
#SBATCH --output={log_dir}/benchmark_%j.out
#SBATCH --error={log_dir}/benchmark_%j.err

# Activate environment
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate {conda_env}

# Run benchmark chunk
python -m sc_tools.bm.cli run-chunk \\
    --catalog {catalog_path} \\
    --config {config_path} \\
    --chunk-id {chunk_id} \\
    --chunk-size {chunk_size} \\
    --output-dir {output_dir} \\
    --strategies {strategies_str}
"""


def generate_benchmark_sbatch(
    catalog: pd.DataFrame,
    config,
    output_dir: str | Path = "benchmark_slurm",
    conda_env: str = "sc_tools_benchmark",
) -> list[Path]:
    """Generate sbatch scripts for distributed benchmark.

    Parameters
    ----------
    catalog
        ROI catalog DataFrame.
    config
        BenchmarkConfig instance.
    output_dir
        Where to write sbatch scripts.
    conda_env
        Conda environment name.

    Returns
    -------
    List of paths to generated sbatch scripts.
    """
    output_dir = Path(output_dir)
    scripts_dir = output_dir / "scripts"
    log_dir = output_dir / "logs"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    # Save catalog and config for reference
    catalog_path = output_dir / "catalog.csv"
    config_path = output_dir / "config.yaml"
    catalog.to_csv(catalog_path, index=False)
    config.to_yaml(config_path)

    # Split into chunks
    n_rois = len(catalog)
    chunk_size = config.slurm_chunk_size
    n_chunks = max(1, (n_rois + chunk_size - 1) // chunk_size)

    # Separate CPU and GPU strategies
    cpu_strategies = [s for s in config.strategies if s in (1, 2)]
    gpu_strategies = [s for s in config.strategies if s in (3, 4)]

    scripts = []

    for chunk_id in range(n_chunks):
        # CPU jobs
        if cpu_strategies:
            script_path = scripts_dir / f"benchmark_cpu_chunk{chunk_id}.sh"
            content = _SBATCH_TEMPLATE.format(
                job_name=f"bm_cpu_{chunk_id}",
                partition=config.slurm_partition_cpu,
                time_limit=config.slurm_time_cpu,
                mem=config.slurm_mem,
                cpus=4,
                gpu_line="",
                log_dir=log_dir,
                conda_env=conda_env,
                catalog_path=catalog_path,
                config_path=config_path,
                chunk_id=chunk_id,
                chunk_size=chunk_size,
                output_dir=output_dir / "results",
                strategies_str=" ".join(str(s) for s in cpu_strategies),
            )
            script_path.write_text(content)
            scripts.append(script_path)

        # GPU jobs
        if gpu_strategies:
            script_path = scripts_dir / f"benchmark_gpu_chunk{chunk_id}.sh"
            content = _SBATCH_TEMPLATE.format(
                job_name=f"bm_gpu_{chunk_id}",
                partition=config.slurm_partition_gpu,
                time_limit=config.slurm_time_gpu,
                mem=config.slurm_mem,
                cpus=4,
                gpu_line="#SBATCH --gres=gpu:1",
                log_dir=log_dir,
                conda_env=conda_env,
                catalog_path=catalog_path,
                config_path=config_path,
                chunk_id=chunk_id,
                chunk_size=chunk_size,
                output_dir=output_dir / "results",
                strategies_str=" ".join(str(s) for s in gpu_strategies),
            )
            script_path.write_text(content)
            scripts.append(script_path)

    logger.info(
        "Generated %d sbatch scripts (%d chunks x %s strategy groups)",
        len(scripts),
        n_chunks,
        "CPU+GPU" if cpu_strategies and gpu_strategies else "single",
    )
    return scripts


def submit_benchmark_jobs(
    scripts: list[Path],
    dry_run: bool = False,
) -> list[str]:
    """Submit sbatch scripts to SLURM.

    Parameters
    ----------
    scripts
        Paths to sbatch scripts.
    dry_run
        If True, print commands without submitting.

    Returns
    -------
    List of SLURM job IDs.
    """
    job_ids = []
    for script in scripts:
        cmd = f"sbatch {script}"
        if dry_run:
            logger.info("[dry-run] %s", cmd)
            job_ids.append("DRY_RUN")
        else:
            result = subprocess.run(
                ["sbatch", str(script)],
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode == 0:
                # "Submitted batch job 12345"
                parts = result.stdout.strip().split()
                job_id = parts[-1] if parts else "unknown"
                job_ids.append(job_id)
                logger.info("Submitted %s -> job %s", script.name, job_id)
            else:
                logger.error("Failed to submit %s: %s", script.name, result.stderr)

    return job_ids


def collect_benchmark_results(
    output_dir: str | Path,
) -> pd.DataFrame:
    """Collect and merge results from all SLURM chunks.

    Parameters
    ----------
    output_dir
        Directory containing chunk results.

    Returns
    -------
    Merged DataFrame of all results.
    """
    output_dir = Path(output_dir)
    results_dir = output_dir / "results"

    if not results_dir.is_dir():
        logger.warning("No results directory found at %s", results_dir)
        return pd.DataFrame()

    csv_files = sorted(results_dir.glob("**/benchmark_results.csv"))
    if not csv_files:
        logger.warning("No result CSV files found")
        return pd.DataFrame()

    dfs = [pd.read_csv(f) for f in csv_files]
    merged = pd.concat(dfs, ignore_index=True)

    # Deduplicate (in case of overlapping resume runs)
    dedup_cols = ["dataset", "roi_id", "strategy", "method"]
    available = [c for c in dedup_cols if c in merged.columns]
    if available:
        merged = merged.drop_duplicates(subset=available, keep="last")

    merged.to_csv(output_dir / "all_benchmark_results.csv", index=False)
    logger.info("Collected %d results from %d files", len(merged), len(csv_files))
    return merged
