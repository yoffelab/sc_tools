"""BenchmarkConfig dataclass with YAML I/O for IMC segmentation benchmarks."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

__all__ = ["BenchmarkConfig"]

logger = logging.getLogger(__name__)


@dataclass
class BenchmarkConfig:
    """Configuration for an IMC segmentation benchmark run.

    Parameters
    ----------
    strategies
        Which strategies to run (1-4). Default: all.
    methods
        Method names within each strategy. Default: all available.
    n_rois
        Max ROIs per dataset (0 = all).
    gpu
        Whether to use GPU for DL methods.
    batch_size
        Batch size for DL inference.
    prob_map_method
        Method for generating probability maps in Strategy 2.
    hf_models
        HuggingFace model names for Strategy 3.
    vit_model_size
        SegFormer model size for Strategy 4.
    vit_n_epochs
        Training epochs for Strategy 4.
    vit_train_split
        Fraction of data for training in Strategy 4.
    slurm_partition_cpu
        SLURM partition for CPU jobs.
    slurm_partition_gpu
        SLURM partition for GPU jobs.
    slurm_time_cpu
        Time limit for CPU jobs (HH:MM:SS).
    slurm_time_gpu
        Time limit for GPU jobs (HH:MM:SS).
    slurm_mem
        Memory per job.
    slurm_chunk_size
        ROIs per SLURM job.
    output_dir
        Where to write benchmark results.
    resume
        Resume from previous run (skip completed ROIs).
    """

    strategies: list[int] = field(default_factory=lambda: [1, 2, 3, 4])
    methods: dict[int, list[str]] = field(
        default_factory=lambda: {
            1: ["cellpose_cyto2", "cellpose_cyto3", "cellpose_nuclei", "stardist", "deepcell"],
            2: ["cellpose_cyto2", "cellpose_nuclei", "stardist", "deepcell"],
            3: ["cellvit_256", "sam_base", "hovernet", "cellpose_cyto3", "stardist"],
            4: ["segformer_b0"],
        }
    )
    n_rois: int = 0
    gpu: bool = True
    batch_size: int = 1
    prob_map_method: str = "gaussian"
    hf_models: list[str] = field(
        default_factory=lambda: [
            "cellvit_256",
            "sam_base",
            "hovernet",
        ]
    )
    vit_model_size: str = "b0"
    vit_n_epochs: int = 50
    vit_train_split: float = 0.8
    slurm_partition_cpu: str = "scu-cpu"
    slurm_partition_gpu: str = "scu-gpu"
    slurm_time_cpu: str = "04:00:00"
    slurm_time_gpu: str = "08:00:00"
    slurm_mem: str = "32G"
    slurm_chunk_size: int = 20
    output_dir: str = "benchmark_results"
    resume: bool = True

    def to_yaml(self, path: str | Path) -> None:
        """Save config to YAML file."""
        try:
            import yaml
        except ImportError as e:
            raise ImportError("pyyaml required: pip install pyyaml") from e

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        data = {
            "strategies": self.strategies,
            "methods": {str(k): v for k, v in self.methods.items()},
            "n_rois": self.n_rois,
            "gpu": self.gpu,
            "batch_size": self.batch_size,
            "prob_map_method": self.prob_map_method,
            "hf_models": self.hf_models,
            "vit_model_size": self.vit_model_size,
            "vit_n_epochs": self.vit_n_epochs,
            "vit_train_split": self.vit_train_split,
            "slurm_partition_cpu": self.slurm_partition_cpu,
            "slurm_partition_gpu": self.slurm_partition_gpu,
            "slurm_time_cpu": self.slurm_time_cpu,
            "slurm_time_gpu": self.slurm_time_gpu,
            "slurm_mem": self.slurm_mem,
            "slurm_chunk_size": self.slurm_chunk_size,
            "output_dir": self.output_dir,
            "resume": self.resume,
        }

        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

        logger.info("Config saved to %s", path)

    @classmethod
    def from_yaml(cls, path: str | Path) -> BenchmarkConfig:
        """Load config from YAML file."""
        try:
            import yaml
        except ImportError as e:
            raise ImportError("pyyaml required: pip install pyyaml") from e

        path = Path(path)
        with open(path) as f:
            data = yaml.safe_load(f)

        # Convert method keys back to int
        if "methods" in data:
            data["methods"] = {int(k): v for k, v in data["methods"].items()}

        return cls(**data)
