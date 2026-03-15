"""Phase 0: Upstream raw data processing and ingestion.

Provides batch manifest management, platform-specific command builders,
modality-aware AnnData loaders, and SLURM sbatch script generation.

Submodules:
- config: Batch manifest parsing and collection
- spaceranger: Space Ranger command builder
- xenium: Xenium Ranger commands and loader
- imc: IMC pipeline commands
- loaders: Modality-specific AnnData loaders and concatenation
- slurm: SLURM sbatch script generation for HPC execution
"""

from __future__ import annotations

from .config import collect_all_batches, load_batch_manifest, validate_manifest
from .imc import IMCPanelMapper, build_imc_composite
from .loaders import (
    concat_samples,
    load_cosmx_sample,
    load_he_image,
    load_imc_sample,
    load_visium_hd_cell_sample,
    load_visium_hd_sample,
    load_visium_sample,
    load_xenium_sample,
)
from .slurm import (
    build_batch_sbatch,
    build_imc_sbatch,
    build_sbatch_header,
    build_spaceranger_sbatch,
    build_xenium_sbatch,
    generate_phase0_inventory,
    load_phase0_status,
    save_phase0_status,
    write_sbatch_script,
)
from .spaceranger import build_batch_commands, build_spaceranger_count_cmd

__all__ = [
    "load_batch_manifest",
    "collect_all_batches",
    "validate_manifest",
    "build_spaceranger_count_cmd",
    "build_batch_commands",
    "load_visium_sample",
    "load_visium_hd_sample",
    "load_visium_hd_cell_sample",
    "load_xenium_sample",
    "load_cosmx_sample",
    "load_imc_sample",
    "load_he_image",
    "concat_samples",
    "IMCPanelMapper",
    "build_imc_composite",
    "build_sbatch_header",
    "build_spaceranger_sbatch",
    "build_xenium_sbatch",
    "build_imc_sbatch",
    "build_batch_sbatch",
    "write_sbatch_script",
    "generate_phase0_inventory",
    "load_phase0_status",
    "save_phase0_status",
]
