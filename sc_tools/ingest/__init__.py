"""Phase 0: Upstream raw data processing and ingestion.

Provides batch manifest management, platform-specific command builders,
and modality-aware AnnData loaders.

Submodules:
- config: Batch manifest parsing and collection
- spaceranger: Space Ranger command builder
- xenium: Xenium Ranger commands and loader
- imc: IMC pipeline commands
- loaders: Modality-specific AnnData loaders and concatenation
"""

from __future__ import annotations

from .config import collect_all_batches, load_batch_manifest, validate_manifest
from .loaders import (
    concat_samples,
    load_imc_sample,
    load_visium_hd_sample,
    load_visium_sample,
    load_xenium_sample,
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
    "load_xenium_sample",
    "load_imc_sample",
    "concat_samples",
]
