"""IMC pipeline command builder for Phase 0a.

Builds shell commands to run the ElementoLab IMC pipeline for segmentation
and single-cell extraction from MCD files. The pipeline outputs to::

    processed/{sample}/tiffs/   # per-channel TIFF images
    processed/{sample}/masks/   # segmentation masks
    processed/{sample}/cells.h5ad  # single-cell quantification (Phase 0b input)

Phase 0b loading is handled by ``sc_tools.ingest.loaders.load_imc_sample()``,
which reads ``processed/{sample}/cells.h5ad``.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def build_imc_pipeline_cmd(
    mcd_file: str,
    panel_csv: str,
    output_dir: str,
    *,
    pipeline_dir: str | None = None,
) -> str:
    """Build shell command to run the ElementoLab IMC pipeline.

    Parameters
    ----------
    mcd_file
        Path to .mcd file.
    panel_csv
        Path to panel CSV defining channels.
    output_dir
        Root output directory. The pipeline writes per-sample results under
        ``output_dir/{sample}/tiffs/``, ``masks/``, and ``cells.h5ad``.
        Typically set to ``processed/`` in the project directory.
    pipeline_dir
        Path to cloned ElementoLab IMC pipeline repo. Defaults to
        ``~/elementolab/imc``.

    Returns
    -------
    Shell command string.

    Notes
    -----
    Assumes the ElementoLab IMC pipeline is cloned and available on the
    system (typically on HPC). The pipeline handles: MCD extraction,
    ilastik segmentation, and single-cell quantification. Phase 0b loading
    then reads ``processed/{sample}/cells.h5ad`` via ``load_imc_sample()``.
    """
    if pipeline_dir is None:
        pipeline_dir = "~/elementolab/imc"

    parts = [
        f"python {pipeline_dir}/run_pipeline.py",
        f"--mcd {mcd_file}",
        f"--panel {panel_csv}",
        f"--output {output_dir}",
    ]
    return " \\\n    ".join(parts)
