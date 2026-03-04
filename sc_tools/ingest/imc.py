"""IMC pipeline command builder for Phase 0.

Builds shell commands to run the ElementoLab IMC pipeline for segmentation
and single-cell extraction from MCD files.
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
        Output directory for processed results.
    pipeline_dir
        Path to cloned IMC pipeline repo. Defaults to ~/elementolab/imc.

    Returns
    -------
    Shell command string.

    Notes
    -----
    Assumes the ElementoLab IMC pipeline is cloned and available on the
    system (typically on HPC). The pipeline handles: MCD extraction,
    ilastik segmentation, and single-cell quantification.
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
