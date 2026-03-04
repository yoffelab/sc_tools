"""Xenium Ranger command builder and loader for Phase 0.

Most Xenium users receive already-decoded output. The primary use case
is loading existing output, not running Xenium Ranger.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def build_xenium_ranger_cmd(
    sample_id: str,
    input_dir: str,
    *,
    output_dir: str = "data",
) -> str:
    """Build a xenium ranger resegment/import command string.

    Parameters
    ----------
    sample_id
        Sample identifier.
    input_dir
        Path to raw Xenium input directory.
    output_dir
        Output directory prefix.

    Returns
    -------
    Shell command string.

    Notes
    -----
    Most Xenium workflows receive pre-decoded output from the instrument.
    This command builder is for cases where reprocessing is needed.
    """
    parts = [
        "xeniumranger resegment",
        f"--id={sample_id}",
        f"--xenium-bundle={input_dir}",
        f"--output-dir={output_dir}/{sample_id}",
    ]
    return " \\\n    ".join(parts)
