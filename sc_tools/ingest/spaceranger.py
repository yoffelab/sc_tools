"""Space Ranger command builder for Phase 0.

Builds spaceranger count CLI commands for Visium and Visium HD samples.
Commands are returned as strings for shell execution or Snakemake rules.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)


def build_spaceranger_count_cmd(
    sample_id: str,
    fastqs: str,
    transcriptome: str,
    *,
    image: str | None = None,
    cytaimage: str | None = None,
    slide: str | None = None,
    area: str | None = None,
    output_dir: str = "data",
    localcores: int = 16,
    localmem: int = 64,
) -> str:
    """Build a spaceranger count CLI command string.

    For Visium: uses --image.
    For Visium HD: uses --cytaimage.

    Parameters
    ----------
    sample_id
        Unique sample identifier (used as --id).
    fastqs
        Path to FASTQ directory.
    transcriptome
        Path to reference transcriptome.
    image
        H&E image path (Visium).
    cytaimage
        CytAssist image path (Visium HD).
    slide
        Slide serial number.
    area
        Capture area (e.g., D1, A1).
    output_dir
        Output directory prefix.
    localcores
        Number of CPU cores.
    localmem
        Memory limit in GB.

    Returns
    -------
    Shell command string.

    Raises
    ------
    ValueError
        If neither image nor cytaimage is provided.
    """
    if not image and not cytaimage:
        raise ValueError("Must provide either 'image' (Visium) or 'cytaimage' (Visium HD)")

    parts = [
        "spaceranger count",
        f"--id={sample_id}",
        f"--transcriptome={transcriptome}",
        f"--fastqs={fastqs}",
    ]

    if cytaimage:
        parts.append(f"--cytaimage={cytaimage}")
    if image:
        parts.append(f"--image={image}")
    if slide:
        parts.append(f"--slide={slide}")
    if area:
        parts.append(f"--area={area}")

    parts.extend(
        [
            f"--output-dir={output_dir}/{sample_id}",
            f"--localcores={localcores}",
            f"--localmem={localmem}",
        ]
    )

    return " \\\n    ".join(parts)


def build_batch_commands(
    manifest: pd.DataFrame,
    transcriptome: str,
    **kwargs,
) -> list[tuple[str, str]]:
    """Build spaceranger count commands for all samples in a manifest.

    Parameters
    ----------
    manifest
        DataFrame with columns: sample_id, fastq_dir, and optionally
        image, cytaimage, slide, area.
    transcriptome
        Path to reference transcriptome.
    **kwargs
        Additional arguments passed to build_spaceranger_count_cmd
        (e.g., localcores, localmem, output_dir).

    Returns
    -------
    List of (sample_id, command_string) tuples.
    """
    commands = []
    for _, row in manifest.iterrows():
        sample_id = row["sample_id"]
        fastqs = row.get("fastq_dir", "")

        cmd_kwargs = {
            "sample_id": sample_id,
            "fastqs": fastqs,
            "transcriptome": transcriptome,
        }

        # Map manifest columns to command arguments
        for col, arg in [
            ("image", "image"),
            ("cytaimage", "cytaimage"),
            ("slide", "slide"),
            ("area", "area"),
        ]:
            val = row.get(col)
            if val is not None and str(val) not in ("", "nan", "None"):
                cmd_kwargs[arg] = str(val)

        cmd_kwargs.update(kwargs)

        try:
            cmd = build_spaceranger_count_cmd(**cmd_kwargs)
            commands.append((sample_id, cmd))
        except ValueError as e:
            logger.warning("Skipping %s: %s", sample_id, e)

    return commands
