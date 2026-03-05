"""SLURM sbatch script generation for Phase 0 HPC execution.

Generates production-grade sbatch scripts from batch manifests and config,
covering SpaceRanger, Xenium Ranger, and IMC pipelines. Scripts include
SLURM headers, input verification, output cleanup, command execution,
and post-run verification -- matching the pattern of hand-written scripts
in the robin project.

Design decisions:
- Shell variables (${SR}, ${CYTAIMAGE}, etc.) in script body for readability
- Per-sample scripts (not array jobs) for easier failure handling
- Pure text generation; no subprocess or job management
"""

from __future__ import annotations

import logging
import stat
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)

# Default SLURM resource settings
_DEFAULT_SLURM = {
    "partition": "scu-cpu",
    "cpus_per_task": 32,
    "mem": "240G",
    "time": "2-00:00:00",
    "nodes": 1,
    "ntasks": 1,
}


def build_sbatch_header(
    job_name: str,
    log_dir: str = "logs",
    *,
    partition: str = "scu-cpu",
    cpus_per_task: int = 32,
    mem: str = "240G",
    time: str = "2-00:00:00",
    nodes: int = 1,
    ntasks: int = 1,
    extra_directives: dict[str, str] | None = None,
) -> str:
    """Generate an ``#SBATCH`` directive block.

    Parameters
    ----------
    job_name
        SLURM job name.
    log_dir
        Directory for stdout/stderr logs (relative to submission dir).
    partition
        SLURM partition name.
    cpus_per_task
        Number of CPU cores per task.
    mem
        Memory allocation (e.g. "240G").
    time
        Wall time limit (e.g. "2-00:00:00").
    nodes
        Number of nodes.
    ntasks
        Number of tasks.
    extra_directives
        Additional ``#SBATCH`` key-value pairs (e.g. ``{"account": "mylab"}``).

    Returns
    -------
    Multi-line string of ``#SBATCH`` directives.
    """
    lines = [
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --output={log_dir}/{job_name}_%j.out",
        f"#SBATCH --error={log_dir}/{job_name}_%j.err",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --nodes={nodes}",
        f"#SBATCH --ntasks={ntasks}",
        f"#SBATCH --cpus-per-task={cpus_per_task}",
        f"#SBATCH --mem={mem}",
        f"#SBATCH --time={time}",
    ]
    if extra_directives:
        for key, val in extra_directives.items():
            lines.append(f"#SBATCH --{key}={val}")
    return "\n".join(lines)


def _merge_slurm_defaults(slurm: dict | None) -> dict:
    """Merge user slurm dict with defaults."""
    merged = dict(_DEFAULT_SLURM)
    if slurm:
        merged.update(slurm)
    return merged


def build_spaceranger_sbatch(
    sample_id: str,
    fastqs: str,
    transcriptome: str,
    output_dir: str,
    *,
    cytaimage: str | None = None,
    image: str | None = None,
    slide: str | None = None,
    area: str | None = None,
    probe_set: str | None = None,
    sample_filter: str | None = None,
    create_bam: bool = False,
    spaceranger_path: str = "spaceranger",
    localcores: int | None = None,
    localmem: int = 220,
    log_dir: str = "logs",
    slurm: dict | None = None,
) -> str:
    """Generate a full sbatch script for SpaceRanger count.

    Parameters
    ----------
    sample_id
        Unique sample identifier.
    fastqs
        Path to FASTQ directory.
    transcriptome
        Path to reference transcriptome.
    output_dir
        Base output directory for SpaceRanger results.
    cytaimage
        CytAssist image path (Visium HD).
    image
        H&E image path (Visium).
    slide
        Slide serial number.
    area
        Capture area (e.g. D1).
    probe_set
        Path to probe set CSV.
    sample_filter
        Value for --sample flag to select FASTQ subset.
    create_bam
        Whether to create BAM output.
    spaceranger_path
        Path to spaceranger binary.
    localcores
        CPU cores for spaceranger (default: use SLURM_CPUS_PER_TASK).
    localmem
        Memory in GB for spaceranger (default: 220).
    log_dir
        Directory for SLURM log files.
    slurm
        SLURM resource overrides (partition, cpus_per_task, mem, time, etc.).

    Returns
    -------
    Complete sbatch script as a string.

    Raises
    ------
    ValueError
        If neither image nor cytaimage is provided.
    """
    if not image and not cytaimage:
        raise ValueError("Must provide either 'image' (Visium) or 'cytaimage' (Visium HD)")

    s = _merge_slurm_defaults(slurm)
    header = build_sbatch_header(
        job_name=f"sr4_{sample_id}",
        log_dir=log_dir,
        partition=s["partition"],
        cpus_per_task=s["cpus_per_task"],
        mem=s["mem"],
        time=s["time"],
        nodes=s.get("nodes", 1),
        ntasks=s.get("ntasks", 1),
    )

    # Build variables section
    var_lines = [
        f'SR="{spaceranger_path}"',
        f'TRANSCRIPTOME="{transcriptome}"',
        f'FASTQS="{fastqs}"',
        f'OUTPUT_DIR="{output_dir}"',
        f'SAMPLE="{sample_id}"',
    ]
    if probe_set:
        var_lines.append(f'PROBE_SET="{probe_set}"')
    if cytaimage:
        var_lines.append(f'CYTAIMAGE="{cytaimage}"')
    if image:
        var_lines.append(f'IMAGE="{image}"')
    if slide:
        var_lines.append(f'SLIDE="{slide}"')
    if area:
        var_lines.append(f'AREA="{area}"')
    if sample_filter:
        var_lines.append(f'SAMPLE_FILTER="{sample_filter}"')

    cores_expr = (
        f'LOCALCORES="{localcores}"' if localcores else 'LOCALCORES="${SLURM_CPUS_PER_TASK}"'
    )
    var_lines.append(cores_expr)
    var_lines.append(f'LOCALMEM="{localmem}"')

    # Build input checks
    checks = [
        'test -x "${SR}" || { echo "[ERROR] SpaceRanger not found: ${SR}"; exit 1; }',
        'test -d "${TRANSCRIPTOME}" || { echo "[ERROR] Transcriptome not found: ${TRANSCRIPTOME}"; exit 1; }',
    ]
    if probe_set:
        checks.append(
            'test -f "${PROBE_SET}" || { echo "[ERROR] Probe set not found: ${PROBE_SET}"; exit 1; }'
        )
    if cytaimage:
        checks.append(
            'test -f "${CYTAIMAGE}" || { echo "[ERROR] CytAssist image not found: ${CYTAIMAGE}"; exit 1; }'
        )
    if image:
        checks.append(
            'test -f "${IMAGE}" || { echo "[ERROR] H&E image not found: ${IMAGE}"; exit 1; }'
        )

    # Build spaceranger command parts
    cmd_parts = [
        '"${SR}" count \\',
        '    --id="${SAMPLE}" \\',
        '    --transcriptome="${TRANSCRIPTOME}" \\',
    ]
    if probe_set:
        cmd_parts.append('    --probe-set="${PROBE_SET}" \\')
    cmd_parts.append('    --fastqs="${FASTQS}" \\')
    if sample_filter:
        cmd_parts.append('    --sample="${SAMPLE_FILTER}" \\')
    if cytaimage:
        cmd_parts.append('    --cytaimage="${CYTAIMAGE}" \\')
    if image:
        cmd_parts.append('    --image="${IMAGE}" \\')
    if slide:
        cmd_parts.append('    --slide="${SLIDE}" \\')
    if area:
        cmd_parts.append('    --area="${AREA}" \\')
    cmd_parts.append(f"    --create-bam={'true' if create_bam else 'false'} \\")
    cmd_parts.append('    --localcores="${LOCALCORES}" \\')
    cmd_parts.append('    --localmem="${LOCALMEM}" \\')
    cmd_parts.append('    --output-dir="${OUTPUT_DIR}/${SAMPLE}"')

    sr_cmd = "\n".join(cmd_parts)

    # Post-run verification
    verification = _spaceranger_post_verification()

    script = f"""#!/bin/bash
{header}

set -euo pipefail

# ---------------------------------------------------------------------------
# Variables
# ---------------------------------------------------------------------------
{chr(10).join(var_lines)}

mkdir -p {log_dir}

# ---------------------------------------------------------------------------
# Banner
# ---------------------------------------------------------------------------
echo "============================================================"
echo "SpaceRanger 4: ${{SAMPLE}}"
echo "Date:      $(date)"
echo "SLURM Job: ${{SLURM_JOB_ID}}"
echo "Output:    ${{OUTPUT_DIR}}/${{SAMPLE}}"
echo "============================================================"
echo ""

# ---------------------------------------------------------------------------
# Input verification
# ---------------------------------------------------------------------------
echo "[CHECK] Verifying inputs..."
{chr(10).join(checks)}
echo "[CHECK] All inputs verified"
echo ""

# ---------------------------------------------------------------------------
# Clean up previous output
# ---------------------------------------------------------------------------
if [[ -d "${{OUTPUT_DIR}}/${{SAMPLE}}" ]]; then
    echo "[CLEANUP] Removing previous output: ${{OUTPUT_DIR}}/${{SAMPLE}}"
    rm -rf "${{OUTPUT_DIR}}/${{SAMPLE}}"
    echo "[CLEANUP] Done"
    echo ""
fi

# ---------------------------------------------------------------------------
# Run SpaceRanger count
# ---------------------------------------------------------------------------
echo "[RUN] Starting SpaceRanger at $(date)"
echo ""

{sr_cmd}

RC=$?

echo ""
echo "============================================================"
echo "SpaceRanger exit code: ${{RC}}"
echo "Finished at: $(date)"
echo "============================================================"

# ---------------------------------------------------------------------------
# Post-run verification
# ---------------------------------------------------------------------------
{verification}

exit ${{RC}}
"""
    return script


def _spaceranger_post_verification() -> str:
    """Return the post-run verification block for SpaceRanger."""
    return r"""if [[ ${RC} -eq 0 ]]; then
    echo ""
    echo "[VERIFY] Checking outputs..."
    OUTS="${OUTPUT_DIR}/${SAMPLE}/outs"

    # Binned outputs
    for bin in square_002um square_008um square_016um; do
        BINDIR="${OUTS}/binned_outputs/${bin}"
        if [[ -d "${BINDIR}" ]]; then
            echo "  [OK] ${bin} present"
        else
            echo "  [MISS] ${bin}"
        fi
    done

    # Cell segmentation (SR4)
    if [[ -d "${OUTS}/cell_segmentation" ]]; then
        echo "  [OK] cell_segmentation/ present"
    elif [[ -d "${OUTS}/segmented_outputs" ]]; then
        echo "  [OK] segmented_outputs/ present"
    else
        echo "  [WARN] No cell segmentation output found"
    fi

    # Cloupe files
    for cf in cloupe_008um.cloupe cloupe_cell.cloupe; do
        if [[ -f "${OUTS}/${cf}" ]]; then
            SIZE=$(du -h "${OUTS}/${cf}" | cut -f1)
            echo "  [OK] ${cf} (${SIZE})"
        fi
    done
else
    echo ""
    echo "[FAIL] SpaceRanger failed for ${SAMPLE}"
    echo "[FAIL] Check error log in: ${OUTPUT_DIR}/${SAMPLE}/"
    ERR=$(find "${OUTPUT_DIR}/${SAMPLE}" -name "_errors" -exec cat {} \; 2>/dev/null | head -20)
    if [[ -n "${ERR}" ]]; then
        echo "[FAIL] Error details:"
        echo "${ERR}"
    fi
fi"""


def build_xenium_sbatch(
    sample_id: str,
    xenium_bundle: str,
    output_dir: str,
    *,
    xenium_ranger_path: str = "xeniumranger",
    log_dir: str = "logs",
    slurm: dict | None = None,
) -> str:
    """Generate a full sbatch script for Xenium Ranger resegment.

    Parameters
    ----------
    sample_id
        Sample identifier.
    xenium_bundle
        Path to raw Xenium bundle directory.
    output_dir
        Output directory prefix.
    xenium_ranger_path
        Path to xeniumranger binary.
    log_dir
        Directory for SLURM log files.
    slurm
        SLURM resource overrides.

    Returns
    -------
    Complete sbatch script as a string.
    """
    s = _merge_slurm_defaults(slurm)
    header = build_sbatch_header(
        job_name=f"xr_{sample_id}",
        log_dir=log_dir,
        partition=s["partition"],
        cpus_per_task=s["cpus_per_task"],
        mem=s["mem"],
        time=s["time"],
        nodes=s.get("nodes", 1),
        ntasks=s.get("ntasks", 1),
    )

    script = f"""#!/bin/bash
{header}

set -euo pipefail

# ---------------------------------------------------------------------------
# Variables
# ---------------------------------------------------------------------------
XR="{xenium_ranger_path}"
XENIUM_BUNDLE="{xenium_bundle}"
OUTPUT_DIR="{output_dir}"
SAMPLE="{sample_id}"

mkdir -p {log_dir}

# ---------------------------------------------------------------------------
# Banner
# ---------------------------------------------------------------------------
echo "============================================================"
echo "Xenium Ranger: ${{SAMPLE}}"
echo "Date:      $(date)"
echo "SLURM Job: ${{SLURM_JOB_ID}}"
echo "Bundle:    ${{XENIUM_BUNDLE}}"
echo "Output:    ${{OUTPUT_DIR}}/${{SAMPLE}}"
echo "============================================================"
echo ""

# ---------------------------------------------------------------------------
# Input verification
# ---------------------------------------------------------------------------
echo "[CHECK] Verifying inputs..."
test -x "${{XR}}" || {{ echo "[ERROR] Xenium Ranger not found: ${{XR}}"; exit 1; }}
test -d "${{XENIUM_BUNDLE}}" || {{ echo "[ERROR] Xenium bundle not found: ${{XENIUM_BUNDLE}}"; exit 1; }}
echo "[CHECK] All inputs verified"
echo ""

# ---------------------------------------------------------------------------
# Clean up previous output
# ---------------------------------------------------------------------------
if [[ -d "${{OUTPUT_DIR}}/${{SAMPLE}}" ]]; then
    echo "[CLEANUP] Removing previous output: ${{OUTPUT_DIR}}/${{SAMPLE}}"
    rm -rf "${{OUTPUT_DIR}}/${{SAMPLE}}"
    echo "[CLEANUP] Done"
    echo ""
fi

# ---------------------------------------------------------------------------
# Run Xenium Ranger
# ---------------------------------------------------------------------------
echo "[RUN] Starting Xenium Ranger at $(date)"
echo ""

"${{XR}}" resegment \\
    --id="${{SAMPLE}}" \\
    --xenium-bundle="${{XENIUM_BUNDLE}}" \\
    --output-dir="${{OUTPUT_DIR}}/${{SAMPLE}}"

RC=$?

echo ""
echo "============================================================"
echo "Xenium Ranger exit code: ${{RC}}"
echo "Finished at: $(date)"
echo "============================================================"

# ---------------------------------------------------------------------------
# Post-run verification
# ---------------------------------------------------------------------------
if [[ ${{RC}} -eq 0 ]]; then
    echo ""
    echo "[VERIFY] Checking outputs..."
    OUTS="${{OUTPUT_DIR}}/${{SAMPLE}}/outs"
    if [[ -d "${{OUTS}}" ]]; then
        echo "  [OK] outs/ directory present"
        ls "${{OUTS}}/" 2>/dev/null | head -10 | sed 's/^/       /'
    else
        echo "  [WARN] outs/ directory NOT found"
    fi
else
    echo ""
    echo "[FAIL] Xenium Ranger failed for ${{SAMPLE}}"
fi

exit ${{RC}}
"""
    return script


def build_imc_sbatch(
    sample_id: str,
    mcd_file: str,
    panel_csv: str,
    output_dir: str,
    *,
    pipeline_dir: str | None = None,
    log_dir: str = "logs",
    slurm: dict | None = None,
) -> str:
    """Generate a full sbatch script for the IMC pipeline.

    Parameters
    ----------
    sample_id
        Sample identifier.
    mcd_file
        Path to MCD file.
    panel_csv
        Path to panel CSV.
    output_dir
        Output directory prefix.
    pipeline_dir
        Path to IMC pipeline installation directory.
    log_dir
        Directory for SLURM log files.
    slurm
        SLURM resource overrides.

    Returns
    -------
    Complete sbatch script as a string.
    """
    from .imc import build_imc_pipeline_cmd

    s = _merge_slurm_defaults(slurm)
    header = build_sbatch_header(
        job_name=f"imc_{sample_id}",
        log_dir=log_dir,
        partition=s["partition"],
        cpus_per_task=s["cpus_per_task"],
        mem=s["mem"],
        time=s["time"],
        nodes=s.get("nodes", 1),
        ntasks=s.get("ntasks", 1),
    )

    imc_cmd = build_imc_pipeline_cmd(mcd_file, panel_csv, output_dir, pipeline_dir=pipeline_dir)

    script = f"""#!/bin/bash
{header}

set -euo pipefail

# ---------------------------------------------------------------------------
# Variables
# ---------------------------------------------------------------------------
MCD_FILE="{mcd_file}"
PANEL_CSV="{panel_csv}"
OUTPUT_DIR="{output_dir}"
SAMPLE="{sample_id}"

mkdir -p {log_dir}

# ---------------------------------------------------------------------------
# Banner
# ---------------------------------------------------------------------------
echo "============================================================"
echo "IMC Pipeline: ${{SAMPLE}}"
echo "Date:      $(date)"
echo "SLURM Job: ${{SLURM_JOB_ID}}"
echo "MCD:       ${{MCD_FILE}}"
echo "Panel:     ${{PANEL_CSV}}"
echo "Output:    ${{OUTPUT_DIR}}"
echo "============================================================"
echo ""

# ---------------------------------------------------------------------------
# Input verification
# ---------------------------------------------------------------------------
echo "[CHECK] Verifying inputs..."
test -f "${{MCD_FILE}}" || {{ echo "[ERROR] MCD file not found: ${{MCD_FILE}}"; exit 1; }}
test -f "${{PANEL_CSV}}" || {{ echo "[ERROR] Panel CSV not found: ${{PANEL_CSV}}"; exit 1; }}
echo "[CHECK] All inputs verified"
echo ""

# ---------------------------------------------------------------------------
# Run IMC pipeline
# ---------------------------------------------------------------------------
echo "[RUN] Starting IMC pipeline at $(date)"
echo ""

{imc_cmd}

RC=$?

echo ""
echo "============================================================"
echo "IMC pipeline exit code: ${{RC}}"
echo "Finished at: $(date)"
echo "============================================================"

# ---------------------------------------------------------------------------
# Post-run verification
# ---------------------------------------------------------------------------
if [[ ${{RC}} -eq 0 ]]; then
    echo ""
    echo "[VERIFY] Checking outputs..."
    if [[ -d "${{OUTPUT_DIR}}/tiffs" ]]; then
        echo "  [OK] tiffs/ directory present"
        NTIFFS=$(ls "${{OUTPUT_DIR}}/tiffs/"*_full.tiff 2>/dev/null | wc -l)
        echo "  [OK] ${{NTIFFS}} TIFF stacks found"
    else
        echo "  [WARN] tiffs/ directory NOT found"
    fi
    if [[ -f "${{OUTPUT_DIR}}/cells.h5ad" ]]; then
        echo "  [OK] cells.h5ad present"
    else
        echo "  [WARN] cells.h5ad NOT found"
    fi
else
    echo ""
    echo "[FAIL] IMC pipeline failed for ${{SAMPLE}}"
fi

exit ${{RC}}
"""
    return script


def build_batch_sbatch(
    manifest: pd.DataFrame,
    modality: str,
    output_dir: str,
    *,
    transcriptome: str | None = None,
    probe_set: str | None = None,
    spaceranger_path: str = "spaceranger",
    xenium_ranger_path: str = "xeniumranger",
    pipeline_dir: str | None = None,
    create_bam: bool = False,
    localcores: int | None = None,
    localmem: int = 220,
    log_dir: str = "logs",
    slurm: dict | None = None,
) -> list[tuple[str, str]]:
    """Generate one sbatch script per manifest row.

    Parameters
    ----------
    manifest
        DataFrame from batch manifest TSV.
    modality
        One of "visium", "visium_hd", "visium_hd_cell", "xenium", "imc".
    output_dir
        Base output directory.
    transcriptome
        Path to reference transcriptome (SpaceRanger modalities).
    probe_set
        Path to probe set CSV (SpaceRanger).
    spaceranger_path
        Path to spaceranger binary.
    xenium_ranger_path
        Path to xeniumranger binary.
    pipeline_dir
        Path to IMC pipeline directory.
    create_bam
        Whether to create BAM output (SpaceRanger).
    localcores
        CPU cores override (default: use SLURM_CPUS_PER_TASK).
    localmem
        Memory in GB for the tool.
    log_dir
        Directory for SLURM log files.
    slurm
        SLURM resource overrides.

    Returns
    -------
    List of (sample_id, script_text) tuples.
    """
    scripts = []

    for _, row in manifest.iterrows():
        sample_id = str(row["sample_id"])

        try:
            if modality in ("visium", "visium_hd", "visium_hd_cell"):
                if not transcriptome:
                    logger.warning(
                        "Skipping %s: transcriptome required for %s", sample_id, modality
                    )
                    continue

                kwargs = {
                    "sample_id": sample_id,
                    "fastqs": str(row.get("fastq_dir", "")),
                    "transcriptome": transcriptome,
                    "output_dir": output_dir,
                    "spaceranger_path": spaceranger_path,
                    "create_bam": create_bam,
                    "localcores": localcores,
                    "localmem": localmem,
                    "log_dir": log_dir,
                    "slurm": slurm,
                }

                # Map optional manifest columns
                for col, arg in [
                    ("cytaimage", "cytaimage"),
                    ("image", "image"),
                    ("slide", "slide"),
                    ("area", "area"),
                    ("sample_filter", "sample_filter"),
                ]:
                    val = row.get(col)
                    if val is not None and str(val) not in ("", "nan", "None"):
                        kwargs[arg] = str(val)

                # probe_set: manifest column overrides function arg
                row_probe = row.get("probe_set")
                if row_probe is not None and str(row_probe) not in ("", "nan", "None"):
                    kwargs["probe_set"] = str(row_probe)
                elif probe_set:
                    kwargs["probe_set"] = probe_set

                script = build_spaceranger_sbatch(**kwargs)

            elif modality == "xenium":
                xenium_dir = str(row.get("xenium_dir", ""))
                script = build_xenium_sbatch(
                    sample_id=sample_id,
                    xenium_bundle=xenium_dir,
                    output_dir=output_dir,
                    xenium_ranger_path=xenium_ranger_path,
                    log_dir=log_dir,
                    slurm=slurm,
                )

            elif modality == "imc":
                mcd_file = str(row.get("mcd_file", ""))
                panel_csv = str(row.get("panel_csv", ""))
                script = build_imc_sbatch(
                    sample_id=sample_id,
                    mcd_file=mcd_file,
                    panel_csv=panel_csv,
                    output_dir=output_dir,
                    pipeline_dir=pipeline_dir,
                    log_dir=log_dir,
                    slurm=slurm,
                )

            else:
                logger.warning("Skipping %s: unsupported modality '%s'", sample_id, modality)
                continue

            scripts.append((sample_id, script))

        except (ValueError, KeyError) as e:
            logger.warning("Skipping %s: %s", sample_id, e)

    return scripts


def write_sbatch_script(
    script_text: str,
    output_path: str | Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Write an sbatch script to disk and make it executable.

    Parameters
    ----------
    script_text
        The sbatch script content.
    output_path
        File path to write to.
    overwrite
        If False, raise FileExistsError when the file already exists.

    Returns
    -------
    Path to the written script.

    Raises
    ------
    FileExistsError
        If the file exists and overwrite is False.
    """
    path = Path(output_path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Script already exists: {path}. Use overwrite=True to replace.")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(script_text)

    # Make executable (owner rwx, group rx, other rx)
    current = path.stat().st_mode
    path.chmod(current | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    logger.info("Wrote sbatch script: %s", path)
    return path
