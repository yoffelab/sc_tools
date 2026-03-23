"""Provenance sidecar file utilities.

Provides write/read/path helpers for .provenance.json sidecar files,
peak memory measurement, provenance record building, and h5ad uns embedding.
"""

from __future__ import annotations

import json
import logging
import os
import platform
import resource
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sc_tools.models.result import ProvenanceRecord

logger = logging.getLogger(__name__)


def sidecar_path_for(artifact_path: str | Path) -> Path:
    """Return the sidecar path for an artifact (D-02).

    Convention: ``{original_filename}.provenance.json``.
    """
    return Path(str(artifact_path) + ".provenance.json")


def get_peak_memory_mb() -> float:
    """Return peak RSS memory in megabytes.

    Uses ``resource.getrusage`` with platform-aware unit handling:
    macOS returns bytes, Linux returns kilobytes.
    """
    ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if platform.system() == "Darwin":
        return ru / (1024 * 1024)  # bytes -> MB
    return ru / 1024  # KB -> MB


def write_sidecar(artifact_path: str | Path, record: ProvenanceRecord) -> Path:
    """Write provenance sidecar atomically (D-01).

    Uses tempfile + os.replace for atomic writes to prevent corruption
    from parallel execution.

    Parameters
    ----------
    artifact_path
        Path to the artifact file.
    record
        ProvenanceRecord to serialize.

    Returns
    -------
    Path to the written sidecar file.
    """
    sp = sidecar_path_for(artifact_path)
    data = record.model_dump(mode="json")
    tmp_fd, tmp_path = tempfile.mkstemp(dir=sp.parent, suffix=".tmp")
    try:
        with os.fdopen(tmp_fd, "w") as f:
            json.dump(data, f, indent=2)
        os.replace(tmp_path, sp)
    except Exception:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise
    return sp


def read_sidecar(artifact_path: str | Path) -> dict | None:
    """Read provenance sidecar if it exists.

    Returns
    -------
    Parsed dict or None if sidecar does not exist.
    """
    sp = sidecar_path_for(artifact_path)
    if not sp.exists():
        return None
    with open(sp) as f:
        return json.load(f)


def build_provenance_record(
    command: str,
    params: dict,
    input_files: list[str],
    runtime_s: float,
    peak_memory_mb: float,
    project_dir: str = ".",
) -> ProvenanceRecord:
    """Build a ProvenanceRecord with InputFile records for each input.

    Computes SHA256 checksums and file sizes for each input file.
    Uses relative paths when possible (D-09).

    Parameters
    ----------
    command
        CLI command name.
    params
        All CLI parameters as-passed (D-07).
    input_files
        List of input file paths.
    runtime_s
        Execution time in seconds.
    peak_memory_mb
        Peak memory usage in MB.
    project_dir
        Project root for relative path computation (D-09).

    Returns
    -------
    ProvenanceRecord instance.
    """
    from sc_tools.models.result import InputFile, ProvenanceRecord

    from sc_tools.provenance.checksum import sha256_file

    inputs = []
    proj = Path(project_dir).resolve()

    for fpath in input_files:
        p = Path(fpath)
        if not p.exists():
            logger.warning("Input file not found for provenance: %s", fpath)
            continue

        # Compute relative path (D-09)
        try:
            rel = str(p.resolve().relative_to(proj))
            path_type = "relative"
        except ValueError:
            rel = str(p.resolve())
            path_type = "absolute"

        inputs.append(
            InputFile(
                path=rel,
                path_type=path_type,
                sha256=sha256_file(p),
                size_bytes=p.stat().st_size,
            )
        )

    return ProvenanceRecord(
        command=command,
        params=params,
        inputs=inputs,
        runtime_s=runtime_s,
        peak_memory_mb=peak_memory_mb,
    )


def embed_provenance_in_adata(
    artifact_path: str | Path,
    record: ProvenanceRecord,
) -> None:
    """Embed provenance in h5ad file's uns/sct_provenance (D-04).

    Opens the h5ad file with h5py in append mode and writes the
    provenance record as a JSON string. Uses h5py directly to avoid
    loading the full AnnData into memory.

    Parameters
    ----------
    artifact_path
        Path to the h5ad file.
    record
        ProvenanceRecord to embed.
    """
    try:
        import h5py
    except ImportError:
        logger.warning("h5py not available, skipping uns provenance embedding")
        return

    try:
        prov_json = json.dumps(record.model_dump(mode="json"))
        with h5py.File(artifact_path, "a") as f:
            # Ensure uns group exists
            if "uns" not in f:
                f.create_group("uns")
            # Delete existing key if present
            if "uns/sct_provenance" in f:
                del f["uns/sct_provenance"]
            f["uns/sct_provenance"] = prov_json
    except Exception:
        logger.warning(
            "Failed to embed provenance in %s", artifact_path, exc_info=True
        )
