"""Lineage trace engine for provenance chain walking.

Provides BFS-based lineage tracing that follows .provenance.json sidecars
recursively, with cycle detection, adata.uns fallback, and SHA256-based
file relocation support.

Per D-11: trace follows sidecar input references recursively.
Per D-12: output is chronological (oldest first).
Per D-10: missing files searched by SHA256.
Per D-05: adata.uns fallback when sidecar missing for h5ad.
"""

from __future__ import annotations

import json
import logging
import os
from collections import deque
from pathlib import Path

logger = logging.getLogger(__name__)

# Maximum files to scan during SHA256 relocation search
_SHA256_SCAN_LIMIT = 1000

# Directories to skip during SHA256 scan
_SKIP_DIRS = {".git", "__pycache__", ".tox", ".mypy_cache", "node_modules", ".eggs"}


def trace_lineage(
    file_path: str | Path,
    project_dir: str | Path = ".",
) -> list[dict]:
    """Walk provenance chains recursively to build full lineage.

    Uses BFS starting from the given file. For each file, reads its
    provenance sidecar (or adata.uns fallback for h5ad files), then
    queues all input files for further tracing.

    Parameters
    ----------
    file_path
        Path to the artifact to trace.
    project_dir
        Project root for relative path resolution and SHA256 relocation search.

    Returns
    -------
    List of lineage steps in chronological order (oldest first).
    Each step is a dict with keys: file, command, timestamp, note.
    """
    from sc_tools.provenance.sidecar import read_sidecar

    file_path = Path(file_path)
    project_dir = Path(project_dir)

    steps: list[dict] = []
    visited: set[str] = set()
    queue: deque[tuple[Path, str | None]] = deque()

    # Start with the target file
    queue.append((file_path, None))

    while queue:
        current_path, relocation_note = queue.popleft()

        # Canonical path for cycle detection
        try:
            canonical = str(current_path.resolve())
        except OSError:
            canonical = str(current_path)

        if canonical in visited:
            continue
        visited.add(canonical)

        # Try to read provenance
        prov = None

        # 1. Try sidecar
        if current_path.exists():
            prov = read_sidecar(str(current_path))

        # 2. Fallback to adata.uns for h5ad files
        if prov is None and str(current_path).endswith(".h5ad") and current_path.exists():
            prov = _read_uns_provenance(current_path)

        if prov is None:
            # No provenance found: origin node
            # Preserve relocation note if file was found via SHA256
            note = (
                "input relocated; origin (no provenance)"
                if relocation_note
                else "origin (no provenance)"
            )
            steps.append(
                {
                    "file": str(current_path),
                    "command": None,
                    "timestamp": None,
                    "note": note,
                }
            )
        else:
            # Provenance found
            note = relocation_note
            steps.append(
                {
                    "file": str(current_path),
                    "command": prov.get("command"),
                    "timestamp": prov.get("timestamp"),
                    "note": note,
                }
            )

            # Queue all inputs for further tracing
            for inp in prov.get("inputs", []):
                inp_path = Path(inp.get("path", ""))
                inp_sha256 = inp.get("sha256")

                if not inp_path.is_absolute():
                    inp_path = project_dir / inp_path

                if inp_path.exists():
                    queue.append((inp_path, None))
                else:
                    # File missing: try adata.uns fallback for h5ad
                    if str(inp_path).endswith(".h5ad"):
                        uns_prov = _read_uns_provenance(inp_path) if inp_path.exists() else None
                        if uns_prov is not None:
                            queue.append((inp_path, None))
                            continue

                    # Try SHA256 relocation (D-10)
                    if inp_sha256:
                        relocated = _find_by_sha256(inp_sha256, project_dir)
                        if relocated is not None:
                            queue.append((relocated, "input relocated"))
                            continue

                    # File truly missing: mark as origin
                    inp_canonical = str(inp_path)
                    if inp_canonical not in visited:
                        visited.add(inp_canonical)
                        steps.append(
                            {
                                "file": str(inp_path),
                                "command": None,
                                "timestamp": None,
                                "note": "origin (no provenance)",
                            }
                        )

    # Reverse for chronological order (oldest first, per D-12)
    steps.reverse()
    return steps


def _read_uns_provenance(path: Path) -> dict | None:
    """Read provenance from h5ad adata.uns['sct_provenance'] without full AnnData load.

    Uses h5py directly to avoid loading the full dataset into memory.

    Parameters
    ----------
    path
        Path to h5ad file.

    Returns
    -------
    Parsed provenance dict, or None if not found or error.
    """
    try:
        import h5py
    except ImportError:
        logger.debug("h5py not available, skipping uns provenance read")
        return None

    try:
        with h5py.File(path, "r") as f:
            if "uns/sct_provenance" not in f:
                return None
            raw = f["uns/sct_provenance"][()]
            if isinstance(raw, bytes):
                raw = raw.decode("utf-8")
            return json.loads(raw)
    except Exception:
        logger.debug("Failed to read uns provenance from %s", path, exc_info=True)
        return None


def _find_by_sha256(sha256: str, project_dir: str | Path) -> Path | None:
    """Search project_dir recursively for a file matching the given SHA256.

    Skips hidden directories, .git, __pycache__, and other non-data dirs.
    Stops after scanning _SHA256_SCAN_LIMIT files to avoid huge directory trees.

    Parameters
    ----------
    sha256
        Expected SHA256 hex digest.
    project_dir
        Root directory to search.

    Returns
    -------
    Path to first matching file, or None.
    """
    from sc_tools.provenance.checksum import sha256_file

    project_dir = Path(project_dir)
    checked = 0

    for dirpath, dirnames, filenames in os.walk(project_dir):
        # Prune hidden dirs and known non-data dirs
        dirnames[:] = [
            d for d in dirnames if not d.startswith(".") and d not in _SKIP_DIRS
        ]

        for fname in filenames:
            if fname.startswith("."):
                continue
            if checked >= _SHA256_SCAN_LIMIT:
                return None

            fpath = Path(dirpath) / fname
            checked += 1

            try:
                if sha256_file(fpath) == sha256:
                    return fpath
            except (OSError, PermissionError):
                continue

    return None
