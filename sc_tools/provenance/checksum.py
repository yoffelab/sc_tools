"""SHA256 streaming file checksum."""

from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_file(path: Path | str, chunk_size: int = 65536) -> str:
    """Compute SHA256 hex digest of a file using streaming reads.

    Parameters
    ----------
    path
        Path to file.
    chunk_size
        Read buffer size in bytes (default 64KB).

    Returns
    -------
    Hex digest string.
    """
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()
