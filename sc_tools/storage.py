"""Storage abstraction: fsspec URI resolution + smart read/write helpers.

Every function accepts either a plain local path (backwards-compatible)
or an fsspec URI.  The correct backend is selected automatically.

URI schemes supported
---------------------
    /path or file://  -- local filesystem (always available)
    sftp://host/path  -- SSH/HPC (requires: pip install sshfs)
    s3://bucket/path  -- AWS S3 (requires: pip install s3fs)
    gs://bucket/path  -- GCS (requires: pip install gcsfs)
    az://cont/path    -- Azure Blob (requires: pip install adlfs)
    box://path        -- Box (requires: pip install boxfs)

Install all remote backends::

    pip install "sc-tools[storage]"
"""

from __future__ import annotations

import io
import logging
import os
import tempfile
from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Any

import pandas as pd

logger = logging.getLogger(__name__)

# Hint messages for missing optional backends
_BACKEND_HINT: dict[str, str] = {
    "sftp": "pip install sshfs",
    "s3": "pip install s3fs",
    "gs": "pip install gcsfs",
    "gcs": "pip install gcsfs",
    "box": "pip install boxfs",
    "az": "pip install adlfs",
    "abfs": "pip install adlfs",
}


def resolve_fs(uri: str | os.PathLike) -> tuple[Any, str]:
    """Return ``(AbstractFileSystem, path)`` for any URI or local path.

    Parameters
    ----------
    uri
        File URI or local path. Examples:
        ``"/local/path"``, ``"s3://bucket/key"``,
        ``"sftp://brb//athena/.../file.h5ad"``.

    Returns
    -------
    tuple[AbstractFileSystem, str]
        Filesystem object and the path within that filesystem.

    Raises
    ------
    ImportError
        When the backend required for the URI scheme is not installed.
    """
    try:
        import fsspec
    except ImportError as e:
        raise ImportError(
            "fsspec is required for URI resolution. Install with: pip install sc-tools[storage]"
        ) from e

    uri = str(uri)
    try:
        fs, path = fsspec.url_to_fs(uri)
    except ImportError as e:
        protocol = uri.split("://")[0] if "://" in uri else "local"
        hint = _BACKEND_HINT.get(protocol, "pip install sc-tools[storage]")
        raise ImportError(
            f"Storage backend for '{protocol}://' not available. Install: {hint}"
        ) from e
    return fs, path


def _is_local(fs: Any) -> bool:
    """Return True when *fs* is a local (POSIX) filesystem."""
    try:
        from fsspec.implementations.local import LocalFileSystem

        return isinstance(fs, LocalFileSystem)
    except ImportError:
        return True  # fsspec not available → assume local


@contextmanager
def open_file(uri: str | os.PathLike, mode: str = "rb") -> Generator[IO, None, None]:
    """Open any URI for reading or writing.

    Parameters
    ----------
    uri
        File URI or local path.
    mode
        Open mode (``"rb"``, ``"wb"``, ``"r"``, ``"w"``).

    Yields
    ------
    IO
        File-like object opened at the given URI.
    """
    fs, path = resolve_fs(uri)
    with fs.open(path, mode=mode) as f:
        yield f


@contextmanager
def with_local_copy(
    uri: str | os.PathLike,
) -> Generator[Path, None, None]:
    """Yield a local *Path*, downloading from remote when necessary.

    For local URIs this is a no-op (returns the original path unchanged).
    For remote URIs the file is downloaded to a temporary file first.

    Useful for libraries (``scanpy``, ``tifffile``) that need a real
    filesystem path rather than a file-like object.

    Parameters
    ----------
    uri
        File URI or local path.

    Yields
    ------
    Path
        Local path to the file.
    """
    fs, path = resolve_fs(str(uri))
    if _is_local(fs):
        yield Path(path)
    else:
        suffix = Path(path).suffix or ".tmp"
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            logger.info("Downloading remote file to local tmp: %s", uri)
            fs.get(path, str(tmp_path))
            yield tmp_path
        finally:
            tmp_path.unlink(missing_ok=True)


def smart_read_h5ad(uri: str | os.PathLike, *, backed: str | None = None):
    """Read an AnnData ``.h5ad`` file from any URI.

    For local paths calls ``anndata.read_h5ad`` directly (backed mode
    supported).  For remote URIs the file is downloaded first (backed
    mode not supported remotely).

    Parameters
    ----------
    uri
        Local path or URI to an ``.h5ad`` file.
    backed
        Backing mode for local files: ``"r"`` (read-only) or ``"r+"``
        (read-write).  Ignored for remote URIs.

    Returns
    -------
    AnnData
    """
    import anndata as ad

    with with_local_copy(uri) as local:
        return ad.read_h5ad(local, backed=backed)


def smart_write_checkpoint(adata: Any, uri: str | os.PathLike, *, fmt: str = "h5ad") -> None:
    """Write AnnData to any URI in h5ad or zarr format.

    For local paths writes directly.  For remote URIs with ``fmt="h5ad"``
    the file is serialised to a buffer and then uploaded.  Zarr format
    delegates to ``adata.write_zarr(uri)`` which uses fsspec natively.

    Parameters
    ----------
    adata
        AnnData object to write.
    uri
        Destination path or URI.
    fmt
        ``"h5ad"`` (default) or ``"zarr"``.  Prefer ``"zarr"`` for cloud
        storage (``s3://``, ``gs://``, ``az://``) — it is chunked and
        cloud-native.
    """
    uri = str(uri)
    fs, path = resolve_fs(uri)
    local = _is_local(fs)

    if fmt == "h5ad":
        if local:
            parent = os.path.dirname(path)
            if parent:
                os.makedirs(parent, exist_ok=True)
            adata.write_h5ad(path)
            logger.info("Wrote h5ad: %s", path)
        else:
            buf = io.BytesIO()
            adata.write_h5ad(buf)
            buf.seek(0)
            with fs.open(path, "wb") as f:
                f.write(buf.read())
            logger.info("Wrote h5ad to remote URI: %s", uri)
    elif fmt == "zarr":
        try:
            import zarr  # noqa: F401
        except ImportError as e:
            raise ImportError("zarr is required for zarr format: pip install zarr") from e
        if local:
            parent = os.path.dirname(path)
            if parent:
                os.makedirs(parent, exist_ok=True)
        adata.write_zarr(uri)
        logger.info("Wrote zarr: %s", uri)
    else:
        raise ValueError(f"Unknown checkpoint format '{fmt}'. Must be 'h5ad' or 'zarr'.")


def smart_read_csv(uri: str | os.PathLike, **kwargs: Any) -> pd.DataFrame:
    """Read a CSV or TSV file from any URI.

    Parameters
    ----------
    uri
        Local path or URI to a CSV or TSV file.
    **kwargs
        Forwarded to ``pd.read_csv``.

    Returns
    -------
    pd.DataFrame
    """
    with open_file(uri, "rb") as f:
        return pd.read_csv(f, **kwargs)


__all__ = [
    "resolve_fs",
    "open_file",
    "with_local_copy",
    "smart_read_h5ad",
    "smart_write_checkpoint",
    "smart_read_csv",
]
