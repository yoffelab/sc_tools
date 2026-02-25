"""
Versioned save utilities: datetime-prefixed filenames for figures and h5ad.

Format: YYDDMM.hh.mm.filename.extension (e.g. 250208.14.30.volcano_faceted.pdf)
"""

from datetime import datetime
from pathlib import Path

# Default DPI for raster figure exports
DEFAULT_FIGURE_DPI = 300


def get_version_prefix(dt: datetime | None = None) -> str:
    """
    Return a version prefix string for filenames: YYDDMM.hh.mm.

    Parameters
    ----------
    dt : datetime, optional
        Use this datetime; if None, uses now().

    Returns
    -------
    str
        e.g. "250208.14.30" for 2025-02-08 14:30
    """
    if dt is None:
        dt = datetime.now()
    return dt.strftime("%y%d%m.%H.%M")


def versioned_filename(basename: str, ext: str, dt: datetime | None = None) -> str:
    """
    Build a versioned filename: YYDDMM.hh.mm.basename.ext (no path).

    Parameters
    ----------
    basename : str
        Base name (no extension), e.g. "volcano_faceted"
    ext : str
        Extension without dot, e.g. "pdf" or "png"
    dt : datetime, optional
        Use this datetime; if None, uses now().

    Returns
    -------
    str
        e.g. "250208.14.30.volcano_faceted.pdf"
    """
    prefix = get_version_prefix(dt)
    ext = ext.lstrip(".")
    return f"{prefix}.{basename}.{ext}"


def versioned_path(
    base_dir: str | Path,
    basename: str,
    ext: str,
    subdir: str | None = None,
    dt: datetime | None = None,
) -> Path:
    """
    Build full path for a versioned file, optionally under a subdir (e.g. pdf/, png/).

    Parameters
    ----------
    base_dir : str or Path
        Root output directory (e.g. figures/process_colocalization)
    basename : str
        Base name (no extension)
    ext : str
        Extension without dot (e.g. "pdf", "png")
    subdir : str, optional
        Subdirectory under base_dir (e.g. "pdf", "png"). Created if missing.
    dt : datetime, optional
        Use this datetime; if None, uses now().

    Returns
    -------
    Path
        e.g. base_dir/pdf/250208.14.30.volcano_faceted.pdf
    """
    base_dir = Path(base_dir)
    if subdir:
        base_dir = base_dir / subdir
    base_dir.mkdir(parents=True, exist_ok=True)
    fname = versioned_filename(basename, ext, dt=dt)
    return base_dir / fname
