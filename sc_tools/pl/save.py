"""
Versioned figure saving: PDF + PNG with datetime prefix, dpi=300, under pdf/ and png/.
"""

from pathlib import Path
from typing import Union, Optional, Tuple
from datetime import datetime

import matplotlib.pyplot as plt

from ..utils.save import (
    versioned_path,
    get_version_prefix,
    DEFAULT_FIGURE_DPI,
)


def save_figure(
    fig: plt.Figure,
    basename: str,
    output_dir: Union[str, Path],
    dpi: int = DEFAULT_FIGURE_DPI,
    dt: Optional[datetime] = None,
    bbox_inches: str = "tight",
    pad_inches: float = 0.1,
    create_pdf_png_folders: bool = True,
) -> Tuple[Path, Path]:
    """
    Save a figure in two formats with a versioned filename under pdf/ and png/.

    Creates:
      - output_dir/pdf/YYDDMM.hh.mm.basename.pdf
      - output_dir/png/YYDDMM.hh.mm.basename.png

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save.
    basename : str
        Base name (no extension), e.g. "volcano_faceted".
    output_dir : str or Path
        Root directory for figures (e.g. "figures/process_colocalization").
        Subdirs "pdf" and "png" are created inside it.
    dpi : int, optional
        DPI for raster (PNG). Default 300.
    dt : datetime, optional
        Timestamp for version prefix; if None, uses now().
    bbox_inches : str, optional
        Passed to savefig. Default "tight".
    pad_inches : float, optional
        Passed to savefig. Default 0.1.
    create_pdf_png_folders : bool, optional
        If True (default), save under output_dir/pdf/ and output_dir/png/.

    Returns
    -------
    tuple of Path
        (path_to_pdf, path_to_png)
    """
    output_dir = Path(output_dir)
    if create_pdf_png_folders:
        pdf_path = versioned_path(output_dir, basename, "pdf", subdir="pdf", dt=dt)
        png_path = versioned_path(output_dir, basename, "png", subdir="png", dt=dt)
    else:
        pdf_path = versioned_path(output_dir, basename, "pdf", subdir=None, dt=dt)
        png_path = versioned_path(output_dir, basename, "png", subdir=None, dt=dt)

    save_kw = {"bbox_inches": bbox_inches, "pad_inches": pad_inches}
    fig.savefig(pdf_path, **save_kw)
    fig.savefig(png_path, dpi=dpi, **save_kw)

    return pdf_path, png_path


__all__ = ["save_figure"]
