"""
Spatial plotting utilities.

Generic helpers for spatial visualization of omics data (H&E image,
categorical and continuous overlays). Built on scanpy.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

__all__ = [
    "plot_spatial_plain_he",
    "plot_spatial_categorical",
    "plot_spatial_continuous",
    "multipage_spatial_pdf",
    "plot_imc_composite",
    "plot_imc_channel",
]


def plot_spatial_plain_he(
    adata,
    library_id: str,
    ax: plt.Axes,
    image_key: str = "hires",
) -> None:
    """
    Plot plain H&E tissue image for a library (no spots overlay).

    Parameters
    ----------
    adata : AnnData
        Full AnnData with adata.uns['spatial'][library_id]['images'][image_key].
    library_id : str
        Key in adata.uns['spatial'].
    ax : Axes
        Matplotlib axes to draw on.
    image_key : str
        Key in spatial['images'] (default 'hires').
    """
    try:
        if library_id not in adata.uns.get("spatial", {}):
            ax.text(
                0.5,
                0.5,
                f"No spatial data for library {library_id}",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title("H&E Tissue", fontsize=12, fontweight="bold")
            return

        spatial_data = adata.uns["spatial"][library_id]
        if "images" not in spatial_data or image_key not in spatial_data["images"]:
            ax.text(
                0.5,
                0.5,
                f"No H&E image found for library {library_id}",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title("H&E Tissue", fontsize=12, fontweight="bold")
            return

        img = spatial_data["images"][image_key]
        ax.imshow(img, aspect="auto")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("H&E Tissue", fontsize=12, fontweight="bold")
    except Exception as e:
        ax.text(
            0.5,
            0.5,
            f"Error loading H&E image:\n{str(e)}",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=10,
        )
        ax.set_title("H&E Tissue", fontsize=12, fontweight="bold")


def plot_spatial_categorical(
    adata,
    library_id: str,
    color: str,
    ax: plt.Axes,
    title: str | None = None,
    palette: dict[str, str] | None = None,
    legend_loc: str = "right margin",
    frameon: bool = False,
    **kwargs: Any,
) -> None:
    """
    Plot spatial overlay of a categorical variable (e.g. annotation, solidity).

    Parameters
    ----------
    adata : AnnData
        Subset AnnData for this library (e.g. adata[adata.obs['library_id'] == library_id]).
    library_id : str
        Key in adata.uns['spatial'].
    color : str
        Column name in adata.obs (categorical).
    ax : Axes
        Matplotlib axes.
    title : str, optional
        Axis title. If None, uses color.
    palette : dict, optional
        Category -> color mapping.
    legend_loc : str
        Passed to scanpy (default 'right margin').
    frameon : bool
        Passed to scanpy (default False).
    **kwargs
        Passed to sc.pl.spatial.
    """
    import scanpy as sc

    if color not in adata.obs.columns:
        ax.text(
            0.5,
            0.5,
            f"{color} not found",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_title(title or color, fontsize=12, fontweight="bold")
        return

    sc.pl.spatial(
        adata,
        color=color,
        library_id=library_id,
        frameon=frameon,
        show=False,
        ax=ax,
        legend_loc=legend_loc,
        palette=palette,
        **kwargs,
    )
    ax.set_title(title or color.replace("_", " ").title(), fontsize=12, fontweight="bold")


def plot_spatial_continuous(
    adata,
    library_id: str,
    color: str,
    ax: plt.Axes,
    title: str | None = None,
    cmap: str = "coolwarm",
    vmin: float | None = None,
    vmax: float | None = None,
    frameon: bool = False,
    values: pd.Series | np.ndarray | None = None,
    **kwargs: Any,
) -> None:
    """
    Plot spatial overlay of a continuous variable (e.g. score).

    Parameters
    ----------
    adata : AnnData
        Subset AnnData for this library.
    library_id : str
        Key in adata.uns['spatial'].
    color : str
        Column name in adata.obs (numeric). Ignored if values is provided.
    ax : Axes
        Matplotlib axes.
    title : str, optional
        Axis title. If None, uses color or "Score".
    cmap : str
        Colormap name (default 'coolwarm').
    vmin, vmax : float, optional
        Color scale limits.
    frameon : bool
        Passed to scanpy (default False).
    values : Series or ndarray, optional
        If provided, use these values for the overlay (length/index must match
        adata.obs_names). Use when scores are in obsm instead of obs.
    **kwargs
        Passed to sc.pl.spatial.
    """
    import scanpy as sc

    if values is not None:
        if isinstance(values, pd.Series):
            plot_values = values.reindex(adata.obs_names).values
        else:
            plot_values = np.asarray(values)
            if len(plot_values) != adata.n_obs:
                ax.text(
                    0.5,
                    0.5,
                    "values length mismatch",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(title or "Score", fontsize=12, fontweight="bold")
                return
        if color not in adata.obs.columns:
            # Temporarily add so scanpy can use it
            adata.obs["_st_continuous_plot"] = plot_values
            color_use = "_st_continuous_plot"
            cleanup = True
        else:
            color_use = color
            cleanup = False
    else:
        if color not in adata.obs.columns:
            ax.text(
                0.5,
                0.5,
                f"{color} not found",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title(title or color, fontsize=12, fontweight="bold")
            return
        color_use = color
        cleanup = False

    try:
        sc.pl.spatial(
            adata,
            color=color_use,
            library_id=library_id,
            frameon=frameon,
            show=False,
            ax=ax,
            cmap=cmap,
            colorbar_loc="right",
            vmin=vmin,
            vmax=vmax,
            **kwargs,
        )
        display_title = (
            (title or (color if color_use == color else "Score")).replace("_", " ").title()
        )
        ax.set_title(display_title, fontsize=12, fontweight="bold")
    finally:
        if cleanup and "_st_continuous_plot" in adata.obs.columns:
            adata.obs.drop(columns=["_st_continuous_plot"], inplace=True)


def plot_imc_composite(
    adata,
    library_id: str,
    ax: plt.Axes,
    image_key: str = "hires",
    title: str | None = None,
) -> None:
    """Plot IMC RGB composite image stored in ``adata.uns['spatial']``.

    Identical API to ``plot_spatial_plain_he`` — reuses the same
    ``adata.uns['spatial'][library_id]['images'][image_key]`` structure so
    that ``sc.pl.spatial(img_key='hires')`` also works.

    Parameters
    ----------
    adata
        AnnData with ``adata.uns['spatial'][library_id]['images'][image_key]``
        holding a ``(H, W, 3)`` uint8 RGB array.
    library_id
        Key in ``adata.uns['spatial']``.
    ax
        Matplotlib axes to draw on.
    image_key
        Key in ``spatial['images']`` (default ``'hires'``).
    title
        Axis title. If ``None``, shows channel info from metadata if available.
    """
    spatial_info = adata.uns.get("spatial", {}).get(library_id)
    if (
        spatial_info is None
        or "images" not in spatial_info
        or image_key not in spatial_info["images"]
    ):
        ax.text(
            0.5,
            0.5,
            f"No IMC composite image for library {library_id}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_title(title or "IMC Composite", fontsize=12, fontweight="bold")
        return

    img = spatial_info["images"][image_key]
    ax.imshow(img, aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])

    if title is None:
        rgb = spatial_info.get("metadata", {}).get("rgb_channels", {})
        if rgb:
            label = f"R={rgb.get('R', '?')} G={rgb.get('G', '?')} B={rgb.get('B', '?')}"
        else:
            label = "IMC Composite"
        ax.set_title(label, fontsize=12, fontweight="bold")
    else:
        ax.set_title(title, fontsize=12, fontweight="bold")


def plot_imc_channel(
    adata,
    library_id: str,
    channel: str,
    ax: plt.Axes,
    *,
    cmap: str = "inferno",
    vmax_percentile: float = 99,
    title: str | None = None,
) -> None:
    """Plot a single IMC channel from the full arcsinh-normalized stack.

    Reads ``adata.uns['spatial'][library_id]['images']['full']`` (shape
    ``(C, H, W)``) and ``metadata['channels']`` to look up the channel index.

    Parameters
    ----------
    adata
        AnnData with IMC image data in ``adata.uns['spatial']``.
    library_id
        Key in ``adata.uns['spatial']``.
    channel
        Marker/channel name (resolved via case-insensitive substring match
        against ``metadata['channels']``).
    ax
        Matplotlib axes to draw on.
    cmap
        Colormap (default ``'inferno'``).
    vmax_percentile
        Percentile used for the upper color scale limit (default 99).
    title
        Axis title. Defaults to the channel name.
    """
    spatial_info = adata.uns.get("spatial", {}).get(library_id)
    if spatial_info is None:
        ax.text(
            0.5,
            0.5,
            f"No spatial data for {library_id}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_title(title or channel, fontsize=12, fontweight="bold")
        return

    full = spatial_info.get("images", {}).get("full")
    channels = spatial_info.get("metadata", {}).get("channels", [])

    if full is None:
        ax.text(
            0.5,
            0.5,
            "No full channel stack (images['full']) found",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_title(title or channel, fontsize=12, fontweight="bold")
        return

    # Resolve channel index
    ch_lower = [c.lower() for c in channels]
    lo = channel.lower()
    idx = None
    if lo in ch_lower:
        idx = ch_lower.index(lo)
    else:
        # Partial match
        matches = [i for i, c in enumerate(ch_lower) if lo in c or c in lo]
        if matches:
            idx = matches[0]

    if idx is None or idx >= full.shape[0]:
        ax.text(
            0.5,
            0.5,
            f"Channel {channel!r} not found",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_title(title or channel, fontsize=12, fontweight="bold")
        return

    img_ch = full[idx]
    vmax = float(np.percentile(img_ch, vmax_percentile))
    im = ax.imshow(img_ch, cmap=cmap, vmin=0, vmax=vmax if vmax > 0 else 1, aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title(title or (channels[idx] if channels else channel), fontsize=12, fontweight="bold")


def multipage_spatial_pdf(
    adata,
    library_id_col: str,
    panels: list[dict],
    output_path: str,
    figsize: tuple[float, float] = (18, 12),
    dpi: int = 300,
) -> None:
    """
    Create a multipage PDF with one page per library and N spatial panels per page.

    Parameters
    ----------
    adata : AnnData
        Full AnnData with obs[library_id_col], uns['spatial'], and any obs columns
        required by the panels.
    library_id_col : str
        Column in adata.obs that identifies the library/sample.
    panels : list of dict
        List of panel specs. Each dict must have a ``"type"`` key
        (``"he"``, ``"categorical"``, or ``"continuous"``).
        ``"he"`` needs no extra keys. ``"categorical"`` needs ``"obs_col"`` and
        ``"title"`` (optional ``"palette"``). ``"continuous"`` needs ``"title"``
        and either ``"obs_col"`` or ``"values"`` (optional ``"cmap"``,
        ``"vmin"``, ``"vmax"``).
    output_path : str
        Path to the output PDF file.
    figsize : tuple
        Figure size per page (default (18, 12)).
    dpi : int
        DPI for saved pages (default 300).
    """
    import os

    library_ids = sorted(adata.obs[library_id_col].dropna().unique())
    n_panels = len(panels)
    n_rows = 2
    n_cols = 3
    if n_panels > n_rows * n_cols:
        n_cols = (n_panels + n_rows - 1) // n_rows

    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with PdfPages(output_path) as pdf:
        for lib_id in library_ids:
            adata_sub = adata[adata.obs[library_id_col] == lib_id].copy()
            if adata_sub.n_obs == 0:
                continue

            fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
            axes = np.atleast_1d(axes).flatten()

            for idx, spec in enumerate(panels):
                if idx >= len(axes):
                    break
                ax = axes[idx]
                ptype = spec.get("type")

                if ptype == "he":
                    plot_spatial_plain_he(
                        adata,
                        lib_id,
                        ax,
                        image_key=spec.get("image_key", "hires"),
                    )
                elif ptype == "categorical":
                    plot_spatial_categorical(
                        adata_sub,
                        lib_id,
                        spec["obs_col"],
                        ax,
                        title=spec.get("title"),
                        palette=spec.get("palette"),
                    )
                elif ptype == "continuous":
                    values = spec.get("values")
                    obs_col = spec.get("obs_col", "")
                    vals_sub = values.reindex(adata_sub.obs_names) if values is not None else None
                    plot_spatial_continuous(
                        adata_sub,
                        lib_id,
                        obs_col or "_",
                        ax,
                        title=spec.get("title"),
                        cmap=spec.get("cmap", "coolwarm"),
                        vmin=spec.get("vmin"),
                        vmax=spec.get("vmax"),
                        values=vals_sub,
                    )
                else:
                    ax.text(
                        0.5,
                        0.5,
                        f"Unknown panel type: {ptype}",
                        ha="center",
                        va="center",
                        transform=ax.transAxes,
                    )

            for j in range(len(panels), len(axes)):
                axes[j].set_visible(False)

            fig.suptitle(f"Library: {lib_id}", fontsize=16, fontweight="bold", y=0.995)
            pdf.savefig(fig, bbox_inches="tight", dpi=dpi)
            plt.close(fig)
