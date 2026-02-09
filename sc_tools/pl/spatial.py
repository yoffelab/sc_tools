"""
Spatial plotting utilities.

Generic helpers for spatial visualization of omics data (H&E image,
categorical and continuous overlays). Built on scanpy.
"""

from typing import Optional, List, Dict, Any
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

__all__ = [
    'plot_spatial_plain_he',
    'plot_spatial_categorical',
    'plot_spatial_continuous',
]


def plot_spatial_plain_he(
    adata,
    library_id: str,
    ax: plt.Axes,
    image_key: str = 'hires',
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
        if library_id not in adata.uns.get('spatial', {}):
            ax.text(
                0.5, 0.5, f'No spatial data for library {library_id}',
                ha='center', va='center', transform=ax.transAxes,
            )
            ax.set_title('H&E Tissue', fontsize=12, fontweight='bold')
            return

        spatial_data = adata.uns['spatial'][library_id]
        if 'images' not in spatial_data or image_key not in spatial_data['images']:
            ax.text(
                0.5, 0.5, f'No H&E image found for library {library_id}',
                ha='center', va='center', transform=ax.transAxes,
            )
            ax.set_title('H&E Tissue', fontsize=12, fontweight='bold')
            return

        img = spatial_data['images'][image_key]
        ax.imshow(img, aspect='auto')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('H&E Tissue', fontsize=12, fontweight='bold')
    except Exception as e:
        ax.text(
            0.5, 0.5, f'Error loading H&E image:\n{str(e)}',
            ha='center', va='center', transform=ax.transAxes, fontsize=10,
        )
        ax.set_title('H&E Tissue', fontsize=12, fontweight='bold')


def plot_spatial_categorical(
    adata,
    library_id: str,
    color: str,
    ax: plt.Axes,
    title: Optional[str] = None,
    palette: Optional[Dict[str, str]] = None,
    legend_loc: str = 'right margin',
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
            0.5, 0.5, f'{color} not found',
            ha='center', va='center', transform=ax.transAxes,
        )
        ax.set_title(title or color, fontsize=12, fontweight='bold')
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
    ax.set_title(title or color.replace('_', ' ').title(), fontsize=12, fontweight='bold')


def plot_spatial_continuous(
    adata,
    library_id: str,
    color: str,
    ax: plt.Axes,
    title: Optional[str] = None,
    cmap: str = 'coolwarm',
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    frameon: bool = False,
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
        Column name in adata.obs (numeric).
    ax : Axes
        Matplotlib axes.
    title : str, optional
        Axis title. If None, uses color.
    cmap : str
        Colormap name (default 'coolwarm').
    vmin, vmax : float, optional
        Color scale limits.
    frameon : bool
        Passed to scanpy (default False).
    **kwargs
        Passed to sc.pl.spatial.
    """
    import scanpy as sc

    if color not in adata.obs.columns:
        ax.text(
            0.5, 0.5, f'{color} not found',
            ha='center', va='center', transform=ax.transAxes,
        )
        ax.set_title(title or color, fontsize=12, fontweight='bold')
        return

    sc.pl.spatial(
        adata,
        color=color,
        library_id=library_id,
        frameon=frameon,
        show=False,
        ax=ax,
        cmap=cmap,
        colorbar_loc='right',
        vmin=vmin,
        vmax=vmax,
        **kwargs,
    )
    display_title = (title or color).replace('_', ' ').title()
    ax.set_title(display_title, fontsize=12, fontweight='bold')
