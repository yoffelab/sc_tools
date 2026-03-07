"""Strategy 2: DNA-only segmentation pipeline.

Extracts DNA channels from multi-channel IMC TIFFs, generates probability
maps (no Ilastik needed), and runs DL segmentation models.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

__all__ = [
    "generate_dna_probability_map",
    "run_strategy2_single",
    "run_all_strategy2",
]

logger = logging.getLogger(__name__)


def generate_dna_probability_map(
    dna_image: np.ndarray,
    method: str = "gaussian",
    sigma: float = 2.0,
) -> np.ndarray:
    """Generate a probability map from DNA channels only (no Ilastik).

    Parameters
    ----------
    dna_image
        2D DNA channel image (H, W), raw ion counts.
    method
        ``"gaussian"``, ``"otsu"``, or ``"multiscale"``.
    sigma
        Gaussian sigma for smoothing.

    Returns
    -------
    Probability map (H, W, 3): [background, nucleus, cytoplasm].
    """
    from sc_tools.data.imc.benchmark.prepare import generate_probability_map

    return generate_probability_map(dna_image, method=method, sigma=sigma)


def run_strategy2_single(
    tiff_path: str | Path,
    method: str,
    channel_csv_path: str | Path | None = None,
    panel_mapper=None,
    prob_map_method: str = "gaussian",
    gpu: bool = False,
    **method_kwargs,
) -> np.ndarray:
    """Run Strategy 2 on a single ROI: extract DNA -> prob map -> DL model.

    Parameters
    ----------
    tiff_path
        Path to ``*_full.tiff``.
    method
        Segmentation method: ``"cellpose_cyto2"``, ``"cellpose_nuclei"``,
        ``"stardist"``, ``"deepcell"``.
    channel_csv_path
        Path to ``*_full.csv``.
    panel_mapper
        Optional ``IMCPanelMapper`` instance.
    prob_map_method
        How to generate probability map from DNA.
    gpu
        Whether to use GPU.
    **method_kwargs
        Additional kwargs for the segmentation method.

    Returns
    -------
    Labeled mask (H, W).
    """
    from sc_tools.data.imc.benchmark.prepare import extract_dna_channels

    # Extract DNA signal
    dna = extract_dna_channels(tiff_path, channel_csv_path, panel_mapper)
    logger.info("DNA image: shape=%s, range=[%.1f, %.1f]", dna.shape, dna.min(), dna.max())

    # Special case: Cellpose nuclei model can run directly on DNA
    if method == "cellpose_nuclei":
        from sc_tools.bm.segment import run_cellpose

        return run_cellpose(
            dna,
            model_type="nuclei",
            gpu=gpu,
            **method_kwargs,
        )

    # Generate probability map for other methods
    prob_map = generate_dna_probability_map(dna, method=prob_map_method)

    return _run_method_on_probmap(prob_map, method, gpu=gpu, **method_kwargs)


def run_all_strategy2(
    tiff_path: str | Path,
    methods: list[str] | None = None,
    channel_csv_path: str | Path | None = None,
    panel_mapper=None,
    prob_map_method: str = "gaussian",
    gpu: bool = False,
) -> dict[str, np.ndarray]:
    """Run all Strategy 2 methods on a single ROI.

    Parameters
    ----------
    tiff_path
        Path to ``*_full.tiff``.
    methods
        Methods to run. Default: all available.
    channel_csv_path
        Path to ``*_full.csv``.
    panel_mapper
        Optional ``IMCPanelMapper``.
    prob_map_method
        How to generate probability map.
    gpu
        Whether to use GPU.

    Returns
    -------
    Dict mapping method name to labeled mask.
    """
    if methods is None:
        methods = ["cellpose_cyto2", "cellpose_nuclei", "stardist", "deepcell"]

    results = {}
    for method in methods:
        try:
            mask = run_strategy2_single(
                tiff_path,
                method,
                channel_csv_path=channel_csv_path,
                panel_mapper=panel_mapper,
                prob_map_method=prob_map_method,
                gpu=gpu,
            )
            results[f"s2_{method}"] = mask
            logger.info("Strategy 2 / %s: %d cells", method, len(np.unique(mask)) - 1)
        except Exception as e:
            logger.error("Strategy 2 / %s failed: %s", method, e)

    return results


def _run_method_on_probmap(
    prob_map: np.ndarray,
    method: str,
    gpu: bool = False,
    **kwargs,
) -> np.ndarray:
    """Dispatch a DL segmentation method on a probability map."""
    if method.startswith("cellpose"):
        from sc_tools.bm.segment import run_cellpose

        model_type = method.replace("cellpose_", "") if "_" in method else "cyto2"
        return run_cellpose(prob_map, model_type=model_type, gpu=gpu, **kwargs)

    elif method == "stardist":
        from sc_tools.bm.segment import run_stardist

        return run_stardist(prob_map, **kwargs)

    elif method == "deepcell":
        from sc_tools.bm.deepcell_runner import run_deepcell

        return run_deepcell(prob_map, **kwargs)

    else:
        raise ValueError(f"Unknown Strategy 2 method: {method!r}")
