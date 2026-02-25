"""
GPU detection and management utilities.
"""

import logging

logger = logging.getLogger(__name__)


def check_gpu_available() -> bool:
    """
    Check if GPU is available for PyTorch.

    Returns
    -------
    bool
        True if GPU is available, False otherwise
    """
    try:
        import torch

        if torch.cuda.is_available():
            device_name = torch.cuda.get_device_name(0)
            device_props = torch.cuda.get_device_properties(0)
            total_memory_gb = device_props.total_memory / 1e9

            logger.info(f"   ✅ GPU available: {device_name}")
            logger.info(f"   GPU memory: {total_memory_gb:.2f} GB")
            return True
        else:
            logger.info("   ⚠️  GPU not available, using CPU")
            return False
    except ImportError:
        logger.info("   ⚠️  PyTorch not available, using CPU")
        return False
    except Exception as e:
        logger.warning(f"   ⚠️  Error checking GPU: {e}, using CPU")
        return False


def get_gpu_setting(use_gpu: bool | None = None) -> bool:
    """
    Determine GPU usage setting with auto-detection.

    Parameters
    ----------
    use_gpu : bool, optional
        Explicit GPU setting:
        - None: Auto-detect GPU availability
        - True: Force GPU (with warning if not available)
        - False: Force CPU

    Returns
    -------
    bool
        True to use GPU, False to use CPU
    """
    if use_gpu is None:
        # Auto-detect
        return check_gpu_available()
    elif use_gpu:
        # User explicitly requested GPU, verify it's available
        if not check_gpu_available():
            logger.warning("   GPU requested but not available, falling back to CPU")
            return False
        return True
    else:
        # User explicitly requested CPU
        logger.info("   Using CPU (use_gpu=False)")
        return False
