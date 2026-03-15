"""Pre-tuned configuration profiles for semi-supervised integration methods.

Provides ``get_scanvi_config()`` and ``get_resolvi_ss_config()`` which return
dictionaries of model and training parameters optimised for spatial omics
data.  Each function accepts a *profile* argument (``"local"`` or ``"hpc"``)
and arbitrary keyword overrides.

Findings are based on cross-platform spatial benchmarking (CosMx + Xenium,
IBD dataset, 2026-03).  The defaults aim to improve bio conservation over
scvi-tools stock values while remaining safe for general use.

Example
-------
::

    from sc_tools.pp.integration_configs import get_scanvi_config

    cfg = get_scanvi_config(profile="hpc", max_epochs_finetune=150)
    # pass cfg entries to run_scanvi or scvi.model.SCANVI

"""

from __future__ import annotations

from copy import deepcopy
from typing import Any

__all__ = [
    "get_scanvi_config",
    "get_resolvi_ss_config",
]

# ---------------------------------------------------------------------------
# scANVI defaults
# ---------------------------------------------------------------------------
# Why each value was chosen (benchmark: CosMx+Xenium IBD, 2026-03):
#
# gene_likelihood="nb"
#   Spatial platforms (CosMx, Xenium, IMC) do not have droplet-based zero
#   inflation.  Negative binomial fits the count distribution better than
#   zero-inflated NB and avoids the extra zero-gate parameters.
#
# classification_ratio=100
#   Doubling the default (50) increases the weight of the classification
#   loss relative to reconstruction, improving celltype separation without
#   over-classifying.
#
# dispersion="gene-batch"
#   Cross-platform data have platform-specific noise profiles.  Estimating
#   a per-gene per-batch dispersion lets the model absorb platform effects
#   that a single global dispersion cannot.
#
# n_latent=30
#   Multi-platform integration with 8+ cell types benefits from a larger
#   latent space.  Values in the 10-20 range lose resolution for rare types.
#
# n_hidden / n_layers
#   Random search (40 configs, 2026-03-13) showed 4 layers is critical for
#   bio conservation.  All top-5 configs used n_layers=4.  2-layer configs
#   clustered below ct_broad_asw=0.18.  HPC profile uses 4 layers + 512
#   hidden (R023: ct_broad_asw=0.396).  Local profile uses 2 layers + 128
#   hidden to stay within 16 GB RAM.
#
# dropout_rate=0.20
#   Higher dropout (0.15-0.20) prevents platform-specific memorization
#   and improves cross-platform generalization.
#
# classification_ratio=192
#   Random search confirmed higher ratios (120-192) push the model to
#   respect celltype structure more aggressively.
#
# kl_warmup_epochs
#   The scvi-tools default warmup (400 steps) is too long for short
#   fine-tuning runs.  We set warmup to ~1/4 of the fine-tuning epochs so
#   the KL term ramps up promptly.

_SCANVI_LOCAL: dict[str, Any] = {
    # --- model ---
    "gene_likelihood": "nb",
    "dispersion": "gene-batch",
    "n_latent": 30,
    "n_hidden": 128,
    "n_layers": 2,
    "dropout_rate": 0.15,
    # --- training (scVI base) ---
    "max_epochs_scvi": 100,
    "early_stopping_scvi": True,
    "batch_size": 256,
    # --- training (scANVI fine-tune) ---
    "max_epochs_finetune": 50,
    "classification_ratio": 100,
    # KL warmup scaled to ~1/4 of fine-tune epochs
    "kl_warmup_epochs": 12,
}

_SCANVI_HPC: dict[str, Any] = {
    # --- model ---
    "gene_likelihood": "nb",
    "dispersion": "gene-batch",
    "n_latent": 53,
    "n_hidden": 512,
    "n_layers": 4,
    "dropout_rate": 0.20,
    # --- training (scVI base) ---
    "max_epochs_scvi": 125,
    "early_stopping_scvi": True,
    "batch_size": 512,
    # --- training (scANVI fine-tune) ---
    "max_epochs_finetune": 30,
    "classification_ratio": 192,
    # KL warmup scaled to ~1/4 of fine-tune epochs
    "kl_warmup_epochs": 8,
}

# ---------------------------------------------------------------------------
# resolVI-SS defaults
# ---------------------------------------------------------------------------
# Why each value was chosen (same benchmark):
#
# prior_diffusion_amount=0.1
#   The stock default (0.3) over-smooths celltype boundaries in tissue.
#   Reducing to 0.1 preserves sharp interfaces while still using spatial
#   context for denoising.
#
# n_hidden=128  (decoder)
#   The stock decoder (32 units) is too small for datasets with many genes
#   and diverse cell types, leading to underfitting.
#
# n_hidden_encoder=256
#   A larger encoder gives the model more capacity to learn platform-
#   specific noise patterns.  The local profile uses 128 for memory.
#
# lr_extra=1e-3
#   The stock value (1e-2) causes diffusion parameters to dominate the
#   optimisation early, producing over-smoothed latent spaces.
#
# n_neighbors=5
#   Fewer spatial neighbors (vs. default 10) reduce over-smoothing,
#   especially in regions with mixed cell types at tissue boundaries.
#
# sparsity_diffusion=5.0
#   Higher sparsity (vs. default 3.0) concentrates diffusion on the
#   nearest neighbours, further limiting unwanted smoothing.

_RESOLVI_SS_LOCAL: dict[str, Any] = {
    # --- model ---
    "n_latent": 10,
    "n_hidden": 64,  # decoder; smaller for CPU
    "n_hidden_encoder": 128,
    "n_layers": 2,
    # --- spatial diffusion ---
    "prior_diffusion_amount": 0.1,
    "n_neighbors": 5,
    "sparsity_diffusion": 5.0,
    # --- training ---
    "max_epochs": 100,
    "batch_size": 256,
    "lr_extra": 1e-3,
}

_RESOLVI_SS_HPC: dict[str, Any] = {
    # --- model ---
    "n_latent": 10,
    "n_hidden": 128,  # decoder
    "n_hidden_encoder": 256,
    "n_layers": 2,
    # --- spatial diffusion ---
    "prior_diffusion_amount": 0.1,
    "n_neighbors": 5,
    "sparsity_diffusion": 5.0,
    # --- training ---
    "max_epochs": 400,
    "batch_size": 512,
    "lr_extra": 1e-3,
}

_PROFILES = {"local", "hpc"}


def _merge(base: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
    """Deep-copy *base* and apply *overrides*."""
    out = deepcopy(base)
    out.update(overrides)
    return out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def get_scanvi_config(profile: str = "local", **overrides: Any) -> dict[str, Any]:
    """Return a scANVI parameter dictionary tuned for spatial omics data.

    Parameters
    ----------
    profile
        ``"local"`` (default) — CPU-friendly, fits in 16 GB RAM.
        ``"hpc"`` — larger architecture and more epochs for GPU training.
    **overrides
        Any key present in the returned dict can be overridden.
        Unknown keys are passed through (useful for forwarding to scvi-tools).

    Returns
    -------
    dict
        Keys map to ``scvi.model.SCVI`` / ``scvi.model.SCANVI`` constructor
        and ``.train()`` arguments.  Callers should split model vs. training
        keys as needed.

    Notes
    -----
    Defaults are based on cross-platform spatial integration benchmarking
    (CosMx + Xenium, IBD dataset).  Key changes from scvi-tools stock:

    * ``gene_likelihood="nb"`` — spatial data lack droplet zero-inflation.
    * ``classification_ratio=192`` (HPC) / ``100`` (local) — stronger
      celltype supervision.
    * ``dispersion="gene-batch"`` — absorbs per-platform noise.
    * ``n_layers=4`` (HPC) — depth is the most impactful architecture
      parameter; all top-5 random search configs used 4 layers.
    * ``dropout_rate=0.20`` (HPC) — prevents platform-specific memorization.
    * Shorter ``kl_warmup_epochs`` scaled to fine-tune length.
    """
    if profile not in _PROFILES:
        msg = f"Unknown profile {profile!r}; choose from {sorted(_PROFILES)}"
        raise ValueError(msg)

    base = _SCANVI_LOCAL if profile == "local" else _SCANVI_HPC
    return _merge(base, overrides)


def get_resolvi_ss_config(profile: str = "local", **overrides: Any) -> dict[str, Any]:
    """Return a resolVI-SS parameter dictionary tuned for spatial omics data.

    Parameters
    ----------
    profile
        ``"local"`` (default) — CPU-friendly, fits in 16 GB RAM.
        ``"hpc"`` — larger architecture and more epochs for GPU training.
    **overrides
        Any key present in the returned dict can be overridden.
        Unknown keys are passed through (useful for forwarding to scvi-tools).

    Returns
    -------
    dict
        Keys map to ``scvi.external.RESOLVI`` constructor and ``.train()``
        arguments.  Callers should split model vs. training keys as needed.

    Notes
    -----
    Defaults are based on cross-platform spatial integration benchmarking
    (CosMx + Xenium, IBD dataset).  Key changes from scvi-tools stock:

    * ``prior_diffusion_amount=0.1`` — less smoothing preserves celltype
      boundaries (stock 0.3).
    * ``n_hidden=128`` decoder — stock 32 is too small for diverse datasets.
    * ``n_hidden_encoder=256`` — larger encoder for cross-platform noise.
    * ``lr_extra=1e-3`` — prevents diffusion params from dominating (stock
      1e-2).
    * ``n_neighbors=5`` — fewer neighbors = less over-smoothing (stock 10).
    * ``sparsity_diffusion=5.0`` — sparser weights (stock 3.0).
    """
    if profile not in _PROFILES:
        msg = f"Unknown profile {profile!r}; choose from {sorted(_PROFILES)}"
        raise ValueError(msg)

    base = _RESOLVI_SS_LOCAL if profile == "local" else _RESOLVI_SS_HPC
    return _merge(base, overrides)
