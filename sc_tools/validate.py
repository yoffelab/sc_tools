"""Checkpoint validation for sc_tools pipeline phases.

Validates AnnData checkpoints against the metadata contracts defined in
Architecture.md Section 2.2. Each phase has required obs/obsm/uns keys
and structural requirements.

Usage:
    from sc_tools.validate import validate_checkpoint, validate_file

    # Validate in-memory AnnData
    issues = validate_checkpoint(adata, phase="p1")

    # Validate from file
    issues = validate_file("results/adata.raw.p1.h5ad", phase="p1")
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np

logger = logging.getLogger(__name__)

# Default obs columns that don't count as "clinical metadata" for p2
_DEFAULT_OBS_COLS = frozenset(
    {
        "sample",
        "raw_data_dir",
        "batch",
        "library_id",
        "n_genes",
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_hb",
        "n_cells",
        "n_counts",
    }
)

# Recognized integration representations for p3
_INTEGRATION_REPS = (
    "X_scvi",
    "X_pca_harmony",
    "X_cytovi",
    "X_pca",
)

# Recognized cluster columns for p3
_CLUSTER_COLS = ("leiden", "cluster", "louvain")


class CheckpointValidationError(Exception):
    """Raised when a checkpoint fails validation in strict mode."""


def validate_checkpoint(
    adata: ad.AnnData,
    phase: str,
    *,
    strict: bool = True,
    fix: bool = False,
) -> list[str]:
    """Validate an AnnData checkpoint against phase requirements.

    Parameters
    ----------
    adata
        AnnData object to validate.
    phase
        Phase identifier: "p1", "p2", "p3", "p35", or "p4".
    strict
        If True, raise CheckpointValidationError on failures.
    fix
        If True, attempt auto-fixes (e.g., rename batch -> raw_data_dir).

    Returns
    -------
    List of validation issue messages. Empty list means valid.

    Raises
    ------
    CheckpointValidationError
        If strict=True and validation fails.
    ValueError
        If phase is not recognized.
    """
    validators = {
        "p1": validate_p1,
        "p2": validate_p2,
        "p3": validate_p3,
        "p35": validate_p35,
        "p4": validate_p4,
    }
    if phase not in validators:
        raise ValueError(f"Unknown phase '{phase}'. Must be one of: {list(validators.keys())}")

    issues = validators[phase](adata, fix=fix)

    if issues:
        for issue in issues:
            logger.warning("Checkpoint %s: %s", phase, issue)
        if strict:
            msg = f"Checkpoint {phase} validation failed with {len(issues)} issue(s):\n"
            msg += "\n".join(f"  - {i}" for i in issues)
            raise CheckpointValidationError(msg)

    return issues


def validate_p1(adata: ad.AnnData, *, fix: bool = False) -> list[str]:
    """Validate Phase 1 checkpoint (raw ingested data).

    Required:
    - obs['sample']
    - obs['raw_data_dir'] (fix: rename from 'batch' if present)
    - obsm['spatial']
    - X contains non-negative values (raw counts)
    """
    issues = []

    # obs['sample']
    if "sample" not in adata.obs.columns:
        issues.append("obs['sample'] is missing")

    # obs['raw_data_dir'] — auto-fix from 'batch'
    if "raw_data_dir" not in adata.obs.columns:
        if fix and "batch" in adata.obs.columns:
            adata.obs["raw_data_dir"] = adata.obs["batch"]
            logger.info("Auto-fix: copied obs['batch'] to obs['raw_data_dir']")
        else:
            hint = " (has 'batch' — use fix=True to rename)" if "batch" in adata.obs.columns else ""
            issues.append(f"obs['raw_data_dir'] is missing{hint}")

    # obsm['spatial']
    if "spatial" not in adata.obsm:
        issues.append("obsm['spatial'] is missing")

    # X non-negative (raw counts)
    if adata.X is not None:
        try:
            from scipy.sparse import issparse

            x_data = adata.X
            if issparse(x_data):
                min_val = x_data.min()
            else:
                min_val = np.nanmin(x_data)
            if min_val < 0:
                issues.append(
                    f"X contains negative values (min={min_val:.4f}); expected raw counts >= 0"
                )
        except Exception as e:
            issues.append(f"Could not check X values: {e}")
    else:
        issues.append("X is None")

    return issues


def validate_p2(adata: ad.AnnData, *, fix: bool = False) -> list[str]:
    """Validate Phase 2 checkpoint (metadata-annotated).

    Required:
    - All of p1
    - At least one clinical metadata column beyond default obs columns
    """
    issues = validate_p1(adata, fix=fix)

    # Check for clinical metadata columns
    clinical_cols = set(adata.obs.columns) - _DEFAULT_OBS_COLS
    if not clinical_cols:
        issues.append(
            "No clinical metadata columns found in obs beyond defaults. "
            "Expected columns from sample_metadata.csv (e.g., diagnosis, age, sex)"
        )

    return issues


def validate_p3(adata: ad.AnnData, *, fix: bool = False) -> list[str]:
    """Validate Phase 3 checkpoint (preprocessed).

    Required:
    - An integration/reduction representation in obsm
      (X_scvi, X_pca_harmony, X_cytovi, or X_pca)
    - A cluster column in obs (leiden, cluster, or louvain)
    - adata.raw is not None (backup of raw counts)
    """
    issues = []

    # Integration representation
    found_reps = [r for r in _INTEGRATION_REPS if r in adata.obsm]
    if not found_reps:
        issues.append(
            f"No integration representation found in obsm. "
            f"Expected one of: {list(_INTEGRATION_REPS)}"
        )

    # Cluster column
    found_clusters = [c for c in _CLUSTER_COLS if c in adata.obs.columns]
    if not found_clusters:
        issues.append(f"No cluster column found in obs. Expected one of: {list(_CLUSTER_COLS)}")

    # adata.raw backup
    if adata.raw is None:
        issues.append("adata.raw is None; expected backup of raw counts")

    return issues


def validate_p35(adata: ad.AnnData, *, fix: bool = False) -> list[str]:
    """Validate Phase 3.5b checkpoint (scored).

    Required:
    - obsm['signature_score']
    - obsm['signature_score_z']
    - uns['signature_score_report']
    """
    issues = []

    if "signature_score" not in adata.obsm:
        issues.append("obsm['signature_score'] is missing")

    if "signature_score_z" not in adata.obsm:
        issues.append("obsm['signature_score_z'] is missing")

    if "signature_score_report" not in adata.uns:
        issues.append("uns['signature_score_report'] is missing")

    return issues


def validate_p4(adata: ad.AnnData, *, fix: bool = False) -> list[str]:
    """Validate Phase 4 checkpoint (cell-typed).

    Required:
    - All of p35
    - obs['celltype']
    - obs['celltype_broad']
    """
    issues = validate_p35(adata, fix=fix)

    if "celltype" not in adata.obs.columns:
        issues.append("obs['celltype'] is missing")

    if "celltype_broad" not in adata.obs.columns:
        issues.append("obs['celltype_broad'] is missing")

    return issues


def validate_file(
    path: str | Path,
    phase: str,
    *,
    strict: bool = True,
    fix: bool = False,
) -> list[str]:
    """Load an h5ad file and validate against phase requirements.

    Parameters
    ----------
    path
        Path to .h5ad file.
    phase
        Phase identifier: "p1", "p2", "p3", "p35", or "p4".
    strict
        If True, raise CheckpointValidationError on failures.
    fix
        If True, attempt auto-fixes and re-save the file.

    Returns
    -------
    List of validation issue messages.
    """
    path = Path(path)
    if not path.exists():
        msg = f"File not found: {path}"
        if strict:
            raise CheckpointValidationError(msg)
        return [msg]

    adata = ad.read_h5ad(path)
    issues = validate_checkpoint(adata, phase, strict=strict, fix=fix)

    if fix and not issues:
        # Re-save only if fixes were applied (no issues remaining)
        adata.write_h5ad(path)
        logger.info("Re-saved fixed checkpoint to %s", path)

    return issues
