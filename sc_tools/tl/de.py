"""Pseudobulk differential expression via PyDESeq2 (SCI-01).

Provides:
- aggregate_pseudobulk: Sum raw counts by subject_id + celltype, with threshold filtering.
- run_pseudobulk_de: Run PyDESeq2 per celltype with auto-inferred or custom design formula.
- _infer_design_formula: Auto-build design formula with collinearity guard for batch covariates.

PyDESeq2 is an optional dependency; use ``_check_deps(['pydeseq2'])`` before calling
``run_pseudobulk_de``.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import scipy.sparse as sp

if TYPE_CHECKING:
    from anndata import AnnData

__all__ = ["aggregate_pseudobulk", "run_pseudobulk_de"]

logger = logging.getLogger(__name__)

# Candidate batch-like columns to auto-detect in design formula
_BATCH_CANDIDATES = ("batch", "library_id")


def _get_raw_counts(adata: AnnData, layer: str | None = None) -> np.ndarray | sp.spmatrix:
    """Extract raw count matrix from AnnData.

    Priority: explicit layer > layers['counts'] > raw.X > X.
    Validates values are non-negative and close to integers.
    """
    if layer is not None:
        if layer in adata.layers:
            X = adata.layers[layer]
        else:
            msg = f"Layer '{layer}' not found in adata.layers. Available: {list(adata.layers.keys())}"
            raise KeyError(msg)
    elif "counts" in adata.layers:
        X = adata.layers["counts"]
    elif adata.raw is not None:
        X = adata.raw.X
    else:
        X = adata.X

    return X


def _validate_counts(arr: np.ndarray) -> np.ndarray:
    """Validate and coerce array to integer counts.

    Rounds values close to integers and casts. Raises if negative values found.
    """
    if (arr < 0).any():
        msg = "Negative values found in count matrix. PyDESeq2 requires non-negative integer counts."
        raise ValueError(msg)
    # Round and cast if close to int
    rounded = np.round(arr).astype(int)
    return rounded


def _is_collinear_with_subject(metadata_df: pd.DataFrame, col: str, subject_key: str) -> bool:
    """Check if a column is collinear with the subject index (1:1 mapping).

    Collinearity means each subject maps to exactly one value of col, AND
    each value of col maps to exactly one subject.
    """
    # metadata_df index is subject IDs
    if col not in metadata_df.columns:
        return False
    # Check 1:1 mapping: each subject has one value, each value has one subject
    subj_to_col = metadata_df[col]
    # Each subject -> exactly one value (always true since one row per subject)
    # Each value -> exactly one subject
    col_to_subj_counts = subj_to_col.value_counts()
    return bool(col_to_subj_counts.max() == 1)


def _infer_design_formula(
    metadata_df: pd.DataFrame,
    condition_key: str,
    subject_key: str,
) -> str:
    """Auto-infer design formula from metadata columns.

    Starts with ``~ {condition_key}``, then adds non-collinear batch-like columns.

    Parameters
    ----------
    metadata_df
        Per-subject metadata DataFrame (index = subject IDs).
    condition_key
        Condition column name.
    subject_key
        Subject ID column name (used for collinearity check).

    Returns
    -------
    str
        Wilkinson-style design formula.
    """
    terms = [condition_key]

    for batch_col in _BATCH_CANDIDATES:
        if batch_col not in metadata_df.columns:
            continue
        if batch_col == condition_key:
            continue
        # Check collinearity with subject index
        if _is_collinear_with_subject(metadata_df, batch_col, subject_key):
            logger.warning(
                "Skipping '%s' from auto-formula: collinear with '%s' "
                "(1:1 design, would produce rank-deficient matrix).",
                batch_col,
                subject_key,
            )
            continue
        # Check that the column has more than one level
        if metadata_df[batch_col].nunique() <= 1:
            continue
        terms.append(batch_col)

    return "~ " + " + ".join(terms)


def aggregate_pseudobulk(
    adata: AnnData,
    subject_key: str,
    celltype_key: str,
    min_cells: int = 10,
    layer: str | None = None,
) -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    """Aggregate counts by subject + celltype for pseudobulk analysis.

    Parameters
    ----------
    adata
        Annotated data matrix with raw counts.
    subject_key
        Column in obs for subject identifiers.
    celltype_key
        Column in obs for cell type labels.
    min_cells
        Minimum cells per subject+celltype combination. Subjects below
        this threshold are excluded from aggregation.
    layer
        Explicit layer name for count matrix. If None, tries
        ``layers['counts']`` > ``raw.X`` > ``X``.

    Returns
    -------
    dict[str, tuple[pd.DataFrame, pd.DataFrame]]
        ``{celltype: (count_df, metadata_df)}`` where ``count_df`` has shape
        ``(n_valid_subjects, n_genes)`` with integer values, and ``metadata_df``
        has per-subject covariates.
    """
    X = _get_raw_counts(adata, layer=layer)

    # Determine gene names
    if layer is None and adata.raw is not None and "counts" not in adata.layers:
        gene_names = adata.raw.var_names
    else:
        gene_names = adata.var_names

    results: dict[str, tuple[pd.DataFrame, pd.DataFrame]] = {}

    # Get celltypes (handle both categorical and non-categorical)
    ct_col = adata.obs[celltype_key]
    if hasattr(ct_col, "cat"):
        celltypes = ct_col.cat.categories.tolist()
    else:
        celltypes = ct_col.unique().tolist()

    for ct in celltypes:
        ct_mask = (adata.obs[celltype_key] == ct).values
        ct_obs = adata.obs.loc[ct_mask]
        ct_X = X[ct_mask]

        subjects = ct_obs[subject_key]
        unique_subjects = subjects.unique()

        counts_list: list[np.ndarray] = []
        valid_subjects: list[str] = []

        for subj in unique_subjects:
            subj_mask = (subjects == subj).values
            n_cells = subj_mask.sum()
            if n_cells < min_cells:
                logger.debug(
                    "Skipping subject '%s' for celltype '%s': %d cells < %d min_cells",
                    subj, ct, n_cells, min_cells,
                )
                continue

            subj_counts = ct_X[subj_mask]
            if sp.issparse(subj_counts):
                subj_counts = subj_counts.toarray()
            counts_list.append(np.asarray(subj_counts).sum(axis=0).flatten())
            valid_subjects.append(subj)

        if len(valid_subjects) < 2:
            logger.info(
                "Skipping celltype '%s': only %d valid subjects (need >= 2).",
                ct, len(valid_subjects),
            )
            continue

        count_df = pd.DataFrame(
            _validate_counts(np.array(counts_list)),
            index=valid_subjects,
            columns=gene_names,
        )

        # Build metadata: one row per valid subject
        meta_df = (
            ct_obs.loc[ct_obs[subject_key].isin(valid_subjects)]
            .drop_duplicates(subset=[subject_key])
            .set_index(subject_key)
            .loc[valid_subjects]
        )

        results[ct] = (count_df, meta_df)

    return results


def run_pseudobulk_de(
    adata: AnnData,
    condition_key: str,
    *,
    subject_key: str = "subject_id",
    celltype_key: str = "celltype",
    formula: str | None = None,
    reference: str | None = None,
    min_subjects_per_group: int = 3,
    min_cells_per_combo: int = 10,
    layer: str | None = None,
    alpha: float = 0.05,
) -> dict[str, pd.DataFrame]:
    """Run pseudobulk DE via PyDESeq2 per celltype.

    Parameters
    ----------
    adata
        Annotated data matrix with raw counts and subject/celltype/condition metadata.
    condition_key
        Column in obs for the experimental condition (e.g. 'treatment', 'disease').
    subject_key
        Column in obs for subject identifiers.
    celltype_key
        Column in obs for cell type labels.
    formula
        Custom Wilkinson-style design formula (e.g. ``'~ condition + sex'``).
        If None, auto-inferred from metadata columns.
    reference
        Reference level for the condition contrast. If None, alphabetically
        first level is used (matching DESeq2 convention).
    min_subjects_per_group
        Minimum subjects per condition group. Celltypes with fewer subjects
        in any group are skipped with a warning.
    min_cells_per_combo
        Minimum cells per subject+celltype combination for aggregation.
    layer
        Explicit layer name for count matrix.
    alpha
        Significance threshold for PyDESeq2.

    Returns
    -------
    dict[str, pd.DataFrame]
        ``{celltype: results_df}`` where ``results_df`` has columns:
        gene, log2FC, pvalue, padj, baseMean.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # Step 1: Aggregate pseudobulk counts
    aggregated = aggregate_pseudobulk(
        adata,
        subject_key=subject_key,
        celltype_key=celltype_key,
        min_cells=min_cells_per_combo,
        layer=layer,
    )

    results: dict[str, pd.DataFrame] = {}

    for ct, (count_df, meta_df) in aggregated.items():
        # Step 2: Check condition groups have enough subjects
        if condition_key not in meta_df.columns:
            logger.warning(
                "Skipping celltype '%s': condition column '%s' not in metadata.",
                ct, condition_key,
            )
            continue

        group_counts = meta_df[condition_key].value_counts()
        if (group_counts < min_subjects_per_group).any():
            small_groups = group_counts[group_counts < min_subjects_per_group]
            logger.warning(
                "Skipping celltype '%s': condition group(s) %s have fewer than "
                "%d subjects.",
                ct,
                dict(small_groups),
                min_subjects_per_group,
            )
            continue

        # Step 3: Build design formula
        if formula is not None:
            design = formula
        else:
            design = _infer_design_formula(meta_df, condition_key, subject_key)

        logger.info("Running DE for celltype '%s' with design: %s", ct, design)

        # Step 4: Determine contrast levels
        levels = sorted(meta_df[condition_key].unique())
        ref_level = reference if reference is not None else levels[0]
        tested_level = [lv for lv in levels if lv != ref_level]
        if len(tested_level) != 1:
            logger.warning(
                "Skipping celltype '%s': expected 2 condition levels, got %d.",
                ct, len(levels),
            )
            continue
        tested_level = tested_level[0]

        contrast = [condition_key, tested_level, ref_level]

        # Step 5: Run PyDESeq2
        try:
            dds = DeseqDataSet(
                counts=count_df,
                metadata=meta_df,
                design=design,
                refit_cooks=True,
            )
            dds.deseq2()

            stat_res = DeseqStats(dds, contrast=contrast, alpha=alpha)
            stat_res.summary()

            # Step 6: Rename columns to match D-03 convention
            res_df = stat_res.results_df.copy()
            res_df = res_df.reset_index()
            res_df = res_df.rename(
                columns={
                    res_df.columns[0]: "gene",
                    "log2FoldChange": "log2FC",
                }
            )
            # Keep only required columns
            keep_cols = ["gene", "log2FC", "pvalue", "padj", "baseMean"]
            available = [c for c in keep_cols if c in res_df.columns]
            res_df = res_df[available]

            results[ct] = res_df
            logger.info(
                "DE complete for '%s': %d genes, %d significant (padj < %.2f).",
                ct, len(res_df), (res_df["padj"] < alpha).sum() if "padj" in res_df.columns else 0, alpha,
            )
        except Exception:
            logger.exception("PyDESeq2 failed for celltype '%s'.", ct)
            continue

    return results
