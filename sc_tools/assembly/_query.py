"""Cross-modal query functions for multi-omic data.

Provides:
- celltype_proportions: Compute cell type proportions per group across modalities.
- aggregate_by_level: N-level abstraction stacking with arbitrary grouping columns.
"""

from __future__ import annotations

import logging

import pandas as pd

__all__ = ["celltype_proportions", "aggregate_by_level"]

logger = logging.getLogger(__name__)


def celltype_proportions(
    mdata,
    *,
    celltype_key: str = "celltype",
    group_by: str = "subject_id",
) -> pd.DataFrame:
    """Compute cell type proportions per group across modalities.

    Parameters
    ----------
    mdata
        MuData object.
    celltype_key
        Column name for cell type annotations in obs.
    group_by
        Column name to group by (e.g. subject_id, sample_id).

    Returns
    -------
    pd.DataFrame
        Columns: [group_by, celltype_key, 'modality', 'count', 'proportion'].
        Modalities missing the celltype_key column are skipped.
    """
    results = []

    for mod_name, mod_adata in mdata.mod.items():
        if celltype_key not in mod_adata.obs.columns:
            logger.info(
                "Skipping modality '%s': '%s' not in obs.columns",
                mod_name,
                celltype_key,
            )
            continue

        if group_by not in mod_adata.obs.columns:
            logger.info(
                "Skipping modality '%s': '%s' not in obs.columns",
                mod_name,
                group_by,
            )
            continue

        # Count cells per (group_by, celltype_key)
        counts = (
            mod_adata.obs
            .groupby([group_by, celltype_key], observed=True)
            .size()
            .reset_index(name="count")
        )

        # Compute proportions within each group
        group_totals = counts.groupby(group_by, observed=True)["count"].transform("sum")
        counts["proportion"] = counts["count"] / group_totals
        counts["modality"] = mod_name

        results.append(counts)

    if not results:
        return pd.DataFrame(
            columns=[group_by, celltype_key, "modality", "count", "proportion"]
        )

    return pd.concat(results, ignore_index=True)


def aggregate_by_level(
    mdata,
    *,
    group_cols: list[str],
    value_col: str = "celltype",
) -> pd.DataFrame:
    """N-level abstraction stacking with arbitrary grouping columns.

    Parameters
    ----------
    mdata
        MuData object.
    group_cols
        List of column names to group by (e.g. ["subject_id"], ["subject_id", "organ"]).
    value_col
        Column name for the value to aggregate (default: "celltype").

    Returns
    -------
    pd.DataFrame
        Columns: [*group_cols, value_col, 'modality', 'count', 'proportion'].
        Modalities missing required columns are skipped.
    """
    results = []

    for mod_name, mod_adata in mdata.mod.items():
        obs = mod_adata.obs

        # Check which group_cols are available in this modality
        available_cols = [c for c in group_cols if c in obs.columns]
        if not available_cols:
            logger.info(
                "Skipping modality '%s': none of %s found in obs.columns",
                mod_name,
                group_cols,
            )
            continue

        if value_col not in obs.columns:
            logger.info(
                "Skipping modality '%s': '%s' not in obs.columns",
                mod_name,
                value_col,
            )
            continue

        # Group by available columns + value_col
        grp_keys = available_cols + [value_col]
        counts = (
            obs
            .groupby(grp_keys, observed=True)
            .size()
            .reset_index(name="count")
        )

        # Compute proportions within each group (excluding value_col)
        if available_cols:
            group_totals = counts.groupby(available_cols, observed=True)["count"].transform("sum")
        else:
            group_totals = counts["count"].sum()
        counts["proportion"] = counts["count"] / group_totals
        counts["modality"] = mod_name

        results.append(counts)

    if not results:
        return pd.DataFrame(
            columns=[*group_cols, value_col, "modality", "count", "proportion"]
        )

    return pd.concat(results, ignore_index=True)
