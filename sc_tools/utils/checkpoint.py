"""Checkpoint path resolution with backwards compatibility.

Maps semantic phase slugs to their canonical checkpoint filenames.  When
an existing project has data written under the old p-code filenames
(e.g. ``adata.raw.p1.h5ad``), ``resolve_checkpoint_path`` will detect the
legacy file and return its path (with a :class:`DeprecationWarning`) so
pipelines do not break.  New runs always write to the current filenames.

Phase slug  -> canonical filename                  (legacy filenames)
----------     --------------------------           ----------------------------------------
qc_filter   -> results/adata.filtered.h5ad         results/adata.raw.h5ad, results/adata.raw.p1.h5ad
metadata_attach -> results/adata.annotated.h5ad     results/adata.annotated.p2.h5ad
preprocess  -> results/adata.normalized.h5ad        results/adata.normalized.p3.h5ad
scoring     -> results/adata.scored.h5ad            results/adata.normalized.scored.p35.h5ad
celltype_manual -> results/adata.celltyped.h5ad     results/adata.celltyped.p4.h5ad
ingest_load -> data/{sample_id}/adata.ingested.h5ad data/{sample_id}/adata.h5ad, data/{sample_id}/adata.p0.h5ad

Usage::

    from sc_tools.utils.checkpoint import resolve_checkpoint_path, read_checkpoint

    # Returns a Path; prefers new name, falls back to legacy if it exists on disk
    path = resolve_checkpoint_path("scoring")
    adata = read_checkpoint("preprocess")

    # Project-level helper (e.g. from a non-cwd working directory)
    path = resolve_checkpoint_path("scoring", project_dir="/athena/.../robin")

    # For ingest_load (per-sample), pass sample_id as kwarg
    path = resolve_checkpoint_path("ingest_load", sample_id="PT01_NAT")
"""

from __future__ import annotations

import warnings
from pathlib import Path

import anndata as ad

# ---------------------------------------------------------------------------
# Canonical (new) and legacy (old) checkpoint filename templates
# ---------------------------------------------------------------------------

_CHECKPOINT_MAP: dict[str, tuple[str, ...]] = {
    "qc_filter": (
        "results/adata.filtered.h5ad",
        "results/adata.raw.h5ad",
        "results/adata.raw.p1.h5ad",
    ),
    "metadata_attach": (
        "results/adata.annotated.h5ad",
        "results/adata.annotated.p2.h5ad",
    ),
    "preprocess": (
        "results/adata.normalized.h5ad",
        "results/adata.normalized.p3.h5ad",
    ),
    "scoring": (
        "results/adata.scored.h5ad",
        "results/adata.normalized.scored.p35.h5ad",
    ),
    "celltype_manual": (
        "results/adata.celltyped.h5ad",
        "results/adata.celltyped.p4.h5ad",
    ),
    "ingest_load": (
        "data/{sample_id}/adata.ingested.h5ad",
        "data/{sample_id}/adata.h5ad",
        "data/{sample_id}/adata.p0.h5ad",
    ),
}


def resolve_checkpoint_path(
    slug: str,
    project_dir: str | Path = ".",
    **kwargs: str,
) -> Path:
    """Return the resolved checkpoint path for a phase slug.

    Prefers the canonical (new) filename.  If only the legacy filename exists
    on disk, returns it and emits a :class:`DeprecationWarning`.  If neither
    exists, returns the canonical path (expected output location for new runs).

    Parameters
    ----------
    slug:
        Phase slug — one of ``qc_filter``, ``metadata_attach``, ``preprocess``,
        ``scoring``, ``celltype_manual``, ``ingest_load``.
    project_dir:
        Root directory of the project (default: current working directory).
    **kwargs:
        Format arguments for the path template, e.g. ``sample_id="PT01_NAT"``
        when ``slug="ingest_load"``.

    Returns
    -------
    Path
        Resolved absolute (or relative) path to the checkpoint file.

    Raises
    ------
    KeyError
        If *slug* is not recognised.
    """
    if slug not in _CHECKPOINT_MAP:
        raise KeyError(
            f"Unknown checkpoint slug '{slug}'. Accepted slugs: {sorted(_CHECKPOINT_MAP)}"
        )

    project_dir = Path(project_dir)
    templates = _CHECKPOINT_MAP[slug]
    canonical_template = templates[0]
    legacy_templates = templates[1:]

    canonical_path = project_dir / canonical_template.format(**kwargs)

    if canonical_path.exists():
        return canonical_path

    for legacy_template in legacy_templates:
        legacy_path = project_dir / legacy_template.format(**kwargs)
        if legacy_path.exists():
            warnings.warn(
                f"Checkpoint found at legacy path '{legacy_path}'. "
                f"Consider renaming to '{canonical_path}' to use the current convention.",
                DeprecationWarning,
                stacklevel=2,
            )
            return legacy_path

    # None exists — return canonical path so callers get a clear FileNotFoundError
    return canonical_path


def read_checkpoint(
    slug: str,
    project_dir: str | Path = ".",
    backed: str | None = None,
    **kwargs: str,
) -> ad.AnnData:
    """Load a checkpoint AnnData, resolving new or legacy path automatically.

    Parameters
    ----------
    slug:
        Phase slug (see :func:`resolve_checkpoint_path`).
    project_dir:
        Root directory of the project.
    backed:
        If not None, open file in backed mode (``'r'``).
    **kwargs:
        Forwarded to :func:`resolve_checkpoint_path` (e.g. ``sample_id=...``).

    Returns
    -------
    AnnData
    """
    path = resolve_checkpoint_path(slug, project_dir, **kwargs)
    return ad.read_h5ad(path, backed=backed)


# ---------------------------------------------------------------------------
# Snakefile helper (importable from within Snakefile Python blocks)
# ---------------------------------------------------------------------------


def snakemake_compat(new: str, old: str) -> str:
    """Resolve a checkpoint path for Snakemake rules.

    Returns *old* if it exists on disk and *new* does not (backwards
    compatibility for projects that still have legacy-named files).
    Otherwise returns *new*.

    Intended to be called at Snakefile parse time to define checkpoint path
    variables used across rules::

        from sc_tools.utils.checkpoint import snakemake_compat
        ADATA_SCORED = snakemake_compat(
            "results/adata.scored.h5ad",
            "results/adata.normalized.scored.p35.h5ad",
        )

    Parameters
    ----------
    new:
        Canonical (new nomenclature) checkpoint path.
    old:
        Legacy (old nomenclature) checkpoint path.

    Returns
    -------
    str
        Path that Snakemake rules should use.
    """
    import os

    return old if os.path.exists(old) and not os.path.exists(new) else new


__all__ = [
    "resolve_checkpoint_path",
    "read_checkpoint",
    "snakemake_compat",
]
