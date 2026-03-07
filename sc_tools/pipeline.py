"""Phase DAG for the sc_tools pipeline.

Defines a semantic, slug-based DAG of pipeline phases with explicit dependencies.
Replaces the p0a/p0b/p1/p2/p3/p35/p4 naming with human-readable slugs.

Quick start
-----------
::

    from sc_tools.pipeline import get_dag, get_available_next, get_phase_checkpoint

    # What can I run from scratch?
    get_available_next([])
    # -> ['ingest_raw']

    # After ingestion completes:
    get_available_next(['ingest_raw', 'ingest_load'])
    # -> ['qc_filter']

    # Where should the qc_filter output go?
    get_phase_checkpoint('qc_filter')
    # -> 'results/adata.raw.h5ad'

    # Register a custom phase for a project:
    from sc_tools.pipeline import extend_dag, PhaseSpec
    extend_dag('spatial_regulon', PhaseSpec(
        label='Spatial Regulon Analysis',
        depends_on=['scoring'],
        branch='regulon',
        checkpoint='results/adata.regulon.h5ad',
    ))

Phase name mapping (new slug -> old code)
-----------------------------------------
    ingest_raw      p0a
    ingest_load     p0b
    qc_filter       p1
    metadata_attach p2
    preprocess      p3
    demographics    p3.5
    scoring         p3.5b
    celltype_manual p4
    biology         p5
    meta_analysis   p6/p7
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class PhaseSpec:
    """Specification for a single pipeline phase.

    Parameters
    ----------
    label:
        Human-readable name shown in reports and UIs.
    depends_on:
        List of phase slugs that must be complete before this phase can start.
        Empty list means this is a root phase (no prerequisites).
    branch:
        Conceptual branch name for grouping parallel tracks.
        Examples: "ingestion", "main", "scoring", "demographics", "meta".
    checkpoint:
        Default output filename template.  May contain ``{sample_id}`` or
        other format fields.  ``None`` means the phase produces no single
        checkpoint file (e.g. figures-only or per-sample outputs).
    optional:
        If True, this phase can be skipped without breaking downstream phases.
    iterative:
        If True, this phase can be re-entered (re-run) after completion without
        it being considered an error.  Used for human-in-loop cycles such as
        manual cell typing.
    """

    label: str
    depends_on: list[str] = field(default_factory=list)
    branch: str = "main"
    checkpoint: str | None = None
    optional: bool = False
    iterative: bool = False


# ---------------------------------------------------------------------------
# Standard pipeline DAG
# ---------------------------------------------------------------------------

STANDARD_PHASES: dict[str, PhaseSpec] = {
    # ── Ingestion ────────────────────────────────────────────────────────────
    "ingest_raw": PhaseSpec(
        label="Raw Data Processing",
        depends_on=[],
        branch="ingestion",
        checkpoint=None,  # platform output dirs (not a single h5ad)
    ),
    "ingest_load": PhaseSpec(
        label="Load into AnnData",
        depends_on=["ingest_raw"],
        branch="ingestion",
        checkpoint="data/{sample_id}/adata.h5ad",  # per-sample
    ),
    # ── QC & Metadata ────────────────────────────────────────────────────────
    "qc_filter": PhaseSpec(
        label="QC Filtering + Concatenation",
        depends_on=["ingest_load"],
        branch="main",
        checkpoint="results/adata.raw.h5ad",
    ),
    "metadata_attach": PhaseSpec(
        label="Metadata Attachment",
        depends_on=["qc_filter"],
        branch="main",
        checkpoint="results/adata.annotated.h5ad",
    ),
    # ── Preprocessing ────────────────────────────────────────────────────────
    "preprocess": PhaseSpec(
        label="Normalize + Integrate + Cluster",
        depends_on=["metadata_attach"],
        branch="main",
        checkpoint="results/adata.normalized.h5ad",
    ),
    # ── Parallel branches from preprocessing ─────────────────────────────────
    "demographics": PhaseSpec(
        label="Cohort Demographics",
        depends_on=["preprocess"],
        branch="demographics",
        checkpoint=None,  # figures only, no checkpoint file
        optional=True,
    ),
    "scoring": PhaseSpec(
        label="Gene Scoring + Auto Cell Typing",
        depends_on=["preprocess"],
        branch="scoring",
        checkpoint="results/adata.scored.h5ad",
    ),
    # ── Cell typing (from scoring) ────────────────────────────────────────────
    "celltype_manual": PhaseSpec(
        label="Manual Cell Typing",
        depends_on=["scoring"],
        branch="celltyping",
        checkpoint="results/adata.celltyped.h5ad",
        optional=True,  # skippable if auto typing is adequate
        iterative=True,  # human-in-loop: annotate → review → repeat
    ),
    # ── Downstream biology ───────────────────────────────────────────────────
    "biology": PhaseSpec(
        label="Downstream Biology",
        depends_on=["scoring"],  # can start from scoring or after celltype_manual
        branch="downstream",
        checkpoint=None,
    ),
    "meta_analysis": PhaseSpec(
        label="Meta Analysis",
        depends_on=["biology"],
        branch="meta",
        checkpoint=None,
        optional=True,
    ),
}

# ---------------------------------------------------------------------------
# Runtime registry (mutable; starts as a copy of STANDARD_PHASES)
# ---------------------------------------------------------------------------

_REGISTRY: dict[str, PhaseSpec] = STANDARD_PHASES.copy()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def extend_dag(slug: str, spec: PhaseSpec) -> None:
    """Register a custom phase (project-specific or experimental).

    Parameters
    ----------
    slug:
        Unique identifier for the phase (e.g. ``"spatial_regulon"``).
    spec:
        :class:`PhaseSpec` describing the phase.

    Raises
    ------
    ValueError
        If any ``depends_on`` slug is not registered.
    """
    missing = [d for d in spec.depends_on if d not in _REGISTRY]
    if missing:
        raise ValueError(
            f"Unknown depends_on slugs for '{slug}': {missing}. Register parent phases first."
        )
    _REGISTRY[slug] = spec


def get_dag() -> dict[str, PhaseSpec]:
    """Return a copy of the current phase DAG (standard + any custom phases)."""
    return _REGISTRY.copy()


def get_phase(slug: str) -> PhaseSpec:
    """Return the :class:`PhaseSpec` for *slug*.

    Raises
    ------
    KeyError
        If the slug is not registered.
    """
    if slug not in _REGISTRY:
        raise KeyError(f"Phase '{slug}' not registered. Use get_dag() to see available phases.")
    return _REGISTRY[slug]


def get_available_next(completed: list[str]) -> list[str]:
    """Return phase slugs whose dependencies are all satisfied.

    Excludes phases that are already complete, except for ``iterative``
    phases which can always be re-entered.

    Parameters
    ----------
    completed:
        List of phase slugs that have already been completed for this project.

    Returns
    -------
    list[str]
        Phase slugs available to run next (order is insertion order of DAG).
    """
    return [
        p
        for p, s in _REGISTRY.items()
        if (p not in completed or s.iterative) and all(d in completed for d in s.depends_on)
    ]


def get_phase_checkpoint(slug: str, **kwargs: str) -> str | None:
    """Return the expected checkpoint path for a phase, with placeholders filled in.

    Parameters
    ----------
    slug:
        Phase identifier (e.g. ``"qc_filter"``).
    **kwargs:
        Format arguments for the checkpoint template (e.g. ``sample_id="s1"``).

    Returns
    -------
    str or None
        Formatted checkpoint path, or ``None`` if the phase has no checkpoint.
    """
    spec = get_phase(slug)
    if spec.checkpoint is None:
        return None
    return spec.checkpoint.format(**kwargs)


def validate_dag() -> list[str]:
    """Check that all ``depends_on`` references point to registered phases.

    Returns
    -------
    list[str]
        Error messages (empty if the DAG is valid).
    """
    errors: list[str] = []
    for slug, spec in _REGISTRY.items():
        for dep in spec.depends_on:
            if dep not in _REGISTRY:
                errors.append(f"Phase '{slug}' depends on unknown phase '{dep}'")
    return errors


__all__ = [
    "PhaseSpec",
    "STANDARD_PHASES",
    "extend_dag",
    "get_dag",
    "get_phase",
    "get_available_next",
    "get_phase_checkpoint",
    "validate_dag",
]
