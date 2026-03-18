"""Phase DAG for the sc_tools pipeline.

Defines a semantic DAG of pipeline phases with explicit dependencies.
DAG nodes are ``(phase_group, subphase)`` tuples.

Quick start
-----------
::

    from sc_tools.pipeline import get_dag, get_available_next, get_phase_checkpoint

    # What can I run from scratch?
    get_available_next([])
    # -> [('data_processing', 'ingest_raw')]

    # After ingestion completes:
    get_available_next([('data_processing', 'ingest_raw'),
                        ('data_processing', 'ingest_load')])
    # -> [('data_processing', 'qc_filter')]

    # Where should the qc_filter output go?
    get_phase_checkpoint('qc_filter')
    # -> 'results/adata.filtered.h5ad'

    # Register a custom phase for a project:
    from sc_tools.pipeline import extend_dag, PhaseSpec
    extend_dag('spatial_regulon', PhaseSpec(
        label='Spatial Regulon Analysis',
        depends_on=[('data_processing', 'scoring')],
        branch='regulon',
        checkpoint='results/adata.regulon.h5ad',
    ))

    # Backward-compat helpers:
    from sc_tools.pipeline import flat_slug_to_tuple, tuple_to_display
    flat_slug_to_tuple('qc_filter')
    # -> ('data_processing', 'qc_filter')
    tuple_to_display(('data_processing', 'qc_filter'))
    # -> 'data_processing/qc_filter'

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

# Type alias for DAG keys
PhaseKey = tuple[str, str]


@dataclass
class PhaseSpec:
    """Specification for a single pipeline phase.

    Parameters
    ----------
    label:
        Human-readable name shown in reports and UIs.
    depends_on:
        List of ``(phase_group, subphase)`` tuples that must be complete
        before this phase can start.  Empty list means this is a root phase.
        For backward compatibility, flat slug strings are also accepted and
        will be auto-converted to ``("data_processing", slug)`` tuples at
        DAG construction time.
    branch:
        Conceptual branch name for grouping parallel tracks.
        Examples: "ingestion", "main", "scoring", "demographics", "meta".
    checkpoint:
        Default output filename template.  May contain ``{sample_id}`` or
        other format fields.  ``None`` means the phase produces no single
        checkpoint file (e.g. figures-only or per-sample outputs).
    phase_group:
        The phase group this phase belongs to.  Standard pipeline phases
        use ``"data_processing"``.  Discovery phases use ``"discovery"``.
    required_obs:
        obs columns that must exist after this phase completes.
    required_obsm:
        obsm keys that must exist after this phase completes.
    x_format:
        Description of what X should contain (e.g. "raw counts", "normalized").
    qc_report:
        Filename template for the QC report produced by this phase (if any).
    old_code:
        Legacy phase code (e.g. "p0a", "p1", "p3.5b") for documentation.
    human_in_loop:
        If True, this phase requires human intervention to complete.
    optional:
        If True, this phase can be skipped without breaking downstream phases.
    iterative:
        If True, this phase can be re-entered (re-run) after completion without
        it being considered an error.  Used for human-in-loop cycles such as
        manual cell typing.
    """

    label: str
    depends_on: list[str | PhaseKey] = field(default_factory=list)
    branch: str = "main"
    checkpoint: str | None = None
    phase_group: str = "data_processing"
    required_obs: list[str] = field(default_factory=list)
    required_obsm: list[str] = field(default_factory=list)
    x_format: str = ""
    qc_report: str | None = None
    old_code: str = ""
    human_in_loop: bool = False
    optional: bool = False
    iterative: bool = False


# ---------------------------------------------------------------------------
# Standard pipeline DAG (flat-slug definitions, converted to tuples below)
# ---------------------------------------------------------------------------

_DP = "data_processing"


def _dp(slug: str) -> PhaseKey:
    """Shorthand for creating a data_processing phase key."""
    return (_DP, slug)


STANDARD_PHASES: dict[str, PhaseSpec] = {
    # -- Ingestion ----------------------------------------------------------
    "ingest_raw": PhaseSpec(
        label="Raw Data Processing",
        depends_on=[],
        branch="ingestion",
        checkpoint=None,
        phase_group=_DP,
        old_code="p0a",
    ),
    "ingest_load": PhaseSpec(
        label="Load into AnnData",
        depends_on=[_dp("ingest_raw")],
        branch="ingestion",
        checkpoint="data/{sample_id}/adata.ingested.h5ad",
        phase_group=_DP,
        required_obs=["sample", "library_id", "raw_data_dir"],
        required_obsm=["spatial"],
        x_format="raw counts",
        old_code="p0b",
    ),
    # -- QC & Metadata ------------------------------------------------------
    "qc_filter": PhaseSpec(
        label="QC Filtering + Concatenation",
        depends_on=[_dp("ingest_load")],
        branch="main",
        checkpoint="results/adata.filtered.h5ad",
        phase_group=_DP,
        required_obs=["sample", "raw_data_dir"],
        required_obsm=["spatial"],
        x_format="raw counts, concatenated",
        qc_report="pre_filter_qc_{date}.html",
        old_code="p1",
    ),
    "metadata_attach": PhaseSpec(
        label="Metadata Attachment",
        depends_on=[_dp("qc_filter")],
        branch="main",
        checkpoint="results/adata.annotated.h5ad",
        phase_group=_DP,
        required_obs=["sample", "raw_data_dir"],
        required_obsm=["spatial"],
        x_format="raw counts, concatenated",
        qc_report="post_filter_qc_{date}.html",
        old_code="p2",
        human_in_loop=True,
    ),
    # -- Preprocessing ------------------------------------------------------
    "preprocess": PhaseSpec(
        label="Normalize + Integrate + Cluster",
        depends_on=[_dp("metadata_attach")],
        branch="main",
        checkpoint="results/adata.normalized.h5ad",
        phase_group=_DP,
        required_obs=["leiden"],
        required_obsm=["X_scvi"],
        x_format="normalized (adata.raw backed up)",
        qc_report="post_integration_qc_{date}.html",
        old_code="p3",
    ),
    # -- Parallel branches from preprocessing -------------------------------
    "demographics": PhaseSpec(
        label="Cohort Demographics",
        depends_on=[_dp("preprocess")],
        branch="demographics",
        checkpoint=None,
        phase_group=_DP,
        old_code="p3.5",
        optional=True,
    ),
    "scoring": PhaseSpec(
        label="Gene Scoring + Auto Cell Typing",
        depends_on=[_dp("preprocess")],
        branch="scoring",
        checkpoint="results/adata.scored.h5ad",
        phase_group=_DP,
        required_obsm=["signature_score", "signature_score_z"],
        x_format="normalized",
        old_code="p3.5b",
    ),
    # -- Cell typing (from scoring) -----------------------------------------
    "celltype_manual": PhaseSpec(
        label="Manual Cell Typing",
        depends_on=[_dp("scoring")],
        branch="celltyping",
        checkpoint="results/adata.celltyped.h5ad",
        phase_group=_DP,
        required_obs=["celltype", "celltype_broad"],
        x_format="normalized",
        qc_report="post_celltyping_qc_{date}.html",
        old_code="p4",
        human_in_loop=True,
        optional=True,
        iterative=True,
    ),
    # -- Downstream biology -------------------------------------------------
    "biology": PhaseSpec(
        label="Downstream Biology",
        depends_on=[_dp("scoring")],
        branch="downstream",
        checkpoint=None,
        phase_group=_DP,
        old_code="p5",
    ),
    "meta_analysis": PhaseSpec(
        label="Meta Analysis",
        depends_on=[_dp("biology")],
        branch="meta",
        checkpoint=None,
        phase_group=_DP,
        old_code="p6/p7",
        optional=True,
    ),
}

# ---------------------------------------------------------------------------
# Runtime registry (mutable; keyed by (phase_group, subphase) tuples)
# ---------------------------------------------------------------------------

_REGISTRY: dict[PhaseKey, PhaseSpec] = {_dp(slug): spec for slug, spec in STANDARD_PHASES.items()}

# Reverse lookup: flat slug -> tuple key (for backward compat)
_SLUG_TO_TUPLE: dict[str, PhaseKey] = {slug: _dp(slug) for slug in STANDARD_PHASES}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _normalize_key(key: str | PhaseKey) -> PhaseKey:
    """Convert a flat slug or tuple to a canonical (phase_group, subphase) tuple."""
    if isinstance(key, tuple):
        return key
    # Flat slug -- look up in the slug-to-tuple map
    if key in _SLUG_TO_TUPLE:
        return _SLUG_TO_TUPLE[key]
    raise KeyError(f"Phase '{key}' not registered. Use get_dag() to see available phases.")


def _normalize_completed(completed: list[str | PhaseKey]) -> set[PhaseKey]:
    """Normalize a list of completed phases (flat slugs or tuples) to a set of tuples."""
    result: set[PhaseKey] = set()
    for item in completed:
        result.add(_normalize_key(item))
    return result


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
        If any ``depends_on`` entry is not registered.
    """
    # Normalize depends_on entries to tuples
    normalized_deps: list[PhaseKey] = []
    for dep in spec.depends_on:
        try:
            normalized_deps.append(_normalize_key(dep))
        except KeyError:
            pass  # Will be caught below

    missing = [d for d in normalized_deps if d not in _REGISTRY]
    # Also check for raw string deps that could not be normalized
    raw_missing = []
    for dep in spec.depends_on:
        if isinstance(dep, str) and dep not in _SLUG_TO_TUPLE:
            raw_missing.append(dep)

    all_missing = raw_missing + [f"{g}/{s}" for g, s in missing]
    if all_missing:
        raise ValueError(
            f"Unknown depends_on slugs for '{slug}': {all_missing}. Register parent phases first."
        )

    # Store with normalized tuple deps
    spec = PhaseSpec(
        label=spec.label,
        depends_on=normalized_deps,
        branch=spec.branch,
        checkpoint=spec.checkpoint,
        phase_group=spec.phase_group,
        required_obs=spec.required_obs,
        required_obsm=spec.required_obsm,
        x_format=spec.x_format,
        qc_report=spec.qc_report,
        old_code=spec.old_code,
        human_in_loop=spec.human_in_loop,
        optional=spec.optional,
        iterative=spec.iterative,
    )

    key = (spec.phase_group, slug)
    _REGISTRY[key] = spec
    _SLUG_TO_TUPLE[slug] = key


def get_dag() -> dict[PhaseKey, PhaseSpec]:
    """Return a copy of the current phase DAG (standard + any custom phases).

    Keys are ``(phase_group, subphase)`` tuples.
    """
    return _REGISTRY.copy()


def get_phase(slug: str | PhaseKey) -> PhaseSpec:
    """Return the :class:`PhaseSpec` for *slug*.

    Accepts either a flat slug (``"qc_filter"``) or a tuple
    (``("data_processing", "qc_filter")``).

    Raises
    ------
    KeyError
        If the slug is not registered.
    """
    key = _normalize_key(slug)
    return _REGISTRY[key]


def get_available_next(
    completed: list[str | PhaseKey],
) -> list[PhaseKey]:
    """Return phase keys whose dependencies are all satisfied.

    Excludes phases that are already complete, except for ``iterative``
    phases which can always be re-entered.

    Parameters
    ----------
    completed:
        List of completed phases.  Accepts flat slugs (``"qc_filter"``)
        or tuples (``("data_processing", "qc_filter")``) or a mix.

    Returns
    -------
    list[tuple[str, str]]
        ``(phase_group, subphase)`` keys available to run next.
    """
    completed_set = _normalize_completed(completed)
    return [
        key
        for key, spec in _REGISTRY.items()
        if (key not in completed_set or spec.iterative)
        and all(_normalize_key(d) in completed_set for d in spec.depends_on)
    ]


def get_phase_checkpoint(slug: str | PhaseKey, **kwargs: str) -> str | None:
    """Return the expected checkpoint path for a phase, with placeholders filled in.

    Accepts either a flat slug or a ``(phase_group, subphase)`` tuple.

    Parameters
    ----------
    slug:
        Phase identifier (e.g. ``"qc_filter"`` or
        ``("data_processing", "qc_filter")``).
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
    for key, spec in _REGISTRY.items():
        for dep in spec.depends_on:
            try:
                dep_key = _normalize_key(dep)
            except KeyError:
                errors.append(f"Phase {key!r} depends on unknown phase {dep!r}")
                continue
            if dep_key not in _REGISTRY:
                errors.append(f"Phase {key!r} depends on unknown phase {dep!r}")
    return errors


# ---------------------------------------------------------------------------
# Backward-compat helpers
# ---------------------------------------------------------------------------


def flat_slug_to_tuple(slug: str) -> PhaseKey:
    """Map a flat phase slug to its ``(phase_group, subphase)`` tuple.

    Parameters
    ----------
    slug:
        Flat phase slug, e.g. ``"qc_filter"``.

    Returns
    -------
    tuple[str, str]
        ``(phase_group, subphase)`` tuple.

    Raises
    ------
    KeyError
        If the slug is not registered.
    """
    if slug in _SLUG_TO_TUPLE:
        return _SLUG_TO_TUPLE[slug]
    raise KeyError(f"Phase '{slug}' not registered. Use get_dag() to see available phases.")


def tuple_to_display(key: PhaseKey) -> str:
    """Format a ``(phase_group, subphase)`` tuple for display.

    Parameters
    ----------
    key:
        ``(phase_group, subphase)`` tuple.

    Returns
    -------
    str
        Display string like ``"data_processing/qc_filter"``.
    """
    return f"{key[0]}/{key[1]}"


__all__ = [
    "PhaseKey",
    "PhaseSpec",
    "STANDARD_PHASES",
    "extend_dag",
    "flat_slug_to_tuple",
    "get_available_next",
    "get_dag",
    "get_phase",
    "get_phase_checkpoint",
    "tuple_to_display",
    "validate_dag",
]
