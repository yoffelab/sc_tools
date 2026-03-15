"""MCP server: sc-registry bookkeeping.

Exposes the sc_tools registry as callable MCP tools so Claude Code can
query project state, checkpoint locations, and data objects.

Tools
-----
    registry_status      -- high-level summary (projects, data, patients)
    list_datasets        -- all checkpoints for a project (queries data table)
    get_checkpoint_uri   -- "where is adata.normalized.h5ad for ggo_visium?"
    register_dataset     -- add a new checkpoint (routes to data table)
    mark_phase_complete  -- mark a pipeline phase as complete for a project
    get_phase_status     -- status and checkpoint details for a specific phase
    set_phase_status     -- update the pipeline status for a phase
    add_subject          -- register a de-identified patient (backward compat)
    list_subjects        -- query patients (backward compat)
    register_biodata     -- register a data object (backward compat)
    list_biodata         -- query data objects (backward compat)
    project_data_summary -- counts by category/platform

Start the server::

    python -m sc_tools.mcp.registry_server

Requires: pip install "sc-tools[mcp,registry]"
"""

from __future__ import annotations

import json
import logging

logger = logging.getLogger(__name__)

try:
    from mcp.server.fastmcp import FastMCP
except ImportError as _mcp_err:
    raise ImportError(
        "mcp is required for the MCP registry server: pip install sc-tools[mcp]"
    ) from _mcp_err

mcp = FastMCP("sc-registry")


def _registry():
    """Return a Registry instance using the configured DB."""
    from sc_tools.registry import Registry

    return Registry()


# ---------------------------------------------------------------------------
# Tool: registry_status
# ---------------------------------------------------------------------------


@mcp.tool()
def registry_status() -> str:
    """Return a high-level status summary of the sc_tools registry.

    Shows total projects, data objects, patients, and per-project status.

    Returns
    -------
    str
        Formatted summary text.
    """
    try:
        reg = _registry()
        s = reg.status()
        lines = [
            "sc_tools Registry Status",
            "=" * 40,
            f"  Projects          : {s['n_projects']}",
            f"  Data objects      : {s['n_data']}",
            f"  Patients          : {s['n_patients']}",
        ]
        if s["active_projects"]:
            lines.append("\n  Active projects:")
            for name in s["active_projects"]:
                lines.append(f"    - {name}")
                phase_summary = s.get("phase_summary", {}).get(name, {})
                if phase_summary:
                    parts = ", ".join(f"{k}={v}" for k, v in sorted(phase_summary.items()))
                    lines.append(f"      data status: {parts}")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_datasets
# ---------------------------------------------------------------------------


@mcp.tool()
def list_datasets(project_name: str, phase: str = "") -> str:
    """List all checkpoint datasets for a project.

    Parameters
    ----------
    project_name
        Project name (e.g. ggo_visium, robin).
    phase
        Optional phase filter (e.g. qc_filter, scoring). Empty means all.

    Returns
    -------
    str
        Table of datasets with id, phase, role, status, and URI.
    """
    try:
        reg = _registry()
        datasets = reg.list_datasets(
            project_name=project_name,
            phase=phase if phase else None,
        )
        if not datasets:
            return (
                f"No datasets found for project '{project_name}'"
                + (f" phase '{phase}'" if phase else "")
                + "."
            )
        lines = [f"Datasets for '{project_name}'" + (f" phase={phase}" if phase else "") + ":"]
        for ds in datasets:
            role = ds.get("file_role", "primary")
            n_obs = ds.get("n_obs")
            obs_str = f" n_obs={n_obs}" if n_obs else ""
            lines.append(
                f"  id={ds['id']} phase={ds['phase']} role={role} "
                f"status={ds['status']}{obs_str}\n"
                f"    uri: {ds['uri']}"
            )
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: get_checkpoint_uri
# ---------------------------------------------------------------------------


@mcp.tool()
def get_checkpoint_uri(project_name: str, phase: str, sample_id: str = "") -> str:
    """Return the URI for a specific checkpoint.

    Parameters
    ----------
    project_name
        Project name (e.g. ggo_visium).
    phase
        Phase slug: ingest_raw, ingest_load, qc_filter, metadata_attach,
        preprocess, scoring, celltype_manual, biology, meta_analysis.
    sample_id
        Optional sample identifier (for per-sample ingest_load checkpoints).

    Returns
    -------
    str
        URI string, or a message indicating the checkpoint is not registered.
    """
    try:
        reg = _registry()
        uri = reg.get_dataset_uri(
            project_name=project_name,
            phase=phase,
            sample_id=sample_id if sample_id else None,
        )
        if uri is None:
            return (
                f"No checkpoint registered for project='{project_name}' "
                f"phase='{phase}'" + (f" sample='{sample_id}'" if sample_id else "") + "."
            )
        return uri
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: register_dataset
# ---------------------------------------------------------------------------


@mcp.tool()
def register_dataset(
    project_name: str,
    phase: str,
    uri: str,
    fmt: str = "h5ad",
    sample_id: str = "",
    status: str = "ready",
    file_role: str = "primary",
    validated: bool = False,
    n_obs: int = 0,
    n_vars: int = 0,
) -> str:
    """Add a new checkpoint to the registry.

    Parameters
    ----------
    project_name
        Project name. Must already exist (call add_project first if needed).
    phase
        Phase slug: ingest_raw, ingest_load, qc_filter, metadata_attach,
        preprocess, scoring, celltype_manual, biology, meta_analysis.
    uri
        Path or URI to the checkpoint file.
    fmt
        File format: h5ad (default), zarr, tiff, tsv.
    sample_id
        Optional sample identifier (for per-sample ingest_load checkpoints).
    status
        Dataset status: ready (default), pending, archived, or error.
    file_role
        File role: primary (default), supplementary, entry_point,
        spatialdata, image, or metadata.
    validated
        True if validate_checkpoint() has passed for this file.
    n_obs
        Number of cells/spots (0 = unknown).
    n_vars
        Number of genes/proteins (0 = unknown).

    Returns
    -------
    str
        Confirmation with assigned dataset id.
    """
    try:
        reg = _registry()
        ds_id = reg.register_dataset(
            project_name=project_name,
            phase=phase,
            uri=uri,
            sample_id=sample_id if sample_id else None,
            fmt=fmt,
            status=status,
            file_role=file_role,
            validated=validated,
            n_obs=n_obs if n_obs else None,
            n_vars=n_vars if n_vars else None,
        )
        return f"Registered dataset id={ds_id}: {project_name}/{phase} [{file_role}] at {uri}"
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_slurm_jobs (deprecated)
# ---------------------------------------------------------------------------


@mcp.tool()
def list_slurm_jobs(active_only: bool = True) -> str:
    """SLURM jobs table has been removed (migration 0008).

    Parameters
    ----------
    active_only
        Ignored (table no longer exists).

    Returns
    -------
    str
        Notice that the table has been removed.
    """
    return "SLURM jobs table has been removed (migration 0008). No jobs to display."


# ---------------------------------------------------------------------------
# Tool: list_agent_tasks (deprecated)
# ---------------------------------------------------------------------------


@mcp.tool()
def list_agent_tasks(running_only: bool = True) -> str:
    """Agent tasks table has been removed.

    Parameters
    ----------
    running_only
        Ignored (table no longer exists).

    Returns
    -------
    str
        Notice that the table has been removed.
    """
    return "Agent tasks table has been removed (migration 0007). No tasks to display."


# ---------------------------------------------------------------------------
# Tool: add_project
# ---------------------------------------------------------------------------


@mcp.tool()
def add_project(
    name: str,
    platform: str = "",
    data_type: str = "",
    domain: str = "",
    imaging_modality: str = "",
    project_type: str = "internal",
    visibility: str = "private",
) -> str:
    """Register a new project in the registry.

    If a project with the same name already exists, its id is returned without
    creating a duplicate (idempotent).

    Parameters
    ----------
    name
        Unique project name (e.g. ggo_visium, robin, lymph_dlbcl).
    platform
        Technology platform string (e.g. visium, imc, xenium).
    domain
        High-level domain: spatial_transcriptomics | spatial_proteomics |
        imaging | single_cell | bulk.

    Returns
    -------
    str
        Confirmation with assigned project id.
    """
    try:
        reg = _registry()
        proj_id = reg.add_project(
            name,
            platform=platform or None,
            domain=domain or None,
        )
        return f"Project '{name}' registered (id={proj_id})."
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: get_available_next_phases
# ---------------------------------------------------------------------------


@mcp.tool()
def get_available_next_phases(project_name: str) -> str:
    """Return which pipeline phases are available to run next for a project.

    Queries the registry for completed phases then applies the pipeline DAG
    to determine what is unblocked.

    Parameters
    ----------
    project_name
        Project name (e.g. ggo_visium).

    Returns
    -------
    str
        List of available phase slugs with labels, checkpoints, and
        dependency information.
    """
    try:
        from sc_tools.pipeline import get_available_next, get_dag

        reg = _registry()
        phase_rows = reg.list_phases(project_name)
        completed = [r["phase"] for r in phase_rows if r["status"] == "complete"]
        in_progress = [r["phase"] for r in phase_rows if r["status"] == "in_progress"]

        dag = get_dag()
        available = get_available_next(completed)

        lines = [
            f"Pipeline status for '{project_name}'",
            "=" * 40,
            f"  Completed  : {', '.join(completed) or 'none'}",
            f"  In progress: {', '.join(in_progress) or 'none'}",
            "",
            "  Available next:",
        ]
        if not available:
            lines.append("    (none -- all phases complete or blocked)")
        else:
            for slug in available:
                spec = dag[slug]
                cp = spec.checkpoint or "(no checkpoint file)"
                flags = []
                if spec.optional:
                    flags.append("optional")
                if spec.iterative:
                    flags.append("iterative")
                flag_str = f"  [{', '.join(flags)}]" if flags else ""
                lines.append(f"    {slug:<20} {spec.label}{flag_str}")
                lines.append(f"      checkpoint: {cp}")
                if spec.depends_on:
                    lines.append(f"      depends on: {', '.join(spec.depends_on)}")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: record_provenance
# ---------------------------------------------------------------------------


@mcp.tool()
def record_provenance(project_name: str, phase: str, env_snapshot_json: str) -> str:
    """Store a runtime environment snapshot alongside a pipeline phase record.

    Parameters
    ----------
    project_name
        Project name (e.g. ggo_visium).
    phase
        Phase slug (e.g. preprocess, scoring).
    env_snapshot_json
        JSON string from environment_info() output, or any provenance note.

    Returns
    -------
    str
        Confirmation message.
    """
    try:
        reg = _registry()
        reg.upsert_phase(project_name, phase, notes=env_snapshot_json)
        return f"Provenance recorded for '{project_name}' phase '{phase}'."
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: update_slurm_job_status (deprecated)
# ---------------------------------------------------------------------------


@mcp.tool()
def update_slurm_job_status(
    slurm_job_id: str,
    status: str,
    error: str = "",
    log_uri: str = "",
) -> str:
    """SLURM jobs table has been removed (migration 0008).

    Parameters
    ----------
    slurm_job_id
        Ignored.
    status
        Ignored.

    Returns
    -------
    str
        Notice that the table has been removed.
    """
    return "SLURM jobs table has been removed (migration 0008). Cannot update."


# ---------------------------------------------------------------------------
# Tool: mark_phase_complete
# ---------------------------------------------------------------------------


@mcp.tool()
def mark_phase_complete(project_name: str, phase: str) -> str:
    """Mark a pipeline phase as complete for a project.

    Updates data rows for this project+phase to status='ready'.

    Parameters
    ----------
    project_name
        Project name (e.g. ggo_visium).
    phase
        Phase slug: ingest_raw, ingest_load, qc_filter, metadata_attach,
        preprocess, scoring, celltype_manual, biology, meta_analysis.

    Returns
    -------
    str
        Confirmation message.
    """
    try:
        reg = _registry()
        reg.mark_phase_complete(project_name, phase)
        return f"Marked phase '{phase}' as complete for project '{project_name}'."
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: get_phase_status
# ---------------------------------------------------------------------------


@mcp.tool()
def get_phase_status(project_name: str, phase: str) -> str:
    """Return the status and checkpoint details for a specific phase of a project.

    Parameters
    ----------
    project_name
        Project name (e.g. ggo_visium).
    phase
        Phase slug: ingest_raw, ingest_load, qc_filter, metadata_attach,
        preprocess, scoring, celltype_manual, biology, meta_analysis.

    Returns
    -------
    str
        Formatted phase details including status, n_obs, n_samples, notes.
    """
    try:
        reg = _registry()
        row = reg.get_phase(project_name, phase)
        if row is None:
            return (
                f"No phase entry for project='{project_name}' phase='{phase}'. "
                "Use set_phase_status to create one."
            )
        lines = [f"Phase '{phase}' for project '{project_name}':"]
        lines.append(f"  status        : {row['status']}")
        lines.append(f"  entry_phase   : {row.get('entry_phase', False)}")
        if row.get("n_obs"):
            lines.append(f"  n_obs         : {row['n_obs']}")
        if row.get("n_vars"):
            lines.append(f"  n_vars        : {row['n_vars']}")
        if row.get("n_samples"):
            lines.append(f"  n_samples     : {row['n_samples']}")
        if row.get("primary_dataset_id"):
            lines.append(f"  dataset_id    : {row['primary_dataset_id']}")
        if row.get("notes"):
            lines.append(f"  notes         : {row['notes']}")
        lines.append(f"  updated_at    : {row.get('updated_at', '')}")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: set_phase_status
# ---------------------------------------------------------------------------


@mcp.tool()
def set_phase_status(
    project_name: str,
    phase: str,
    status: str,
    notes: str = "",
    n_obs: int = 0,
    n_vars: int = 0,
    n_samples: int = 0,
) -> str:
    """Update the pipeline status for a phase.

    Creates a data row if none exists for this phase (upsert).

    Parameters
    ----------
    project_name
        Project name (e.g. ggo_visium).
    phase
        Phase slug: ingest_raw, ingest_load, qc_filter, metadata_attach,
        preprocess, scoring, celltype_manual, biology, meta_analysis.
    status
        New status: not_started | in_progress | complete | failed | skipped.
    notes
        Optional free-text notes (e.g. method selection rationale).
    n_obs
        Number of cells/spots at this phase (0 to skip).
    n_vars
        Number of genes/proteins at this phase (0 to skip).
    n_samples
        Number of samples at this phase (0 to skip).

    Returns
    -------
    str
        Confirmation message.
    """
    try:
        reg = _registry()
        reg.upsert_phase(
            project_name,
            phase,
            status=status,
            notes=notes if notes else None,
            n_obs=n_obs if n_obs else None,
            n_vars=n_vars if n_vars else None,
            n_samples=n_samples if n_samples else None,
        )
        parts = [f"Set phase '{phase}' to status='{status}' for project '{project_name}'."]
        if n_obs:
            parts.append(f"n_obs={n_obs:,}")
        if n_vars:
            parts.append(f"n_vars={n_vars:,}")
        if n_samples:
            parts.append(f"n_samples={n_samples}")
        return " ".join(parts)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: register_data_source
# ---------------------------------------------------------------------------


@mcp.tool()
def register_data_source(
    name: str,
    uri: str,
    description: str = "",
    platform: str = "",
    domain: str = "",
    imaging_modality: str = "",
    source_type: str = "",
    organism: str = "",
    tissue: str = "",
    disease: str = "",
    n_samples: int = 0,
    n_cells: int = 0,
    publication: str = "",
    access_notes: str = "",
    status: str = "available",
) -> str:
    """Register a raw data source in the catalog.

    Parameters
    ----------
    name
        Unique identifier (e.g. saha_ibd_cosmx, 10x_visium_human_breast).
    uri
        Primary location: HPC path, URL, GEO accession, DOI, etc.
    source_type
        hpc_lab | hpc_collaborator | public_10x | public_geo |
        public_zenodo | public_portal.
    status
        available | restricted | pending_download | archived.

    Returns
    -------
    str
        Confirmation with assigned data source id.
    """
    try:
        reg = _registry()
        src_id = reg.register_data_source(
            name,
            uri,
            description=description or None,
            platform=platform or None,
            domain=domain or None,
            imaging_modality=imaging_modality or None,
            source_type=source_type or None,
            organism=organism or None,
            tissue=tissue or None,
            disease=disease or None,
            n_samples=n_samples or None,
            n_cells=n_cells or None,
            publication=publication or None,
            access_notes=access_notes or None,
            status=status,
        )
        return f"Data source '{name}' registered (id={src_id})."
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_data_sources
# ---------------------------------------------------------------------------


@mcp.tool()
def list_data_sources(
    platform: str = "",
    source_type: str = "",
    disease: str = "",
    tissue: str = "",
) -> str:
    """List data sources in the catalog, with optional filters.

    Parameters
    ----------
    platform
        Filter by platform (e.g. cosmx, xenium, imc, visium_hd).
    source_type
        Filter by source type (hpc_collaborator, public_10x, public_geo, ...).
    disease
        Filter by disease keyword (partial match, e.g. IBD, cancer).
    tissue
        Filter by tissue keyword (partial match, e.g. colon, brain).

    Returns
    -------
    str
        Formatted table of matching data sources.
    """
    try:
        reg = _registry()
        sources = reg.list_data_sources(
            platform=platform or None,
            source_type=source_type or None,
            disease=disease or None,
            tissue=tissue or None,
        )
        if not sources:
            return "No data sources found matching the given filters."
        lines = [f"Data sources ({len(sources)} total):"]
        lines.append(f"  {'ID':<4}  {'Name':<40}  {'Platform':<14}  {'Source':<20}  {'Status'}")
        lines.append("  " + "-" * 96)
        for s in sources:
            lines.append(
                f"  {s['id']:<4}  {(s['name'] or ''):<40}  "
                f"{(s['platform'] or ''):<14}  {(s['source_type'] or ''):<20}  "
                f"{s['status']}"
            )
            if s.get("description"):
                lines.append(f"       {s['description'][:80]}")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: link_project_data_source
# ---------------------------------------------------------------------------


@mcp.tool()
def link_project_data_source(
    project_name: str,
    data_source_name: str,
    role: str = "input",
    notes: str = "",
) -> str:
    """Link a project to a data source.

    Parameters
    ----------
    project_name
        Project name (must already exist in registry).
    data_source_name
        Data source name (must already exist in registry).
    role
        input (primary data used) | reference (e.g. scRNA ref for deconvolution) |
        supplementary (additional context).
    notes
        Optional free-text notes about the link.

    Returns
    -------
    str
        Confirmation message.
    """
    try:
        reg = _registry()
        link_id = reg.link_project_data_source(
            project_name, data_source_name, role=role, notes=notes or None
        )
        return (
            f"Linked project '{project_name}' -> data source '{data_source_name}' "
            f"(role={role}, link_id={link_id})."
        )
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: add_subject (backward compat -> add_patient)
# ---------------------------------------------------------------------------


@mcp.tool()
def add_subject(
    subject_id: str,
    organism: str = "human",
    sex: str = "",
    age_at_collection: float = 0,
    diagnosis: str = "",
    diagnosis_code: str = "",
    disease_stage: str = "",
    treatment_status: str = "",
    tissue_of_origin: str = "",
    vital_status: str = "",
    survival_days: float = 0,
) -> str:
    """Register a de-identified patient in the registry.

    Patients are cross-project. No PHI (names, DOB, MRNs) should be stored.

    Parameters
    ----------
    subject_id
        Unique de-identified identifier (e.g. PT001).
    organism
        human (default), mouse, rat, etc.
    sex
        M, F, or unknown.
    diagnosis
        Free text diagnosis (e.g. DLBCL, UC, GGO).

    Returns
    -------
    str
        Confirmation with assigned patient DB id.
    """
    try:
        reg = _registry()
        kwargs: dict = {"organism": organism}
        if sex:
            kwargs["sex"] = sex
        if age_at_collection:
            kwargs["age_at_collection"] = age_at_collection
        if diagnosis:
            kwargs["diagnosis"] = diagnosis
        if diagnosis_code:
            kwargs["diagnosis_code"] = diagnosis_code
        if disease_stage:
            kwargs["disease_stage"] = disease_stage
        if treatment_status:
            kwargs["treatment_status"] = treatment_status
        if tissue_of_origin:
            kwargs["tissue_of_origin"] = tissue_of_origin
        if vital_status:
            kwargs["vital_status"] = vital_status
        if survival_days:
            kwargs["survival_days"] = survival_days
        sid = reg.add_patient(subject_id, metadata=kwargs)
        return f"Patient '{subject_id}' registered (id={sid})."
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_subjects (backward compat -> list_patients)
# ---------------------------------------------------------------------------


@mcp.tool()
def list_subjects(
    project_name: str = "",
    diagnosis: str = "",
    tissue: str = "",
) -> str:
    """List patients in the registry, with optional filters.

    Parameters
    ----------
    project_name
        Filter by linked project name.
    diagnosis
        Filter by diagnosis keyword (partial match).
    tissue
        Filter by tissue_of_origin keyword (partial match).

    Returns
    -------
    str
        Formatted list of patients.
    """
    try:
        reg = _registry()
        subjects = reg.list_subjects(
            project_name=project_name or None,
            diagnosis=diagnosis or None,
            tissue=tissue or None,
        )
        if not subjects:
            return "No patients found matching the given filters."
        lines = [f"Patients ({len(subjects)} total):"]
        for s in subjects:
            parts = [f"  {s['patient_id']}"]
            if s.get("organism"):
                parts.append(f"organism={s['organism']}")
            if s.get("sex"):
                parts.append(f"sex={s['sex']}")
            if s.get("diagnosis"):
                parts.append(f"dx={s['diagnosis']}")
            if s.get("tissue_of_origin"):
                parts.append(f"tissue={s['tissue_of_origin']}")
            lines.append(" ".join(parts))
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: add_sample (deprecated)
# ---------------------------------------------------------------------------


@mcp.tool()
def add_sample(
    sample_id: str,
    subject_id: str = "",
    project_name: str = "",
    tissue: str = "",
    tissue_region: str = "",
    fixation_method: str = "",
    sample_type: str = "",
    batch: str = "",
    notes: str = "",
) -> str:
    """Samples table has been merged into patients (migration 0008).

    Sample info is now stored in patients.metadata JSONB.

    Returns
    -------
    str
        Deprecation notice.
    """
    return (
        "Samples table has been removed (migration 0008). "
        "Sample info is now stored in patients.metadata JSONB. "
        "Use add_subject to register a patient with sample metadata."
    )


# ---------------------------------------------------------------------------
# Tool: list_samples (deprecated)
# ---------------------------------------------------------------------------


@mcp.tool()
def list_samples(
    project_name: str = "",
    subject_id: str = "",
    batch: str = "",
) -> str:
    """Samples table has been merged into patients (migration 0008).

    Returns
    -------
    str
        Deprecation notice.
    """
    return (
        "Samples table has been removed (migration 0008). "
        "Sample info is now in patients.metadata JSONB."
    )


# ---------------------------------------------------------------------------
# Tool: register_biodata (backward compat -> register_data)
# ---------------------------------------------------------------------------


@mcp.tool()
def register_biodata(
    project_name: str,
    category: str,
    platform: str,
    uri: str,
    fmt: str = "",
    status: str = "pending",
    file_role: str = "primary",
    phase: str = "",
    n_obs: int = 0,
    n_vars: int = 0,
) -> str:
    """Register a data object in the registry.

    Parameters
    ----------
    project_name
        Project name (must already exist).
    category
        Data category: spatial_seq | image | rnaseq | epigenomics | genome_seq.
    platform
        Platform slug (e.g. visium, imc, xenium, chromium_3p).
    uri
        Path or URI to the data file.
    fmt
        File format: h5ad, zarr, tiff, tsv, fastq, bam, bed.
    phase
        Pipeline phase slug (optional).

    Returns
    -------
    str
        Confirmation with assigned data id.
    """
    try:
        reg = _registry()
        bd_id = reg.register_data(
            project_name,
            phase=phase or "unknown",
            uri=uri,
            fmt=fmt or None,
            platform=platform,
            category=category,
            status=status,
            file_role=file_role,
            n_obs=n_obs or None,
            n_vars=n_vars or None,
        )
        return f"Data[{category}] registered (id={bd_id}) for '{project_name}' at {uri}"
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_biodata (backward compat -> list_datasets)
# ---------------------------------------------------------------------------


@mcp.tool()
def list_biodata(
    project_name: str = "",
    category: str = "",
    platform: str = "",
    modality: str = "",
) -> str:
    """List data objects in the registry, with optional filters.

    Parameters
    ----------
    project_name
        Filter by project name.
    category
        Filter by category: spatial_seq | image | rnaseq | epigenomics | genome_seq.
    platform
        Filter by platform slug (e.g. visium, imc).
    modality
        Filter by modality string.

    Returns
    -------
    str
        Formatted list of data objects.
    """
    try:
        reg = _registry()
        items = reg.list_biodata(
            project_name=project_name or None,
            category=category or None,
            platform=platform or None,
        )
        if modality:
            items = [bd for bd in items if bd.get("modality") == modality]
        if not items:
            return "No data objects found matching the given filters."
        lines = [f"Data objects ({len(items)} total):"]
        for bd in items:
            role = bd.get("file_role", "primary")
            mod = bd.get("modality", "")
            mod_str = f" modality={mod}" if mod else ""
            lines.append(
                f"  id={bd['id']} category={bd.get('category', '')} "
                f"platform={bd.get('platform', '')} "
                f"role={role} status={bd.get('status', '')}{mod_str}"
            )
            lines.append(f"    uri: {bd['uri']}")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_modalities
# ---------------------------------------------------------------------------


@mcp.tool()
def list_modalities(biodata_type: str = "") -> str:
    """List all known BioData modalities from the platform registry.

    Parameters
    ----------
    biodata_type
        Optional filter: spatial_seq | image | rnaseq | epigenomics | genome_seq.
        Empty means all.

    Returns
    -------
    str
        Formatted list of modalities with platform counts.
    """
    try:
        from sc_tools.biodata import list_modalities as _list_mods
        from sc_tools.biodata import list_platforms_by_modality

        mods = _list_mods(biodata_type=biodata_type or None)
        if not mods:
            return "No modalities found."
        lines = [f"BioData Modalities ({len(mods)} total):"]
        for mod in mods:
            platforms = list_platforms_by_modality(mod)
            plat_names = ", ".join(p.name for p in platforms[:8])
            if len(platforms) > 8:
                plat_names += f", ... (+{len(platforms) - 8} more)"
            lines.append(f"  {mod} ({len(platforms)} platforms)")
            lines.append(f"    {plat_names}")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: project_data_summary
# ---------------------------------------------------------------------------


@mcp.tool()
def project_data_summary(project_name: str) -> str:
    """Return a summary of data objects for a project, grouped by category and platform.

    Parameters
    ----------
    project_name
        Project name.

    Returns
    -------
    str
        Formatted summary with counts.
    """
    try:
        reg = _registry()
        summary = reg.project_data_summary(project_name)
        lines = [
            f"Data summary for '{project_name}':",
            f"  Total data objects: {summary['total']}",
        ]
        if summary["by_category"]:
            lines.append("  By category:")
            for cat, count in sorted(summary["by_category"].items()):
                lines.append(f"    {cat}: {count}")
        if summary.get("by_modality"):
            lines.append("  By modality:")
            for mod, count in sorted(summary["by_modality"].items()):
                lines.append(f"    {mod}: {count}")
        if summary["by_platform"]:
            lines.append("  By platform:")
            for plat, count in sorted(summary["by_platform"].items()):
                lines.append(f"    {plat}: {count}")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    mcp.run()
