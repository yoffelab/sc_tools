"""MCP server: sc-registry bookkeeping.

Exposes the sc_tools registry as callable MCP tools so Claude Code can
query project state, checkpoint locations, SLURM jobs, and agent tasks.

Tools
-----
    registry_status      -- high-level summary (projects, jobs, tasks, phase counts)
    list_datasets        -- all checkpoints for a project with URIs
    get_checkpoint_uri   -- "where is adata.normalized.h5ad for ggo_visium?"
    register_dataset     -- add a new checkpoint to the registry
    list_slurm_jobs      -- running/recent SLURM jobs with statuses
    list_agent_tasks     -- running agent task log
    mark_phase_complete  -- mark a pipeline phase as complete for a project
    get_phase_status     -- status and checkpoint details for a specific phase
    set_phase_status     -- update the pipeline status for a phase
    add_subject          -- register a de-identified subject
    list_subjects        -- query subjects
    add_sample           -- register a sample
    list_samples         -- query samples
    register_biodata     -- register a typed BioData object
    list_biodata         -- query BioData objects
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

    Shows total projects, datasets, active SLURM jobs, running agent tasks,
    the list of active project names, and per-project phase completion counts.

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
            f"  Datasets          : {s['n_datasets']}",
            f"  Subjects          : {s['n_subjects']}",
            f"  Samples           : {s['n_samples']}",
            f"  BioData objects   : {s['n_biodata']}",
            f"  Active SLURM jobs : {s['active_slurm_jobs']}",
            f"  Running tasks     : {s['running_agent_tasks']}",
        ]
        if s["active_projects"]:
            lines.append("\n  Active projects:")
            for name in s["active_projects"]:
                proj = reg.get_project(name)
                phases = json.loads(proj.get("phases_complete") or "[]") if proj else []
                label = f" [phases complete: {', '.join(phases)}]" if phases else ""
                lines.append(f"    - {name}{label}")
                # Per-phase summary from project_phases table
                phase_summary = s.get("phase_summary", {}).get(name, {})
                if phase_summary:
                    parts = ", ".join(f"{k}={v}" for k, v in sorted(phase_summary.items()))
                    lines.append(f"      phase counts: {parts}")
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
            validated = ds.get("validated", False)
            n_obs = ds.get("n_obs")
            obs_str = f" n_obs={n_obs}" if n_obs else ""
            lines.append(
                f"  id={ds['id']} phase={ds['phase']} role={role} "
                f"status={ds['status']} validated={validated}{obs_str}\n"
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
# Tool: list_slurm_jobs
# ---------------------------------------------------------------------------


@mcp.tool()
def list_slurm_jobs(active_only: bool = True) -> str:
    """List SLURM jobs tracked in the registry.

    Parameters
    ----------
    active_only
        If True (default), only show submitted/running jobs.
        If False, show all jobs.

    Returns
    -------
    str
        Formatted job list.
    """
    try:
        reg = _registry()
        if active_only:
            jobs = reg.list_active_jobs()
        else:
            with reg._session() as sess:
                jobs = [reg._to_dict(r) for r in sess.query(reg._SlurmJob).all()]
        if not jobs:
            return "No SLURM jobs found."
        lines = ["SLURM jobs:"]
        for j in jobs:
            lines.append(
                f"  id={j['id']} slurm_id={j['slurm_job_id']} "
                f"cluster={j['cluster']} phase={j['phase']} "
                f"status={j['status']}"
            )
            if j.get("sample_id"):
                lines[-1] += f" sample={j['sample_id']}"
            if j.get("error_msg"):
                lines.append(f"    error: {j['error_msg']}")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_agent_tasks
# ---------------------------------------------------------------------------


@mcp.tool()
def list_agent_tasks(running_only: bool = True) -> str:
    """List agent tasks tracked in the registry.

    Parameters
    ----------
    running_only
        If True (default), only show running tasks.
        If False, show all tasks.

    Returns
    -------
    str
        Formatted task list.
    """
    try:
        reg = _registry()
        if running_only:
            tasks = reg.list_running_tasks()
        else:
            with reg._session() as sess:
                tasks = [reg._to_dict(r) for r in sess.query(reg._AgentTask).all()]
        if not tasks:
            return "No agent tasks found."
        lines = ["Agent tasks:"]
        for t in tasks:
            lines.append(
                f"  id={t['id']} type={t['task_type']} "
                f"status={t['status']} started={t.get('started_at', '')}"
            )
            if t.get("error"):
                lines.append(f"    error: {t['error']}")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


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
    data_type
        Data type (may match platform for single-modality projects).
    domain
        High-level domain: spatial_transcriptomics | spatial_proteomics |
        imaging | single_cell | bulk.
    imaging_modality
        Imaging modality: brightfield | fluorescence | multiplexed_fluorescence |
        probe_based | mass_spec_imaging | sequencing_based.
    project_type
        internal (lab-led, default) or external (collaboration / contract).
    visibility
        private (restricted access, default) or public (shareable / published).

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
            data_type=data_type or None,
            domain=domain or None,
            imaging_modality=imaging_modality or None,
            project_type=project_type,
            visibility=visibility,
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
            lines.append("    (none — all phases complete or blocked)")
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

    Call this after environment_info() to link the compute environment to
    the phase that produced a checkpoint.  The snapshot is stored in the
    notes field of the project_phases row.

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
# Tool: update_slurm_job_status
# ---------------------------------------------------------------------------


@mcp.tool()
def update_slurm_job_status(
    slurm_job_id: str,
    status: str,
    error: str = "",
    log_uri: str = "",
) -> str:
    """Update the status of a tracked SLURM job.

    Parameters
    ----------
    slurm_job_id
        The SLURM job ID string (as returned by sbatch, e.g. "12345678").
    status
        New status: submitted | running | completed | failed.
    error
        Optional error message (for failed jobs).
    log_uri
        Optional URI to the job log file.

    Returns
    -------
    str
        Confirmation message.
    """
    try:
        reg = _registry()
        reg.update_job_status(
            slurm_job_id,
            status,
            error=error or None,
            log_uri=log_uri or None,
        )
        return f"SLURM job {slurm_job_id} status updated to '{status}'."
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: mark_phase_complete
# ---------------------------------------------------------------------------


@mcp.tool()
def mark_phase_complete(project_name: str, phase: str) -> str:
    """Mark a pipeline phase as complete for a project.

    Updates both the legacy ``phases_complete`` JSON list and the
    ``project_phases`` table.

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

    Creates the row if it does not exist (upsert).

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

    Data sources are distinct from projects (lab work) and datasets (processed
    checkpoints). They represent raw input data: HPC directories, public dataset
    portals, GEO accessions, Zenodo records, etc.

    Parameters
    ----------
    name
        Unique identifier (e.g. saha_ibd_cosmx, 10x_visium_human_breast).
    uri
        Primary location: HPC path (cayuga:/athena/project-saha/data_IBD),
        URL, GEO accession, DOI, etc.
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
            f"Linked project '{project_name}' → data source '{data_source_name}' "
            f"(role={role}, link_id={link_id})."
        )
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: add_subject
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
    """Register a de-identified subject in the registry.

    Subjects are cross-project and can be linked to multiple projects.
    No PHI (names, DOB, MRNs) should be stored.

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
        Confirmation with assigned subject DB id.
    """
    try:
        reg = _registry()
        kwargs: dict = {}
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
        sid = reg.add_subject(subject_id, organism=organism, **kwargs)
        return f"Subject '{subject_id}' registered (id={sid})."
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_subjects
# ---------------------------------------------------------------------------


@mcp.tool()
def list_subjects(
    project_name: str = "",
    diagnosis: str = "",
    tissue: str = "",
) -> str:
    """List subjects in the registry, with optional filters.

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
        Formatted list of subjects.
    """
    try:
        reg = _registry()
        subjects = reg.list_subjects(
            project_name=project_name or None,
            diagnosis=diagnosis or None,
            tissue=tissue or None,
        )
        if not subjects:
            return "No subjects found matching the given filters."
        lines = [f"Subjects ({len(subjects)} total):"]
        for s in subjects:
            parts = [f"  {s['subject_id']}"]
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
# Tool: add_sample
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
    """Register a sample in the registry.

    Samples link subjects to physical specimens within projects.

    Parameters
    ----------
    sample_id
        Sample identifier (e.g. S001_A1).
    subject_id
        De-identified subject ID (optional, must already exist).
    project_name
        Project name (optional, must already exist).
    tissue
        Tissue type (e.g. colon, lymph_node, lung).
    fixation_method
        FFPE, fresh_frozen, OCT, PFA.
    sample_type
        biopsy, resection, TMA, organoid, xenograft.
    batch
        Experimental batch identifier.

    Returns
    -------
    str
        Confirmation with assigned sample DB id.
    """
    try:
        reg = _registry()
        kwargs: dict = {}
        if tissue:
            kwargs["tissue"] = tissue
        if tissue_region:
            kwargs["tissue_region"] = tissue_region
        if fixation_method:
            kwargs["fixation_method"] = fixation_method
        if sample_type:
            kwargs["sample_type"] = sample_type
        if batch:
            kwargs["batch"] = batch
        if notes:
            kwargs["notes"] = notes
        sid = reg.add_sample(
            sample_id,
            subject_id=subject_id or None,
            project_name=project_name or None,
            **kwargs,
        )
        return f"Sample '{sample_id}' registered (id={sid})."
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_samples
# ---------------------------------------------------------------------------


@mcp.tool()
def list_samples(
    project_name: str = "",
    subject_id: str = "",
    batch: str = "",
) -> str:
    """List samples in the registry, with optional filters.

    Parameters
    ----------
    project_name
        Filter by project name.
    subject_id
        Filter by subject de-identified ID.
    batch
        Filter by batch identifier.

    Returns
    -------
    str
        Formatted list of samples.
    """
    try:
        reg = _registry()
        samples = reg.list_samples(
            project_name=project_name or None,
            subject_id=subject_id or None,
            batch=batch or None,
        )
        if not samples:
            return "No samples found matching the given filters."
        lines = [f"Samples ({len(samples)} total):"]
        for s in samples:
            parts = [f"  id={s['id']} {s['sample_id']}"]
            if s.get("tissue"):
                parts.append(f"tissue={s['tissue']}")
            if s.get("batch"):
                parts.append(f"batch={s['batch']}")
            if s.get("fixation_method"):
                parts.append(f"fix={s['fixation_method']}")
            lines.append(" ".join(parts))
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: register_biodata
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
    """Register a typed BioData object in the registry.

    The category determines the data type: spatial_seq, image, rnaseq,
    epigenomics, or genome_seq. Default field values are auto-filled from
    the platform registry when available.

    Parameters
    ----------
    project_name
        Project name (must already exist).
    category
        BioData type: spatial_seq | image | rnaseq | epigenomics | genome_seq.
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
        Confirmation with assigned BioData id.
    """
    try:
        reg = _registry()
        bd_id = reg.register_biodata(
            project_name,
            category,
            platform,
            uri,
            fmt=fmt or None,
            status=status,
            file_role=file_role,
            phase=phase or None,
            n_obs=n_obs or None,
            n_vars=n_vars or None,
        )
        return f"BioData[{category}] registered (id={bd_id}) for '{project_name}' at {uri}"
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_biodata
# ---------------------------------------------------------------------------


@mcp.tool()
def list_biodata(
    project_name: str = "",
    category: str = "",
    platform: str = "",
    modality: str = "",
) -> str:
    """List BioData objects in the registry, with optional filters.

    Parameters
    ----------
    project_name
        Filter by project name.
    category
        Filter by category: spatial_seq | image | rnaseq | epigenomics | genome_seq.
    platform
        Filter by platform slug (e.g. visium, imc).
    modality
        Filter by modality string (e.g. "Spatial Proteomics - Mass Spec").

    Returns
    -------
    str
        Formatted list of BioData objects.
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
            return "No BioData objects found matching the given filters."
        lines = [f"BioData objects ({len(items)} total):"]
        for bd in items:
            role = bd.get("file_role", "primary")
            mod = bd.get("modality", "")
            mod_str = f" modality={mod}" if mod else ""
            lines.append(
                f"  id={bd['id']} type={bd['type']} platform={bd.get('platform', '')} "
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

    The modality tier sits between BioDataType (5 JTI types) and platform slug,
    grouping related platforms (e.g. "Spatial Proteomics - Mass Spec" groups
    imc, mibi, maldi_ims).

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
    """Return a summary of BioData objects for a project, grouped by category and platform.

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
            f"  Total BioData objects: {summary['total']}",
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
