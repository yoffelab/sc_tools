"""MCP server: sc-registry bookkeeping.

Exposes the sc_tools registry as callable MCP tools so Claude Code can
query project state, checkpoint locations, and data objects.

Tools
-----
    registry_status         -- high-level summary (projects, data sources, inventory, datasets)
    register_inventory_item -- add an inventory item (Layer 1)
    list_inventory_items    -- query inventory items with filters
    create_dataset          -- create a named dataset (Layer 2)
    add_dataset_member      -- add an inventory item to a dataset
    link_project_dataset    -- link a project to a dataset
    list_datasets           -- datasets linked to a project
    record_provenance       -- record transformation provenance
    get_provenance          -- query provenance records
    set_phase_status        -- upsert a pipeline phase status
    get_phase_status        -- status and details for a specific phase
    mark_phase_complete     -- mark a pipeline phase as complete
    get_available_next_phases -- what phases can run next
    project_data_summary    -- dataset/inventory summary for a project
    add_project             -- register a new project
    register_data_source    -- register a raw data source
    list_data_sources       -- query data sources
    add_subject / list_subjects -- patient management
    list_modalities         -- BioData modality catalog
    register_dataset        -- backward compat wrapper (deprecated)
    register_biodata        -- backward compat wrapper (deprecated)
    list_biodata            -- backward compat wrapper (deprecated)
    get_checkpoint_uri      -- deprecated, use get_phase_status
    link_project_data_source -- deprecated, use link_project_dataset

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

    Shows total projects, data sources, inventory items, datasets, and patients,
    plus per-project phase status.

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
            f"  Data sources      : {s['n_data_sources']}",
            f"  Inventory items   : {s['n_inventory_items']}",
            f"  Datasets          : {s['n_datasets']}",
            f"  Patients          : {s['n_patients']}",
        ]
        if s["active_projects"]:
            lines.append("\n  Active projects:")
            for name in s["active_projects"]:
                lines.append(f"    - {name}")
                phase_summary = s.get("phase_summary", {}).get(name, {})
                if phase_summary:
                    parts = ", ".join(f"{k}={v}" for k, v in sorted(phase_summary.items()))
                    lines.append(f"      phases: {parts}")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: register_inventory_item
# ---------------------------------------------------------------------------


@mcp.tool()
def register_inventory_item(
    name: str,
    uri: str,
    modality: str,
    data_source_name: str = "",
    platform: str = "",
    fmt: str = "h5ad",
    n_obs: int = 0,
    n_vars: int = 0,
    size_mb: float = 0.0,
    organism: str = "",
    tissue: str = "",
) -> str:
    """Register an inventory item (Layer 1 data object) in the registry.

    Inventory items represent individual data files (AnnData, images, etc.)
    that can later be assembled into multi-modal datasets.

    Parameters
    ----------
    name
        Unique identifier for this item (e.g. ibd_cosmx_rna, robin_visium_raw).
    uri
        Path or URI to the data file.
    modality
        Data modality (e.g. spatial_transcriptomics, imaging_mass_cytometry).
    data_source_name
        Name of the data source this item comes from (must exist). Empty to skip.
    platform
        Technology platform (e.g. cosmx, visium, xenium).
    fmt
        File format: h5ad (default), zarr, tiff, tsv, mudata.
    n_obs
        Number of observations/cells/spots (0 to skip).
    n_vars
        Number of variables/genes/proteins (0 to skip).
    size_mb
        File size in megabytes (0.0 to skip).
    organism
        Organism (e.g. human, mouse).
    tissue
        Tissue type (e.g. colon, brain, lung).

    Returns
    -------
    str
        Confirmation with assigned inventory item id.
    """
    try:
        reg = _registry()
        item_id = reg.register_inventory_item(
            name,
            uri,
            modality,
            data_source_name=data_source_name or None,
            platform=platform or None,
            fmt=fmt,
            n_obs=n_obs or None,
            n_vars=n_vars or None,
            size_mb=size_mb or None,
            organism=organism or None,
            tissue=tissue or None,
        )
        return f"Inventory item '{name}' registered (id={item_id})."
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_inventory_items
# ---------------------------------------------------------------------------


@mcp.tool()
def list_inventory_items(modality: str = "", platform: str = "") -> str:
    """List inventory items in the registry, with optional filters.

    Parameters
    ----------
    modality
        Filter by modality string (e.g. spatial_transcriptomics). Empty means all.
    platform
        Filter by platform string (e.g. cosmx, visium). Empty means all.

    Returns
    -------
    str
        Formatted table of inventory items.
    """
    try:
        reg = _registry()
        items = reg.list_inventory_items(
            modality=modality or None,
            platform=platform or None,
        )
        if not items:
            return "No inventory items found matching the given filters."
        lines = [f"Inventory items ({len(items)} total):"]
        lines.append(
            f"  {'ID':<4}  {'Name':<40}  {'Modality':<28}  {'Platform':<14}  {'Format':<8}  {'n_obs'}"
        )
        lines.append("  " + "-" * 110)
        for it in items:
            n_obs_str = str(it.get("n_obs", "")) if it.get("n_obs") else ""
            lines.append(
                f"  {it['id']:<4}  {(it.get('name') or ''):<40}  "
                f"{(it.get('modality') or ''):<28}  "
                f"{(it.get('platform') or ''):<14}  "
                f"{(it.get('format') or ''):<8}  "
                f"{n_obs_str}"
            )
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: create_dataset
# ---------------------------------------------------------------------------


@mcp.tool()
def create_dataset(
    name: str,
    description: str = "",
    fmt: str = "mudata",
) -> str:
    """Create a named dataset (Layer 2) that groups inventory items.

    Datasets represent multi-modal collections of inventory items (e.g.
    RNA + protein from the same experiment). Use add_dataset_member to
    populate the dataset after creation.

    Parameters
    ----------
    name
        Unique dataset name (e.g. ibd_cosmx_multimodal, robin_discovery).
    description
        Optional human-readable description.
    fmt
        Format: mudata (default), h5ad, zarr.

    Returns
    -------
    str
        Confirmation with assigned dataset id.
    """
    try:
        reg = _registry()
        ds_id = reg.create_dataset(
            name,
            description=description or None,
            fmt=fmt,
        )
        return f"Dataset '{name}' created (id={ds_id})."
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: add_dataset_member
# ---------------------------------------------------------------------------


@mcp.tool()
def add_dataset_member(
    dataset_name: str,
    inventory_name: str,
    modality_key: str,
) -> str:
    """Add an inventory item to a dataset as a named modality slot.

    Each modality_key must be unique within the dataset (e.g. 'rna',
    'protein', 'morphology').

    Parameters
    ----------
    dataset_name
        Name of the target dataset (must exist, uses current version).
    inventory_name
        Name of the inventory item to add (must exist).
    modality_key
        Key for this modality within the dataset (e.g. 'rna', 'protein').

    Returns
    -------
    str
        Confirmation with assigned member id.
    """
    try:
        reg = _registry()
        member_id = reg.add_dataset_member(dataset_name, inventory_name, modality_key)
        return (
            f"Added '{inventory_name}' as '{modality_key}' to dataset "
            f"'{dataset_name}' (member_id={member_id})."
        )
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: link_project_dataset
# ---------------------------------------------------------------------------


@mcp.tool()
def link_project_dataset(
    project_name: str,
    dataset_name: str,
    role: str = "primary",
    notes: str = "",
) -> str:
    """Link a project to a dataset.

    A project can have multiple datasets (e.g. primary analysis dataset,
    reference scRNA dataset for deconvolution).

    Parameters
    ----------
    project_name
        Project name (must already exist in registry).
    dataset_name
        Dataset name (must already exist, uses current version).
    role
        primary | reference | supplementary.
    notes
        Optional free-text notes about the link.

    Returns
    -------
    str
        Confirmation with assigned link id.
    """
    try:
        reg = _registry()
        link_id = reg.link_project_dataset(
            project_name, dataset_name, role=role, notes=notes or None
        )
        return (
            f"Linked project '{project_name}' -> dataset '{dataset_name}' "
            f"(role={role}, link_id={link_id})."
        )
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_datasets
# ---------------------------------------------------------------------------


@mcp.tool()
def list_datasets(project_name: str = "") -> str:
    """List datasets, optionally filtered by project.

    Parameters
    ----------
    project_name
        Filter by project name. Empty means all datasets.

    Returns
    -------
    str
        Table of datasets with name, version, format, and is_current flag.
    """
    try:
        reg = _registry()
        datasets = reg.list_datasets(
            project_name=project_name or None,
        )
        if not datasets:
            msg = "No datasets found"
            if project_name:
                msg += f" for project '{project_name}'"
            return msg + "."
        lines = [f"Datasets ({len(datasets)} total):"]
        lines.append(f"  {'ID':<4}  {'Name':<40}  {'Ver':<4}  {'Format':<10}  {'Current'}")
        lines.append("  " + "-" * 72)
        for ds in datasets:
            lines.append(
                f"  {ds['id']:<4}  {(ds.get('name') or ''):<40}  "
                f"v{ds.get('version', 1):<3}  "
                f"{(ds.get('format') or ''):<10}  "
                f"{'yes' if ds.get('is_current') else 'no'}"
            )
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: record_provenance
# ---------------------------------------------------------------------------


@mcp.tool()
def record_provenance(
    tool: str,
    tool_version: str = "",
    phase_id: int = 0,
    dataset_id: int = 0,
    reference_genome: str = "",
    reference_dataset: str = "",
    signature_source: str = "",
    params_json: str = "",
    environment_json: str = "",
    script_uri: str = "",
    agent: str = "",
) -> str:
    """Record provenance metadata for a transformation or analysis step.

    Exactly one of phase_id or dataset_id must be non-zero
    to identify the target of this provenance record.

    Parameters
    ----------
    tool
        Name of the tool or function (e.g. scanpy.pp.normalize_total, scvi.SCVI).
    tool_version
        Version string of the tool.
    phase_id
        Target phase id (0 to skip).
    dataset_id
        Target dataset id (0 to skip).
    reference_genome
        Reference genome used (e.g. GRCh38).
    reference_dataset
        Reference dataset used (e.g. HuBMAP gut atlas).
    signature_source
        Source of gene signatures (e.g. MSigDB, CellTypist).
    params_json
        JSON string of parameters used (e.g. '{"n_top_genes": 2000}').
    environment_json
        JSON string of environment snapshot (e.g. '{"python": "3.11", ...}').
    script_uri
        URI of the script that performed the transformation.
    agent
        Name of the agent that performed the transformation.

    Returns
    -------
    str
        Confirmation with assigned provenance id.
    """
    try:
        # Parse JSON strings
        params = json.loads(params_json) if params_json else None
        environment = json.loads(environment_json) if environment_json else None

        reg = _registry()
        prov_id = reg.record_provenance(
            tool,
            tool_version=tool_version or None,
            phase_id=phase_id or None,
            dataset_id=dataset_id or None,
            reference_genome=reference_genome or None,
            reference_dataset=reference_dataset or None,
            signature_source=signature_source or None,
            params=params,
            environment=environment,
            script_uri=script_uri or None,
            agent=agent or None,
        )
        return f"Provenance recorded (id={prov_id}) for tool='{tool}'."
    except json.JSONDecodeError as exc:
        return f"ERROR: Invalid JSON: {exc}"
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: get_provenance
# ---------------------------------------------------------------------------


@mcp.tool()
def get_provenance(
    phase_id: int = 0,
    dataset_id: int = 0,
) -> str:
    """Query provenance records for a specific target.

    Provide one of phase_id or dataset_id to filter.
    If both are 0, returns all provenance records.

    Parameters
    ----------
    phase_id
        Filter by phase id (0 to skip).
    dataset_id
        Filter by dataset id (0 to skip).

    Returns
    -------
    str
        Formatted list of provenance records.
    """
    try:
        reg = _registry()
        records = reg.get_provenance(
            phase_id=phase_id or None,
            dataset_id=dataset_id or None,
        )
        if not records:
            return "No provenance records found."
        lines = [f"Provenance records ({len(records)} total):"]
        for r in records:
            lines.append(
                f"  id={r['id']} tool={r.get('tool', '')} "
                f"target={r.get('target_type', '')} "
                f"version={r.get('tool_version', '') or ''}"
            )
            if r.get("script_uri"):
                lines.append(f"    script: {r['script_uri']}")
            if r.get("agent"):
                lines.append(f"    agent: {r['agent']}")
            if r.get("params"):
                lines.append(f"    params: {json.dumps(r['params'])}")
            if r.get("created_at"):
                lines.append(f"    created: {r['created_at']}")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: set_phase_status
# ---------------------------------------------------------------------------


@mcp.tool()
def set_phase_status(
    project_name: str,
    dataset_name: str,
    phase_group: str,
    subphase: str,
    status: str,
    notes: str = "",
    n_obs: int = 0,
    n_vars: int = 0,
) -> str:
    """Update the pipeline status for a phase (upsert).

    Creates a phase row if none exists for this project + dataset +
    phase_group + subphase combination.

    Parameters
    ----------
    project_name
        Project name (e.g. ggo_visium).
    dataset_name
        Dataset name (must exist, uses current version).
    phase_group
        Phase group: data_processing | discovery.
    subphase
        Subphase slug (e.g. qc_filter, clustering_v1, scoring).
    status
        New status: not_started | in_progress | complete | failed | skipped | ready.
    notes
        Optional free-text notes (e.g. method selection rationale).
    n_obs
        Number of cells/spots at this phase (0 to skip).
    n_vars
        Number of genes/proteins at this phase (0 to skip).

    Returns
    -------
    str
        Confirmation message.
    """
    try:
        reg = _registry()
        reg.upsert_phase(
            project_name,
            dataset_name,
            phase_group,
            subphase,
            status=status,
            notes=notes or None,
            n_obs=n_obs or None,
            n_vars=n_vars or None,
        )
        parts = [
            f"Set phase '{phase_group}/{subphase}' to status='{status}' "
            f"for project '{project_name}' dataset '{dataset_name}'."
        ]
        if n_obs:
            parts.append(f"n_obs={n_obs:,}")
        if n_vars:
            parts.append(f"n_vars={n_vars:,}")
        return " ".join(parts)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: get_phase_status
# ---------------------------------------------------------------------------


@mcp.tool()
def get_phase_status(
    project_name: str,
    phase_group: str,
    subphase: str,
) -> str:
    """Return the status and details for a specific phase of a project.

    Parameters
    ----------
    project_name
        Project name (e.g. ggo_visium).
    phase_group
        Phase group: data_processing | discovery.
    subphase
        Subphase slug (e.g. qc_filter, clustering_v1).

    Returns
    -------
    str
        Formatted phase details including status, n_obs, n_vars, notes.
    """
    try:
        reg = _registry()
        row = reg.get_phase(project_name, phase_group, subphase)
        if row is None:
            return (
                f"No phase entry for project='{project_name}' "
                f"phase='{phase_group}/{subphase}'. "
                "Use set_phase_status to create one."
            )
        lines = [f"Phase '{phase_group}/{subphase}' for project '{project_name}':"]
        lines.append(f"  status        : {row['status']}")
        if row.get("uri"):
            lines.append(f"  uri           : {row['uri']}")
        if row.get("n_obs"):
            lines.append(f"  n_obs         : {row['n_obs']}")
        if row.get("n_vars"):
            lines.append(f"  n_vars        : {row['n_vars']}")
        if row.get("dataset_id"):
            lines.append(f"  dataset_id    : {row['dataset_id']}")
        if row.get("notes"):
            lines.append(f"  notes         : {row['notes']}")
        lines.append(f"  updated_at    : {row.get('updated_at', '')}")
        return "\n".join(lines)
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: mark_phase_complete
# ---------------------------------------------------------------------------


@mcp.tool()
def mark_phase_complete(
    project_name: str,
    phase_group: str,
    subphase: str,
) -> str:
    """Mark a pipeline phase as complete (status='ready') for a project.

    Parameters
    ----------
    project_name
        Project name (e.g. ggo_visium).
    phase_group
        Phase group: data_processing | discovery.
    subphase
        Subphase slug (e.g. qc_filter, clustering_v1).

    Returns
    -------
    str
        Confirmation message.
    """
    try:
        reg = _registry()
        reg.mark_phase_complete(project_name, phase_group, subphase)
        return f"Marked phase '{phase_group}/{subphase}' as complete for project '{project_name}'."
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
        completed = [r["subphase"] for r in phase_rows if r["status"] in ("complete", "ready")]
        in_progress = [r["subphase"] for r in phase_rows if r["status"] == "in_progress"]

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
# Tool: project_data_summary
# ---------------------------------------------------------------------------


@mcp.tool()
def project_data_summary(project_name: str) -> str:
    """Return a summary of datasets and inventory items for a project.

    Parameters
    ----------
    project_name
        Project name.

    Returns
    -------
    str
        Formatted summary with dataset details and member counts.
    """
    try:
        reg = _registry()
        summary = reg.project_data_summary(project_name)
        lines = [
            f"Data summary for '{project_name}':",
            f"  Datasets        : {summary['n_datasets']}",
            f"  Inventory items : {summary['n_inventory_items']}",
        ]
        for ds in summary.get("datasets", []):
            current = " (current)" if ds.get("is_current") else ""
            lines.append(
                f"\n  Dataset: {ds['name']} v{ds.get('version', 1)}{current} "
                f"[{ds.get('format', '')}] role={ds.get('role', '')}"
            )
            lines.append(f"    Members ({ds.get('n_members', 0)}):")
            for m in ds.get("members", []):
                n_obs_str = f" n_obs={m['n_obs']}" if m.get("n_obs") else ""
                lines.append(
                    f"      {m['modality_key']}: {m.get('inventory_name', '')} "
                    f"({m.get('modality', '')}){n_obs_str}"
                )
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
# Tool: register_dataset (backward compat wrapper -- DEPRECATED)
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
    """[DEPRECATED] Register a checkpoint dataset. Use register_inventory_item instead.

    This is a backward-compatibility wrapper. It calls the deprecated
    Registry.register_dataset() method which internally creates an
    inventory item.

    Parameters
    ----------
    project_name
        Project name.
    phase
        Phase slug.
    uri
        Path or URI to the checkpoint file.
    fmt
        File format (default h5ad).
    sample_id
        Optional sample identifier.
    status
        Dataset status (default ready).
    file_role
        File role (default primary).
    validated
        Whether checkpoint has been validated.
    n_obs
        Number of observations (0 to skip).
    n_vars
        Number of variables (0 to skip).

    Returns
    -------
    str
        Confirmation with assigned id. Includes deprecation notice.
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
        return (
            f"[DEPRECATED] Registered dataset id={ds_id}: {project_name}/{phase} "
            f"[{file_role}] at {uri}. "
            "Please migrate to register_inventory_item + create_dataset."
        )
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: register_biodata (backward compat wrapper -- DEPRECATED)
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
    """[DEPRECATED] Register a data object. Use register_inventory_item instead.

    This is a backward-compatibility wrapper.

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
        Confirmation with assigned data id. Includes deprecation notice.
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
        return (
            f"[DEPRECATED] Data[{category}] registered (id={bd_id}) "
            f"for '{project_name}' at {uri}. "
            "Please migrate to register_inventory_item."
        )
    except Exception as exc:
        return f"ERROR: {exc}"


# ---------------------------------------------------------------------------
# Tool: list_biodata (backward compat wrapper -- DEPRECATED)
# ---------------------------------------------------------------------------


@mcp.tool()
def list_biodata(
    project_name: str = "",
    category: str = "",
    platform: str = "",
    modality: str = "",
) -> str:
    """[DEPRECATED] List data objects. Use list_inventory_items instead.

    This is a backward-compatibility wrapper.

    Parameters
    ----------
    project_name
        Filter by project name.
    category
        Filter by category.
    platform
        Filter by platform slug.
    modality
        Filter by modality string.

    Returns
    -------
    str
        Formatted list of data objects. Includes deprecation notice.
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
            return "[DEPRECATED] No data objects found. Use list_inventory_items instead."
        lines = [f"[DEPRECATED] Data objects ({len(items)} total) -- use list_inventory_items:"]
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
# Tool: link_project_data_source (DEPRECATED stub)
# ---------------------------------------------------------------------------


@mcp.tool()
def link_project_data_source(
    project_name: str,
    data_source_name: str,
    role: str = "input",
    notes: str = "",
) -> str:
    """[DEPRECATED] Link a project to a data source.

    This tool is deprecated. Use link_project_dataset instead.
    Data sources are now linked to projects through datasets:
    register_inventory_item -> create_dataset -> add_dataset_member -> link_project_dataset.

    Parameters
    ----------
    project_name
        Project name.
    data_source_name
        Data source name.
    role
        Link role.
    notes
        Optional notes.

    Returns
    -------
    str
        Deprecation notice.
    """
    return (
        "Deprecated. Use link_project_dataset instead. "
        "Data sources are now linked to projects through the inventory/dataset layer: "
        "register_inventory_item -> create_dataset -> add_dataset_member -> link_project_dataset."
    )


# ---------------------------------------------------------------------------
# Tool: get_checkpoint_uri (DEPRECATED stub)
# ---------------------------------------------------------------------------


@mcp.tool()
def get_checkpoint_uri(project_name: str, phase: str = "", sample_id: str = "") -> str:
    """[DEPRECATED] Get checkpoint URI. Use get_phase_status instead.

    Checkpoint URIs are now stored in the project_phases table. Use
    get_phase_status to retrieve phase details including URI.

    Parameters
    ----------
    project_name
        Project name.
    phase
        Phase slug.
    sample_id
        Sample identifier (ignored in new schema).

    Returns
    -------
    str
        Deprecation notice with migration guidance.
    """
    return (
        "Deprecated. Checkpoint URIs are now stored in the project_phases table. "
        "Use get_phase_status(project_name, phase_group, subphase) to retrieve "
        "phase details including the URI."
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    mcp.run()
