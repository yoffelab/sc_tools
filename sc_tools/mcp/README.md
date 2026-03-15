# sc_tools MCP Servers

Two Model Context Protocol (MCP) servers expose sc_tools capabilities to Claude Code
and other MCP clients via stdio JSON-RPC (FastMCP default transport).

| Server | Name | Purpose |
|--------|------|---------|
| `tools_server.py` | `sc-tools` | Analysis tools: validate checkpoints, run QC, score signatures, generate HPC scripts |
| `registry_server.py` | `sc-registry` | Bookkeeping: query projects, datasets, phases, subjects, samples, SLURM jobs |

---

## Installation

```bash
# sc-tools server only
pip install "sc-tools[mcp]"

# sc-registry server (includes registry dependencies)
pip install "sc-tools[mcp,registry]"
```

---

## Starting the servers

```bash
# Analysis tools server
python -m sc_tools.mcp.tools_server

# Registry bookkeeping server
python -m sc_tools.mcp.registry_server
```

Both servers communicate over stdio JSON-RPC.  Add them to your MCP client
configuration (e.g. Claude Desktop `claude_desktop_config.json`) as stdio
transport servers.

---

## sc-tools server tools

Analysis tools for running pipeline steps and inspecting data.

| Tool | Description |
|------|-------------|
| `validate_checkpoint` | Validate an AnnData checkpoint against the sc_tools metadata contract for a given phase. |
| `generate_qc_report` | Generate a shell command to produce a date-versioned QC HTML report (pre_filter, post_filter, post_integration, post_celltyping). |
| `score_gene_signatures` | Generate a shell command to score gene signatures on an AnnData checkpoint and write obsm results. |
| `run_deconvolution` | Generate a shell command to run cell-type deconvolution (cell2location, tangram, destvi) on a spatial AnnData. |
| `generate_sbatch_spaceranger` | Generate a Space Ranger SLURM sbatch script for a Visium or Visium HD sample. |
| `generate_sbatch_imc` | Generate an IMC pipeline SLURM sbatch script for a given sample and MCD file. |
| `generate_sbatch_xenium` | Generate a Xenium Ranger SLURM sbatch script for a raw Xenium bundle. |
| `collect_batch_manifests` | Collect per-batch TSV files in a metadata/phase0/ directory into all_samples.tsv. |
| `load_sample` | Load one sample into AnnData via the appropriate modality loader and report its shape and keys. |
| `inspect_checkpoint` | Inspect an AnnData checkpoint structure (shape, obs columns, obsm/uns/layers keys) without loading arrays into memory. |
| `list_signatures` | List gene signatures available in a JSON gene sets file, grouped by pathway/category with gene counts. |
| `environment_info` | Report the current runtime environment: Python version, key package versions, GPU devices, and rapids/dask availability. |

---

## sc-registry server tools

Bookkeeping tools for querying and updating the sc_tools registry database.

| Tool | Description |
|------|-------------|
| `registry_status` | Return a high-level status summary of the registry: project counts, dataset counts, active SLURM jobs, running agent tasks, and per-project phase completion. |
| `list_datasets` | List all checkpoint datasets for a project, optionally filtered by phase, with URIs and validation status. |
| `get_checkpoint_uri` | Return the URI for a specific project-phase checkpoint (or per-sample ingest_load checkpoint). |
| `register_dataset` | Add a new checkpoint to the registry with phase, URI, format, role, and validation metadata. |
| `list_slurm_jobs` | List SLURM jobs tracked in the registry, optionally filtering to active (submitted/running) jobs only. |
| `list_agent_tasks` | List agent tasks tracked in the registry, optionally filtering to running tasks only. |
| `add_project` | Register a new project in the registry (idempotent: returns existing id if already present). |
| `get_available_next_phases` | Return which pipeline phases are available to run next for a project, based on completed phases and the pipeline DAG. |
| `record_provenance` | Store a runtime environment snapshot alongside a pipeline phase record to track what compute environment produced a checkpoint. |
| `update_slurm_job_status` | Update the status of a tracked SLURM job (submitted, running, completed, failed). |
| `mark_phase_complete` | Mark a pipeline phase as complete for a project, updating both the legacy phases_complete list and the project_phases table. |
| `get_phase_status` | Return the status and checkpoint details (n_obs, n_vars, n_samples, notes) for a specific phase of a project. |
| `set_phase_status` | Update the pipeline status for a phase (upsert), optionally recording n_obs, n_vars, n_samples, and notes. |
| `register_data_source` | Register a raw data source in the catalog (HPC directory, GEO accession, Zenodo record, etc.). |
| `list_data_sources` | List data sources in the catalog, with optional filters by platform, source type, disease, or tissue. |
| `link_project_data_source` | Link a project to a data source with a role (input, reference, supplementary). |
| `add_subject` | Register a de-identified subject in the registry with clinical metadata (no PHI). |
| `list_subjects` | List subjects in the registry, with optional filters by project, diagnosis, or tissue. |
| `add_sample` | Register a sample linking a subject to a physical specimen within a project. |
| `list_samples` | List samples in the registry, with optional filters by project, subject, or batch. |
| `register_biodata` | Register a typed BioData object (spatial_seq, image, rnaseq, epigenomics, genome_seq) for a project. |
| `list_biodata` | List BioData objects in the registry, with optional filters by project, category, platform, or modality. |
| `list_modalities` | List all known BioData modalities from the platform registry, grouped by BioDataType. |
| `project_data_summary` | Return a summary of BioData objects for a project, grouped by category, modality, and platform with counts. |

---

## How the two servers work together

The two servers are designed to be used together within the same Claude Code session.
A typical analysis workflow uses both:

1. **sc-registry** provides checkpoint locations.  Call `get_checkpoint_uri` to
   retrieve the path to `adata.normalized.h5ad` for a project before passing it
   to an analysis tool.

2. **sc-tools** consumes those URIs.  Call `validate_checkpoint`, `inspect_checkpoint`,
   or `score_gene_signatures` with the URI returned by the registry.

3. After analysis completes, **sc-registry** records the result.  Call
   `register_dataset` to store the new checkpoint URI, then `set_phase_status` to
   advance the phase from `in_progress` to `complete` with cell and gene counts.

4. **sc-registry** then reports what can run next.  Call `get_available_next_phases`
   to read the pipeline DAG and determine which phases are now unblocked.

This separation of concerns keeps analysis logic (sc-tools) independent of
bookkeeping state (sc-registry), while allowing agents to coordinate both in a
single session.

---

## Transport

Both servers use FastMCP with the default **stdio JSON-RPC** transport.
No HTTP port is opened.  Clients communicate by launching the server process
and writing JSON-RPC messages to its stdin.

For Claude Desktop, add each server to `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "sc-tools": {
      "command": "python",
      "args": ["-m", "sc_tools.mcp.tools_server"]
    },
    "sc-registry": {
      "command": "python",
      "args": ["-m", "sc_tools.mcp.registry_server"]
    }
  }
}
```
