# Registry Guide

## 1. Overview

The sc_tools registry is a lightweight database that tracks what data you have, where it lives on disk, and what processing has been done to it. It does not store biological data itself -- every record is a pointer to a file (h5ad, h5mu, zarr) on an HPC scratch directory or cloud bucket. Think of it as a lab notebook for your computational projects: it records what exists, how it got there, and what state it is in.

The registry organizes information into four layers:

```
Layer 0: Data Sources          -- raw data dumps (GEO, HPC dirs, collaborator hand-offs)
    |
    v  (agent ingest)
Layer 1: Inventory Items       -- clean AnnData files, one per source per modality
    |
    v  (assemble)
Layer 2: Datasets              -- named, versioned collections of inventory items
    |
    v  (link to project)
Layer 3: Projects & Phases     -- analysis projects with tracked processing steps
```

Each layer builds on the one below it. Raw data enters at Layer 0, gets cleaned into Layer 1, gets assembled into Layer 2, and gets analyzed in Layer 3.


## 2. The 4 Layers

### Layer 0: Data Sources

A data source is a raw, uncurated data dump. It could be a GEO accession, an HPC directory from a collaborator, a 10x Genomics public dataset, or a download from a portal. You register these when you find or receive data, before any processing happens.

Data sources are registered with metadata about where they came from (platform, organism, tissue, disease) and how to access them (URI, source type, access notes).

**When things get added:** When you discover or receive a new dataset. This is the first thing you do.

**Example entries:**

| name | uri | platform | source_type |
|------|-----|----------|-------------|
| `10x_human_brain_xenium` | `https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-for-xenium-explorer-demo-1-standard` | xenium | public_10x |
| `GSE176171_ulcerative_colitis` | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176171` | chromium_3p | public_geo |
| `hubmap_intestine_codex` | `https://portal.hubmapconsortium.org/browse/dataset/abc123` | codex | public_portal |


### Layer 1: Inventory Items

An inventory item is a single clean AnnData file (h5ad) produced by processing a data source. Each item covers one modality from one source. If a data source contains both RNA and protein data (e.g., CITE-seq), it produces two inventory items.

Inventory items are the building blocks for datasets. They can be shared across multiple projects without duplication.

**When things get added:** When an agent processes raw data from a data source into a usable h5ad file.

**Example entries:**

| name | modality | platform | n_obs | n_vars |
|------|----------|----------|-------|--------|
| `10x_brain_xenium_rna` | rna | xenium | 150,000 | 313 |
| `GSE176171_uc_scrnaseq` | rna | chromium_3p | 300,000 | 20,000 |
| `hubmap_intestine_protein` | protein | codex | 50,000 | 200 |


### Layer 2: Datasets

A dataset is a named, versioned assembly of inventory items. It can be single-modality (one AnnData) or multi-modal (multiple inventory items composed into a MuData). This is what you actually work with in a project.

Datasets are logical assemblies. The `dataset_members` table defines which inventory items belong to a dataset, keyed by modality. Physical MuData files are materialized only when needed -- individual modalities can be loaded independently via their inventory URIs.

Datasets are versionable. When you need to add a modality or change composition, you bump the version. Old versions remain linked to any project phases that used them.

**When things get added:** When you decide which inventory items to combine for an analysis.

**Example entries:**

| name | version | format | members |
|------|---------|--------|---------|
| `uc_multimodal` | 1 | mudata | rna (10x_brain_xenium_rna), rna_ref (GSE176171_uc_scrnaseq) |
| `uc_multimodal` | 2 | mudata | rna, rna_ref, protein (hubmap_intestine_protein) |
| `brain_xenium` | 1 | h5ad | rna (10x_brain_xenium_rna) |


### Layer 3: Projects and Phases

A project is your analysis unit. It links to one or more datasets (with roles like `primary`, `reference`, `validation`) and tracks processing progress through phases.

Phases are organized into two groups:

- **`data_processing`** -- structured pipeline steps following the DAG: `ingest_raw` → `qc_filter` → `preprocess` → `scoring` → `celltype_manual` → `biology`, etc. These have dependencies defined in `pipeline.py`.
- **`discovery`** -- free-form iteration cycles with descriptive names (`clustering_v1`, `deconv_attempt_2`, `final_annotation`). No ordering enforced -- anything exploratory.

Each phase points to a checkpoint file (h5ad) on disk. The checkpoint contains the actual transformed data at that stage.

**When things get added:** When you create a project and start processing.

**Example entries:**

| project | dataset | phase_group | subphase | status | n_obs |
|---------|---------|-------------|----------|--------|-------|
| uc_atlas | uc_multimodal | data_processing | qc_filter | ready | 300,000 |
| uc_atlas | uc_multimodal | data_processing | preprocess | ready | 250,000 |
| uc_atlas | uc_multimodal | discovery | clustering_v1 | ready | 250,000 |
| uc_atlas | uc_multimodal | discovery | deconv_attempt_2 | in_progress | 250,000 |


## 3. Provenance

Provenance records track how each transformation was performed: what tool was used, which version, what parameters were passed, and what environment it ran in.

**What it tracks:**
- `tool` and `tool_version` -- e.g., scanpy 1.11.5, scvi-tools 1.4.2
- `params` -- the exact parameter dict (min_genes, n_top_genes, batch_key, etc.)
- `environment` -- container SHA, conda env, HPC node
- `reference_genome` -- if applicable
- `script_uri` -- path to the script that ran the transformation
- `agent` -- which Claude Code agent performed the work

**Why it matters:** When you write a methods section six months later, you need exact tool versions and parameters. When a result looks wrong, you need to know what produced it. Provenance gives you both.

**Where it attaches:**
- **Inventory items** -- how was this raw data ingested and cleaned?
- **Datasets** -- how was this assembly created?
- **Project phases** -- what tool and parameters produced this checkpoint?

Each provenance record uses a `target_type` discriminator (`inventory`, `phase`, or `dataset`) and links to exactly one target via the corresponding foreign key.


## 4. How Files and DB Interact

The registry is a ledger of pointers, not a data store.

```
Registry DB                              Filesystem
-----------                              ----------
inventory_items.uri  ───────────────►  /scratch/user/inventory/GSE176171_uc_scrnaseq.h5ad
datasets.uri         ───────────────►  /scratch/user/datasets/uc_multimodal.h5mu  (optional)
project_phases.uri   ───────────────►  /scratch/user/projects/uc_atlas/checkpoints/qc_filter.h5ad
```

**Key principles:**

- Checkpoint files (h5ad/h5mu) contain the actual biological data -- expression matrices, cell metadata, embeddings, spatial coordinates.
- Per-cell detail (e.g., `adata.obs["qc_pass"]`, `adata.obs["cell_type"]`) lives inside the files, not in the DB. The DB only stores aggregate counts like `n_obs` and `n_vars`.
- The pipeline reads a checkpoint, transforms it, writes the next checkpoint, and registers the result in the DB.
- Provenance records what produced each file, so you can trace any checkpoint back to its inputs and parameters.


## 5. Common Workflows

### Registering a new data source

```python
from sc_tools.registry import Registry

reg = Registry()

# Layer 0: register the raw data location
reg.register_data_source(
    name="GSE176171_ulcerative_colitis",
    uri="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176171",
    platform="chromium_3p",
    organism="human",
    tissue="colon",
    disease="ulcerative_colitis",
    source_type="public_geo",
)
```

### Ingesting into inventory

```python
# Layer 1: after an agent processes raw data into a clean h5ad
reg.register_inventory_item(
    data_source_name="GSE176171_ulcerative_colitis",
    name="GSE176171_uc_scrnaseq",
    uri="/scratch/user/inventory/GSE176171_uc_scrnaseq.h5ad",
    modality="rna",
    platform="chromium_3p",
    n_obs=300000,
    n_vars=20000,
)
```

### Assembling a dataset

```python
# Layer 2: combine inventory items into a named dataset
reg.create_dataset(name="uc_multimodal", format="mudata")
reg.add_dataset_member("uc_multimodal", "10x_brain_xenium_rna", modality_key="rna")
reg.add_dataset_member("uc_multimodal", "GSE176171_uc_scrnaseq", modality_key="rna_ref")
```

### Linking to a project and tracking phases

```python
# Layer 3: link dataset to project, then track processing
reg.link_project_dataset("uc_atlas", "uc_multimodal", role="primary")

# After QC completes:
reg.upsert_phase(
    project_name="uc_atlas",
    dataset_name="uc_multimodal",
    phase_group="data_processing",
    subphase="qc_filter",
    status="ready",
    uri="/scratch/user/projects/uc_atlas/checkpoints/qc_filter.h5ad",
    n_obs=250000,
    n_vars=18000,
)
```

### Recording provenance

```python
reg.record_provenance(
    tool="scanpy",
    tool_version="1.11.5",
    params={"min_genes": 200, "min_cells": 10, "max_pct_mito": 20},
    environment={"container": "sc_tools.sif", "sha256": "abc123..."},
    phase_id=phase_id,  # returned from upsert_phase
)
```

### Versioning a dataset

```python
# Need to add CODEX protein data for a validation experiment
reg.bump_dataset_version("uc_multimodal")  # creates v2, copies existing members
reg.add_dataset_member("uc_multimodal", "hubmap_intestine_protein", modality_key="protein")
# v2 now has: rna + rna_ref + protein
# v1 still exists with: rna + rna_ref (any phases linked to v1 stay on v1)
```

### Querying

```python
# What datasets does this project use?
reg.list_project_datasets("uc_atlas")

# What phases have been completed?
reg.list_phases("uc_atlas")
# --> [{"phase_group": "data_processing", "subphase": "qc_filter", "status": "ready", ...}, ...]

# Get details for a specific phase
reg.get_phase("uc_atlas", "data_processing", "qc_filter")
# --> {"status": "ready", "uri": "/scratch/.../qc_filter.h5ad", "n_obs": 250000, ...}

# What provenance does a phase have?
reg.get_provenance(phase_id=42)
# --> [{"tool": "scanpy", "tool_version": "1.11.5", "params": {...}, ...}]

# List all inventory items for a modality
reg.list_inventory_items(modality="rna")

# List all data sources for a platform
reg.list_data_sources(platform="chromium_3p")
```


## 6. MCP Tools (for Claude Code agents)

Agents interact with the registry through MCP tools, not Python directly. The MCP server (`sc_tools.mcp.registry_server`) wraps the `Registry` class methods.

| MCP Tool | Registry Method | Purpose |
|----------|----------------|---------|
| `registry_status` | `reg.status()` | High-level summary: counts, active projects |
| `add_project` | `reg.add_project()` | Register a new project |
| `register_data_source` | `reg.register_data_source()` | Register a raw data source (Layer 0) |
| `list_data_sources` | `reg.list_data_sources()` | Query data sources with filters |
| `register_inventory_item` | `reg.register_inventory_item()` | Register a clean h5ad (Layer 1) |
| `list_inventory_items` | `reg.list_inventory_items()` | Query inventory items |
| `create_dataset` | `reg.create_dataset()` | Create a named dataset (Layer 2) |
| `add_dataset_member` | `reg.add_dataset_member()` | Add an inventory item to a dataset |
| `list_datasets` | `reg.list_datasets()` | List datasets, optionally filtered by project |
| `link_project_dataset` | `reg.link_project_dataset()` | Link a dataset to a project (Layer 3) |
| `set_phase_status` | `reg.upsert_phase()` | Create or update a processing phase |
| `get_phase_status` | `reg.get_phase()` | Get status and checkpoint for a phase |
| `mark_phase_complete` | `reg.mark_phase_complete()` | Mark a phase as complete |
| `get_available_next_phases` | (uses `pipeline.get_available_next()`) | What pipeline phases are unblocked |
| `record_provenance` | `reg.record_provenance()` | Record tool/version/params for a transformation |
| `get_provenance` | `reg.get_provenance()` | Query provenance for an item/phase/dataset |
| `project_data_summary` | `reg.project_data_summary()` | Counts by category and platform |
| `add_subject` | `reg.add_patient()` | Register a de-identified patient |
| `list_subjects` | `reg.list_subjects()` | Query patients |
| `list_modalities` | (uses `biodata` module) | List known modalities and platforms |

**Backward-compatible tools** (deprecated, will route to new methods):

| MCP Tool | Routes to |
|----------|-----------|
| `register_dataset` | `reg.register_data()` (auto-creates inventory + dataset) |
| `register_biodata` | `reg.register_data()` (auto-creates inventory + dataset) |
| `list_biodata` | `reg.list_biodata()` |
| `link_project_data_source` | `reg.link_project_data_source()` |
| `get_checkpoint_uri` | `reg.get_dataset_uri()` |


## 7. Modality Conventions

The `modality` field on inventory items follows a controlled vocabulary. The rule is: one inventory item per modality per data source.

| Technology | Inventory items produced | Modality values |
|------------|------------------------|----------------|
| **Visium / Visium HD** | 1 | `rna` (spatial coords live in `adata.obsm`, not a separate modality) |
| **CosMx** | 1 | `rna` |
| **Xenium** | 1 | `rna` |
| **MERFISH** | 1 | `rna` |
| **scRNA-seq** (10x Chromium, etc.) | 1 | `rna` |
| **CITE-seq** | 2 | `rna` and `protein` |
| **Multiome** (RNA + ATAC) | 2 | `rna` and `atac` |
| **IMC** | 1 | `imc` (protein + spatial bundled together) |
| **Spatial proteomics** (CODEX, etc.) | 1 | `protein` |

When assembling a multi-modal dataset, the `modality_key` in `dataset_members` can be more descriptive than the raw modality. For example, a deconvolution project might have:

```python
reg.add_dataset_member("uc_multimodal", "10x_brain_xenium_rna", modality_key="rna")
reg.add_dataset_member("uc_multimodal", "GSE176171_uc_scrnaseq", modality_key="rna_ref")
```

Both inventory items have `modality="rna"`, but the dataset distinguishes them with different keys.


## 8. FAQ

**Q: Where does per-cell QC status live?**
A: In `adata.obs["qc_pass"]` inside the checkpoint h5ad file, not in the DB. The registry only stores aggregate counts (`n_obs`, `n_vars`).

**Q: Can the same inventory item be in multiple datasets?**
A: Yes -- that is the point. A public scRNA-seq atlas can serve as a reference in five different project datasets. The inventory item is registered once; each dataset adds it as a member.

**Q: Can a project use multiple datasets?**
A: Yes. A project can link to multiple datasets with different roles. For example, a primary spatial dataset and a validation scRNA-seq dataset.

**Q: What happens when I bump a dataset version?**
A: A new version is created with the same members copied over. You can then add or remove members on the new version. Project phases linked to the old version stay linked to it -- they do not automatically move to the new version.

**Q: What is the difference between `data_processing` and `discovery` phases?**
A: `data_processing` phases follow the pipeline DAG (`ingest_raw` → `qc_filter` → `preprocess` → `scoring` → ...) and have dependency constraints. `discovery` phases are free-form with descriptive names (`clustering_v1`, `deconv_attempt_2`) -- no ordering enforced. Use `data_processing` for the standard pipeline and `discovery` for exploratory analysis.

**Q: Can I delete an inventory item that is part of a dataset?**
A: No. The `dataset_members` table uses `ON DELETE RESTRICT` on the inventory foreign key. You must remove the member from all datasets before deleting the inventory item.

**Q: How do I load a multi-modal dataset?**
A: Query the dataset members, then load each inventory item individually and combine:

```python
import scanpy as sc
from mudata import MuData

members = reg.get_dataset_members("uc_multimodal")
if len(members) == 1:
    adata = sc.read_h5ad(members[0]["uri"])
else:
    mods = {m["modality_key"]: sc.read_h5ad(m["uri"]) for m in members}
    mdata = MuData(mods)
```

**Q: Where is the database file?**
A: By default, `~/.sc_tools/registry.db` (SQLite). Set the `SC_TOOLS_REGISTRY_URL` environment variable to use PostgreSQL instead.
