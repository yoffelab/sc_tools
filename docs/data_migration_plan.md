# Data Migration Plan: Four-Layer Registry Schema Redesign

**Date:** 2026-03-16
**Status:** Draft (Rev 3 — post-review rounds 1-3)
**Migration:** 0013_four_layer_schema, 0014_drop_legacy_tables

## 1. Motivation

The current registry uses 6 tables with a flat data model. Processing phases, raw data
inventory, and project-specific checkpoints are all stored in one `data_processing_phase`
table, making it difficult to:

- Share raw data across projects without duplication
- Track data provenance from source through processing
- Version multi-modal dataset assemblies (MuData)
- Separate "what data exists" from "what processing has been done"

The new schema introduces a 4-layer hierarchy that cleanly separates concerns.

## 2. Current Schema (6 tables)

| Table | ORM Model | Row count source | Purpose |
|-------|-----------|-----------------|---------|
| `projects` | `Project` | `n_projects` | Project metadata |
| `data_processing_phase` | `Data` | `n_data` | Checkpoints + phase markers (mixed) |
| `data_inventory` | `DataSource` | -- | Raw data sources (HPC, GEO, etc.) |
| `data_project_map` | `ProjectDataSource` | -- | Project-to-data-source links |
| `patients` | `Patient` | `n_patients` | De-identified patients |
| `patient_data_map` | `PatientDataMap` | -- | Patient-to-data links (FK to data_processing_phase) |

### Current FK Relationships

```
data_processing_phase.project_id --> projects.id (CASCADE)
data_project_map.project_id --> projects.id (CASCADE)
data_project_map.data_source_id --> data_inventory.id (CASCADE)
patient_data_map.patient_id --> patients.id (CASCADE)
patient_data_map.data_id --> data_processing_phase.id (CASCADE)
```

### Current ORM (registry.py)

`_build_models()` returns a 6-tuple:
`(Project, Data, Patient, DataSource, ProjectDataSource, PatientDataMap)`

The `Data` model maps to `data_processing_phase` and serves dual duty: it holds both
real file checkpoints (with URIs like `sftp://...`) and phase-tracking markers (with
`phase://` URIs and `file_role='phase_marker'`).

### Current MCP Tools (registry_server.py)

Active tools: `registry_status`, `list_datasets`, `get_checkpoint_uri`,
`register_dataset`, `add_project`, `get_available_next_phases`, `record_provenance`,
`mark_phase_complete`, `get_phase_status`, `set_phase_status`, `register_data_source`,
`list_data_sources`, `link_project_data_source`, `add_subject`, `list_subjects`,
`register_biodata`, `list_biodata`, `list_modalities`, `project_data_summary`.

Deprecated stubs: `list_slurm_jobs`, `list_agent_tasks`, `update_slurm_job_status`,
`add_sample`, `list_samples`.

## 3. New Schema (10 tables)

**Note on DDL syntax:** All DDL below uses PostgreSQL syntax for illustration
(SERIAL, JSONB, partial indexes). Actual table creation uses SQLAlchemy ORM types
(`Column(Integer, primary_key=True, autoincrement=True)`, `Column(JSON)`) which are
dialect-neutral and work on both PostgreSQL and SQLite.

### Layer 0: data_sources (renamed from current data_inventory)

Raw, uncurated data dumps. Agents clean these into inventory items.

```sql
data_sources (
    id SERIAL PRIMARY KEY,
    name VARCHAR UNIQUE NOT NULL,
    description TEXT,
    uri TEXT NOT NULL,
    platform VARCHAR,
    domain VARCHAR,
    imaging_modality VARCHAR,
    source_type VARCHAR,         -- "geo", "hpc_directory", "arrayexpress", "manual"
    organism VARCHAR,
    tissue VARCHAR,
    disease VARCHAR,
    n_samples INT,
    n_cells INT,
    publication VARCHAR,
    access_notes TEXT,
    status VARCHAR DEFAULT 'discovered',
    metadata JSONB DEFAULT '{}',
    created_at VARCHAR
);
```

### Layer 1: inventory_items (NEW -- clean ingested raw data)

Each row is one AnnData (h5ad) produced by processing a data_source.
Named `inventory_items` (not `data_inventory`) to avoid naming collision with the
renamed Layer 0 table. ORM model: `InventoryItem`.

```sql
inventory_items (
    id SERIAL PRIMARY KEY,
    data_source_id INT REFERENCES data_sources(id) ON DELETE SET NULL,
    name VARCHAR UNIQUE NOT NULL,
    uri TEXT NOT NULL,
    modality VARCHAR NOT NULL,   -- "rna", "protein", "spatial", "morphology", "imc"
    platform VARCHAR,
    format VARCHAR DEFAULT 'h5ad',
    n_obs INT,
    n_vars INT,
    size_mb FLOAT,
    organism VARCHAR,
    tissue VARCHAR,
    metadata JSONB DEFAULT '{}',
    created_at VARCHAR,
    updated_at VARCHAR
);
```

### Layer 2: datasets (versioned assemblies)

Named, versioned collections of inventory items. Can be MuData (multi-modal) or AnnData.

```sql
datasets (
    id SERIAL PRIMARY KEY,
    name VARCHAR NOT NULL,
    version INT NOT NULL DEFAULT 1,
    description TEXT,
    uri TEXT,                    -- path to .h5mu or .h5ad (may be NULL)
    format VARCHAR DEFAULT 'mudata',
    n_obs INT,
    size_mb FLOAT,
    is_current BOOLEAN DEFAULT TRUE,
    metadata JSONB DEFAULT '{}',
    created_at VARCHAR,
    updated_at VARCHAR,
    UNIQUE(name, version)
);
```

### dataset_members (composition of a dataset)

```sql
dataset_members (
    id SERIAL PRIMARY KEY,
    dataset_id INT REFERENCES datasets(id) ON DELETE CASCADE,
    inventory_id INT REFERENCES inventory_items(id) ON DELETE RESTRICT,
    modality_key VARCHAR NOT NULL,
    UNIQUE(dataset_id, modality_key)
);
```

Note: `ON DELETE RESTRICT` prevents accidental deletion of inventory items that are
part of a dataset. You must remove the member first.

### Layer 3: projects (UNCHANGED)

```sql
projects (
    id SERIAL PRIMARY KEY,
    name VARCHAR UNIQUE NOT NULL,
    platform VARCHAR,
    domain VARCHAR,
    status VARCHAR DEFAULT 'active',
    created_at VARCHAR
);
```

### project_datasets (replaces data_project_map)

```sql
project_datasets (
    id SERIAL PRIMARY KEY,
    project_id INT REFERENCES projects(id) ON DELETE CASCADE,
    dataset_id INT REFERENCES datasets(id) ON DELETE CASCADE,
    role VARCHAR DEFAULT 'primary',
    notes TEXT,
    UNIQUE(project_id, dataset_id)
);
```

### project_phases (replaces data_processing_phase for project-specific processing)

```sql
project_phases (
    id SERIAL PRIMARY KEY,
    project_id INT REFERENCES projects(id) ON DELETE CASCADE,
    dataset_id INT REFERENCES datasets(id) ON DELETE RESTRICT NOT NULL,
    phase_group VARCHAR NOT NULL,  -- 'data_processing' | 'discovery'
    subphase VARCHAR NOT NULL,
    status VARCHAR DEFAULT 'pending',
    uri TEXT,
    n_obs INT,
    n_vars INT,
    metadata JSONB DEFAULT '{}',
    created_at VARCHAR,
    updated_at VARCHAR,
    UNIQUE(project_id, dataset_id, phase_group, subphase),
    CHECK (phase_group IN ('data_processing', 'discovery'))
);
```

Phase groups:
- **data_processing** subphases: `ingest_raw`, `ingest_load`, `qc_filter`,
  `metadata_attach`, `preprocess`, `demographics`, `scoring`, `celltype_manual`,
  `biology`, `meta_analysis`. These are the full set of existing pipeline phases,
  kept as-is. No CHECK constraint on subphase values -- `pipeline.py` DAG is the
  source of truth for valid data_processing subphases.
- **discovery** subphases: free-form strings, no constraints. Users can use
  whatever naming they want (`clustering_v1`, `deconv_attempt_3`, `final_final`).
  Discovery is unstructured by design.

**Design decision:** `dataset_id` is NOT NULL with `ON DELETE RESTRICT`. A dataset
must be assembled and linked to the project before any phase tracking begins.
The flow is: data_sources -> inventory_items -> dataset -> link to project -> start phases.
Processing history should not silently vanish when a dataset is deleted.

**Design decision:** The CHECK constraint enforces only valid `phase_group` values
(`data_processing` or `discovery`). Subphase validation for `data_processing` is
handled by the `pipeline.py` DAG, not a DB constraint. Discovery subphases are
entirely free-form.

### patients (UNCHANGED)

```sql
patients (
    id SERIAL PRIMARY KEY,
    patient_id VARCHAR UNIQUE NOT NULL,
    metadata JSONB DEFAULT '{}',
    created_at VARCHAR
);
```

### patient_data_map (MODIFIED -- FK retargeted)

```sql
patient_data_map (
    id SERIAL PRIMARY KEY,
    patient_id INT REFERENCES patients(id) ON DELETE CASCADE,
    inventory_id INT REFERENCES inventory_items(id) ON DELETE CASCADE,
    role VARCHAR DEFAULT 'source',
    notes TEXT,
    UNIQUE(patient_id, inventory_id)
);
```

### provenance (NEW)

Tracks tool versions, params, and environment for all transformation points.

```sql
provenance (
    id SERIAL PRIMARY KEY,
    target_type VARCHAR NOT NULL,  -- 'inventory' | 'phase' | 'dataset'
    inventory_id INT REFERENCES inventory_items(id) ON DELETE SET NULL,
    phase_id INT REFERENCES project_phases(id) ON DELETE SET NULL,
    dataset_id INT REFERENCES datasets(id) ON DELETE SET NULL,
    tool VARCHAR NOT NULL,
    tool_version VARCHAR,
    reference_genome VARCHAR,
    reference_dataset TEXT,
    signature_source TEXT,
    n_input_obs INT,
    n_output_obs INT,
    params JSONB DEFAULT '{}',
    environment JSONB DEFAULT '{}',
    script_uri TEXT,
    agent VARCHAR,
    created_at VARCHAR,
    CHECK (
        (target_type = 'inventory' AND inventory_id IS NOT NULL AND phase_id IS NULL AND dataset_id IS NULL) OR
        (target_type = 'phase' AND inventory_id IS NULL AND phase_id IS NOT NULL AND dataset_id IS NULL) OR
        (target_type = 'dataset' AND inventory_id IS NULL AND phase_id IS NULL AND dataset_id IS NOT NULL)
    )
);
```

The `target_type` discriminator column enables efficient querying without checking
all three nullable FKs. Index on `(target_type, inventory_id)` etc. for fast lookups.

## 4. Relationship Diagram

```
data_sources ──(agent ingest)──→ inventory_items ──→ dataset_members ──→ datasets
                                       │                                     │
                                 patient_data_map                     project_datasets
                                       │                                     │
                                   patients                              projects
                                                                            │
                                                                      project_phases
                                                                       (phase_group +
                                                                        subphase with
                                                                        CHECK constraint)
                                                                            │
                                                                       provenance
                                                              (target_type discriminator,
                                                               also → inventory, dataset)
```

**Full FK map (10 tables):**
```
inventory_items.data_source_id     → data_sources.id       (SET NULL)
dataset_members.dataset_id         → datasets.id           (CASCADE)
dataset_members.inventory_id       → inventory_items.id     (RESTRICT)
project_datasets.project_id        → projects.id           (CASCADE)
project_datasets.dataset_id        → datasets.id           (CASCADE)
project_phases.project_id          → projects.id           (CASCADE)
project_phases.dataset_id          → datasets.id           (RESTRICT)
patient_data_map.patient_id        → patients.id           (CASCADE)
patient_data_map.inventory_id      → inventory_items.id     (CASCADE)
provenance.inventory_id            → inventory_items.id     (SET NULL)
provenance.phase_id                → project_phases.id     (SET NULL)
provenance.dataset_id              → datasets.id           (SET NULL)
```

### How the DB and Files Interact

The DB is a **ledger**, not a data store. It tracks what exists and whether it is done.

- Each `project_phases` row points to a concrete checkpoint file (h5ad/h5mu) via
  the `uri` column. The full data lives in the files; the DB records metadata about them.
- Provenance records what tool, version, and parameters produced each file.
- The pipeline DAG (`pipeline.py`) reads the current checkpoint, applies transforms,
  writes the next checkpoint, and records the result in the DB.
- **Per-sample tracking is NOT in the DB.** The DB tracks phase completion at the
  dataset level. Per-sample detail (e.g., `adata.obs['qc_pass']`, cell-level QC flags)
  lives in the h5ad checkpoint files. The registry is a ledger pointing to files,
  not a data store.

## 5. Implementation Tasks

### Task 1a: Migration 0013 -- Add New Tables (additive only)

**File:** `sc_tools/migrations/versions/0013_four_layer_schema.py`

This migration is **additive only** — no drops, no renames. Old and new tables coexist
for a burn-in period. This de-risks the big-bang approach.

Steps (in order):

1. Create `data_sources` table (copy schema from current `data_inventory` + add `metadata` JSONB)
2. Create `inventory_items` table (clean ingested data, FK to data_sources)
3. Create `datasets` table
4. Create `dataset_members` table
5. Create `project_datasets` table
6. Create `project_phases` table (with CHECK constraint)
7. Create `provenance` table (with discriminator + CHECK)
8. Create all indexes listed in the Indexes section
9. Copy existing `data_inventory` rows → `data_sources`
10. Migrate `data_processing_phase` rows:
    - Rows with real file URIs (not `phase://` markers) → `inventory_items`
    - Phase marker rows (`file_role='phase_marker'`) → `project_phases`
      (map slugs per table above, store original slug in metadata JSONB)
    - Re-map `patient_data_map.data_id` → add new `inventory_id` column (nullable),
      populate from mapping, keep old `data_id` column intact for now
11. Copy `data_project_map` links → `project_datasets` (auto-create placeholder
    datasets where needed)

### Task 1b: Migration 0014 -- Drop Legacy Tables (after burn-in)

**File:** `sc_tools/migrations/versions/0014_drop_legacy_tables.py`

Run after confirming 0013 data is correct and all code paths use new tables.

Steps:
1. Drop `patient_data_map.data_id` column (old FK)
2. Drop `data_processing_phase` table
3. Drop `data_project_map` table
4. Rename old `data_inventory` table to `_data_inventory_legacy` (keep for safety)
   or drop if confirmed safe

**Phase slug mapping** (current --> new):

All existing phase slugs are kept as-is. The mapping is 1:1 into `data_processing`:

| Current phase slug | New phase_group | New subphase |
|-------------------|----------------|-------------|
| `ingest_raw` | data_processing | ingest_raw |
| `ingest_load` | data_processing | ingest_load |
| `qc_filter` | data_processing | qc_filter |
| `metadata_attach` | data_processing | metadata_attach |
| `preprocess` | data_processing | preprocess |
| `demographics` | data_processing | demographics |
| `scoring` | data_processing | scoring |
| `celltype_manual` | data_processing | celltype_manual |
| `biology` | data_processing | biology |
| `meta_analysis` | data_processing | meta_analysis |

All phases are structured pipeline steps and belong in `data_processing`.
The `discovery` phase_group is for NEW unstructured iteration work only
(e.g., `clustering_v1`, `deconv_attempt_3`), not for migrating existing phases.

### Indexes

Beyond PK and UNIQUE constraints, add explicit indexes for common query patterns:

```sql
CREATE INDEX ix_inventory_items_data_source_id ON inventory_items(data_source_id);
CREATE INDEX ix_dataset_members_dataset_id ON dataset_members(dataset_id);
CREATE INDEX ix_project_datasets_project_id ON project_datasets(project_id);
CREATE INDEX ix_project_phases_project_id ON project_phases(project_id);
CREATE INDEX ix_project_phases_project_dataset ON project_phases(project_id, dataset_id);
CREATE INDEX ix_patient_data_map_inventory_id ON patient_data_map(inventory_id);
CREATE INDEX ix_provenance_inventory ON provenance(inventory_id) WHERE inventory_id IS NOT NULL;
CREATE INDEX ix_provenance_phase ON provenance(phase_id) WHERE phase_id IS NOT NULL;
CREATE INDEX ix_provenance_dataset ON provenance(dataset_id) WHERE dataset_id IS NOT NULL;
CREATE UNIQUE INDEX ix_datasets_current ON datasets(name) WHERE is_current = TRUE;
```

Note: Partial indexes (`WHERE ... IS NOT NULL` and `WHERE is_current = TRUE`) work
in PostgreSQL. For SQLite, partial WHERE clauses are not enforced -- the
`ix_datasets_current` uniqueness gap on SQLite must be handled in application code
(`bump_dataset_version()` must unset `is_current` before inserting the new version).
Regular (non-partial) indexes still provide query performance benefits on SQLite.

The `ix_provenance_target` index (on `target_type` alone) was removed -- too low
selectivity to be useful. The three partial indexes on the FK columns provide
efficient lookups by target type.

### Task 2: ORM Models Update

**File:** `sc_tools/registry.py` -- `_build_models()` function

Change return from 6-tuple to `SimpleNamespace` with named attributes:

```python
from types import SimpleNamespace

def _build_models(Base):
    # ... define all 10 model classes ...
    return SimpleNamespace(
        Project=Project,
        DataSource=DataSource,
        InventoryItem=InventoryItem,
        Dataset=Dataset,
        DatasetMember=DatasetMember,
        ProjectDataset=ProjectDataset,
        ProjectPhase=ProjectPhase,
        Patient=Patient,
        PatientDataMap=PatientDataMap,
        Provenance=Provenance,
    )
```

Callers use `models.Project` instead of positional unpacking.

| New ORM Model | Table | Notes |
|--------------|-------|-------|
| `Project` | `projects` | Unchanged |
| `DataSource` | `data_sources` | Renamed table, add `metadata` column |
| `InventoryItem` | `inventory_items` | NEW: clean ingested data |
| `Dataset` | `datasets` | NEW: versioned assemblies |
| `DatasetMember` | `dataset_members` | NEW: composition |
| `ProjectDataset` | `project_datasets` | Replaces `ProjectDataSource` |
| `ProjectPhase` | `project_phases` | Replaces `Data` for phase tracking |
| `Patient` | `patients` | Unchanged |
| `PatientDataMap` | `patient_data_map` | FK retargeted: `inventory_id` |
| `Provenance` | `provenance` | NEW: tool/version/params |

Key relationship changes:
- `Project.data` relationship removed (was to `Data`/`data_processing_phase`)
- `Project.phase_links` added (to `ProjectPhase`)
- `Project.dataset_links` added (to `ProjectDataset`)
- `PatientDataMap.data` renamed to `PatientDataMap.inventory_item`
- `PatientDataMap.data_id` renamed to `PatientDataMap.inventory_id`
- `DataSource.project_links` removed (was via `data_project_map`)

Inverse relationships to add:
- `InventoryItem.data_source` -> DataSource
- `InventoryItem.dataset_memberships` -> DatasetMember
- `Dataset.members` -> DatasetMember
- `Dataset.project_links` -> ProjectDataset
- `Dataset.phases` -> ProjectPhase (via project_phases.dataset_id)
- `Provenance.inventory_item` -> InventoryItem
- `Provenance.phase` -> ProjectPhase
- `Provenance.dataset` -> Dataset

### Task 3: Registry Methods Rewrite

**File:** `sc_tools/registry.py` -- `Registry` class

#### New methods

| Method | Purpose |
|--------|---------|
| `register_inventory_item(data_source_name, name, uri, modality, ...)` | Create data_inventory row |
| `get_inventory_item(name)` | Return inventory item dict |
| `list_inventory_items(modality=None, platform=None)` | List inventory items |
| `create_dataset(name, description=None, format='mudata')` | Create datasets row |
| `get_dataset(name, version=None)` | Return dataset dict (latest if no version) |
| `add_dataset_member(dataset_name, inventory_name, modality_key)` | Add member |
| `remove_dataset_member(dataset_name, modality_key)` | Remove member |
| `get_dataset_members(dataset_name, version=None)` | List members |
| `bump_dataset_version(dataset_name)` | Create new version, copy members |
| `link_project_dataset(project_name, dataset_name, role='primary')` | Link |
| `list_project_datasets(project_name)` | List datasets for project |
| `record_provenance(tool, tool_version=None, ..., inventory_id=None, phase_id=None, dataset_id=None)` | Record provenance |
| `get_provenance(inventory_id=None, phase_id=None, dataset_id=None)` | Query provenance |

#### Rewritten methods (change FK paths)

| Method | Change |
|--------|--------|
| `upsert_phase(project_name, dataset_name, phase_group, subphase, ...)` | Queries `project_phases` instead of `Data`; new signature requires `dataset_name`, `phase_group`, `subphase` |
| `get_phase(project_name, phase_group, subphase)` | Queries `project_phases` |
| `list_phases(project_name, phase_group=None)` | Queries `project_phases` |
| `mark_phase_complete(project_name, phase_group, subphase)` | Updates `project_phases` |
| `list_datasets(project_name=None)` | Queries new `datasets` table via `project_datasets` |
| `get_dataset_uri(project_name, dataset_name)` | Via `project_datasets` --> `datasets` |
| `project_data_summary(project_name)` | Via `project_datasets` --> `datasets` --> `dataset_members` --> `data_inventory` |
| `status()` | Updated counts: data_sources, inventory_items, datasets, phases |

#### Renamed/updated methods

| Method | Change |
|--------|--------|
| `register_data_source()` | Same signature; table now `data_sources` |
| `list_data_sources()` | Same signature; table now `data_sources` |

#### Backward compat wrappers (deprecated, emit DeprecationWarning)

| Wrapper | Routes to |
|---------|-----------|
| `register_data(project_name, phase, uri, ...)` | `register_inventory_item()` |
| `register_dataset(project_name, phase, uri, ...)` | `register_inventory_item()` |
| `register_biodata(project_name, category, platform, uri, ...)` | `register_inventory_item()` |
| `list_biodata(project_name=None, ...)` | `list_inventory_items()` |
| `link_project_data_source(project_name, data_source_name, ...)` | `link_project_dataset()` |
| `list_project_data_sources(project_name)` | `list_project_datasets()` |

Thin wrappers: emit `DeprecationWarning`, delegate to the closest new method with
minimal auto-inference. Do NOT auto-create datasets or inventory items. Callers must
migrate to the new API.

#### Error handling

`add_dataset_member()`, `bump_dataset_version()`, and `link_project_dataset()` must
do pre-checks and raise `ValueError` with descriptive messages (e.g., "Dataset
'foo' not found", "Inventory item 'bar' is already a member of dataset 'foo'").
Raw `IntegrityError` from the DB must not propagate to callers.

#### Transaction safety

Multi-step operations (any method that touches >1 table) must use a single session
scope. If any step fails, the entire operation rolls back -- no partial data.

### Task 4: MCP Server Update

**File:** `sc_tools/mcp/registry_server.py`

#### New MCP tools

| Tool | Calls |
|------|-------|
| `register_inventory_item(name, uri, modality, ...)` | `reg.register_inventory_item()` |
| `create_dataset(name, description, format)` | `reg.create_dataset()` |
| `add_dataset_member(dataset_name, inventory_name, modality_key)` | `reg.add_dataset_member()` |
| `link_project_dataset(project_name, dataset_name, role)` | `reg.link_project_dataset()` |
| `record_provenance(tool, tool_version, ...)` | `reg.record_provenance()` |
| `get_provenance(inventory_name=None, dataset_name=None, ...)` | `reg.get_provenance()` |
| `list_inventory_items(modality, platform)` | `reg.list_inventory_items()` |

#### Updated MCP tools

| Tool | Change |
|------|--------|
| `set_phase_status` | New params: `phase_group`, `subphase` (replaces flat `phase`) |
| `get_phase_status` | New params: `phase_group`, `subphase` |
| `mark_phase_complete` | New params: `phase_group`, `subphase` |
| `list_datasets` | Queries new `datasets` table |
| `project_data_summary` | Updated query path |
| `registry_status` | Updated counts (data_sources, inventory, datasets, phases) |
| `register_dataset` | Backward compat wrapper |
| `register_biodata` | Backward compat wrapper |

#### Deprecated MCP tools (keep as stubs)

| Tool | Message |
|------|---------|
| `link_project_data_source` | "Use link_project_dataset instead" |
| `get_checkpoint_uri` | Updated to query datasets table |

#### Remove deprecated stubs

These can finally be removed (they were stubs since migration 0008):
- `list_slurm_jobs`
- `list_agent_tasks`
- `update_slurm_job_status`
- `add_sample`
- `list_samples`

### Task 4.5: Pipeline DAG Update

**File:** `sc_tools/pipeline.py`

The current DAG uses flat phase slugs. Update to use `(phase_group, subphase)` tuples:

- DAG nodes become `("data_processing", "ingest_raw")`, `("data_processing", "qc_filter")`, etc.
- `get_available_next_phases()` returns these tuples
- `get_available_next_phases` MCP tool formats them as `"data_processing/ingest_raw"` for display
- Discovery phases are not in the DAG (they are unstructured iterations)
- Dependency chain within data_processing:
  `ingest_raw -> ingest_load -> qc_filter -> metadata_attach -> preprocess ->
   demographics -> scoring -> celltype_manual -> biology -> meta_analysis`
- `pipeline.py` is the source of truth for valid `data_processing` subphases and
  their dependency ordering. No DB CHECK constraint duplicates this.

### Task 5: Test Rewrite

**File:** `sc_tools/tests/test_registry.py`

Complete rewrite. Test classes:

| Class | Tests |
|-------|-------|
| `TestProjects` | add, get, list, delete (unchanged) |
| `TestDataSources` | register, get, list, filter by platform/disease/tissue |
| `TestInventoryItems` | register linked to data_source, register standalone, get, list, filter |
| `TestDatasets` | create, get, list, add_member, remove_member, get_members |
| `TestDatasetVersioning` | bump_version copies members, is_current flag, get latest |
| `TestProjectDatasets` | link, list, role filtering |
| `TestProjectPhases` | upsert data_processing, upsert discovery, get, list, mark_complete |
| `TestPatients` | add, get, list (unchanged) |
| `TestPatientDataMap` | link patient to inventory item |
| `TestProvenance` | record for inventory, record for phase, record for dataset, get, CHECK constraint |
| `TestBackwardCompat` | register_data, register_dataset, register_biodata, list_biodata wrappers |
| `TestStatus` | status() returns updated counts |
| `TestEdgeCases` | cross-project data sharing, dataset mutability, empty datasets |
| `TestTransactionRollback` | failure mid-way through multi-step ops leaves no partial data |
| `TestMigration0013` | loads 6-table fixture, runs upgrade, verifies row counts and FK integrity |

### Task 6: Update docs/Plan.md

Add reference to this migration plan. Update the "Technical Decisions" section.

## 6. MuData Scalability Notes

Datasets are **logical assemblies** by default. The `dataset_members` table IS the
dataset structure. Physical `.h5mu` files are materialized only when needed.

- `format` column distinguishes `mudata` vs `h5ad` (single modality)
- Individual modalities can be loaded independently via inventory URIs
- Recommended utility pattern:

```python
def load_dataset(dataset_id) -> MuData | AnnData:
    members = registry.get_dataset_members(dataset_id)
    if len(members) == 1:
        return sc.read_h5ad(members[0].uri)
    mods = {m.modality_key: sc.read_h5ad(m.inventory_uri) for m in members}
    return MuData(mods)
```

## 7. Execution Order

```
Task 1a (migration 0013 — additive)
    |
    v
Task 2 (ORM models)
    |
    v
Task 3 (registry methods)
    |
    +---> Task 4 (MCP server)      [parallel]
    +---> Task 4.5 (pipeline DAG)  [parallel]
    +---> Task 5 (tests)           [parallel]
    |
    v
Task 6 (docs)
    |
    v
[burn-in period — verify data integrity]
    |
    v
Task 1b (migration 0014 — drop legacy)
```

Estimated effort: ~4 agent dispatches (migration+ORM, methods, MCP+pipeline+tests, cleanup).

## 8. Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Data migration splits rows across two tables | Data loss if FK re-mapping is wrong | Two-phase migration (0013 additive, 0014 drops); test on prod DB copy first |
| Backward compat breakage | Agents and MCP clients fail | Keep all wrapper methods with DeprecationWarning; run full test suite |
| Pipeline DAG breakage | `get_available_next_phases` returns wrong data | Task 4.5 explicitly updates pipeline.py to use (phase_group, subphase) tuples |
| Dual-DB migration | Must work on both Supabase (PostgreSQL) and dev SQLite | Test migration on both engines; CHECK constraints written in expanded OR form |
| Phase slug mapping ambiguity | Old slugs could lose meaning if compressed | All 10 existing slugs kept as-is with 1:1 mapping into data_processing; no information loss |
| Dataset version bump mid-processing | project_phases reference old dataset_id | project_phases stays linked to the dataset version it was created with; new version = new phases |
| Concurrent version bump | Two rows with is_current=TRUE for the same dataset name | Partial unique index `ix_datasets_current` prevents this on PostgreSQL; application-level check on SQLite |
| Transaction rollback on multi-step ops | Partial data left in DB | All multi-step methods use a single session scope; failure rolls back everything |
| inventory_items deletion while in dataset | Silent dataset corruption | `ON DELETE RESTRICT` on dataset_members.inventory_id prevents this |

### Migration Test Strategy

Before running 0014 (drop legacy), verify 0013 data integrity:

1. Dump current registry.db row counts per table
2. Run migration 0013 on a copy
3. Verify: `count(inventory_items)` = `count(data_processing_phase WHERE uri NOT LIKE 'phase://%')`
4. Verify: `count(project_phases)` = `count(data_processing_phase WHERE file_role = 'phase_marker')`
5. Verify: all `patient_data_map.inventory_id` values exist in `inventory_items`
6. Verify: `PRAGMA foreign_key_check` (SQLite) or equivalent FK validation (Postgres)
7. Add a test class `TestMigration0013` with a fixture that loads a representative
   6-table snapshot and verifies upgrade + downgrade preserves all data

Additionally, remove dead test classes: `TestSlurmJobs`, `TestAgentTasks` (stubs
since migration 0008, no longer exercised).

## 9. Rollback Plan

If the migration fails or introduces regressions:

1. The Alembic downgrade path in 0013 must restore:
   - `data_processing_phase` table with original data
   - `data_project_map` table
   - `patient_data_map.data_id` column (rename from `inventory_id`)
   - Rename `data_sources` back to `data_inventory`
2. Drop new tables: `provenance`, `project_phases`, `project_datasets`,
   `dataset_members`, `datasets`, new `data_inventory`
3. Git revert of ORM/Registry/MCP changes

## 10. Files Modified

| File | Action |
|------|--------|
| `sc_tools/migrations/versions/0013_four_layer_schema.py` | CREATE |
| `sc_tools/registry.py` | REWRITE (`_build_models`, `Registry` class) |
| `sc_tools/mcp/registry_server.py` | UPDATE (new tools, updated signatures) |
| `sc_tools/tests/test_registry.py` | REWRITE |
| `docs/Plan.md` | UPDATE (reference to this plan) |
