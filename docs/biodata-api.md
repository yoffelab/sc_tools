# BioData API Reference

API for registering, querying, and managing typed biological data objects, subjects, and samples in the sc_tools registry.

**Install:** `pip install "sc-tools[registry]"`

```python
from sc_tools.registry import Registry
reg = Registry()
```

---

## 1. BioData Operations

### register_biodata

Register a typed BioData object. Returns the `bio_data` DB id.

```python
reg.register_biodata(
    project_name: str,       # must already exist in registry
    category: str,           # "image" | "rnaseq" | "spatial_seq" | "epigenomics" | "genome_seq"
    platform: str,           # slug from KNOWN_PLATFORMS (e.g. "xenium", "imc", "chromium_3p")
    uri: str,                # file URI (local path, s3://, sftp://, etc.)
    subcategory: str = None, # auto-filled from platform registry
    measurement: str = None, # auto-filled from platform registry
    resolution: str = None,  # auto-filled from platform registry
    spatial: bool = None,    # auto-filled from platform registry
    fmt: str = None,         # auto-inferred from URI extension
    status: str = "pending", # "pending" | "complete" | "failed"
    file_role: str = "primary",
    validated: bool = False,
    n_obs: int = None,
    n_vars: int = None,
    phase: str = None,       # pipeline phase slug (e.g. "qc_filter")
    size_mb: float = None,
    md5: str = None,
    subject_id: str = None,  # de-identified subject ID
    sample_db_id: int = None,
    **type_kwargs,           # subtype-specific columns (see below)
) -> int
```

**Subtype-specific kwargs by category:**

| Category | Model | Extra columns |
|----------|-------|---------------|
| `image` | BioImage | `n_channels`, `pixel_size_um`, `staining_protocol`, `image_type` |
| `rnaseq` | RNASeqData | `library_type`, `chemistry`, `reference_genome` |
| `spatial_seq` | SpatialSeqData | `spatial_resolution`, `panel_size`, `bin_size_um`, `fov_count` |
| `epigenomics` | EpigenomicsData | `assay_type`, `target_protein`, `n_peaks` |
| `genome_seq` | GenomeSeqData | `sequencing_type`, `panel_name`, `coverage` |

**Example:**

```python
bd_id = reg.register_biodata(
    "lymph_dlbcl", "image", "imc",
    "sftp://brb//athena/.../roi001_full.tiff",
    phase="ingest_load",
    n_channels=40,
    pixel_size_um=1.0,
)
```

### get_biodata

Retrieve a single BioData record by its DB id. Returns a dict with all columns (including subtype-specific ones), or `None`.

```python
reg.get_biodata(biodata_id: int) -> dict | None
```

### list_biodata

Query BioData objects with optional filters. All filters are optional.

```python
reg.list_biodata(
    project_name: str = None,
    category: str = None,     # filters by JTI type discriminator
    platform: str = None,
    sample_db_id: int = None,
    phase: str = None,
) -> list[dict]
```

### project_data_summary

Return aggregate counts grouped by category, modality, and platform for a project.

```python
reg.project_data_summary(project_name: str) -> dict
# Returns: {"project": str, "total": int,
#           "by_category": {str: int},
#           "by_modality": {str: int},
#           "by_platform": {str: int}}
```

### migrate_datasets_to_biodata

Migrate all legacy `datasets` rows to `bio_data`. Sets bidirectional FK links (`datasets.bio_data_id` and `bio_data.legacy_dataset_id`). Skips already-migrated rows. Returns count migrated.

```python
reg.migrate_datasets_to_biodata() -> int
```

---

## 2. Subject / Sample API

### Subjects

```python
# Register a de-identified subject (idempotent)
reg.add_subject(
    subject_id: str,            # e.g. "PT001"
    organism: str = "human",    # "human" | "mouse" | "rat" | ...
    # Standard clinical columns (all optional):
    sex: str = None,
    age_at_collection: str = None,
    diagnosis: str = None,
    diagnosis_code: str = None,
    disease_stage: str = None,
    treatment_status: str = None,
    tissue_of_origin: str = None,
    cause_of_death: str = None,
    survival_days: int = None,
    vital_status: str = None,
    clinical_metadata_json: str = None,  # JSON overflow for extra fields
) -> int  # subject DB id

# Retrieve by de-identified ID
reg.get_subject(subject_id: str) -> dict | None

# Query with filters (all optional)
reg.list_subjects(
    project_name: str = None,   # filter to subjects linked to this project
    diagnosis: str = None,      # case-insensitive substring match
    tissue: str = None,         # case-insensitive substring match on tissue_of_origin
) -> list[dict]

# Link subject to project (idempotent)
reg.link_subject_to_project(
    subject_id: str,
    project_name: str,
    role: str = "enrolled",     # "enrolled" | "reference" | "control"
) -> int  # link DB id
```

### Samples

```python
# Register a sample
reg.add_sample(
    sample_id: str,             # e.g. "S001_A1"
    subject_id: str = None,     # de-identified ID (must exist if provided)
    project_name: str = None,   # must exist if provided
    tissue: str = None,
    tissue_region: str = None,
    collection_date: str = None,
    fixation_method: str = None,
    sample_type: str = None,
    batch: str = None,
    notes: str = None,
) -> int  # sample DB id

# Retrieve by sample_id, optionally scoped to project
reg.get_sample(sample_id: str, project_name: str = None) -> dict | None

# Query with filters (all optional)
reg.list_samples(
    project_name: str = None,
    subject_id: str = None,
    batch: str = None,
) -> list[dict]
```

**Relationship model:**

```
Subject (global, de-identified) 1--* Sample (project-scoped)
Sample 1--* BioData (multiple data types per sample)
Subject *--* Project (via subject_project_links, with role)
```

---

## 3. MCP Tool Mapping

The `sc-registry` MCP server exposes equivalent tools for use by Claude Code agents.

| Python API | MCP tool name | Notes |
|------------|---------------|-------|
| `reg.register_biodata(...)` | `register_biodata` | Same params; strings only (no None, use empty string) |
| `reg.list_biodata(...)` | `list_biodata` | Filters via string params |
| `reg.project_data_summary(...)` | `project_data_summary` | Returns JSON string |
| `reg.add_subject(...)` | `add_subject` | Clinical fields as individual string params |
| `reg.list_subjects(...)` | `list_subjects` | Filters via string params |
| `reg.add_sample(...)` | `add_sample` | All params as strings |
| `reg.list_samples(...)` | `list_samples` | Filters via string params |
| `reg.get_biodata(id)` | -- | No direct MCP equivalent; use `list_biodata` |
| `reg.get_subject(id)` | -- | No direct MCP equivalent; use `list_subjects` |
| `reg.get_sample(id)` | -- | No direct MCP equivalent; use `list_samples` |
| `reg.migrate_datasets_to_biodata()` | -- | Run via Python directly |
| `reg.link_subject_to_project(...)` | -- | Run via Python directly |
| `biodata.list_modalities()` | `list_modalities` | Platform registry query (no DB) |

---

## 4. Auto-Fill Behavior

When `platform` is a known slug in `KNOWN_PLATFORMS` (95+ platforms pre-seeded), `register_biodata` auto-fills fields at two levels:

### Level 1: PlatformSpec base fields

These are set from the `PlatformSpec` when not explicitly provided:

| Field | Source |
|-------|--------|
| `subcategory` | `PlatformSpec.subcategory` (e.g. `"sequencing_based"`, `"multiplexed"`) |
| `measurement` | `PlatformSpec.measurement` (e.g. `"rna"`, `"protein"`, `"dna"`) |
| `resolution` | `PlatformSpec.resolution` (e.g. `"single_cell"`, `"spot"`, `"bulk"`) |
| `spatial` | `PlatformSpec.spatial` (bool) |
| `modality` | `PlatformSpec.modality` (e.g. `"Spatial Proteomics - Mass Spec"`) |

### Level 2: PlatformSpec.defaults (child-table columns)

The `defaults` dict on each `PlatformSpec` auto-fills subtype-specific columns:

| Platform | Auto-filled defaults |
|----------|---------------------|
| `imc` | `staining_protocol="IMC"`, `image_type="multiplexed"` |
| `visium` | `spatial_resolution="spot"` |
| `visium_hd` | `spatial_resolution="spot"`, `bin_size_um=8.0` |
| `cosmx_1k` | `panel_size=1000` |
| `chromium_3p` | `chemistry="chromium_v3"`, `library_type="single_cell"` |
| `atac_seq` | `assay_type="atac_seq"` |
| `illumina_wgs` | `sequencing_type="wgs"` |

Explicit values always take precedence over auto-filled defaults.

### Level 3: Format from URI

When `fmt` is not provided, it is inferred from the URI file extension:

| Extension | Inferred format |
|-----------|----------------|
| `.h5ad` | `h5ad` |
| `.zarr` | `zarr` |
| `.tiff` | `tiff` |
| `.csv` | `csv` |
| `.tsv` | `tsv` |
| `.fastq` | `fastq` |
| `.bam` | `bam` |
| `.bed` | `bed` |

---

## 5. Deprecation Timeline (register_dataset)

| Phase | Version | Behavior |
|-------|---------|----------|
| **A** | v1.0 (current) | `register_dataset()` emits `DeprecationWarning`. All new code should use `register_biodata()`. |
| **B** | v1.1 | `register_dataset()` dual-writes: creates both a `datasets` row and a `bio_data` row, linking them via `datasets.bio_data_id`. Existing callers continue to work with no code changes. |
| **C** | future | `register_dataset()` removed. All callers must use `register_biodata()`. The `datasets` table remains read-only for backward compatibility. |

---

## 6. Migration Guide

### Upgrade the database schema

```bash
# From the repo root:
cd sc_tools
alembic upgrade head
```

This applies all migrations in `sc_tools/migrations/versions/`, including `0005_biodata_hierarchy.py` which creates the `bio_data`, `bio_images`, `rnaseq_data`, `spatial_seq_data`, `epigenomics_data`, `genome_seq_data`, `subjects`, `samples`, and `subject_project_links` tables.

For new installs, `Registry()` calls `create_all()` automatically -- no manual migration needed.

### Backfill existing data

```python
from sc_tools.registry import Registry

reg = Registry()

# 1. Migrate legacy datasets -> bio_data
count = reg.migrate_datasets_to_biodata()
print(f"Migrated {count} datasets")

# 2. Verify bidirectional links
for ds in reg.list_datasets("my_project"):
    if ds.get("bio_data_id"):
        bd = reg.get_biodata(ds["bio_data_id"])
        assert bd["legacy_dataset_id"] == ds["id"]

# 3. Register subjects and samples (manual step)
subj_id = reg.add_subject("PT001", sex="F", diagnosis="DLBCL")
reg.link_subject_to_project("PT001", "lymph_dlbcl")
sample_id = reg.add_sample("S001_A1", subject_id="PT001",
                           project_name="lymph_dlbcl",
                           tissue="lymph_node", batch="batch1")
```

### Verify migration

```bash
python -m sc_tools registry status
```

The status output includes `bio_data` counts alongside legacy `datasets` counts. Both systems coexist until Phase C of the deprecation timeline.
