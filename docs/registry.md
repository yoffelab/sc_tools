# sc_tools Registry and Storage

Reference for the storage layer, registry database schema, BioData hierarchy, and MCP servers.
For pipeline phases and checkpoint conventions, see [[Architecture]].

---

## 1. Storage Layer (`sc_tools/storage.py`)

All data loading and writing goes through the storage layer so code works identically against
local paths, HPC scratch (SFTP), S3, GCS, Azure, or Box.

| URI scheme | Backend | Install |
|------------|---------|---------|
| `/path` or `file://` | Local filesystem | stdlib (always available) |
| `sftp://brb//athena/...` | SSH/HPC | `pip install sshfs` |
| `s3://bucket/key` | AWS S3 | `pip install s3fs` |
| `gs://bucket/key` | GCS | `pip install gcsfs` |
| `az://container/blob` | Azure Blob | `pip install adlfs` |
| `box://folder/path` | Box | `pip install boxfs` |

**Install all remote backends:** `pip install "sc-tools[storage]"`

Key functions:

| Function | Purpose |
|----------|---------|
| `resolve_fs(uri)` | Returns `(AbstractFileSystem, path)` for any URI |
| `open_file(uri, mode)` | Context manager; opens any URI for read/write |
| `with_local_copy(uri)` | Downloads remote file to tmp; yields local Path |
| `smart_read_h5ad(uri)` | Reads .h5ad from local or remote URI |
| `smart_write_checkpoint(adata, uri, fmt)` | Writes h5ad or zarr to any URI |
| `smart_read_csv(uri)` | Reads CSV/TSV from any URI |

**Data placement by storage tier:**

| Data type | Location | Format | URI example |
|-----------|----------|--------|-------------|
| Raw FASTQs, TIFFs, MCDs | HPC scratch | native | `sftp://brb//athena/.../sample1/` |
| Active checkpoints (p0-p4) | HPC scratch | h5ad | `sftp://brb//athena/.../adata.filtered.h5ad` |
| Archived checkpoints | AWS S3 | zarr | `s3://yoffelab-sc/projects/ggo_visium/results/adata.celltyped.zarr` |
| Manuscript figures | Box | PNG/PDF | `box://YoffeLab/sc_tools/ggo_visium/figures/` |
| Registry DB | Local or lab server | SQLite/PG | `~/.sc_tools/registry.db` |

---

## 2. Registry Database (`sc_tools/registry.py`)

SQLite by default (zero-config). PostgreSQL via `SC_TOOLS_REGISTRY_URL` env var.

```bash
# Default (local SQLite)
python -m sc_tools registry status

# Postgres (lab server)
export SC_TOOLS_REGISTRY_URL="postgresql://user:pass@host:5432/sc_tools"
```

**Install:** `pip install "sc-tools[registry]"` (adds SQLAlchemy + Alembic)

Schema migrations are managed via Alembic in `sc_tools/migrations/`. For new installs,
`Registry()` calls `create_all()` automatically.

### Schema tables

| Table | Tracks |
|-------|--------|
| `projects` | project name, platform, domain, imaging_modality, status |
| `datasets` | AnnData checkpoints with URI, phase, status, md5 (legacy; use `bio_data` for new data) |
| `project_phases` | per-phase pipeline status (status, n_obs, n_vars, n_samples, notes) |
| `slurm_jobs` | submitted SLURM jobs with cluster, status, log URI |
| `agent_tasks` | running Claude Code agent tasks with inputs/outputs |
| `data_sources` | raw data sources (HPC dirs, public datasets, GEO accessions) |
| `project_data_sources` | many-to-many join between projects and data sources |
| **`subjects`** | cross-project de-identified patient/subject records |
| **`samples`** | physical specimens linking subjects to data within projects |
| **`subject_project_links`** | many-to-many join between subjects and projects (role: enrolled/reference/control) |
| **`bio_data`** | JTI base for typed biological data objects (replaces `datasets` for new data) |
| **`bio_images`** | image-specific fields (IMC, CODEX, H&E): n_channels, pixel_size, staining_protocol |
| **`rnaseq_data`** | RNA-seq-specific fields: library_type, chemistry, reference_genome |
| **`spatial_seq_data`** | spatial seq-specific fields: spatial_resolution, panel_size, bin_size, fov_count |
| **`epigenomics_data`** | epigenomics-specific fields: assay_type, target_protein, n_peaks |
| **`genome_seq_data`** | genome seq-specific fields: sequencing_type, panel_name, coverage |

### BioData type hierarchy (Joined Table Inheritance)

```
BioData (bio_data) ── base table, shared columns (incl. modality)
├── BioImage (bio_images)             ── IMC, CODEX, H&E, IF, IHC
├── RNASeqData (rnaseq_data)          ── bulk + single-cell RNA-seq
├── SpatialSeqData (spatial_seq_data) ── Visium, Xenium, CosMx, MERFISH
├── EpigenomicsData (epigenomics_data)── ATAC-seq, ChIP-seq, CUT&Tag
└── GenomeSeqData (genome_seq_data)   ── WGS, WES, targeted panels
```

Each `bio_data` row has a `type` discriminator column and a `modality` column.
Full taxonomy: [[biodata-hierarchy]].

### Subject/Sample model

```
Subject (global, de-identified) 1──* Sample (project-scoped)
Sample 1──* BioData (multiple data types per sample)
Subject can appear across multiple projects (via subject_project_links)
```

775 patients from 4 cohorts (DRAKE, ROS, GGO-IMC, PLM).
Field definitions and per-cohort mappings: [[clinical-data-schema]].

### Three-tier classification

`BioDataType` (5 JTI types) → `BioDataModality` (19 human-readable families) → `BioDataPlatform` (95+ slugs).
The modality tier groups related platforms (e.g. "Spatial Proteomics - Mass Spec" groups imc, mibi, maldi_ims).
Full list: [[biodata-hierarchy]].

### Platform registry

`sc_tools.biodata` provides a pure-Python registry of 95+ platforms. Each `PlatformSpec` maps to the correct
BioData subtype, modality, and default field values. Auto-fill at registration time populates modality, child
columns, and format from URI extension. Extensible via `register_platform()` — no schema migration needed.

### Migration path

`register_dataset()` now dual-writes to both `datasets` and `bio_data` (Phase B deprecation).
`datasets.bio_data_id` FK links forward; `bio_data.legacy_dataset_id` links backward.
Both systems coexist; `register_dataset()` emits a deprecation warning.

---

## 3. MCP Servers (`sc_tools/mcp/`)

Two FastMCP servers registered in `.mcp.json` (auto-discovered by Claude Code):

| Server | Module | Tools |
|--------|--------|-------|
| `sc-tools` | `sc_tools.mcp.tools_server` | validate_checkpoint, generate_qc_report, score_gene_signatures, run_deconvolution, generate_sbatch_spaceranger, generate_sbatch_imc, collect_batch_manifests, load_sample |
| `sc-registry` | `sc_tools.mcp.registry_server` | registry_status, list_datasets, get_checkpoint_uri, register_dataset, list_slurm_jobs, list_agent_tasks, mark_phase_complete, add_subject, list_subjects, add_sample, list_samples, register_biodata, list_biodata, list_modalities, project_data_summary |

**Install:** `pip install "sc-tools[mcp,registry]"`

**CLI status check:**

```bash
python -m sc_tools registry status
```
