# External Integrations

**Analysis Date:** 2026-03-20

## APIs & External Services

**Public Data Downloads:**
- IMC benchmark datasets via `urllib.request` (in `sc_tools/data/imc/benchmark/public.py`)
  - Downloads zip archives to local cache
  - Used for public benchmark reference datasets

**Sequencing Platforms (Data Ingestion Only):**
- 10x Genomics Space Ranger (Visium) - reads output directory structure, no API
- 10x Genomics Space Ranger HD (Visium HD) - reads output directory structure, no API
- 10x Genomics Xenium - reads output directory structure, no API
- Nanostring CosMx - reads output directory structure, no API
- Imaging Mass Cytometry (IMC) - reads output directory structure, no API

**No REST APIs or third-party cloud service integrations detected** — all data ingestion is file-based from platform outputs.

## Data Storage

**Databases:**

**SQLite (Default):**
- Location: `~/.sc_tools/registry.db`
- Purpose: Project metadata, pipeline phase tracking, data provenance, inventory management
- Client: SQLAlchemy ORM (`sqlalchemy >= 2.0`)
- Schema version: 18 migrations (as of 2026-03-20)
- Migrations managed by: Alembic (`alembic >= 1.13`)
- Migration directory: `sc_tools/migrations/versions/`
- Config: `alembic.ini` at repo root

**PostgreSQL (Optional):**
- Connection: via `SC_TOOLS_REGISTRY_URL` environment variable
- Format: `postgresql://user:password@host:port/database`
- Driver: `psycopg2-binary >= 2.9`
- Full compatibility with SQLite schema (same ORM models)
- Use cases: Multi-user deployments, shared registry across HPC clusters

**Primary Data Format:**
- AnnData (`.h5ad` files) - HDF5-based single-cell matrix storage
- Stored in: HPC scratch filesystems or remote storage backends
- No SQL tables for expression matrices (SQLite registry tracks metadata only)

**File Storage Backends:**

**Local Filesystem:**
- Primary: HPC cluster scratch directories (`/athena` on brb, similar on cayuga)
- Supported: Native Python pathlib access

**Remote Storage (fsspec-abstracted):**
- **AWS S3** - `s3://bucket/path` URIs via `s3fs >= 2024.1`
  - Authentication: `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY` environment variables
  - Use case: Cloud-based data archives, checkpoint sharing

- **Google Cloud Storage (GCS)** - `gs://bucket/path` URIs via `gcsfs >= 2024.1`
  - Authentication: `GOOGLE_APPLICATION_CREDENTIALS` environment variable (JSON keyfile path)

- **Azure Blob Storage** - `az://container/path` URIs via `adlfs >= 2024.1`
  - Authentication: `AZURE_STORAGE_ACCOUNT_NAME`, `AZURE_STORAGE_ACCOUNT_KEY` environment variables

- **SFTP/SSH** - `sftp://host/path` URIs via `sshfs >= 0.13`
  - Use case: Remote HPC cluster access (e.g., submitting jobs on brb from cayuga)
  - SSH keys: Uses system SSH config (~/.ssh/)

- **Box.com** - `box://path` URIs via `boxfs >= 0.2`
  - Authentication: OAuth flow at runtime (interactive)
  - Use case: Enterprise file sharing (separate `[storage-box]` extra)

- **Zarr Format** - `zarr >= 2.18`, `ome-zarr >= 0.9`
  - Supported on all backends (local, S3, GCS, Azure)
  - Use case: Large-scale distributed array storage, HPC-friendly chunking

**URI Resolution:**
- Unified abstraction in `sc_tools/storage.py`
- Function `resolve_fs(uri)` returns `(AbstractFileSystem, path)` for any URI
- Automatic backend selection based on scheme
- Missing backend raises ImportError with installation hint

## Authentication & Identity

**Auth Provider:**
- Custom (none) - sc_tools is single-user scientific software
- No user identity/login system
- Assumes: Single user running analysis on HPC cluster or local machine
- Permission model: Unix file permissions on HPC scratch

**Database Access Control:**
- SQLite: File-system permissions only
- PostgreSQL: Database user credentials via `SC_TOOLS_REGISTRY_URL` connection string

**Remote Storage Auth:**
- **S3:** IAM credentials (env vars or AWS credential file at `~/.aws/credentials`)
- **GCS:** Service account JSON (path via `GOOGLE_APPLICATION_CREDENTIALS`)
- **Azure:** Storage account credentials (env vars)
- **SFTP:** SSH public key authentication (system SSH config)
- **Box.com:** OAuth 2.0 (interactive flow at runtime)

## Monitoring & Observability

**Error Tracking:**
- None - errors logged to stdout/stderr via Python logging module

**Logging:**
- Framework: Python standard `logging` module
- Configured per module (e.g., `logger = logging.getLogger(__name__)`)
- Log levels: DEBUG, INFO, WARNING, ERROR
- Output: Console/stderr (no external log aggregation)

**Progress Tracking:**
- `tqdm >= 4.65` - Progress bars for long-running operations
- No webhooks, external progress tracking, or status callbacks

**Benchmarking & Metrics:**
- `scib-metrics >= 0.4` - Single-cell integration quality metrics (in-process)
- Plotly HTML reports for visualization of results (saved locally)

**No external observability:**
- No APM (application performance monitoring)
- No metrics aggregation (Prometheus, Datadog, etc.)
- No real-time alerting

## CI/CD & Deployment

**Hosting:**
- GitHub (repository: https://github.com/yoffelab/sc_tools)
- ReadTheDocs (documentation: https://sc-tools.readthedocs.io)
- PyPI (package distribution - assumed, not yet published)

**CI Pipeline:**
- GitHub Actions (`.github/workflows/ci.yml`)
- Stages: lint, test, docs, docker
- No deployment pipeline (manual docker build / push)

**Artifact Storage:**
- GitHub Actions artifacts: HTML docs (14-day retention)
- Docker images: Built locally, pushed manually (not auto-published)

**HPC Job Submission:**
- SLURM scheduler (brb and cayuga clusters)
- Job scripts generated by sc_tools (see `sc_tools/ingest/slurm.py`, `sc_tools/bm/slurm.py`)
- No integration with job monitoring services

## Environment Configuration

**Required env vars for full deployment:**
- `SC_TOOLS_REGISTRY_URL` - PostgreSQL connection string (optional; defaults to SQLite)

**Optional env vars for remote storage:**
- `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY` - S3 credentials
- `GOOGLE_APPLICATION_CREDENTIALS` - GCS service account JSON path
- `AZURE_STORAGE_ACCOUNT_NAME`, `AZURE_STORAGE_ACCOUNT_KEY` - Azure credentials

**Secrets location:**
- Environment variables (set by HPC cluster, local shell, GitHub Actions secrets)
- AWS credentials: `~/.aws/credentials` or environment variables
- SSH keys: `~/.ssh/` (system default)
- No `.env` file in repository

**Build-time configuration:**
- `pyproject.toml` - Dependency versions, installation extras
- `Dockerfile` - RAPIDS base image, system deps, PyTorch/rapids-singlecell installation order

## Webhooks & Callbacks

**Incoming Webhooks:**
- None

**Outgoing Webhooks/Callbacks:**
- None - sc_tools is pull-based (analyzes data when invoked)

**Integration Entry Points:**
- MCP servers (`sc_tools.mcp.registry_server`, `sc_tools.mcp.tools_server`) - invoked by Claude Code
- CLI (`python -m sc_tools registry status`) - manual invocation
- Python API (direct import of `sc_tools` modules) - programmatic use
- SLURM sbatch scripts (generated, not bidirectional callbacks)

## Data Format Standards

**Input Formats:**
- AnnData (`.h5ad`) - Single-cell expression matrices with metadata
- Parquet (`.parquet`) - Spatial coordinates, sample manifests
- CSV/TSV (`.csv`, `.tsv`) - Metadata, sample lists, gene lists
- TIFF (`.tif`, `.tiff`) - Microscopy images (IMC, Xenium)
- JSON (`.json`) - Configuration, platform specs

**Output Formats:**
- AnnData (`.h5ad`) - Checkpoints at each pipeline phase
- Zarr (`.zarr/`) - Distributed arrays, alternative checkpoint format
- OME-Zarr - NGFF standard (spatial image + expression)
- HTML reports (`.html`) - QC reports, analysis summaries (Jinja2 + Plotly)
- CSV/TSV - Results tables, metrics

**Metadata Standards:**
- AnnData obs/var/obsm keys standardized per phase
- BioData platform specs in `sc_tools/biodata.py` (registry-independent)
- Phase slugs (p0a, p0b, p1, p2, p3, p3.5, p3.5b, p4, p5, p6/p7)

## Third-Party Package Integrations

**Biological Data Tools (via Python imports):**
- `scanpy` - Single-cell analysis (wrapped/extended by sc_tools)
- `squidy` - Spatial analysis tools
- `scvi-tools` - VAE models for integration/deconvolution
- `tangram-sc` - Cell-type deconvolution
- `harmonypy` - Batch integration
- `bbknn`, `scanorama` - Alternative batch correction methods
- `leidenalg` + `igraph` - Graph-based clustering
- `celltypist` - Automated cell typing
- `gseapy` - Gene set enrichment
- `decoupler` - Pathway/TF activity inference
- `cellpose`, `stardist`, `deepcell` - Image segmentation

**No direct API calls to these packages** - all are imported as Python libraries, not remote services.

## Schema / Data Model

**Four-Layer Registry Schema:**

1. **Layer 1 - Inventory (InventoryItems):**
   - Raw files from sequencing/imaging platforms
   - Tracks: file URI, format, size, modality, platform
   - Referenced by: DataSources and Datasets

2. **Layer 2 - Datasets:**
   - Named collections of inventory items or AnnData checkpoints
   - Tracks: dataset name, URI, phase, format, size
   - Links projects to data

3. **Layer 3 - Projects:**
   - Metadata for analysis projects
   - Tracks: project name, platform, data_type, status, phases_complete

4. **Layer 4 - Provenance:**
   - Transformation records (tool versions, parameters, environment)
   - Tracks: input URIs, output URIs, tool name, parameters, environment

**Subjects/Patients (optional):**
- Patient/donor metadata (name, age, condition, etc.)
- Linked to inventory items for tracking subject samples

**BioData Modality Registry (in-memory):**
- Pure Python (no DB) registry of known platforms (Visium, Xenium, CosMx, IMC, etc.)
- Auto-populated at module import
- Extensible via `register_platform()` API

## Data Lifecycle

**Checkpoint Management:**
- Phase checkpoints saved as AnnData (`.h5ad`)
- URI tracked in registry (local or remote)
- Validation available via `validate_checkpoint()` MCP tool
- Recovery from failed phases via stored checkpoint URIs

**Provenance Tracking:**
- `record_provenance()` MCP tool logs transformation steps
- Captures: input data, tool, version, parameters, environment snapshot
- Queryable via `get_provenance()` for audit trails

---

*Integration audit: 2026-03-20*
