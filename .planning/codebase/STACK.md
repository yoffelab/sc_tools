# Technology Stack

**Analysis Date:** 2026-03-20

## Languages

**Primary:**
- Python 3.11+ (3.11, 3.12, 3.13, 3.14 supported) - Core analysis and pipeline infrastructure

## Runtime

**Environment:**
- Python 3.11 minimum (specified in `pyproject.toml` line 11)
- CUDA 13.0+ compatible (GPU support via RAPIDS)
- Linux/macOS/Windows support (tested on ubuntu-latest, macos-latest via CI)

**Package Manager:**
- pip (primary)
- uv (fast pip replacement, used in Docker)
- setuptools >= 61 (build system)

**Lockfile:**
- No pip.lock; depends on direct `pyproject.toml` constraints for reproducibility

## Frameworks & Core Libraries

**Scientific Computing:**
- `numpy >= 1.24` - Array operations
- `scipy >= 1.10` - Scientific algorithms
- `pandas >= 2.0` - Data manipulation (TSV/CSV handling)
- `scikit-learn >= 1.3` - Machine learning algorithms

**Single-Cell / Spatial Analysis:**
- `scanpy >= 1.9` - Single-cell analysis toolkit (wrapped by sc_tools)
- `squidpy >= 1.3` - Spatial transcriptomics analysis
- `anndata >= 0.10` - AnnData format standard (primary data container)

**Visualization:**
- `matplotlib >= 3.7` - Low-level plotting
- `seaborn >= 0.12` - Statistical plotting
- `plotly >= 5.18` - Interactive HTML reports (optional, in [benchmark] extra)
- `statannotations >= 0.6` - Statistical annotation on plots
- `pingouin >= 0.5` - Statistical tests for plotting
- `marsilea >= 0.4` - Publication-quality composite figures (optional, in [viz] extra)

**Data Processing & Performance:**
- `dask >= 2024.1` - Distributed computing / large data handling
- `dask-expr >= 1.0` - Dataframe expressions for dask
- `tqdm >= 4.65` - Progress bars
- `tifffile >= 2023.1.0` - TIFF image I/O
- `Pillow >= 10.0` - Image processing
- `statsmodels >= 0.14` - Statistical modeling

**Deconvolution & Integration (Optional):**
- `scvi-tools >= 1.0` - Variational inference for single-cell (included in [pipeline] extra)
- `tangram-sc >= 1.0` - Cell-type deconvolution (included in [pipeline] extra)
- `harmonypy` - Batch integration (included in [pipeline] extra)
- `bbknn` - Batch correction via balanced k-NN (included in [pipeline] extra)
- `scanorama` - Batch integration (included in [pipeline] extra)
- `leidenalg` - Leiden clustering algorithm (included in [pipeline] extra)
- `igraph` - Graph algorithms (included in [pipeline] extra)
- `lifelines >= 0.27` - Survival analysis (included in [pipeline] extra)

**Gene Set & Pathway Analysis (Optional):**
- `gseapy >= 1.0` - Gene set enrichment analysis (included in [geneset] extra)
- `pyucell >= 0.1` - uCell gene set scoring (included in [geneset] extra)
- `decoupler >= 1.4` - Transcription factor/pathway activity (included in [decoupler] extra)

**Segmentation & Benchmarking (Optional):**
- `scikit-image >= 0.21` - Image analysis (included in [benchmark] extra)
- `scib-metrics >= 0.4` - Single-cell integration metrics (included in [pipeline] and [benchmark] extras)
- `cellpose >= 3.0,<4.0` - Cell segmentation (included in [benchmark] extra)
- `stardist >= 0.9` - Instance segmentation (included in [benchmark] extra)
- `deepcell >= 0.12` - Deep learning segmentation (included in [benchmark-extended] extra)

**Cell Typing (Optional):**
- `celltypist >= 1.3` - Automated cell type annotation (included in [celltyping] extra)
- `scarches >= 0.5` - scArches for transfer learning (included in [celltyping] extra)

**Foundation Models (Optional):**
- `transformers >= 4.36` - Hugging Face transformers (included in [foundation] extra)
- `datasets >= 2.16` - Dataset library (included in [foundation] extra)
- `accelerate >= 0.26` - Distributed training (included in [foundation] extra)

**GPU Acceleration (Optional):**
- `torch >= 2.0` - PyTorch (GPU compute, installed first in Docker before rapids)
- `rapids-singlecell` - GPU-accelerated scanpy drop-in (auto-detected in `sc_tools.pp._gpu`)
- Uses cuML, cuGraph, cuDF from RAPIDS base image for GPU-accelerated neighbors, leiden, umap, pca

## Database & ORM

**Database Support:**
- `sqlalchemy >= 2.0` - SQL ORM (included in [registry] extra)
- `alembic >= 1.13` - Database migrations (included in [registry] extra)
- `psycopg2-binary >= 2.9` - PostgreSQL driver (included in [registry] extra)

**Default Backend:**
- SQLite at `~/.sc_tools/registry.db` (zero-config, no extra deps)
- PostgreSQL via `SC_TOOLS_REGISTRY_URL` environment variable (optional)

**Schema:**
- Four-layer schema: Projects → Datasets → BioData (inventory) → DataSources
- Alembic migrations in `sc_tools/migrations/versions/` (18 versions as of 2026-03-20)
- Migration entry point: `sc_tools/migrations/env.py`

## Storage & Remote Backends

**Local Filesystem:**
- Native Python pathlib (always available)

**Remote Storage (Optional):**
- `fsspec >= 2024.1` - Unified filesystem abstraction (included in [storage] extra)
- `s3fs >= 2024.1` - AWS S3 access
- `sshfs >= 0.13` - SFTP/HPC access (for `sftp://host/path` URIs)
- `gcsfs >= 2024.1` - Google Cloud Storage
- `adlfs >= 2024.1` - Azure Blob Storage
- `zarr >= 2.18` - Zarr format I/O (included in [storage] extra)
- `ome-zarr >= 0.9` - OME-Zarr format support
- `boxfs >= 0.2` - Box.com cloud storage (separate [storage-box] extra, requires OAuth)

**URI Schemes Supported:**
- `file://` or `/path` - local filesystem
- `sftp://host/path` - SSH/HPC clusters
- `s3://bucket/path` - AWS S3
- `gs://bucket/path` - Google Cloud Storage
- `az://container/path` - Azure Blob
- `box://path` - Box.com (requires boxfs OAuth)

## Model Context Protocol (MCP)

**Framework:**
- `mcp >= 1.0` - Model Context Protocol (included in [mcp] extra)

**Servers:**
- `sc_tools.mcp.registry_server` - Exposes sc-registry as MCP tools
- `sc_tools.mcp.tools_server` - Exposes analysis tools as MCP tools

**MCP Tools Exposed:**
- Registry tools: `registry_status`, `set_phase_status`, `get_phase_status`, `mark_phase_complete`, etc.
- Analysis tools: `validate_checkpoint`, `generate_qc_report`, `score_gene_signatures`, `run_deconvolution`, etc.

## Development & Testing

**Testing Framework:**
- `pytest >= 7.0` - Test runner
- `pytest-cov >= 4.0` - Coverage reporting

**Linting & Formatting:**
- `ruff >= 0.1.0` - Fast Python linter and formatter
- Target version: Python 3.11 (specified in `pyproject.toml` tool.ruff section)
- Line length: 100 characters
- Rules: E (pycodestyle errors), W (warnings), F (Pyflakes), I (isort), B (flake8-bugbear), C4 (comprehensions), UP (pyupgrade)

**Documentation:**
- `sphinx >= 7.2` - Documentation generator
- `pydata-sphinx-theme >= 0.15` - PyData theme
- `sphinx-autodoc-typehints >= 2.0` - Type hints in docs
- `sphinx-copybutton >= 0.5` - Copy code button
- `sphinx-design >= 0.5` - Design components
- `myst-nb >= 1.0` - Jupyter notebook support
- `ipykernel` - Jupyter kernel

## Containerization

**Docker:**
- Base image: `nvcr.io/nvidia/rapidsai/base:25.12-cuda13-py3.12` (RAPIDS 25.02 with CUDA 13.0, Python 3.12)
- Built with GPU support (CUDA 13.0 compatible with brb and cayuga HPC systems)
- System deps: build-essential, gfortran, libhdf5-dev (HDF5 for h5ad/AnnData), git, curl

**Containerization Tools:**
- Docker build (`docker build -t sc_tools:latest .`)
- Apptainer/Singularity conversion (`apptainer build containers/sc_tools.sif docker-daemon://sc_tools:latest`)

**Container Extras Baked In:**
- [pipeline] extra (full runtime deps): scvi-tools, tangram, harmony, bbknn, scanorama, leidenalg, igraph, lifelines, scib-metrics, plotly, sqlalchemy, marsilea, fsspec, s3fs, zarr

## CI/CD Pipeline

**Version Control:**
- Git (GitHub)
- Repository: https://github.com/yoffelab/sc_tools

**CI/CD Framework:**
- GitHub Actions (`.github/workflows/ci.yml`)
- Stages: lint (ruff), test (multi-Python 3.11-3.14, multi-OS), docs, docker build
- Concurrency control with auto-cancel on new push to same ref
- Test matrix: Python 3.11, 3.12, 3.13, 3.14 × (Ubuntu, macOS)

**Documentation Hosting:**
- ReadTheDocs (https://sc-tools.readthedocs.io)
- Artifacts: HTML docs uploaded after main branch builds (14-day retention)

## Configuration Files

**Python Configuration:**
- `pyproject.toml` - Source of truth for all deps, metadata, tool configs
- `setup.py` / `setup.cfg` - Not used (setuptools.build_meta via pyproject.toml)

**Linting & Formatting:**
- Ruff config in `pyproject.toml` [tool.ruff] section
- No separate `.eslintrc`, `.prettierrc`, or `black` config (uses ruff for all)

**Testing:**
- Pytest config in `pyproject.toml` [tool.pytest.ini_options]
- Test discovery path: `sc_tools/tests/`
- Additional pytest plugins: none specified

**Build Tools:**
- Dockerfile at repo root
- No Makefile or custom build scripts
- Installation via `pip install -e ".[pipeline]"` (editable mode with pipeline extras)

## Environment Variables

**Critical:**
- `SC_TOOLS_REGISTRY_URL` - Override default SQLite registry with PostgreSQL URL (format: `postgresql://user:pass@host:port/dbname`)

**Optional Storage Auth:**
- AWS S3: `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY` (via s3fs)
- GCS: `GOOGLE_APPLICATION_CREDENTIALS` (service account JSON path)
- Azure: `AZURE_STORAGE_ACCOUNT_NAME`, `AZURE_STORAGE_ACCOUNT_KEY`
- Box.com: OAuth flow at runtime (boxfs)

**HPC Cluster Environment:**
- SLURM job scheduling (brb/cayuga clusters use SLURM)
- CUDA paths for GPU acceleration
- Conda environment activation for SLURM jobs

## Python Environment

**Editable Installation:**
```bash
pip install -e ".[pipeline,registry,mcp]"  # Full stack for pipeline runs
pip install -e ".[dev]"                    # For development + testing
pip install -e ".[storage,deconvolution]"  # A la carte extras
```

**Conda Support:**
- Not the primary distribution method
- `.nvmrc` / `.python-version` files not present; relies on GitHub Actions matrix for version testing

---

*Stack analysis: 2026-03-20*
