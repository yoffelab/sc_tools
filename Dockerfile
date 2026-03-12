# sc_tools pipeline container
# Build:  docker build -t sc_tools:latest .
# To SIF: apptainer build containers/sc_tools.sif docker-daemon://sc_tools:latest
# Run:    ./scripts/run_container.sh projects/visium/ggo_visium [command]
#         (auto-detects runtime; uses Docker on macOS, Apptainer on Linux)
#
# GPU (HPC): apptainer exec --nv --bind /athena containers/sc_tools.sif python script.py

FROM continuumio/miniconda3:latest

LABEL maintainer="Junbum Kim"
LABEL description="sc_tools: spatial and single-cell multiomics analysis"
LABEL org.opencontainers.image.source="https://github.com/yoffelab/sc_tools"

# System deps: HDF5 (for h5py/AnnData), build tools, curl (for uv)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    libhdf5-dev \
    git \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install uv (fast pip replacement)
ADD --chmod=755 https://astral.sh/uv/install.sh /uv-installer.sh
RUN /uv-installer.sh && rm /uv-installer.sh
ENV PATH="/root/.local/bin:$PATH"

# Create conda env with Python 3.11
RUN conda create -n sc_tools python=3.11 -y
ENV PATH="/opt/conda/envs/sc_tools/bin:$PATH"

WORKDIR /workspace

# Install sc_tools package (pyproject.toml defines all deps)
# The [pipeline] extra bundles all HPC runtime dependencies:
#   scvi-tools, tangram, harmony, bbknn, scanorama, leidenalg, igraph,
#   lifelines, scib-metrics, plotly, sqlalchemy, marsilea, fsspec/s3fs/zarr
COPY pyproject.toml README.md ./
COPY sc_tools/ ./sc_tools/

ENV UV_LINK_MODE=copy

# Step 1: Install PyTorch with CUDA 12.1 support FIRST so scvi-tools uses it.
# Both brb (L40S/A40/RTX6000) and cayuga (A100/A40/H100) run CUDA 12.x.
# CPU fallback still works if no GPU is present.
RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install --python /opt/conda/envs/sc_tools/bin/python \
    torch --index-url https://download.pytorch.org/whl/cu121

# Step 2: Install sc_tools with the full pipeline extra
RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install --python /opt/conda/envs/sc_tools/bin/python -e ".[pipeline]"

# Copy shared scripts (project files are bind-mounted at runtime, not baked in)
COPY scripts/ ./scripts/

# Verify install
RUN python -c "import sc_tools; import scanpy; import squidpy; import scvi; import harmonypy; import lifelines; print('sc_tools OK')"

# Default workdir; overridden by run_container.sh via -w flag
WORKDIR /workspace

ENTRYPOINT []
CMD ["bash"]
