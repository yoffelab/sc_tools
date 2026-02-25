# sc_tools pipeline container
# Based on pyproject.toml. Use for Nextflow (Docker local / Singularity HPC).
# Uses conda env sc_tools + uv for fast package installation.
#
# Build:   docker build -t sc_tools:latest .
# Run:     docker run -v $(pwd):/workspace sc_tools:latest bash -c "cd projects/visium/ggo_visium && python scripts/loupe2adata.py"
# Or:      ./scripts/run_docker.sh projects/visium/ggo_visium

FROM continuumio/miniconda3:latest

LABEL maintainer="Junbum Kim"
LABEL description="sc_tools: spatial and single-cell multiomics analysis"
LABEL org.opencontainers.image.source="https://github.com/yoffelab/sc_tools"

# Install system deps (scientific stack + build tools + curl for uv)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libhdf5-dev \
    git \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install uv for fast package installation
ADD --chmod=755 https://astral.sh/uv/install.sh /uv-installer.sh
RUN /uv-installer.sh && rm /uv-installer.sh
ENV PATH="/root/.local/bin:$PATH"

# Create conda env sc_tools with Python 3.10
RUN conda create -n sc_tools python=3.10 -y
ENV PATH="/opt/conda/envs/sc_tools/bin:$PATH"

# Dask: enable new query-planning backend (legacy removed in 2024+)
ENV DASK_DATAFRAME__QUERY_PLANNING=true

WORKDIR /workspace

# Copy package definition and install sc_tools with deconvolution (uv pip into conda env)
COPY pyproject.toml README.md ./
COPY sc_tools/ ./sc_tools/

ENV UV_LINK_MODE=copy
RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install --python /opt/conda/envs/sc_tools/bin/python -e ".[deconvolution]"

# Copy shared scripts
COPY scripts/ ./scripts/

# Copy ggo_visium project (metadata, Makefile; data/results/figures excluded via .dockerignore)
RUN mkdir -p projects/visium
COPY projects/visium/ggo_visium/ ./projects/visium/ggo_visium/

# Ensure project dirs exist for volume mounts
RUN mkdir -p projects/visium/ggo_visium/data \
             projects/visium/ggo_visium/results \
             projects/visium/ggo_visium/figures \
             projects/visium/ggo_visium/outputs

# Default: run from project dir
WORKDIR /workspace/projects/visium/ggo_visium

# Verify install (runs in conda env via PATH)
RUN python -c "import sc_tools; import scanpy; import squidpy; print('sc_tools OK')"

# Default: bash (override to run pipeline)
ENTRYPOINT []
CMD ["bash"]
