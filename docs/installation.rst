Installation
============

Requirements
------------

- Python 3.11 or later
- `conda <https://conda.io>`_ (recommended for environment management)

Environment setup
-----------------

Create the conda environment from the repository root:

.. code-block:: bash

   conda env create -f environment.yml
   conda activate sc_tools

Core install
------------

Install sc_tools in editable mode (core dependencies only):

.. code-block:: bash

   pip install -e .

Optional extras
---------------

sc_tools uses optional dependencies for heavy ML and GPU workflows.
Install only what you need:

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Extra
     - Contents
     - Install command
   * - ``deconvolution``
     - scvi-tools, tangram-sc
     - ``pip install -e ".[deconvolution]"``
   * - ``geneset``
     - gseapy, pyucell
     - ``pip install -e ".[geneset]"``
   * - ``integration``
     - harmonypy, bbknn, scanorama
     - ``pip install -e ".[integration]"``
   * - ``spatial``
     - utag (spatial-aware clustering)
     - ``pip install -e ".[spatial]"``
   * - ``gpu``
     - torch, rapids-singlecell
     - ``pip install -e ".[gpu]"``
   * - ``storage``
     - fsspec, s3fs, sshfs, gcsfs, adlfs, zarr, ome-zarr
     - ``pip install -e ".[storage]"``
   * - ``storage-box``
     - boxfs (Box storage via OAuth)
     - ``pip install -e ".[storage-box]"``
   * - ``registry``
     - SQLAlchemy, Alembic (project/dataset tracking)
     - ``pip install -e ".[registry]"``
   * - ``mcp``
     - Model Context Protocol servers
     - ``pip install -e ".[mcp]"``
   * - ``benchmark``
     - scikit-image, scib-metrics, plotly, cellpose, stardist
     - ``pip install -e ".[benchmark]"``
   * - ``benchmark-extended``
     - All of benchmark + deepcell, torch, transformers, segmentation-models-pytorch
     - ``pip install -e ".[benchmark-extended]"``
   * - ``viz``
     - marsilea (publication composite figures)
     - ``pip install -e ".[viz]"``
   * - ``decoupler``
     - decoupleR (TF/pathway activity)
     - ``pip install -e ".[decoupler]"``
   * - ``dev``
     - pytest, ruff, igraph, leidenalg
     - ``pip install -e ".[dev]"``
   * - ``docs``
     - sphinx, pydata-sphinx-theme, myst-nb
     - ``pip install -e ".[docs]"``

Full pipeline install (deconvolution + geneset + integration):

.. code-block:: bash

   pip install -e ".[deconvolution,geneset,integration]"

Container usage
---------------

All pipeline runs are containerized. On Linux/HPC use Apptainer (primary);
on macOS/Windows use Docker (fallback).

.. code-block:: bash

   # Run a script inside the container (auto-detects Apptainer vs Docker)
   ./scripts/run_container.sh projects/visium/ggo_visium python scripts/foo.py

   # Interactive shell
   ./scripts/run_container.sh projects/visium/ggo_visium

   # Force Docker on any platform
   SC_TOOLS_RUNTIME=docker ./scripts/run_container.sh projects/visium/ggo_visium

Build the Apptainer SIF from the Docker image:

.. code-block:: bash

   docker build -t sc_tools:latest .
   apptainer build containers/sc_tools.sif docker-daemon://sc_tools:latest

Build the documentation
-----------------------

.. code-block:: bash

   pip install -e ".[docs]"
   make docs           # build HTML
   make docs-open      # build and open in browser (macOS)
   make docs-clean     # remove build artifacts
