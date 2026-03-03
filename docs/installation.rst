Installation
============

Requirements
------------

- Python 3.10 or later
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
     - harmonypy
     - ``pip install -e ".[integration]"``
   * - ``spatial``
     - utag (spatial-aware clustering)
     - ``pip install -e ".[spatial]"``
   * - ``gpu``
     - torch, rapids-singlecell
     - ``pip install -e ".[gpu]"``
   * - ``dev``
     - pytest, ruff
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
