sc_tools.registry — Project Registry
======================================

SQLAlchemy-based registry for tracking projects, datasets, SLURM jobs,
agent tasks, and pipeline phase status. SQLite by default; PostgreSQL
via ``SC_TOOLS_REGISTRY_URL`` environment variable.

.. code-block:: python

   from sc_tools.registry import Registry

   reg = Registry()
   reg.add_project("ggo_visium", platform="visium")
   reg.register_dataset("ggo_visium", "results/adata.raw.h5ad", phase="qc_filter")
   reg.mark_phase_complete("ggo_visium", "qc_filter")
   print(reg.status())

Install: ``pip install sc-tools[registry]``

CLI: ``python -m sc_tools registry status``

.. currentmodule:: sc_tools.registry

.. autoclass:: Registry
   :members:
   :undoc-members: False
