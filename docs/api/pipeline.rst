sc_tools.pipeline — Phase DAG
==============================

Defines the pipeline phase graph as semantic slugs with explicit dependencies.
Use this module to query available next steps and expected checkpoint paths.

.. code-block:: python

   from sc_tools.pipeline import get_available_next, get_phase_checkpoint

   # What can run after QC and metadata?
   completed = {"ingest_raw", "ingest_load", "qc_filter", "metadata_attach"}
   next_phases = get_available_next(completed)
   # -> {"preprocess"}

   # Expected checkpoint path for a phase
   path = get_phase_checkpoint("preprocess")
   # -> "results/adata.normalized.h5ad"

.. currentmodule:: sc_tools.pipeline

.. autoclass:: PhaseSpec

.. autofunction:: get_phase

.. autofunction:: get_available_next

.. autofunction:: get_phase_checkpoint

.. autofunction:: get_dag

.. autofunction:: extend_dag

.. autofunction:: validate_dag
