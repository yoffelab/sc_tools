sc_tools.validate — Checkpoint Validation
==========================================

Validates AnnData checkpoint files against the metadata contracts
defined in Architecture.md Section 2.2.

.. code-block:: python

   from sc_tools.validate import validate_checkpoint

   # Validate a Phase 1 checkpoint
   validate_checkpoint("results/adata.raw.h5ad", phase="qc_filter")

   # Validate with auto-fix (renames obs['batch'] -> obs['raw_data_dir'])
   validate_checkpoint("results/adata.raw.h5ad", phase="qc_filter", fix=True)

.. currentmodule:: sc_tools.validate

.. autofunction:: validate_checkpoint

.. autofunction:: validate_file

.. autofunction:: validate_p1

.. autofunction:: validate_p2

.. autofunction:: validate_p3

.. autofunction:: validate_p35

.. autofunction:: validate_p4

.. autoclass:: CheckpointValidationError
