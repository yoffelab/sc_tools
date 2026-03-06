sc_tools.ingest — Data Ingestion
=================================

Phase 0 data loading: batch manifests, platform command builders,
modality-specific AnnData loaders, and sample concatenation.

.. code-block:: python

   import sc_tools.ingest as ingest

   # Load batch manifest
   manifest = ingest.load_batch_manifest("metadata/phase0/batch1_samples.tsv")

   # Load a single IMC sample
   adata = ingest.load_imc_sample("processed/sample_01", sample_id="sample_01")

   # Concatenate all samples
   adata_all = ingest.concat_samples([adata1, adata2, adata3])

Batch Manifests
---------------

.. currentmodule:: sc_tools.ingest.config

.. autofunction:: load_batch_manifest

.. autofunction:: collect_all_batches

.. autofunction:: validate_manifest

Modality Loaders
----------------

Each loader sets ``obs['sample']``, ``obs['library_id']``,
``obs['raw_data_dir']``, ``obsm['spatial']``, and ensures ``X`` contains raw counts.

.. currentmodule:: sc_tools.ingest.loaders

.. autofunction:: load_visium_sample

.. autofunction:: load_visium_hd_sample

.. autofunction:: load_visium_hd_cell_sample

.. autofunction:: load_xenium_sample

.. autofunction:: load_imc_sample

.. autofunction:: load_he_image

.. autofunction:: concat_samples

SLURM / HPC Helpers
--------------------

.. currentmodule:: sc_tools.ingest.slurm

.. autofunction:: build_sbatch_header

.. autofunction:: write_sbatch_script
