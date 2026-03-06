sc_tools.bm — Benchmarking
===========================

Integration and segmentation benchmarking tools.

Integration Benchmark
---------------------

Compare batch-correction methods and select the best by batch score.

.. code-block:: python

   from sc_tools.bm import run_full_integration_workflow

   adata, comparison_df, best_method = run_full_integration_workflow(
       adata,
       modality="visium",
       batch_key="library_id",
       methods=["harmony", "combat", "scvi"],
       output_dir="results",
   )

.. currentmodule:: sc_tools.bm.integration

.. autofunction:: run_full_integration_workflow

.. autofunction:: run_integration_benchmark

.. autofunction:: compute_integration_metrics

.. autofunction:: compute_composite_score

.. autofunction:: compare_integrations

Segmentation Benchmark
----------------------

Compare cell segmentation methods with panoptic quality, boundary F1,
and cell-type preservation metrics.

.. currentmodule:: sc_tools.bm.segmentation

.. autofunction:: compute_segmentation_accuracy

.. autofunction:: score_segmentation

.. currentmodule:: sc_tools.bm.segment

.. autofunction:: run_cellpose

.. autofunction:: run_stardist

.. autofunction:: run_deepcell

Mask I/O
--------

.. currentmodule:: sc_tools.bm.mask_io

.. autofunction:: load_mask

.. autofunction:: load_tiff_mask

.. autofunction:: load_cellpose_mask

.. autofunction:: load_stardist_mask

.. autofunction:: load_deepcell_mask

.. autofunction:: load_cellprofiler_mask
