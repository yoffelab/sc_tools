sc_tools.qc — Quality Control
==============================

QC metrics, spot filtering, sample classification, and HTML report generation.

.. code-block:: python

   import sc_tools.qc as qc

   qc.calculate_qc_metrics(adata, mt_pattern="^MT-")
   qc.filter_cells(adata, min_genes=200)
   qc.filter_genes(adata, min_cells=3)
   fig = qc.qc_2x2_grid(adata)

.. currentmodule:: sc_tools.qc

Metrics and Filtering
---------------------

Thin wrappers around scanpy ``pp.calculate_qc_metrics``, ``pp.filter_cells``,
``pp.filter_genes``, and ``pp.highly_variable_genes``.

.. autofunction:: calculate_qc_metrics

.. autofunction:: filter_cells

.. autofunction:: filter_genes

.. autofunction:: highly_variable_genes

Spatially Variable Genes
------------------------

Wraps squidpy ``gr.spatial_autocorr`` (Moran's I).

.. autofunction:: spatially_variable_genes

.. autofunction:: spatially_variable_genes_per_library

Sample-Level QC
---------------

Per-sample metrics, adaptive MAD outlier detection, pass/fail classification.

.. autofunction:: filter_spots

.. autofunction:: compute_sample_metrics

.. autofunction:: classify_samples

.. autofunction:: save_pass_fail_lists

.. autofunction:: apply_qc_filter

HTML Report
-----------

.. autofunction:: generate_qc_report

QC Plots
--------

.. currentmodule:: sc_tools.qc.plots

.. autofunction:: qc_2x2_grid

.. autofunction:: qc_2x4_pre_post

.. autofunction:: qc_spatial_multipage

.. autofunction:: qc_violin_metrics

.. autofunction:: qc_scatter_counts_genes

.. autofunction:: plot_highly_variable_genes

.. autofunction:: plot_spatially_variable_genes

.. autofunction:: qc_sample_comparison_bar

.. autofunction:: qc_sample_violin_grouped

.. autofunction:: qc_sample_scatter_matrix
