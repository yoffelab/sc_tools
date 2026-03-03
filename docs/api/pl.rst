sc_tools.pl — Plotting
=======================

Plotting utilities for spatial omics data following the scanpy API pattern.

.. code-block:: python

   import sc_tools.pl as pl

   pl.multipage_spatial_pdf(adata, keys=["Hallmark/HYPOXIA"], output_path="out.pdf")
   pl.qc_2x2_grid(adata)
   pl.plot_gsea_dotplot(results_df)

Spatial Plots
-------------

.. currentmodule:: sc_tools.pl.spatial

.. autofunction:: plot_spatial_plain_he

.. autofunction:: plot_spatial_categorical

.. autofunction:: plot_spatial_continuous

.. autofunction:: multipage_spatial_pdf

Heatmaps
--------

.. currentmodule:: sc_tools.pl.heatmaps

.. autofunction:: signature_score_heatmap

.. autofunction:: cluster_within_groups

.. autofunction:: annotation_colors_from_categories

.. autofunction:: get_obs_category_colors

Enrichment
----------

.. currentmodule:: sc_tools.pl

.. autofunction:: plot_gsea_dotplot

QC Plots
--------

Re-exported from :mod:`sc_tools.qc.plots` for ``st.pl.qc_*`` usage.
See :doc:`qc` for full documentation.

.. currentmodule:: sc_tools.pl

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

Utilities
---------

.. autofunction:: save_figure
