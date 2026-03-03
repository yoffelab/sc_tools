API Reference
=============

sc_tools follows the `scanpy <https://scanpy.readthedocs.io>`_ API pattern.
Import the package as ``import sc_tools as st`` and access modules via
``st.pp``, ``st.pl``, ``st.tl``, ``st.qc``, ``st.memory``, ``st.utils``.

.. list-table:: Module overview
   :header-rows: 1
   :widths: 15 85

   * - Module
     - Description
   * - :doc:`pp`
     - **Preprocessing** — normalization, batch integration, dimensionality reduction, clustering.
   * - :doc:`pl`
     - **Plotting** — spatial plots, heatmaps, QC plots, enrichment dot plots.
   * - :doc:`tl`
     - **Tools** — signature scoring, gene set loading, ORA/GSEA, deconvolution, colocalization.
   * - :doc:`qc`
     - **Quality control** — QC metrics, spot filtering, sample classification, HTML report.
   * - :doc:`memory`
     - **Memory / GPU** — GPU detection, memory profiling, aggressive cleanup.
   * - :doc:`utils`
     - **Utilities** — signature helpers, versioned save paths.

.. toctree::
   :maxdepth: 2
   :hidden:

   pp
   pl
   tl
   qc
   memory
   utils
