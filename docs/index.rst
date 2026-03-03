sc_tools
========

**sc_tools** is a reusable Python toolkit for spatial and single-cell multiomics analysis.
It follows the `scanpy <https://scanpy.readthedocs.io>`_ API pattern
(``sc_tools.pp``, ``sc_tools.pl``, ``sc_tools.tl``, ``sc_tools.qc``)
and covers the full pipeline from raw data ingestion through downstream biology.

.. code-block:: python

   import sc_tools as st
   import anndata as ad

   # Load preprocessed AnnData (Phase 3 checkpoint)
   adata = ad.read_h5ad("results/adata.normalized.p3.h5ad")

   # Score gene signatures (Phase 3.5b)
   signatures = st.tl.merge_gene_signatures(
       project_sigs, st.tl.load_hallmark()
   )
   st.tl.score_signature(adata, signatures)

   # Spatial visualization
   st.pl.multipage_spatial_pdf(adata, keys=["Hallmark/HYPOXIA"], output="figures/spatial.pdf")

Where to start
--------------

New to sc_tools? Start with the :doc:`installation` guide, then follow the
:doc:`tutorials/index`.

.. grid:: 1 2 2 3
   :gutter: 2

   .. grid-item-card:: Installation
      :link: installation
      :link-type: doc

      Get sc_tools installed with the right optional extras for your workflow.

   .. grid-item-card:: Tutorials
      :link: tutorials/index
      :link-type: doc

      Step-by-step notebooks covering each pipeline phase with synthetic data.

   .. grid-item-card:: API Reference
      :link: api/index
      :link-type: doc

      Complete function reference organized by module (pp, pl, tl, qc, memory, utils).

.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   tutorials/index
   api/index
