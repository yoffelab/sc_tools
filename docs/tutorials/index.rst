Tutorials
=========

The tutorials below use synthetic AnnData objects so no project data is required.
They follow the pipeline phases described in :doc:`../api/index`.

.. toctree::
   :maxdepth: 1
   :hidden:

   01_getting_started
   02_qc
   03_preprocessing
   04_signature_scoring
   05_deconvolution
   06_spatial_analysis

.. grid:: 1 2 2 3
   :gutter: 2

   .. grid-item-card:: 1. Getting Started
      :link: 01_getting_started
      :link-type: doc

      End-to-end minimal workflow: load data, preprocess, score signatures, plot.

   .. grid-item-card:: 2. QC and Filtering
      :link: 02_qc
      :link-type: doc

      Phase 1 — QC metrics, spot filtering, cross-sample comparison, HTML report.

   .. grid-item-card:: 3. Preprocessing
      :link: 03_preprocessing
      :link-type: doc

      Phase 3 — normalize, integrate (scVI / Harmony), cluster, UMAP.

   .. grid-item-card:: 4. Gene Signature Scoring
      :link: 04_signature_scoring
      :link-type: doc

      Phase 3.5b — load Hallmark sets, score_signature, ORA, GSEA dot plot.

   .. grid-item-card:: 5. Cell-type Deconvolution
      :link: 05_deconvolution
      :link-type: doc

      Phase 3.5b — extract reference profiles, run deconvolution, inspect proportions.

   .. grid-item-card:: 6. Spatial Analysis
      :link: 06_spatial_analysis
      :link-type: doc

      Phase 5 — colocalization, Morans I, neighborhood enrichment, spatial PDF.
