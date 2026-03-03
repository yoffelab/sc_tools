sc_tools.tl — Analysis Tools
==============================

Analysis tools for spatial omics data following the scanpy API pattern.

.. code-block:: python

   import sc_tools.tl as tl

   hallmark = tl.load_hallmark()
   combined = tl.merge_gene_signatures(project_sigs, hallmark)
   tl.score_signature(adata, combined)

.. currentmodule:: sc_tools.tl

Gene Signature Scoring
----------------------

Scores are stored in ``adata.obsm['signature_score']`` (raw) and
``adata.obsm['signature_score_z']`` (z-scored). Column names are full paths,
e.g. ``Hallmark/HYPOXIA``, ``Myeloid/Macrophage_Core``.

.. autofunction:: score_signature

Gene Set Loaders and Curation
------------------------------

.. autofunction:: load_hallmark

.. autofunction:: load_msigdb_json

.. autofunction:: load_gmt

.. autofunction:: list_gene_sets

.. autofunction:: validate_gene_signatures

.. autofunction:: merge_gene_signatures

.. autofunction:: update_gene_symbols

.. autofunction:: save_gene_signatures

Enrichment Testing
------------------

Group-level enrichment. Optional dependency: ``pip install sc-tools[geneset]``.

.. autofunction:: run_ora

.. autofunction:: run_gsea_pseudobulk

Cell-type Deconvolution
-----------------------

Backend registry: cell2location, tangram, destvi.
Optional dependency: ``pip install sc-tools[deconvolution]``.

.. autofunction:: deconvolution

.. autofunction:: extract_reference_profiles

.. autofunction:: select_signature_genes

Colocalization
--------------

.. currentmodule:: sc_tools.tl.colocalization

.. autofunction:: truncated_similarity

.. autofunction:: pearson_correlation

.. autofunction:: morans_i

.. autofunction:: morans_i_batch

.. autofunction:: neighborhood_enrichment

.. autofunction:: neighborhood_enrichment_batch
