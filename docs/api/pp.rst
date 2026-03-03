sc_tools.pp — Preprocessing
============================

Modality-aware preprocessing with GPU auto-detection
(rapids-singlecell, falls back to scanpy).

.. code-block:: python

   import sc_tools.pp as pp

   # One-call recipe
   pp.preprocess(adata, modality="visium", integration="scvi", batch_key="library_id")

   # Or step-by-step
   pp.backup_raw(adata)
   pp.normalize_total(adata)
   pp.log_transform(adata)
   pp.pca(adata)
   pp.cluster(adata, resolution=0.8)

.. currentmodule:: sc_tools.pp

Recipe
------

.. autofunction:: preprocess

Normalization
-------------

.. autofunction:: backup_raw

.. autofunction:: normalize_total

.. autofunction:: log_transform

.. autofunction:: scale

.. autofunction:: arcsinh_transform

.. autofunction:: filter_genes_by_pattern

Integration
-----------

All integration functions are **soft dependencies**; install hints are printed
if the required package is missing.

.. autofunction:: run_scvi

.. autofunction:: run_harmony

.. autofunction:: run_cytovi

Dimensionality Reduction and Clustering
----------------------------------------

.. autofunction:: pca

.. autofunction:: neighbors

.. autofunction:: umap

.. autofunction:: leiden

.. autofunction:: cluster

.. autofunction:: run_utag
