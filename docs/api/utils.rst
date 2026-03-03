sc_tools.utils — Utilities
===========================

General-purpose helpers for signature score retrieval and versioned file paths.

.. code-block:: python

   from sc_tools.utils.signatures import get_signature_columns, get_signature_df

   cols = get_signature_columns(adata)
   score_df = get_signature_df(adata)

.. currentmodule:: sc_tools.utils

Signature Helpers
-----------------

.. autofunction:: get_signature_columns

.. autofunction:: filter_signatures

.. autofunction:: clean_sig_name

Versioned Save Paths
--------------------

.. autofunction:: get_version_prefix

.. autofunction:: versioned_filename

.. autofunction:: versioned_path
