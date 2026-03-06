sc_tools.storage — Storage Abstraction
========================================

fsspec-based URI resolution for reading and writing data across
local filesystem, SFTP/HPC, S3, GCS, Azure, and Box.

.. code-block:: python

   from sc_tools.storage import smart_read_h5ad, smart_write_checkpoint

   # Read from HPC scratch via SFTP
   adata = smart_read_h5ad("sftp://brb//athena/.../adata.raw.h5ad")

   # Write to S3
   smart_write_checkpoint(adata, "s3://bucket/results/adata.normalized.h5ad")

Install remote backends: ``pip install sc-tools[storage]``

.. currentmodule:: sc_tools.storage

.. autofunction:: resolve_fs

.. autofunction:: open_file

.. autofunction:: with_local_copy

.. autofunction:: smart_read_h5ad

.. autofunction:: smart_write_checkpoint

.. autofunction:: smart_read_csv
