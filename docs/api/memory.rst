sc_tools.memory — GPU and Memory
==================================

GPU detection and memory profiling utilities.

.. code-block:: python

   import sc_tools.memory as memory

   if memory.check_gpu_available():
       print("GPU ready")

   memory.log_memory("before deconvolution")

.. currentmodule:: sc_tools.memory

GPU
---

.. autofunction:: check_gpu_available

Memory Profiling
----------------

.. autofunction:: get_memory_usage

.. autofunction:: log_memory

.. autofunction:: aggressive_cleanup

.. autofunction:: estimate_adata_memory

.. autofunction:: check_memory_threshold
