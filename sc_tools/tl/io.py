"""
Result saving for tl (e.g. h5ad). Versioned by default.

Re-exports from sc_tools.data.io for use as st.tl.io_save.write_h5ad.
"""

from ..data.io import write_h5ad

__all__ = ["write_h5ad"]
