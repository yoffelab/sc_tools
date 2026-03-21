"""
sc-tools: Generic reusable utilities for spatial omics analysis.

This package provides modular, reusable tools following scanpy's API pattern:
- sc_tools.pp: Preprocessing (normalization, integration, clustering, modality recipes)
- sc_tools.pl: Plotting utilities (spatial, heatmaps, statistical annotations)
- sc_tools.tl: Analysis tools (statistical testing, colocalization, etc.)
- sc_tools.qc: QC metrics, filters, spatially variable genes, QC report plotting
- sc_tools.data: Data loading, preprocessing, and I/O
- sc_tools.memory: Memory management and GPU utilities

All code in this package should be generic and not project-specific.
"""

__version__ = "0.1.0"

# Lazy imports: heavy submodules (scanpy, dask, anndata) are only loaded on
# first attribute access.  This keeps ``sct --help`` fast (CLI-08).
_SUBMODULES = {
    "bm", "gr", "ingest", "memory", "pl", "pp", "qc",
    "storage", "tl", "utils", "validate", "registry",
}

__all__ = sorted(_SUBMODULES)


def __getattr__(name: str):
    if name in _SUBMODULES:
        import importlib

        try:
            mod = importlib.import_module(f".{name}", __name__)
        except ImportError:
            if name == "registry":
                raise AttributeError(name) from None
            raise
        globals()[name] = mod
        return mod
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
