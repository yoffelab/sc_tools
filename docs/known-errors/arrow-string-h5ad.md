# ArrowStringArray h5ad Serialization Error

## Error

```
anndata._io.specs.registry.IORegistryError: No method registered for writing
<class 'pandas.arrays.ArrowStringArray'> into <class 'h5py._hl.group.Group'>
Error raised while writing key '_index' of <class 'h5py._hl.group.Group'> to /obs
```

## Cause

Pandas 2.x with pyarrow installed silently converts string columns to `ArrowStringArray` (backed by `pyarrow.StringType`) during:
- Reading CSVs / parquet files
- Concatenating DataFrames
- h5ad round-trips on newer anndata versions

The `h5py` backend used by `anndata.write_h5ad()` has no registered writer for `ArrowStringArray`. This affects `obs`, `var`, and their indices.

Categorical columns can also hide Arrow-backed categories underneath, causing the same error.

## Fix

**Use `sc_tools.data.io.write_h5ad()`** instead of calling `adata.write_h5ad()` directly. The centralized function automatically coerces Arrow strings before writing.

For robin project scripts, use `robin_utils.save_checkpoint()` which calls `_coerce_arrow_strings()`.

If you must call `adata.write_h5ad()` directly (e.g., in tests or one-off scripts):

```python
from sc_tools.data.io import _coerce_arrow_strings

_coerce_arrow_strings(adata)
adata.write_h5ad(path)
```

## Prevention

Add this at the top of scripts that read external data:

```python
import pandas as pd
pd.set_option("future.infer_string", False)
try:
    pd.options.mode.string_storage = "python"
except Exception:
    pass
```

This prevents pandas from creating Arrow-backed strings in the first place.

## History

This error has been independently fixed 5+ times across the codebase:

| Location | Fix | Date |
|----------|-----|------|
| `sc_tools/scripts/run_qc_report.py` | `_coerce_arrow_strings()` on read | 2026 |
| `sc_tools/scripts/run_preprocessing.py` | `_coerce_arrow_strings()` on read + `pd.set_option` | 2026 |
| `sc_tools/scripts/run_integration_benchmark.py` | `_coerce_arrow_strings_df()` + `pd.set_option` | 2026 |
| `robin/scripts/robin_utils.py` | `_coerce_arrow_strings()` in `save_checkpoint()` | 2026-03 |
| `robin/scripts/run_signature_scoring_v2.py` | Inline coercion before `write_h5ad()` | 2026-03-23 |
| `robin/scripts/bin_categorization.py` | Inline coercion before `write_h5ad()` | 2026-03-23 |
| **`sc_tools/sc_tools/data/io.py`** | **Centralized fix in `write_h5ad()`** | 2026-03-23 |

## Rule

**Never call `adata.write_h5ad()` directly in production scripts.** Use `sc_tools.data.io.write_h5ad()` or `robin_utils.save_checkpoint()`.
