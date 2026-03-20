# Coding Conventions

**Analysis Date:** 2026-03-20

## Naming Patterns

**Files:**
- All lowercase with underscores: `normalize.py`, `sample_qc.py`, `report.py`
- Test files: `test_<module>.py` (e.g., `test_pp.py`, `test_qc.py`)
- Package directories: lowercase, no underscores (e.g., `sc_tools/pp`, `sc_tools/qc`)

**Functions:**
- snake_case throughout: `backup_raw()`, `normalize_total()`, `filter_genes_by_pattern()`, `calculate_qc_metrics()`
- Private helper functions prefixed with underscore: `_make_data_source()`, `_minimal_adata()`, `_protein_modalities`
- Constants: UPPER_SNAKE_CASE (e.g., `DEFAULT_FILTER_PATTERNS`, `_SPOT_FILTER_DEFAULTS`)

**Variables:**
- snake_case: `n_obs`, `n_vars`, `adata`, `var_names`, `mt_pattern`
- Loop variables: single letters acceptable (e.g., `i`, `j`, `p`)
- Constants in module scope: UPPER_SNAKE_CASE (prefixed with `_` if private): `_SPOT_FILTER_DEFAULTS`, `DEFAULT_FILTER_PATTERNS`

**Types:**
- Python 3.11+ type hints: `AnnData`, `list[str]`, `dict[str, Any]`, `float | None`
- Use `from __future__ import annotations` in all files for forward compatibility
- Parameters documented with type hints in docstring (see section below)

**Imports:**
```python
from __future__ import annotations  # Always first

import logging  # Standard library
import re
from typing import Any  # typing module last in stdlib

import numpy as np  # Third-party
import pandas as pd
from anndata import AnnData
from scipy import sparse

from ._gpu import get_backend  # Relative imports last
```

Order:
1. `from __future__ import annotations`
2. Standard library (sorted alphabetically)
3. Third-party packages (numpy, pandas, scipy, anndata, etc.)
4. Relative imports from same package

## Code Style

**Formatting:**
- Tool: `ruff format` (enforced by CI)
- Line length: 100 characters (configured in `pyproject.toml`)
- Quote style: Double quotes (`"string"` not `'string'`)

**Linting:**
- Tool: `ruff check` (configured in `pyproject.toml`)
- Active rules: E, W, F, I (isort), B (flake8-bugbear), C4, UP
- Ignored: E501 (line length handled by formatter), B008, E741 (Moran's I), E402 (lazy loading in scripts)

**Line Length:**
- Target: 100 characters per line
- Docstrings can exceed but should remain readable

## Import Organization

**Order (strictly enforced by ruff lint rule I):**
1. `from __future__ import annotations`
2. Standard library: `logging`, `re`, `from typing import Any`
3. Data science stack: `numpy`, `pandas`, `scipy`, `anndata`
4. Relative imports: `from ._gpu import get_backend`

**Path Aliases:**
- None configured in current codebase
- Use absolute imports from package root: `from sc_tools.pp import backup_raw`

## Error Handling

**Patterns:**
- Use specific exception types: `ValueError`, `NotImplementedError`
- Raise with descriptive messages formatted with f-strings:
```python
raise ValueError(
    f"Unsupported IMC normalization method '{method}'. Choose from: 'arcsinh', 'clr'."
)
raise NotImplementedError(
    "CLR normalization for IMC is planned for Phase 2. Use method='arcsinh' for now."
)
```
- Check for required columns/attributes:
```python
if sample_col not in adata.obs.columns:
    raise ValueError(f"sample_col={sample_col!r} not in adata.obs.columns")
```
- Use try/except for optional operations (e.g., parsing Space Ranger metrics):
```python
try:
    # Optional operation
    metrics = parse_space_ranger_metrics(...)
except Exception:
    logger.warning("Could not parse Space Ranger metrics for %s", sample)
```

**No custom exception classes** — use built-in types with clear messages.

## Logging

**Framework:** Python built-in `logging` module

**Pattern (all modules):**
```python
import logging
logger = logging.getLogger(__name__)
```

**Usage:**
- Info level for major steps: `logger.info("Normalizing adata (target_sum=%s)", target_sum)`
- Warning level for recoverable issues: `logger.warning("Could not parse metrics for %s", sample)`
- Debug level rarely used; info is preferred for transparency
- Include relevant metrics: counts, dimensions, parameters

**Examples from codebase:**
```python
logger.info("normalize_total (target_sum=%s, backend=%s)", target_sum, name)
logger.info("Removing %d/%d genes matching patterns %s", n_matching, adata.n_vars, patterns)
logger.info("Filtered AnnData saved: %s (%d obs)", op, adata.n_obs)
```

## Docstrings

**Format:** NumPy/SciPy docstring style (widely used in bioinformatics)

**Structure:**
```python
def function_name(
    adata: AnnData,
    *,
    mt_pattern: str | re.Pattern | None = "^(MT-|mt-|Mt)",
    inplace: bool = True,
    **kwargs: Any,
) -> pd.DataFrame | None:
    """One-line summary of what function does.

    Optional longer description explaining behavior,
    context, or important notes.

    Parameters
    ----------
    adata : AnnData
        Description of this parameter.
    mt_pattern : str or compiled regex or None
        Description. Default value and example.
    inplace : bool
        If True, modify in place; otherwise return copy (default True).
    **kwargs
        Passed to underlying function.

    Returns
    -------
    DataFrame or None
        What is returned. If inplace=True, return None; if False, return object.

    Notes
    -----
    Optional section for implementation details or warnings.
    """
```

**Key patterns:**
- First line is concise one-liner
- Parameters section always present
- Type hints in docstring match function signature
- Return section specifies None when inplace=True
- Defaults documented in parameter descriptions
- Use "if True/False" language for boolean params with effects
- Use "or None" to indicate optional parameters

**Example:**
```python
def backup_raw(adata: AnnData) -> None:
    """Save a copy of the current adata to adata.raw (no-op if already set).

    Parameters
    ----------
    adata
        Annotated data matrix. Modified in place.
    """
```

## Comments

**When to Comment:**
- Explain WHY, not WHAT (code should be self-documenting)
- Clarify non-obvious algorithms or parameter choices
- Document intentional design decisions (e.g., lenient defaults)
- Comment out large blocks of old code (don't delete if uncertain)

**Examples:**
```python
# Default gene patterns to exclude: mitochondrial, ribosomal, hemoglobin
DEFAULT_FILTER_PATTERNS = [
    r"^MT-",  # mitochondrial
    r"^RP[SL]",  # ribosomal
    r"^HB[^(P)]",  # hemoglobin (but not HBEGF, HBP1, etc.)
]

# Default thresholds — intentionally very lenient
_SPOT_FILTER_DEFAULTS: dict[str, dict[str, Any]] = {
    "visium": {"min_counts": 50, "min_genes": 20, "max_pct_mt": None},
    ...
}
```

## Function Design

**Size:**
- Keep functions focused on single responsibility
- 20-60 lines typical for public functions
- Longer functions should break into helpers

**Parameters:**
- Use keyword-only arguments for optional params: `def func(adata, *, option=default):`
- Group related boolean flags
- Defaults should be safe/lenient (especially for QC thresholds)
- Accept `**kwargs` to pass through to backend functions

**Return Values:**
- Inplace operations return None (pattern: `inplace=True` by default)
- Copy operations return AnnData or result object
- Functions that return None should document this in return section

**Inplace Pattern (consistent throughout codebase):**
```python
def transform(
    adata: AnnData,
    inplace: bool = True,
) -> AnnData | None:
    """Transform adata.

    Parameters
    ----------
    adata
        Annotated data.
    inplace
        If True, modify in place. Otherwise return copy.

    Returns
    -------
    AnnData or None
        Modified adata if inplace=False, else None.
    """
    if not inplace:
        adata = adata.copy()

    # Modify adata
    adata.X = np.arcsinh(adata.X / cofactor)

    if not inplace:
        return adata
    return None
```

## Module Design

**Exports:**
- All public functions listed in `__all__`:
```python
__all__ = [
    "normalize_total",
    "log_transform",
    "scale",
    "arcsinh_transform",
    "filter_genes_by_pattern",
    "backup_raw",
    "normalize_imc",
]
```
- Private functions/constants prefixed with `_` (not in `__all__`)

**Barrel Files:**
- Package `__init__.py` files re-export public functions for convenience:
```python
# sc_tools/pp/__init__.py
from .normalize import (
    arcsinh_transform,
    backup_raw,
    filter_genes_by_pattern,
    ...
)
```
- Users can `from sc_tools.pp import backup_raw` directly

**Module Docstring:**
- First line: brief description
- Optional: usage examples with code blocks
- Example from `sc_tools/pp/__init__.py`:
```python
"""sc_tools.pp: Modality-aware preprocessing for spatial omics data.

Provides normalization, batch integration, dimensionality reduction,
clustering, and full preprocessing recipes for Visium, Visium HD,
Xenium, CosMx, and IMC data.

GPU acceleration via rapids-singlecell is auto-detected; falls back to scanpy.

Usage (recipe)::

    import sc_tools.pp as pp
    adata = pp.preprocess(adata, modality="visium", batch_key="library_id")

Usage (individual steps)::

    import sc_tools.pp as pp
    pp.backup_raw(adata)
    pp.normalize_total(adata)
    pp.log_transform(adata)
"""
```

## Type Hints

**Style:**
- Use `from __future__ import annotations` for forward compatibility
- Prefer pipe operator: `float | None` over `Union[float, None]`
- Explicit types: `list[str]`, `dict[str, Any]`, not `List`, `Dict`

**Examples:**
```python
def func(
    adata: AnnData,
    patterns: list[str] | None = None,
    exclude: bool = True,
) -> AnnData | None:
```

**Generic types:**
```python
from typing import Any
def func(adata: AnnData, **kwargs: Any) -> None:
```

---

*Convention analysis: 2026-03-20*
