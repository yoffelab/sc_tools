# Testing Patterns

**Analysis Date:** 2026-03-20

## Test Framework

**Runner:**
- pytest >= 7.0
- Config: `pyproject.toml [tool.pytest.ini_options]`

**Assertion Library:**
- Built-in `assert` statements
- NumPy testing: `numpy.testing.assert_allclose()`, `assert_array_equal()`

**Run Commands:**
```bash
pytest sc_tools/tests                    # Run all tests
pytest sc_tools/tests -v                 # Verbose output
pytest sc_tools/tests --tb=short         # Short traceback
pytest sc_tools/tests -k test_backup     # Run tests matching pattern
pytest sc_tools/tests --cov=sc_tools     # With coverage
```

**Configuration (pyproject.toml):**
```toml
[tool.pytest.ini_options]
testpaths = ["sc_tools/tests"]
addopts = "-v --tb=short"
filterwarnings = [
    "ignore::DeprecationWarning",
    "ignore::FutureWarning",
]
```

## Test File Organization

**Location:**
- Separate directory: `sc_tools/tests/` (not co-located with source)
- One test file per module: `test_pp.py` tests `sc_tools.pp`, `test_qc.py` tests `sc_tools.qc`

**Naming:**
- Test files: `test_<module>.py`
- Test classes: `Test<Feature>` (e.g., `TestNormalizeTotal`, `TestFilterGenesByPattern`)
- Test functions: `test_<behavior>()` (e.g., `test_backup_creates_raw`, `test_backup_noop_if_exists`)

**File Count:** 35 test files covering core functionality

**Structure:**
```
sc_tools/tests/
├── test_pp.py                 # Preprocessing tests
├── test_qc.py                 # Quality control tests
├── test_registry.py           # Database/registry tests
├── test_benchmark_*.py        # Benchmarking tests
├── test_integration_*.py      # Integration tests
├── test_ingest*.py            # Data ingestion tests
└── test_*_real_data.py        # Real data validation tests
```

## Test Structure

**Suite Organization:**
```python
"""Unit tests for sc_tools.pp preprocessing module.

Tests use synthetic AnnData fixtures. Integration tests (scVI, Harmony, CytoVI)
are skipped if the required packages are not installed.
"""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData
from scipy import sparse

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def adata_counts():
    """Synthetic count matrix with MT/RP/HB genes and library_id."""
    np.random.seed(42)
    n_obs, n_vars = 200, 500
    X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars)).astype("float32")
    gene_names = [f"GENE{i}" for i in range(n_vars - 8)]
    gene_names += ["MT-CO1", "MT-CO2", "MT-ND1"]
    gene_names += ["RPS6", "RPL11"]
    gene_names += ["HBA1", "HBA2", "HBB"]

    adata = AnnData(
        X,
        obs={"library_id": (["L1"] * 100) + (["L2"] * 100)},
    )
    adata.var_names = gene_names
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    return adata

# ---------------------------------------------------------------------------
# Test Class: Feature
# ---------------------------------------------------------------------------

class TestBackupRaw:
    def test_backup_creates_raw(self, adata_counts):
        from sc_tools.pp import backup_raw

        assert adata_counts.raw is None
        backup_raw(adata_counts)
        assert adata_counts.raw is not None
        assert adata_counts.raw.n_vars == 500

    def test_backup_noop_if_exists(self, adata_counts):
        from sc_tools.pp import backup_raw

        adata_counts.raw = adata_counts.copy()
        original_raw_shape = adata_counts.raw.shape
        backup_raw(adata_counts)
        assert adata_counts.raw.shape == original_raw_shape
```

**Patterns:**
- Module docstring explains scope and dependencies
- Fixtures section marked with `# -----------` dividers
- Test classes group related tests (one feature per class)
- Imports inline in test functions (ensures fresh state, avoids circular imports)
- Assertions use simple `assert` statements

## Test Structure

**Fixture Patterns:**
```python
@pytest.fixture()
def adata_counts():
    """Synthetic count matrix with MT/RP/HB genes and library_id."""
    # ... setup code
    return adata

@pytest.fixture()
def adata_sparse(adata_counts):
    """Same as adata_counts but with sparse X."""
    adata = adata_counts.copy()
    adata.X = sparse.csr_matrix(adata.X)
    return adata
```

**Fixture scope:** Function (default) — fresh fixture per test
- No setup/teardown needed for synthetic data
- Fixtures can depend on other fixtures

**Test Methods in Classes:**
```python
class TestNormalizeTotal:
    def test_normalize_inplace(self, adata_counts):
        from sc_tools.pp import normalize_total

        original_sums = adata_counts.X.sum(axis=1)
        normalize_total(adata_counts, target_sum=1e4)
        new_sums = adata_counts.X.sum(axis=1)
        np.testing.assert_allclose(new_sums, 1e4, rtol=1e-5)
        assert not np.allclose(original_sums, new_sums)

    def test_normalize_sparse(self, adata_sparse):
        from sc_tools.pp import normalize_total

        normalize_total(adata_sparse, target_sum=1e4)
        sums = np.asarray(adata_sparse.X.sum(axis=1)).flatten()
        np.testing.assert_allclose(sums, 1e4, rtol=1e-5)
```

## Mocking

**Framework:** None (minimal mocking philosophy)

**Approach:**
- Use synthetic AnnData fixtures instead of mocks (lightweight, deterministic)
- Import tested functions inside test methods to avoid import-time side effects
- No monkeypatch of external libraries

**When NOT to Mock:**
- Data loading functions — use fixtures or tmp_path
- QC/preprocessing functions — use synthetic AnnData
- Database operations — use file-backed SQLite in tmp_path (see `test_registry.py`)

**Example: Registry Tests (Database Isolation)**
```python
@pytest.fixture()
def reg(tmp_path):
    """Return a file-backed SQLite Registry for testing."""
    from sc_tools.registry import Registry

    db_url = f"sqlite:///{tmp_path / 'test_registry.db'}"
    return Registry(db_url=db_url)

def test_add_project(self, reg):
    pid = reg.add_project("proj_a", platform="visium", data_type="visium")
    assert isinstance(pid, int)
    assert pid > 0
```

## Fixtures and Factories

**Test Data:**
```python
def _minimal_adata(n_obs=80, n_vars=100, mt_genes=5):
    """Synthetic adata with MT- genes for QC metrics."""
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
    var_names = [f"MT-{i}" for i in range(mt_genes)] + [f"g{i}" for i in range(mt_genes, n_vars)]
    adata = sc.AnnData(
        X,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=var_names),
    )
    return adata

def _multi_sample_adata(n_samples=6, n_obs_per=50, n_vars=100, mt_genes=5):
    """Synthetic multi-sample adata for sample-level QC."""
    # ... setup multiple samples, concatenate
```

**Location:**
- Module-level helper functions (prefixed with `_`)
- Defined at top of test file after imports
- Called by test functions/classes

**Naming convention:**
- Fixture functions (pytest): `@pytest.fixture() def adata_counts():`
- Factory functions (manual): `def _minimal_adata(...):`

## Coverage

**Requirements:** No explicit coverage targets enforced

**View Coverage:**
```bash
pytest sc_tools/tests --cov=sc_tools --cov-report=html
open htmlcov/index.html  # View HTML report
pytest sc_tools/tests --cov=sc_tools --cov-report=term-missing
```

**Strategy:**
- Core functionality (QC, preprocessing) has comprehensive unit tests
- Real data tests optional (skipped if data unavailable)
- Edge cases tested (sparse matrices, small AnnData, missing columns)

## Test Types

**Unit Tests (primary):**
- Scope: Single function or class method
- Data: Synthetic AnnData fixtures with known shapes/values
- Location: `test_pp.py`, `test_qc.py`, etc.
- Examples:
  - `test_backup_creates_raw()` — validates backup creates adata.raw
  - `test_normalize_inplace()` — validates normalization math
  - `test_filter_genes_by_pattern()` — validates regex matching

**Integration Tests:**
- Scope: Multiple functions in sequence (e.g., preprocessing pipeline)
- Data: Synthetic or real project data
- Location: `test_pp_real_data.py`, `test_pipeline.py`
- Examples:
  - `test_preprocess_visium()` — full pipeline: backup → normalize → log → scale
  - `test_apply_qc_filter()` — backup → filter spots → filter genes → save

**Real Data Tests:**
- Scope: Validate against actual project data (IMC, Visium, CosMx)
- Data: Stored as checkpoints in `.sc_tools/checkpoints/`
- Location: `test_gr_real_data.py`, `test_pp_real_data.py`
- Status: Skipped if checkpoint not found
- Example:
  ```python
  @pytest.mark.skipif(not _IMC_AVAILABLE, reason=f"IMC checkpoint not found: {_IMC_PATH}")
  def test_preprocess_imc(imc_adata):
      # Run full QC/preprocessing on real IMC data
      result = sc_tools.pp.preprocess(imc_adata, modality="imc")
      assert result.shape[0] > 0
  ```

**E2E Tests:**
- Not formally separated; integration tests serve this purpose
- Some MCP tests execute full phase pipelines: `test_mcp_run_full_phase.py`

## Common Patterns

**Async Testing:**
- Not applicable (no async code in codebase)

**Error Testing:**
```python
def test_register_invalid_data_source_raises(self, reg):
    """ValueError raised when DataSource not found."""
    with pytest.raises(ValueError, match="DataSource.*not found"):
        reg.register_inventory_item(
            "item_bad",
            uri="/data/bad.h5ad",
            modality="rna",
            data_source_name="nonexistent",  # Invalid reference
        )

def test_compute_sample_metrics_missing_col(self):
    """Raises ValueError when sample_col not in obs."""
    adata = _minimal_adata()
    with pytest.raises(ValueError, match="not in adata.obs.columns"):
        compute_sample_metrics(adata, sample_col="nonexistent")
```

**Parametrized Tests (Limited):**
- Not extensively used; each scenario tested separately in class
- Example:
  ```python
  class TestFilterSpots:
      def test_filter_spots_visium_hd(self):
          """Spot filtering with Visium HD defaults."""
          # Test with visium_hd modality

      def test_filter_spots_xenium(self):
          """Spot filtering with Xenium defaults."""
          # Test with xenium modality
  ```

**Temporary Files (tmp_path):**
```python
def test_qc_2x2_grid_save(tmp_path):
    """Test saving figures to disk."""
    adata = _minimal_adata()
    calculate_qc_metrics(adata, inplace=True, percent_top=(10, 20))
    qc_2x2_grid(adata, output_dir=tmp_path, basename="qc_test")
    assert (tmp_path / "qc_test.pdf").exists()
    assert (tmp_path / "qc_test.png").exists()
```

**Conditional Skipping:**
```python
@pytest.mark.skipif(not _IMC_AVAILABLE, reason=f"IMC checkpoint not found: {_IMC_PATH}")
def test_preprocess_imc(imc_adata):
    # Only runs if checkpoint exists

pytest.importorskip("squidpy", reason="squidy required for gr tests")
# Skip entire test if import fails

@pytest.mark.skip(reason="squidpy optional; run when squidpy installed")
def test_spatially_variable_genes():
    # Explicitly marked as manual-run only
```

**Assertions on AnnData:**
```python
# Shape checks
assert adata.n_obs == expected_n_obs
assert adata.n_vars == expected_n_vars

# Column existence
assert "highly_variable" in adata.var.columns
assert "pct_counts_mt" in adata.obs.columns

# Array values
np.testing.assert_allclose(adata.X.sum(axis=1), 1e4, rtol=1e-5)
assert adata.X.max() <= 10.0 + 1e-6

# Sparse matrix handling
if sparse.issparse(adata.X):
    sums = np.asarray(adata.X.sum(axis=1)).flatten()
else:
    sums = adata.X.sum(axis=1)
```

## Optional Dependency Handling

**Pattern: pytest.importorskip()**
```python
# At module top level
sqlalchemy = pytest.importorskip("sqlalchemy", reason="sqlalchemy not installed")

# In individual tests
def test_jinja_report(self, tmp_path):
    pytest.importorskip("jinja2")
    # ... test code
```

**Pattern: @pytest.mark.skipif()**
```python
from pathlib import Path
_ROBIN_HALLMARK = Path("~/.sc_tools/data/h.all.v7.2.json").expanduser()

@pytest.mark.skipif(not _ROBIN_HALLMARK.exists(), reason="Robin h.all JSON not found")
def test_gene_set_enrichment(self):
    # Only runs if data file exists
```

**Test Result:**
- Skipped tests show as 's' in pytest output
- No failure recorded — optional dependencies don't block CI

## File Paths in Tests

**Import functions inside tests:**
```python
def test_backup_creates_raw(self, adata_counts):
    from sc_tools.pp import backup_raw  # Fresh import per test

    assert adata_counts.raw is None
    backup_raw(adata_counts)
    assert adata_counts.raw is not None
```

**Access via module paths:**
- `sc_tools/pp/normalize.py` — implements `backup_raw()`
- `sc_tools/tests/test_pp.py` — tests via `from sc_tools.pp import backup_raw`
- `sc_tools/qc/metrics.py` — implements `calculate_qc_metrics()`
- `sc_tools/tests/test_qc.py` — tests via `from sc_tools.qc import calculate_qc_metrics`

---

*Testing analysis: 2026-03-20*
