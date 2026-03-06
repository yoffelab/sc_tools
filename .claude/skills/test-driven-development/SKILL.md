---
name: test-driven-development
tier: 1
description: "[Tier 1 — Core] Write the failing test first, then implement. Use for all new sc_tools functions, bug fixes, and behavior changes. No production code without a failing test first."
allowed-tools: Bash Read Write Edit Glob Grep
---

# Test-Driven Development (TDD)

**Core principle:** If you didn't watch the test fail, you don't know if it tests the right thing.

```
NO PRODUCTION CODE WITHOUT A FAILING TEST FIRST
```

## When to Use

**Always:** new features, bug fixes, refactoring, behavior changes.

**Exceptions (require explicit approval):** throwaway prototypes, generated config files.

## Red-Green-Refactor Cycle

### RED — Write Failing Test
Write one minimal test showing what should happen. One behavior, clear name, real code (no mocks unless unavoidable).

```python
# sc_tools example
def test_filter_spots_removes_low_count_spots():
    adata = make_test_adata(n_obs=100)
    adata.obs['total_counts'] = np.random.randint(0, 1000, 100)
    result = filter_spots(adata, min_counts=500)
    assert (result.obs['total_counts'] >= 500).all()
```

### Verify RED — Watch It Fail (MANDATORY, never skip)
```bash
pytest sc_tools/tests/test_qc.py::test_filter_spots_removes_low_count_spots -v
```
Confirm: test **fails** with expected message (feature missing, not a typo).

### GREEN — Minimal Code
Write the simplest code that passes the test. Do not add features, refactor other code, or "improve" beyond what the test requires.

### Verify GREEN — Watch It Pass (MANDATORY)
```bash
pytest sc_tools/tests/ -v --tb=short -q
```
All tests pass. No warnings.

### REFACTOR — Clean Up
After green only: remove duplication, improve names, extract helpers. Keep tests green. Do not add behavior.

## sc_tools-specific standards

- Test files: `sc_tools/tests/test_{module}.py`
- Fixtures: use `make_test_adata()` and other helpers in `conftest.py`
- Mock external deps (scVI, GPU) with `pytest.importorskip` or `unittest.mock`
- Run full suite before committing: `pytest sc_tools/tests/ -v --tb=short -q`
- Lint before committing: `make lint`

## Red Flags — Stop and Start Over

- Code written before test
- Test passes immediately without implementation
- Cannot explain why the test failed
- "I'll write tests after" rationalization
- Mocks that test the mock, not the code

## Verification Checklist

Before marking work complete:
- [ ] Every new function has a test
- [ ] Watched each test fail before implementing
- [ ] Test failed for the right reason (feature missing)
- [ ] Wrote minimal code to pass each test
- [ ] All tests pass with fresh run
- [ ] Lint clean
