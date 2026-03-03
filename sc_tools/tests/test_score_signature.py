"""
Unit tests for sc_tools.tl.score_signature.

- Minimal: synthetic adata + dict (no project files).
- Robin/ggo (synthetic): project gene_signatures.json + synthetic adata; skip if JSON missing.
- Robin/ggo (real): project adata + project gene_signatures.json; skip if either missing.
  Subset to max_obs to keep runtime low.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from sc_tools.tl import score_signature

# Repo root: assume tests run from repo root (e.g. pytest sc_tools/tests/)
REPO_ROOT = Path(__file__).resolve().parents[2]
ROBIN_JSON = REPO_ROOT / "projects" / "visium_hd" / "robin" / "metadata" / "gene_signatures.json"
ROBIN_P3 = REPO_ROOT / "projects" / "visium_hd" / "robin" / "results" / "adata.normalized.p3.h5ad"
GGO_JSON = REPO_ROOT / "projects" / "visium" / "ggo_visium" / "metadata" / "gene_signatures.json"
GGO_P2 = REPO_ROOT / "projects" / "visium" / "ggo_visium" / "results" / "adata.annotated.p2.h5ad"

# Cap obs for real-adata tests to keep test time reasonable
MAX_OBS_REAL = 500


def _minimal_adata(n_obs: int = 50, n_vars: int = 100):
    """Build minimal adata with var_names G1..Gn_vars and backup in raw (enough for control sampling)."""
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
    var_names = [f"G{i}" for i in range(1, n_vars + 1)]
    adata = sc.AnnData(
        X,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=var_names),
    )
    adata.raw = adata.copy()
    return adata


def test_score_signature_minimal():
    """Minimal test: synthetic adata + two-level dict; no project files."""
    adata = _minimal_adata(n_obs=50, n_vars=100)
    signatures = {"GroupA": {"Sig1": ["G1", "G2", "G3"]}}

    result = score_signature(adata, signatures, use_raw=True, save_path=None, copy=False)

    assert result is adata
    assert "signature_score" in adata.obsm
    assert "signature_score_z" in adata.obsm
    df_raw = adata.obsm["signature_score"]
    df_z = adata.obsm["signature_score_z"]
    assert df_raw.shape[0] == 50 and df_raw.shape[1] == 1
    assert df_z.shape[0] == 50 and df_z.shape[1] == 1
    assert "GroupA/Sig1" in df_raw.columns
    assert "GroupA/Sig1" in df_z.columns
    assert "signature_score_report" in adata.uns
    report = adata.uns["signature_score_report"]
    assert "signature" in report.columns and "n_present" in report.columns
    assert "status" in report.columns
    assert list(report["signature"]) == ["GroupA/Sig1"]
    # No score columns left in obs (temp column removed)
    assert not any(c.startswith("_tmp_sig_") for c in adata.obs.columns)
    assert "GroupA/Sig1" not in adata.obs.columns


def _synthetic_adata_from_genes(gene_list: list[str], n_obs: int = 100, extra_genes: int = 10):
    """Build adata whose var_names include gene_list plus extra dummy genes."""
    all_genes = list(dict.fromkeys(gene_list + [f"DUMMY_{i}" for i in range(extra_genes)]))
    n_vars = len(all_genes)
    np.random.seed(43)
    X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
    adata = sc.AnnData(
        X,
        obs=pd.DataFrame(index=[f"spot_{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=all_genes),
    )
    adata.raw = adata.copy()
    return adata


def _collect_first_genes_from_json(json_path: Path, max_genes: int = 25):
    """Flatten first few signatures and collect up to max_genes."""
    import json

    with open(json_path) as f:
        data = json.load(f)
    genes = []
    for k, v in data.items():
        if k == "_meta" or not isinstance(v, dict):
            continue
        for _sub, val in v.items():
            if isinstance(val, list):
                for g in val:
                    if isinstance(g, str) and g not in genes:
                        genes.append(g)
                        if len(genes) >= max_genes:
                            return genes
    return genes


@pytest.mark.skipif(not ROBIN_JSON.exists(), reason="Robin gene_signatures.json not found")
def test_score_signature_robin_signatures():
    """Run score_signature with robin metadata/gene_signatures.json."""
    genes = _collect_first_genes_from_json(ROBIN_JSON, max_genes=25)
    assert len(genes) >= 3, "Need at least 3 genes from robin JSON"
    adata = _synthetic_adata_from_genes(genes, n_obs=100)

    score_signature(adata, str(ROBIN_JSON), use_raw=True, save_path=None, copy=False)

    assert "signature_score" in adata.obsm
    assert "signature_score_z" in adata.obsm
    assert adata.obsm["signature_score"].shape[0] == 100
    assert adata.obsm["signature_score"].shape[1] >= 1
    assert all("/" in c for c in adata.obsm["signature_score"].columns)
    assert "signature_score_report" in adata.uns
    assert not any(c.startswith("_tmp_sig_") for c in adata.obs.columns)


@pytest.mark.skipif(not GGO_JSON.exists(), reason="ggo_visium gene_signatures.json not found")
def test_score_signature_ggo_visium_signatures():
    """Run score_signature with ggo_visium metadata/gene_signatures.json (synthetic adata)."""
    genes = _collect_first_genes_from_json(GGO_JSON, max_genes=25)
    assert len(genes) >= 3, "Need at least 3 genes from ggo_visium JSON"
    adata = _synthetic_adata_from_genes(genes, n_obs=100)

    score_signature(adata, str(GGO_JSON), use_raw=True, save_path=None, copy=False)

    assert "signature_score" in adata.obsm
    assert "signature_score_z" in adata.obsm
    assert adata.obsm["signature_score"].shape[0] == 100
    assert adata.obsm["signature_score"].shape[1] >= 1
    assert all("/" in c for c in adata.obsm["signature_score"].columns)
    assert "signature_score_report" in adata.uns
    assert not any(c.startswith("_tmp_sig_") for c in adata.obs.columns)


def _load_and_subset_adata(path: Path, max_obs: int = MAX_OBS_REAL):
    """Load adata from path and subset to at most max_obs to keep test fast."""
    adata = sc.read_h5ad(path)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    if adata.n_obs > max_obs:
        adata = adata[:max_obs].copy()
    return adata


@pytest.mark.skipif(
    not (ROBIN_P3.exists() and ROBIN_JSON.exists()),
    reason="Robin adata.normalized.p3.h5ad and/or gene_signatures.json not found",
)
def test_score_signature_robin_real():
    """Run score_signature on real robin p3 adata and metadata/gene_signatures.json."""
    adata = _load_and_subset_adata(ROBIN_P3)
    n_obs = adata.n_obs

    score_signature(adata, str(ROBIN_JSON), use_raw=True, save_path=None, copy=False)

    assert "signature_score" in adata.obsm
    assert "signature_score_z" in adata.obsm
    assert adata.obsm["signature_score"].shape[0] == n_obs
    assert adata.obsm["signature_score"].shape[1] >= 1
    assert all("/" in c for c in adata.obsm["signature_score"].columns)
    assert "signature_score_report" in adata.uns
    report = adata.uns["signature_score_report"]
    assert len(report) >= 1
    assert not any(c.startswith("_tmp_sig_") for c in adata.obs.columns)


@pytest.mark.skipif(
    not (GGO_P2.exists() and GGO_JSON.exists()),
    reason="ggo_visium adata.annotated.p2.h5ad and/or gene_signatures.json not found",
)
def test_score_signature_ggo_visium_real():
    """Run score_signature on real ggo_visium p2 adata and metadata/gene_signatures.json."""
    adata = _load_and_subset_adata(GGO_P2)
    n_obs = adata.n_obs

    score_signature(adata, str(GGO_JSON), use_raw=True, save_path=None, copy=False)

    assert "signature_score" in adata.obsm
    assert "signature_score_z" in adata.obsm
    assert adata.obsm["signature_score"].shape[0] == n_obs
    assert adata.obsm["signature_score"].shape[1] >= 1
    assert all("/" in c for c in adata.obsm["signature_score"].columns)
    assert "signature_score_report" in adata.uns
    report = adata.uns["signature_score_report"]
    assert len(report) >= 1
    assert not any(c.startswith("_tmp_sig_") for c in adata.obs.columns)


# ---------------------------------------------------------------------------
# Tests for method= parameter
# ---------------------------------------------------------------------------


def test_score_signature_method_scanpy_explicit():
    """method='scanpy' is the default; verify uns["scoring_method"] is recorded."""
    adata = _minimal_adata(n_obs=50, n_vars=100)
    signatures = {"GroupA": {"Sig1": ["G1", "G2", "G3"]}}
    score_signature(adata, signatures, method="scanpy", use_raw=True, copy=False)
    assert adata.uns.get("scoring_method") == "scanpy"


def test_score_signature_invalid_method():
    adata = _minimal_adata(n_obs=20, n_vars=50)
    with pytest.raises(ValueError, match="method="):
        score_signature(adata, {"A": {"Sig": ["G1", "G2", "G3"]}}, method="invalid_method")


def test_score_signature_method_ucell():
    """method='ucell': scores written to obsm; scoring_method recorded."""
    pytest.importorskip("pyucell")
    adata = _minimal_adata(n_obs=50, n_vars=100)
    # Remove raw since ucell does not use it
    signatures = {"GroupA": {"Sig1": ["G1", "G2", "G3", "G4", "G5"]}}
    score_signature(adata, signatures, method="ucell", copy=False)
    assert "signature_score" in adata.obsm
    assert "signature_score_z" in adata.obsm
    assert adata.uns.get("scoring_method") == "ucell"
    df = adata.obsm["signature_score"]
    assert df.shape == (50, 1)
    assert "GroupA/Sig1" in df.columns


def test_score_signature_method_ssgsea():
    """method='ssgsea': scores written to obsm; scoring_method recorded."""
    pytest.importorskip("gseapy")
    adata = _minimal_adata(n_obs=20, n_vars=100)
    signatures = {"GroupA": {"Sig1": ["G1", "G2", "G3", "G4", "G5"]}}
    score_signature(adata, signatures, method="ssgsea", use_raw=False, copy=False)
    assert "signature_score" in adata.obsm
    assert "signature_score_z" in adata.obsm
    assert adata.uns.get("scoring_method") == "ssgsea"


def test_score_signature_ucell_import_error(monkeypatch):
    """method='ucell' raises ImportError with install hint when pyucell is missing."""
    import sys

    monkeypatch.setitem(sys.modules, "pyucell", None)
    adata = _minimal_adata(n_obs=20, n_vars=50)
    with pytest.raises((ImportError, TypeError)):
        score_signature(adata, {"A": {"Sig": ["G1", "G2", "G3"]}}, method="ucell")


def test_score_signature_ssgsea_import_error(monkeypatch):
    """method='ssgsea' raises ImportError with install hint when gseapy is missing."""
    import sys

    monkeypatch.setitem(sys.modules, "gseapy", None)
    adata = _minimal_adata(n_obs=20, n_vars=50)
    with pytest.raises((ImportError, TypeError)):
        score_signature(adata, {"A": {"Sig": ["G1", "G2", "G3"]}}, method="ssgsea")
