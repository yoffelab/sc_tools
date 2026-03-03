"""
Unit tests for sc_tools.tl.gene_sets.

Tests all loaders and curation utilities:
- load_hallmark: bundled data
- load_msigdb_json: MSigDB flat JSON
- load_gmt: GMT file format
- list_gene_sets: bundled collection names
- validate_gene_signatures: validation report
- merge_gene_signatures: combining dicts
- update_gene_symbols: alias replacement
- save_gene_signatures: JSON roundtrip
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pandas as pd
import pytest

from sc_tools.tl.gene_sets import (
    list_gene_sets,
    load_gmt,
    load_hallmark,
    load_msigdb_json,
    merge_gene_signatures,
    save_gene_signatures,
    update_gene_symbols,
    validate_gene_signatures,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
ROBIN_HALLMARK = (
    REPO_ROOT / "projects" / "visium_hd" / "robin" / "metadata" / "h.all.v2025.1.Hs.json"
)


# ---------------------------------------------------------------------------
# load_hallmark
# ---------------------------------------------------------------------------


def test_load_hallmark_returns_two_level_dict():
    result = load_hallmark()
    assert isinstance(result, dict)
    assert "Hallmark" in result
    inner = result["Hallmark"]
    assert isinstance(inner, dict)
    assert len(inner) == 50


def test_load_hallmark_strips_prefix():
    result = load_hallmark()
    inner = result["Hallmark"]
    # No key should start with HALLMARK_
    for k in inner:
        assert not k.startswith("HALLMARK_"), f"Key still has prefix: {k}"


def test_load_hallmark_gene_lists_are_nonempty():
    result = load_hallmark()
    for name, genes in result["Hallmark"].items():
        assert isinstance(genes, list), f"{name}: expected list"
        assert len(genes) >= 1, f"{name}: empty gene list"


def test_load_hallmark_unsupported_organism():
    with pytest.raises(NotImplementedError, match="mouse"):
        load_hallmark(organism="mouse")


def test_load_hallmark_known_set_present():
    result = load_hallmark()
    assert "TNFA_SIGNALING_VIA_NFKB" in result["Hallmark"]
    assert "HYPOXIA" in result["Hallmark"]


# ---------------------------------------------------------------------------
# load_msigdb_json
# ---------------------------------------------------------------------------


def test_load_msigdb_json_richer_format(tmp_path):
    """Test loading MSigDB richer format (geneSymbols key)."""
    data = {
        "MY_SET_A": {"geneSymbols": ["GeneA", "GeneB"], "collection": "test"},
        "MY_SET_B": {"geneSymbols": ["GeneC"]},
    }
    fp = tmp_path / "my.collection.v1.json"
    fp.write_text(json.dumps(data))

    result = load_msigdb_json(fp)
    assert "my.collection" in result
    inner = result["my.collection"]
    assert inner["MY_SET_A"] == ["GeneA", "GeneB"]
    assert inner["MY_SET_B"] == ["GeneC"]


def test_load_msigdb_json_flat_format(tmp_path):
    """Test loading flat MSigDB format (list values)."""
    data = {"SET_X": ["G1", "G2", "G3"]}
    fp = tmp_path / "custom.json"
    fp.write_text(json.dumps(data))

    result = load_msigdb_json(fp, category_name="Custom")
    assert "Custom" in result
    assert result["Custom"]["SET_X"] == ["G1", "G2", "G3"]


def test_load_msigdb_json_explicit_category(tmp_path):
    data = {"SET_A": ["G1", "G2"]}
    fp = tmp_path / "data.json"
    fp.write_text(json.dumps(data))

    result = load_msigdb_json(fp, category_name="MyCategory")
    assert "MyCategory" in result


def test_load_msigdb_json_missing_file():
    with pytest.raises(FileNotFoundError):
        load_msigdb_json("/nonexistent/path.json")


@pytest.mark.skipif(not ROBIN_HALLMARK.exists(), reason="Robin h.all JSON not found")
def test_load_msigdb_json_robin_hallmark():
    result = load_msigdb_json(ROBIN_HALLMARK, category_name="Hallmark")
    assert "Hallmark" in result
    assert len(result["Hallmark"]) == 50


# ---------------------------------------------------------------------------
# load_gmt
# ---------------------------------------------------------------------------


def test_load_gmt_basic(tmp_path):
    content = textwrap.dedent("""\
        SET_A\tdescription_a\tGENE1\tGENE2\tGENE3
        SET_B\tdescription_b\tGENEX\tGENEY
    """)
    fp = tmp_path / "test.gmt"
    fp.write_text(content)

    result = load_gmt(fp)
    assert "test" in result
    inner = result["test"]
    assert inner["SET_A"] == ["GENE1", "GENE2", "GENE3"]
    assert inner["SET_B"] == ["GENEX", "GENEY"]


def test_load_gmt_explicit_category(tmp_path):
    content = "SET_A\tdesc\tGENE1\tGENE2\n"
    fp = tmp_path / "my_gmt.gmt"
    fp.write_text(content)

    result = load_gmt(fp, category_name="MyGMT")
    assert "MyGMT" in result


def test_load_gmt_empty_lines_skipped(tmp_path):
    content = "SET_A\tdesc\tGENE1\tGENE2\n\n\nSET_B\tdesc\tGENE3\n"
    fp = tmp_path / "test.gmt"
    fp.write_text(content)

    result = load_gmt(fp)
    assert len(result["test"]) == 2


def test_load_gmt_missing_file():
    with pytest.raises(FileNotFoundError):
        load_gmt("/nonexistent/file.gmt")


# ---------------------------------------------------------------------------
# list_gene_sets
# ---------------------------------------------------------------------------


def test_list_gene_sets_returns_list():
    result = list_gene_sets()
    assert isinstance(result, list)
    assert "hallmark_human" in result


# ---------------------------------------------------------------------------
# validate_gene_signatures
# ---------------------------------------------------------------------------


def _sample_signatures():
    return {
        "CatA": {
            "SigGood": ["GENE1", "GENE2", "GENE3", "GENE4"],
            "SigEmpty": [],
            "SigSmall": ["ONLY_ONE"],
            "SigDup": ["G1", "G1", "G2", "G3"],
        }
    }


def test_validate_basic_structure():
    sigs = _sample_signatures()
    df = validate_gene_signatures(sigs)
    assert isinstance(df, pd.DataFrame)
    assert {"signature", "n_genes", "n_unique", "n_duplicates", "status"}.issubset(df.columns)
    assert len(df) == 4


def test_validate_detects_empty():
    sigs = {"Cat": {"Empty": []}}
    df = validate_gene_signatures(sigs)
    row = df[df["signature"] == "Cat/Empty"].iloc[0]
    assert row["status"] == "empty"


def test_validate_detects_below_min():
    sigs = {"Cat": {"Small": ["G1", "G2"]}}
    df = validate_gene_signatures(sigs, min_genes=3)
    row = df[df["signature"] == "Cat/Small"].iloc[0]
    assert "below_min" in row["status"]


def test_validate_detects_duplicates():
    sigs = {"Cat": {"Dup": ["G1", "G1", "G2", "G3", "G4"]}}
    df = validate_gene_signatures(sigs)
    row = df[df["signature"] == "Cat/Dup"].iloc[0]
    assert row["n_duplicates"] == 1
    assert row["n_unique"] == 4


def test_validate_with_var_names():
    sigs = {"Cat": {"Sig": ["GENE1", "GENE2", "GENE3", "GENE4", "UNKNOWN"]}}
    var_names = ["GENE1", "GENE2", "GENE3", "GENE4"]
    df = validate_gene_signatures(sigs, var_names=var_names)
    assert "n_present" in df.columns
    row = df.iloc[0]
    assert row["n_present"] == 4
    assert row["n_missing"] == 1
    assert abs(row["pct_coverage"] - 0.8) < 1e-6


def test_validate_from_json_path(tmp_path):
    sigs = {"Cat": {"Sig": ["G1", "G2", "G3"]}}
    fp = tmp_path / "sigs.json"
    fp.write_text(json.dumps(sigs))
    df = validate_gene_signatures(fp)
    assert len(df) == 1
    assert df.iloc[0]["status"] == "ok"


def test_validate_skips_meta_key():
    sigs = {"_meta": {"updated": "2025-01-01"}, "Cat": {"Sig": ["G1", "G2", "G3"]}}
    df = validate_gene_signatures(sigs)
    # _meta should be skipped
    assert all("_meta" not in s for s in df["signature"])


# ---------------------------------------------------------------------------
# merge_gene_signatures
# ---------------------------------------------------------------------------


def test_merge_two_dicts():
    d1 = {"CatA": {"SigA": ["G1", "G2"]}}
    d2 = {"CatB": {"SigB": ["G3", "G4"]}}
    result = merge_gene_signatures(d1, d2)
    assert "CatA" in result
    assert "CatB" in result
    assert result["CatA"]["SigA"] == ["G1", "G2"]
    assert result["CatB"]["SigB"] == ["G3", "G4"]


def test_merge_overwrites_on_collision():
    d1 = {"CatA": {"Sig": ["OLD"]}}
    d2 = {"CatA": {"Sig": ["NEW"]}}
    result = merge_gene_signatures(d1, d2)
    assert result["CatA"]["Sig"] == ["NEW"]


def test_merge_skips_meta():
    d1 = {"_meta": {"ts": "2025"}, "CatA": {"Sig": ["G1"]}}
    d2 = {"CatA": {"Sig2": ["G2"]}}
    result = merge_gene_signatures(d1, d2)
    assert "_meta" not in result
    assert "Sig" in result["CatA"]
    assert "Sig2" in result["CatA"]


def test_merge_three_dicts():
    d1 = {"A": {"s1": ["G1"]}}
    d2 = {"B": {"s2": ["G2"]}}
    d3 = {"C": {"s3": ["G3"]}}
    result = merge_gene_signatures(d1, d2, d3)
    assert set(result.keys()) == {"A", "B", "C"}


def test_merge_with_hallmark():
    project = {"Myeloid": {"Macrophage": ["CD68", "CSF1R", "CD163"]}}
    hallmark = load_hallmark()
    combined = merge_gene_signatures(project, hallmark)
    assert "Myeloid" in combined
    assert "Hallmark" in combined
    assert len(combined["Hallmark"]) == 50


# ---------------------------------------------------------------------------
# update_gene_symbols
# ---------------------------------------------------------------------------


def test_update_gene_symbols_basic():
    sigs = {"Cat": {"Sig": ["OLD_SYMBOL", "VALID_GENE"]}}
    alias_map = {"OLD_SYMBOL": "NEW_SYMBOL"}
    result = update_gene_symbols(sigs, alias_map)
    assert result["Cat"]["Sig"] == ["NEW_SYMBOL", "VALID_GENE"]


def test_update_gene_symbols_no_mutation():
    sigs = {"Cat": {"Sig": ["OLD"]}}
    alias_map = {"OLD": "NEW"}
    _ = update_gene_symbols(sigs, alias_map)
    # Original unchanged
    assert sigs["Cat"]["Sig"] == ["OLD"]


def test_update_gene_symbols_nested():
    sigs = {"A": {"B": {"Sig": ["OLD", "KEEP"]}}}
    alias_map = {"OLD": "REPLACED"}
    result = update_gene_symbols(sigs, alias_map)
    assert result["A"]["B"]["Sig"] == ["REPLACED", "KEEP"]


def test_update_gene_symbols_preserves_meta():
    sigs = {"_meta": {"ts": "2025"}, "Cat": {"Sig": ["G1"]}}
    alias_map = {"G1": "G1_new"}
    result = update_gene_symbols(sigs, alias_map)
    assert result["_meta"] == {"ts": "2025"}
    assert result["Cat"]["Sig"] == ["G1_new"]


# ---------------------------------------------------------------------------
# save_gene_signatures
# ---------------------------------------------------------------------------


def test_save_and_reload(tmp_path):
    sigs = {"CatA": {"Sig1": ["G1", "G2", "G3"]}}
    out = tmp_path / "saved.json"
    save_gene_signatures(sigs, out)

    assert out.exists()
    with open(out) as fh:
        loaded = json.load(fh)

    assert "CatA" in loaded
    assert loaded["CatA"]["Sig1"] == ["G1", "G2", "G3"]
    # _meta datestamp added
    assert "_meta" in loaded
    assert "updated" in loaded["_meta"]


def test_save_creates_parent_dirs(tmp_path):
    sigs = {"Cat": {"Sig": ["G1"]}}
    out = tmp_path / "subdir" / "nested" / "sigs.json"
    save_gene_signatures(sigs, out)
    assert out.exists()
