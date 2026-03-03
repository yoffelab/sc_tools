"""
Gene set loaders and curation utilities.

Standardised loading of gene sets into the two-level nested dict format
{category: {name: [genes]}} consumed by score_signature.

Functions
---------
load_hallmark          Load bundled MSigDB Hallmark (50 sets, human).
load_msigdb_json       Load any MSigDB-format JSON file.
load_gmt               Load standard GMT file.
list_gene_sets         List available bundled collections.
validate_gene_signatures Validate a nested signature dict or JSON file.
merge_gene_signatures  Combine multiple signature dicts.
update_gene_symbols    Replace deprecated gene symbols using an alias map.
save_gene_signatures   Write a signature dict to JSON with a datestamp.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd

# Path to bundled data shipped with the package
_DATA_DIR = Path(__file__).resolve().parent.parent / "data"
_HALLMARK_HUMAN = _DATA_DIR / "hallmark_human.json"

_BUNDLED = {
    "hallmark_human": _HALLMARK_HUMAN,
}


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


def load_hallmark(organism: str = "human") -> dict:
    """
    Load bundled MSigDB Hallmark gene sets.

    Returns the 50 Hallmark gene sets as a two-level dict::

        {"Hallmark": {"TNFA_SIGNALING_VIA_NFKB": ["ABCA1", ...], ...}}

    The ``HALLMARK_`` prefix is stripped from set names so that column names
    produced by ``score_signature`` are concise (e.g. ``Hallmark/HYPOXIA``
    instead of ``Hallmark/HALLMARK_HYPOXIA``).

    Parameters
    ----------
    organism : str
        Only ``"human"`` is currently supported.

    Returns
    -------
    dict
        Two-level nested dict ``{category: {name: [genes]}}``.
    """
    if organism != "human":
        raise NotImplementedError(
            f"organism={organism!r} is not yet supported. Only 'human' is available."
        )
    if not _HALLMARK_HUMAN.exists():
        raise FileNotFoundError(
            f"Bundled Hallmark data not found at {_HALLMARK_HUMAN}. "
            "Re-install sc_tools or run the data-bundling script."
        )
    with open(_HALLMARK_HUMAN) as fh:
        raw = json.load(fh)

    # Strip HALLMARK_ prefix for cleaner column names
    sets = {name.removeprefix("HALLMARK_"): genes for name, genes in raw.items()}
    return {"Hallmark": sets}


def load_msigdb_json(
    path: str | Path,
    category_name: str | None = None,
) -> dict:
    """
    Load an MSigDB-format JSON file.

    MSigDB JSON format: a flat dict where each key is a set name and the value
    is either a plain list of gene symbols or an object with a ``geneSymbols``
    key (the richer export format).

    Parameters
    ----------
    path : str or Path
        Path to the MSigDB JSON file.
    category_name : str or None
        Category name for the outer dict key. If None, inferred from the
        filename stem (e.g. ``"h.all.v2025.1.Hs"`` becomes ``"h.all"``).

    Returns
    -------
    dict
        Two-level nested dict ``{category_name: {set_name: [genes]}}``.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"MSigDB JSON not found: {path}")

    if category_name is None:
        # Use first two dot-separated parts of stem as category, fall back to stem
        stem = path.stem
        parts = stem.split(".")
        category_name = ".".join(parts[:2]) if len(parts) >= 2 else stem

    with open(path) as fh:
        raw = json.load(fh)

    if not isinstance(raw, dict):
        raise ValueError(f"Expected a JSON object (dict) at top level in {path}")

    sets: dict[str, list[str]] = {}
    for set_name, value in raw.items():
        if isinstance(value, list):
            sets[set_name] = [g for g in value if isinstance(g, str)]
        elif isinstance(value, dict):
            # Richer MSigDB format: look for geneSymbols key
            genes = value.get("geneSymbols", value.get("genes", []))
            sets[set_name] = [g for g in genes if isinstance(g, str)]
        else:
            continue  # Skip unexpected types

    return {category_name: sets}


def load_gmt(
    path: str | Path,
    category_name: str | None = None,
) -> dict:
    """
    Load a GMT (Gene Matrix Transposed) file.

    GMT format: tab-separated; each line is::

        SET_NAME\\tDESCRIPTION\\tGENE1\\tGENE2\\t...

    Parameters
    ----------
    path : str or Path
        Path to the GMT file.
    category_name : str or None
        Category name for the outer dict key. If None, the file stem is used.

    Returns
    -------
    dict
        Two-level nested dict ``{category_name: {set_name: [genes]}}``.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"GMT file not found: {path}")

    if category_name is None:
        category_name = path.stem

    sets: dict[str, list[str]] = {}
    with open(path) as fh:
        for _line_no, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue  # Need at least name, description, one gene
            set_name = parts[0]
            genes = [g for g in parts[2:] if g]
            sets[set_name] = genes

    return {category_name: sets}


def list_gene_sets() -> list[str]:
    """
    List names of all bundled gene set collections.

    Returns
    -------
    list[str]
        Names of bundled collections (usable as ``organism`` hints).
    """
    return list(_BUNDLED.keys())


# ---------------------------------------------------------------------------
# Curation utilities
# ---------------------------------------------------------------------------


def validate_gene_signatures(
    signatures: dict | str | Path,
    var_names: list[str] | None = None,
    min_genes: int = 3,
) -> pd.DataFrame:
    """
    Validate a nested gene signature dict or JSON file.

    Checks every leaf gene list and reports coverage, duplicates, and empty
    sets. Optionally checks presence against a provided gene universe
    (e.g. ``adata.var_names``).

    Parameters
    ----------
    signatures : dict, str, or Path
        Two-level nested dict or path to a JSON file containing signatures.
    var_names : list[str] or None
        Gene universe (e.g. ``list(adata.var_names)``). If provided, reports
        ``n_present`` and ``pct_coverage``.
    min_genes : int
        Minimum number of genes required per set (flags sets below threshold).

    Returns
    -------
    pd.DataFrame
        One row per signature with columns: ``signature``, ``n_genes``,
        ``n_unique``, ``n_duplicates``, ``n_present`` (if var_names given),
        ``n_missing`` (if var_names given), ``pct_coverage`` (if var_names
        given), ``status``.
    """
    if isinstance(signatures, (str, Path)):
        path = Path(signatures)
        if not path.exists():
            raise FileNotFoundError(f"Signatures file not found: {path}")
        with open(path) as fh:
            signatures = json.load(fh)
    if not isinstance(signatures, dict):
        raise TypeError("signatures must be a dict or path to a JSON file")

    var_set = {g.upper(): g for g in (var_names or [])} if var_names else None

    rows = []
    leaves = _flatten_to_leaves(signatures)
    for col_name, genes in leaves:
        n_genes = len(genes)
        unique_genes = list(dict.fromkeys(genes))
        n_unique = len(unique_genes)
        n_duplicates = n_genes - n_unique

        row: dict = {
            "signature": col_name,
            "n_genes": n_genes,
            "n_unique": n_unique,
            "n_duplicates": n_duplicates,
        }

        if var_set is not None:
            present = [g for g in unique_genes if g.upper() in var_set]
            missing = [g for g in unique_genes if g.upper() not in var_set]
            row["n_present"] = len(present)
            row["n_missing"] = len(missing)
            row["pct_coverage"] = len(present) / n_unique if n_unique > 0 else 0.0

        # Determine status
        if n_genes == 0:
            status = "empty"
        elif n_unique < min_genes:
            status = f"below_min ({n_unique}<{min_genes})"
        elif var_set is not None and row.get("n_present", n_unique) < min_genes:
            status = f"low_coverage ({row['n_present']} present)"
        else:
            status = "ok"
        row["status"] = status
        rows.append(row)

    return pd.DataFrame(rows)


def merge_gene_signatures(*dicts: dict) -> dict:
    """
    Combine multiple two-level signature dicts.

    Later dicts overwrite earlier ones on key collision at both the category
    and set-name level. The special ``_meta`` key is skipped.

    Parameters
    ----------
    *dicts : dict
        Two-level nested dicts to merge.

    Returns
    -------
    dict
        Merged two-level nested dict.

    Examples
    --------
    >>> project = {"Myeloid": {"Macrophage": ["CD68", "CSF1R"]}}
    >>> hallmark = load_hallmark()
    >>> combined = merge_gene_signatures(project, hallmark)
    """
    result: dict = {}
    for d in dicts:
        if not isinstance(d, dict):
            raise TypeError(f"Expected dict, got {type(d)}")
        for category, value in d.items():
            if category == "_meta":
                continue
            if isinstance(value, dict):
                existing = result.setdefault(category, {})
                for name, genes in value.items():
                    if name == "_meta":
                        continue
                    existing[name] = genes
            elif isinstance(value, list):
                # Flat one-level dict: treat category as both category and set name
                result.setdefault("_flat", {})[category] = value
    return result


def update_gene_symbols(
    signatures: dict,
    alias_map: dict[str, str],
) -> dict:
    """
    Replace deprecated or alias gene symbols in a signature dict.

    Does NOT fetch from the internet; the caller provides the alias map
    (e.g. from an HGNC download or a per-project correction list).

    Parameters
    ----------
    signatures : dict
        Two-level nested dict of gene signatures.
    alias_map : dict[str, str]
        Mapping from old symbol to new symbol, e.g. ``{"FAM19A5": "TAFA5"}``.
        Case-sensitive. Symbols not in the map are left unchanged.

    Returns
    -------
    dict
        New nested dict with symbols updated. Original dict is not mutated.
    """
    if not isinstance(signatures, dict):
        raise TypeError("signatures must be a dict")

    def _update_list(genes: list) -> list:
        return [alias_map.get(g, g) for g in genes]

    def _recurse(d: dict) -> dict:
        out = {}
        for k, v in d.items():
            if k == "_meta":
                out[k] = v
            elif isinstance(v, dict):
                out[k] = _recurse(v)
            elif isinstance(v, list):
                out[k] = _update_list(v)
            else:
                out[k] = v
        return out

    return _recurse(signatures)


def save_gene_signatures(signatures: dict, path: str | Path) -> None:
    """
    Write a gene signature dict to JSON with consistent formatting.

    Adds (or updates) a ``_meta`` key at the top level with a datestamp.
    Keys are sorted; indent is 2.

    Parameters
    ----------
    signatures : dict
        Two-level nested dict of gene signatures.
    path : str or Path
        Output JSON path. Parent directories are created if needed.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    out = dict(signatures)  # Shallow copy
    out["_meta"] = {
        "updated": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
    }

    with open(path, "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------


def _flatten_to_leaves(d: dict, prefix: tuple = ()) -> list[tuple[str, list]]:
    """Flatten nested dict to (full_path_string, genes) pairs. Skips _meta."""
    out = []
    for k, v in d.items():
        if k == "_meta":
            continue
        if isinstance(v, dict):
            out.extend(_flatten_to_leaves(v, prefix + (k,)))
        elif isinstance(v, list):
            genes = [g for g in v if isinstance(g, str)]
            col = "/".join(prefix + (k,))
            out.append((col, genes))
    return out
