"""
Shared utilities for QC HTML report generation.

Provides figure encoding, template rendering, date stamping, metric table
builders, embedding auto-detection, and optional segmentation/integration
section computation.
"""

from __future__ import annotations

import base64
import io
import logging
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

__all__ = [
    "fig_to_base64",
    "plotly_to_html",
    "render_template",
    "get_date_stamp",
    "get_modality_terms",
    "build_metrics_table_rows",
    "auto_detect_embeddings",
    "compute_segmentation_section",
    "compute_integration_section",
    "_find_latest_report",
    "_extract_body_content",
    "_extract_head_css",
    "_wrap_with_tabs",
]

logger = logging.getLogger(__name__)

# Modalities that use protein markers instead of gene expression
_PROTEIN_MODALITIES = {"imc"}

# Modalities that produce single-cell (not spot) observations
_CELL_MODALITIES = {
    "imc",
    "xenium",
    "cosmx",
    "cosmx_1k",
    "cosmx_6k",
    "cosmx_full_library",
    "visium_hd_cell",
}


def get_modality_terms(modality: str) -> dict[str, str | bool]:
    """Return display-label terminology adapted to the data modality.

    Internal column names (``n_genes_by_counts``, ``total_counts``) are
    unchanged -- only **display labels** in plots and HTML templates differ.

    Parameters
    ----------
    modality
        One of ``visium``, ``visium_hd``, ``xenium``, ``cosmx``, ``imc``, etc.

    Returns
    -------
    dict
        Keys: ``feature``, ``features``, ``feature_lower``, ``features_lower``,
        ``observation``, ``observations``, ``observation_lower``,
        ``observations_lower``, ``intensity_label``, ``feature_count_label``,
        ``has_mt``.
    """
    is_protein = modality in _PROTEIN_MODALITIES
    is_cell = modality in _CELL_MODALITIES
    return {
        "feature": "Protein" if is_protein else "Gene",
        "features": "Proteins" if is_protein else "Genes",
        "feature_lower": "protein" if is_protein else "gene",
        "features_lower": "proteins" if is_protein else "genes",
        "observation": "Cell" if is_cell else "Spot",
        "observations": "Cells" if is_cell else "Spots",
        "observation_lower": "cell" if is_cell else "spot",
        "observations_lower": "cells" if is_cell else "spots",
        "intensity_label": "Total intensity" if is_protein else "Total counts",
        "feature_count_label": "Proteins detected" if is_protein else "Genes detected",
        "has_mt": not is_protein,
    }


# Known integration embedding keys in priority order
_KNOWN_EMBEDDINGS = {
    "scVI": "X_scVI",
    "scVI (raw)": "X_scVI_raw",
    "scANVI": "X_scANVI",
    "scANVI (raw)": "X_scANVI_raw",
    "Harmony": "X_pca_harmony",
    "CytoVI": "X_cytovi",
    "CytoVI (arcsinh)": "X_cytovi",
    "BBKNN": "X_umap_bbknn",
    "ComBat": "X_pca_combat",
    "Scanorama": "X_scanorama",
    "DestVI": "X_destvi",
    "IMC Phenotyping": "X_pca_imc_phenotyping",
    "IMC Pheno + Harmony": "X_imc_pheno_harmony",
    "Z-score + Harmony": "X_zscore_harmony",
    "Unintegrated (PCA)": "X_pca",
}


def fig_to_base64(fig: plt.Figure, dpi: int = 150) -> str:
    """Render a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def plotly_to_html(fig) -> str:
    """Convert a plotly figure to an embeddable HTML div string.

    Parameters
    ----------
    fig
        A ``plotly.graph_objects.Figure``.

    Returns
    -------
    str
        HTML ``<div>`` with embedded plotly chart (no full page wrapper).
    """
    return fig.to_html(full_html=False, include_plotlyjs=False)


def render_template(
    template_name: str | Path,
    context: dict,
    assets_dir: Path | None = None,
) -> str:
    """Render a Jinja2 template by name from the assets directory.

    Uses ``Environment(loader=FileSystemLoader(...))`` so that
    ``{% extends %}`` and ``{% include %}`` directives work correctly.

    Parameters
    ----------
    template_name
        Template filename (e.g. ``"pre_filter_qc_template.html"``) looked up
        inside *assets_dir*, **or** an absolute/relative ``Path`` object for
        backward compatibility (the parent becomes *assets_dir* and the stem
        becomes *template_name*).
    context
        Dict of variables passed to ``Template.render()``.
    assets_dir
        Directory containing Jinja2 templates.  Defaults to
        ``<package_root>/assets/``.  Ignored when *template_name* is a
        ``Path`` object with a parent directory.

    Returns
    -------
    str
        Rendered HTML string.
    """
    try:
        from jinja2 import Environment, FileSystemLoader
    except ImportError as e:
        raise ImportError(
            "jinja2 is required for HTML report generation. Install with: pip install jinja2"
        ) from e

    # Backward-compat: if caller passes a Path, extract dir + filename.
    p = Path(template_name)
    if p.is_absolute() or (p.parent != Path(".")):
        # Caller passed a full/relative path — honour it directly.
        assets_dir = p.parent
        template_name = p.name
    else:
        template_name = str(template_name)

    if assets_dir is None:
        assets_dir = Path(__file__).parent.parent / "assets"

    env = Environment(loader=FileSystemLoader(str(assets_dir)), autoescape=False)
    template = env.get_template(template_name)
    return template.render(**context)


def get_date_stamp(date_stamp: str | None = None) -> str:
    """Return a YYYYMMDD date stamp (today if *date_stamp* is ``None``)."""
    if date_stamp is not None:
        return date_stamp
    return datetime.now().strftime("%Y%m%d")


def build_metrics_table_rows(classified: pd.DataFrame) -> list[dict]:
    """Build per-sample table rows from a ``classify_samples`` DataFrame.

    Parameters
    ----------
    classified
        Output of ``sc_tools.qc.classify_samples``.

    Returns
    -------
    list[dict]
        Each dict has keys: ``sample``, ``qc_pass``, ``qc_fail_reasons``,
        ``n_spots``, ``total_counts_median``, ``n_genes_median``,
        ``n_genes_detected``, ``pct_mt_median``, ``pct_mt_gt5``.
    """
    table_rows: list[dict] = []
    for sample in classified.index:
        row: dict = {"sample": str(sample)}
        row["qc_pass"] = bool(classified.loc[sample, "qc_pass"])
        row["qc_fail_reasons"] = (
            str(classified.loc[sample, "qc_fail_reasons"])
            if classified.loc[sample, "qc_fail_reasons"]
            else ""
        )
        for col in [
            "n_spots",
            "total_counts_median",
            "n_genes_median",
            "n_genes_detected",
            "pct_mt_median",
            "pct_mt_gt5",
        ]:
            val = classified.loc[sample, col] if col in classified.columns else None
            if val is not None and not (isinstance(val, float) and np.isnan(val)):
                row[col] = val
            else:
                row[col] = None
        table_rows.append(row)
    return table_rows


def auto_detect_embeddings(adata: AnnData) -> dict[str, str]:
    """Scan ``adata.obsm`` for known integration embedding keys.

    Returns
    -------
    dict[str, str]
        Mapping of human-readable name to ``obsm`` key, for keys that exist.
    """
    found: dict[str, str] = {}
    for name, key in _KNOWN_EMBEDDINGS.items():
        if key in adata.obsm:
            found[name] = key
    return found


def compute_segmentation_section(
    adata: AnnData,
    masks_dir: str | Path,
    sample_col: str = "library_id",
) -> dict | None:
    """Load masks and compute per-ROI segmentation scores.

    Parameters
    ----------
    adata
        AnnData with sample annotations.
    masks_dir
        Directory containing mask TIFF files (one per ROI/sample).
    sample_col
        Column in ``adata.obs`` identifying samples.

    Returns
    -------
    dict or None
        Dict with ``"df"`` (DataFrame of per-ROI scores) and ``"plots"``
        (dict of plot name to base64 string). Returns ``None`` if no masks
        are found.
    """
    masks_dir = Path(masks_dir)
    if not masks_dir.exists():
        logger.info("Segmentation masks dir not found: %s", masks_dir)
        return None

    mask_files = sorted(masks_dir.glob("*.tif")) + sorted(masks_dir.glob("*.tiff"))
    if not mask_files:
        logger.info("No mask TIFFs found in %s", masks_dir)
        return None

    try:
        from sc_tools.bm.mask_io import load_mask
        from sc_tools.bm.segmentation import score_segmentation
    except ImportError:
        logger.warning("sc_tools.bm not available; skipping segmentation scoring")
        return None

    rows = []
    for mf in mask_files:
        try:
            mask = load_mask(str(mf))
            scores = score_segmentation(mask)
            scores["roi"] = mf.stem
            rows.append(scores)
        except Exception:
            logger.warning("Failed to score mask %s", mf.name, exc_info=True)
    if not rows:
        return None

    df = pd.DataFrame(rows).set_index("roi")

    # Composite score bar chart
    plots: dict[str, str] = {}
    if "composite_score" in df.columns:
        fig, ax = plt.subplots(figsize=(max(6, len(df) * 0.5), 4))
        ax.barh(df.index, df["composite_score"], color="#3498db")
        ax.set_xlabel("Composite Score")
        ax.set_title("Segmentation Quality per ROI")
        fig.tight_layout()
        plots["seg_composite_bar"] = fig_to_base64(fig)

    return {"df": df, "plots": plots}


def compute_integration_section(
    adata: AnnData,
    embedding_keys: dict[str, str],
    batch_key: str,
    celltype_key: str | None = None,
    comparison_df: pd.DataFrame | None = None,
) -> dict | None:
    """Compute integration comparison and generate plots.

    Parameters
    ----------
    adata
        AnnData with embeddings.
    embedding_keys
        Mapping of method name to ``obsm`` key.
    batch_key
        Column in ``adata.obs`` with batch labels.
    celltype_key
        Column in ``adata.obs`` with cell type labels. ``None`` to skip bio
        metrics.
    comparison_df
        Pre-computed benchmark DataFrame. When provided, skips calling
        ``compare_integrations()`` and uses this directly for plotting.

    Returns
    -------
    dict or None
        Dict with ``"comparison_df"`` and ``"plots"`` (dict of name to HTML
        or base64 strings). Returns ``None`` on failure.
    """
    if not embedding_keys:
        logger.info("No embeddings provided for integration comparison")
        return None

    if comparison_df is None:
        try:
            from sc_tools.bm.integration import compare_integrations

            comparison_df = compare_integrations(
                adata,
                embedding_keys,
                batch_key=batch_key,
                celltype_key=celltype_key,
                include_unintegrated=True,
                use_scib="sklearn",
            )
        except ImportError:
            logger.warning("sc_tools.bm.integration not available")
            return None
        except Exception:
            logger.warning("Integration comparison failed", exc_info=True)
            return None

    if comparison_df.empty:
        return None

    # Preserve scib_fallback flag from DataFrame.attrs (set by compare_integrations).
    scib_fallback: bool = bool(comparison_df.attrs.get("scib_fallback", False))

    plots: dict[str, str] = {}

    # Radar chart
    try:
        from sc_tools.pl.benchmarking import plot_integration_radar

        fig_radar = plot_integration_radar(comparison_df)
        plots["integration_radar"] = plotly_to_html(fig_radar)
    except ImportError:
        logger.warning("Integration benchmarking plot modules not available")
    except Exception:
        logger.debug("Integration radar plot failed", exc_info=True)

    # Batch vs bio scatter (supplementary)
    if "batch_score" in comparison_df.columns and "bio_score" in comparison_df.columns:
        try:
            from sc_tools.pl.benchmarking import plot_batch_vs_bio

            fig_bvb = plot_batch_vs_bio(comparison_df)
            plots["batch_vs_bio"] = plotly_to_html(fig_bvb)
        except ImportError:
            pass  # already warned above
        except Exception:
            logger.debug("Batch vs bio plot failed", exc_info=True)

    return {"comparison_df": comparison_df, "plots": plots, "scib_fallback": scib_fallback}


# ---------------------------------------------------------------------------
# Tab-navigation helpers (Plan B)
# ---------------------------------------------------------------------------


def _find_latest_report(output_dir: Path, report_type: str) -> Path | None:
    """Find the most recent ``{report_type}_qc_YYYYMMDD.html`` in *output_dir*.

    Globs for ``{report_type}_qc_????????.html`` and returns the file with
    the lexicographically largest 8-digit date suffix. Returns ``None`` when
    no matching files exist.

    Parameters
    ----------
    output_dir
        Directory to search.
    report_type
        Report type prefix, e.g. ``"pre_filter"``, ``"post_filter"``.

    Returns
    -------
    Path or None
    """
    candidates = sorted(output_dir.glob(f"{report_type}_qc_????????.html"))
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stem[-8:])


def _extract_body_content(html_str: str) -> str:
    """Extract inner content between ``<body>`` and ``</body>``.

    Falls back to returning the full string when no ``<body>`` tag is found.

    Parameters
    ----------
    html_str
        Full HTML string.

    Returns
    -------
    str
        Content inside the ``<body>`` element.
    """
    import re

    match = re.search(r"<body[^>]*>(.*)</body>", html_str, re.DOTALL | re.IGNORECASE)
    return match.group(1) if match else html_str


def _extract_head_css(html_str: str) -> str:
    """Extract the CSS text from the first ``<style>`` block in ``<head>``.

    Returns an empty string when no ``<style>`` element is found.

    Parameters
    ----------
    html_str
        Full HTML string.

    Returns
    -------
    str
        Raw CSS text, or ``""`` if absent.
    """
    import re

    match = re.search(r"<style[^>]*>(.*?)</style>", html_str, re.DOTALL | re.IGNORECASE)
    return match.group(1) if match else ""


def _wrap_with_tabs(
    current_label: str,
    current_content: str,
    previous_tabs: list[tuple[str, str, str]],
    current_css: str,
    current_head_extras: str = "",
    title: str = "",
) -> str:
    """Assemble a complete tabbed HTML document.

    The active tab is always the *current* report; previous reports are shown
    in reverse chronological order (most recent first).

    Parameters
    ----------
    current_label
        Display label for the active tab, e.g. ``"Report 2: Post-filter"``.
    current_content
        Inner HTML body of the active report (between ``<body>`` tags).
    previous_tabs
        List of ``(label, body_html, head_css)`` tuples for older reports, in
        the order they should appear (leftmost = first in list).
    current_css
        CSS text extracted from the current report's ``<style>`` block.
    current_head_extras
        Extra ``<head>`` content (e.g. CDN ``<script>`` tags) needed by the
        current report. Plotly CDN is always included when this string contains
        ``plotly`` or when any previous tab content references Plotly.
    title
        ``<title>`` tag text (defaults to *current_label*).

    Returns
    -------
    str
        A complete, self-contained HTML document with tab navigation.
    """
    if not title:
        title = current_label

    # Always include Plotly CDN (some reports embed it, unified wrapper needs it)
    plotly_cdn = '<script src="https://cdn.plot.ly/plotly-3.4.0.min.js"></script>'

    # Build list of all tabs: current first, then previous in order
    all_tabs: list[tuple[str, str, str]] = [(current_label, current_content, current_css)]
    all_tabs.extend(previous_tabs)

    # Assign stable HTML IDs
    def _tab_id(idx: int) -> str:
        return "tab-current" if idx == 0 else f"tab-prev-{idx}"

    btn_id = lambda idx: "btn-" + _tab_id(idx)  # noqa: E731

    # --- Tab navigation bar ---
    btn_html_parts: list[str] = []
    for i, (label, _, _) in enumerate(all_tabs):
        active_cls = " active" if i == 0 else ""
        btn_html_parts.append(
            f'<button class="tab-btn{active_cls}" id="{btn_id(i)}" '
            f"onclick=\"showTab('{_tab_id(i)}')\">{label}</button>"
        )
    nav_html = "\n    ".join(btn_html_parts)

    # --- Panel divs ---
    panel_parts: list[str] = []
    for i, (_, body_html, css_text) in enumerate(all_tabs):
        display = "block" if i == 0 else "none"
        style_block = f"<style>{css_text}</style>" if css_text.strip() else ""
        panel_parts.append(
            f'<div class="tab-panel" id="{_tab_id(i)}" style="display:{display};">\n'
            f"{style_block}\n"
            f"{body_html}\n"
            f"</div>"
        )
    panels_html = "\n".join(panel_parts)

    bootstrap_css = '<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootswatch@5.3.3/dist/flatly/bootstrap.min.css">'
    bootstrap_js = '<script src="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.8/js/bootstrap.bundle.min.js"></script>'

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
{bootstrap_css}
{plotly_cdn}
{current_head_extras}
<style>
  *, *::before, *::after {{ box-sizing: border-box; }}
  body {{ margin: 0; padding: 0; }}
  .tab-nav {{
    position: sticky; top: 0; z-index: 1000;
    background: #1a1a2e; padding: 0 8px;
    display: flex; flex-wrap: wrap; gap: 2px;
  }}
  .tab-btn {{
    background: transparent; border: none; border-bottom: 3px solid transparent;
    color: #ccc; padding: 10px 18px; cursor: pointer; font-size: 0.9rem;
    white-space: nowrap; transition: color 0.15s, border-color 0.15s;
  }}
  .tab-btn:hover {{ color: #fff; border-bottom-color: #7f8c8d; }}
  .tab-btn.active {{ color: #fff; border-bottom-color: #3498db; font-weight: 600; }}
  .tab-panel {{ min-height: 100vh; }}
</style>
<script>
  function showTab(id) {{
    document.querySelectorAll('.tab-panel').forEach(function(p) {{
      p.style.display = 'none';
    }});
    document.querySelectorAll('.tab-btn').forEach(function(b) {{
      b.classList.remove('active');
    }});
    document.getElementById(id).style.display = 'block';
    document.getElementById('btn-' + id).classList.add('active');
  }}
  document.addEventListener('DOMContentLoaded', function() {{
    showTab('tab-current');
  }});
</script>
</head>
<body>
<div class="tab-nav">
  {nav_html}
</div>
{panels_html}
{bootstrap_js}
</body>
</html>"""
    return html
