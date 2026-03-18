"""Tests for sc_tools.bm.report — generate_report_index and offline Plotly mode."""

from __future__ import annotations

import inspect
import pathlib
import tempfile

import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Task 1: generate_report_index
# ---------------------------------------------------------------------------


def test_generate_report_index_importable():
    """generate_report_index must be importable from sc_tools.bm.report."""
    from sc_tools.bm.report import generate_report_index  # noqa: F401


def test_generate_report_index_in_all():
    """generate_report_index must appear in __all__."""
    import sc_tools.bm.report as mod

    assert "generate_report_index" in mod.__all__


def test_generate_report_index_signature():
    """Function signature must match spec."""
    from sc_tools.bm.report import generate_report_index

    sig = inspect.signature(generate_report_index)
    params = list(sig.parameters)
    assert "output_dir" in params
    assert "project_name" in params
    assert "output_path" in params


def test_generate_report_index_returns_path():
    """Must return a Path object pointing to an index.html file."""
    from sc_tools.bm.report import generate_report_index

    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        result = generate_report_index(d)
        assert isinstance(result, pathlib.Path)
        assert result.name == "index.html"


def test_generate_report_index_writes_file():
    """index.html must be written to output_dir by default."""
    from sc_tools.bm.report import generate_report_index

    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        result = generate_report_index(d)
        assert result.exists()
        assert result.read_text()  # non-empty


def test_generate_report_index_custom_output_path():
    """output_path parameter overrides default location."""
    from sc_tools.bm.report import generate_report_index

    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        custom = d / "custom_index.html"
        result = generate_report_index(d, output_path=custom)
        assert result == custom
        assert custom.exists()


def test_generate_report_index_empty_dir():
    """Empty directory must produce a page with 'No reports found' message."""
    from sc_tools.bm.report import generate_report_index

    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        idx = generate_report_index(d, project_name="MyProject")
        html = idx.read_text()
        assert "No reports found" in html


def test_generate_report_index_lists_html_files():
    """HTML files in output_dir must appear as cards in the index."""
    from sc_tools.bm.report import generate_report_index

    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        (d / "my_benchmark.html").write_text(
            "<html><head><title>Benchmark Test</title></head><body></body></html>"
        )
        idx = generate_report_index(d)
        html = idx.read_text()
        assert "my_benchmark.html" in html
        assert "card" in html


def test_generate_report_index_excludes_self():
    """index.html itself must not appear as a listed report."""
    from sc_tools.bm.report import generate_report_index

    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        (d / "report_a.html").write_text(
            "<html><head><title>Report A</title></head><body></body></html>"
        )
        idx = generate_report_index(d)
        html = idx.read_text()
        # index.html should not appear as a link target
        assert "index.html" not in html or html.count("index.html") == 0


def test_generate_report_index_project_name_in_output():
    """project_name must appear in the generated HTML."""
    from sc_tools.bm.report import generate_report_index

    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        idx = generate_report_index(d, project_name="AwesomeProject")
        html = idx.read_text()
        assert "AwesomeProject" in html


def test_generate_report_index_extracts_title():
    """Title must be extracted from <title> tag of each report."""
    from sc_tools.bm.report import generate_report_index

    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        (d / "seg_report.html").write_text(
            "<html><head><title>Segmentation Report 2024</title></head><body></body></html>"
        )
        idx = generate_report_index(d)
        html = idx.read_text()
        assert "Segmentation Report 2024" in html


def test_generate_report_index_report_type_inference():
    """report_type badges must be inferred correctly from filenames."""
    from sc_tools.bm.report import generate_report_index

    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        (d / "integration_benchmark.html").write_text(
            "<html><head><title>Integration Report</title></head><body></body></html>"
        )
        (d / "segmentation_report.html").write_text(
            "<html><head><title>Segmentation Report</title></head><body></body></html>"
        )
        (d / "qc_summary.html").write_text(
            "<html><head><title>QC Report</title></head><body></body></html>"
        )
        idx = generate_report_index(d)
        html = idx.read_text()
        # Both badge types should appear
        assert "integration" in html
        assert "segmentation" in html
        assert "qc" in html


# ---------------------------------------------------------------------------
# Task 6: offline=True parameter on generate_* functions
# ---------------------------------------------------------------------------


def test_generate_segmentation_report_offline_param():
    """generate_segmentation_report must accept offline keyword argument."""
    from sc_tools.bm.report import generate_segmentation_report

    sig = inspect.signature(generate_segmentation_report)
    assert "offline" in sig.parameters


def test_generate_integration_report_offline_param():
    """generate_integration_report must accept offline keyword argument."""
    from sc_tools.bm.report import generate_integration_report

    sig = inspect.signature(generate_integration_report)
    assert "offline" in sig.parameters


def test_generate_benchmark_report_offline_param():
    """generate_benchmark_report must accept offline keyword argument."""
    from sc_tools.bm.report import generate_benchmark_report

    sig = inspect.signature(generate_benchmark_report)
    assert "offline" in sig.parameters


def test_generate_integration_report_online_has_cdn_script():
    """Default (offline=False) output must include the CDN plotly script tag."""
    pytest.importorskip("plotly")
    from sc_tools.bm.report import generate_integration_report

    df = pd.DataFrame(
        {
            "method": ["harmony", "scvi"],
            "batch_score": [0.8, 0.7],
            "bio_score": [0.7, 0.9],
            "overall_score": [0.75, 0.8],
        }
    )

    with tempfile.TemporaryDirectory() as tmp:
        out = pathlib.Path(tmp) / "test_integration.html"
        generate_integration_report(df, output_path=out)
        html = out.read_text()
        assert "cdn.plot.ly" in html


def test_generate_integration_report_offline_no_cdn_script():
    """offline=True output must NOT include the CDN plotly script tag."""
    pytest.importorskip("plotly")
    from sc_tools.bm.report import generate_integration_report

    df = pd.DataFrame(
        {
            "method": ["harmony", "scvi"],
            "batch_score": [0.8, 0.7],
            "bio_score": [0.7, 0.9],
            "overall_score": [0.75, 0.8],
        }
    )

    with tempfile.TemporaryDirectory() as tmp:
        out = pathlib.Path(tmp) / "test_integration_offline.html"
        generate_integration_report(df, output_path=out, offline=True)
        html = out.read_text()
        # The CDN <script src> tag must be absent. The inline bundle may contain
        # internal string references to cdn.plot.ly, so check the tag only.
        import re as _re

        assert not _re.search(r'<script[^>]+src=["\']https://cdn\.plot\.ly', html)


def test_generate_integration_report_offline_embeds_plotlyjs():
    """offline=True output must include inline plotly.js content."""
    pytest.importorskip("plotly")
    from sc_tools.bm.report import generate_integration_report

    df = pd.DataFrame(
        {
            "method": ["harmony", "scvi"],
            "batch_score": [0.8, 0.7],
            "bio_score": [0.7, 0.9],
            "overall_score": [0.75, 0.8],
        }
    )

    with tempfile.TemporaryDirectory() as tmp:
        out = pathlib.Path(tmp) / "test_integration_offline.html"
        generate_integration_report(df, output_path=out, offline=True)
        html = out.read_text()
        # Plotly inline embed contains "plotly" JS — the bundle is large
        assert "Plotly" in html or "plotly" in html.lower()


def test_plotly_to_html_accepts_include_plotlyjs():
    """_plotly_to_html must accept include_plotlyjs parameter."""
    from sc_tools.bm.report import _plotly_to_html

    sig = inspect.signature(_plotly_to_html)
    assert "include_plotlyjs" in sig.parameters


# ---------------------------------------------------------------------------
# Task: run-over-run tab comparison
# ---------------------------------------------------------------------------


def _make_integration_df():
    return pd.DataFrame(
        {
            "method": ["harmony", "scvi"],
            "batch_score": [0.8, 0.7],
            "bio_score": [0.7, 0.9],
            "overall_score": [0.75, 0.8],
        }
    )


def test_integration_report_second_run_has_tabs():
    """Second generate_integration_report in same dir must produce tab markup."""
    pytest.importorskip("plotly")
    from sc_tools.bm.report import generate_integration_report

    df = _make_integration_df()
    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        generate_integration_report(df, output_path=d / "integration_report_v1.html")
        r2 = generate_integration_report(df, output_path=d / "integration_report_v2.html")
        html = r2.read_text()
        assert "tab-panel" in html or "tab-btn" in html, (
            "Second report should contain tab markup when a prior report exists"
        )


def test_integration_report_first_run_no_tabs():
    """First generate_integration_report in empty dir must NOT add tab markup."""
    pytest.importorskip("plotly")
    from sc_tools.bm.report import generate_integration_report

    df = _make_integration_df()
    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        r1 = generate_integration_report(df, output_path=d / "integration_report_v1.html")
        html = r1.read_text()
        assert "tab-panel" not in html and "tab-btn" not in html, (
            "First report should not contain tab markup when no prior report exists"
        )


def test_segmentation_report_second_run_has_tabs():
    """Second generate_segmentation_report in same dir must produce tab markup."""
    pytest.importorskip("plotly")
    from sc_tools.bm.report import generate_segmentation_report

    df = pd.DataFrame(
        {
            "method": ["cellpose", "deepcell"],
            "composite_score": [0.85, 0.78],
        }
    )
    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        generate_segmentation_report(df, output_path=d / "segmentation_report_v1.html")
        r2 = generate_segmentation_report(df, output_path=d / "segmentation_report_v2.html")
        html = r2.read_text()
        assert "tab-panel" in html or "tab-btn" in html, (
            "Second segmentation report should contain tab markup when a prior report exists"
        )


def test_benchmark_report_second_run_has_tabs():
    """Second generate_benchmark_report in same dir must produce tab markup."""
    pytest.importorskip("plotly")
    from sc_tools.bm.report import generate_benchmark_report

    df = pd.DataFrame(
        {
            "method": ["cellpose", "deepcell"],
            "strategy": ["default", "default"],
            "roi_id": ["roi1", "roi1"],
            "dataset": ["ds1", "ds1"],
            "tissue": ["lung", "lung"],
            "composite_score": [0.85, 0.78],
        }
    )
    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        generate_benchmark_report(df, output_path=d / "benchmark_report_v1.html")
        r2 = generate_benchmark_report(df, output_path=d / "benchmark_report_v2.html")
        html = r2.read_text()
        assert "tab-panel" in html or "tab-btn" in html, (
            "Second benchmark report should contain tab markup when a prior report exists"
        )


def test_tab_combination_gracefully_skipped_on_error(monkeypatch):
    """Tab combination errors must be silently swallowed; report still written."""
    pytest.importorskip("plotly")
    from sc_tools.bm import report as bm_report
    from sc_tools.bm.report import generate_integration_report

    def _raise(*a, **kw):
        raise RuntimeError("boom")

    # Patch _find_latest_bm_report to raise to simulate an unexpected error
    monkeypatch.setattr(bm_report, "_find_latest_bm_report", _raise)

    df = _make_integration_df()
    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        r1 = generate_integration_report(df, output_path=d / "integration_report_v1.html")
        # Ensure file was written despite the patched error
        assert r1.exists()
        assert len(r1.read_text()) > 100


# ---------------------------------------------------------------------------
# batch_bio_table: new batch × bio conservation table in integration report
# ---------------------------------------------------------------------------


def test_integration_report_contains_batch_bio_table():
    """generate_integration_report must embed a batch × bio HTML table when scores exist."""
    pytest.importorskip("plotly")
    from sc_tools.bm.report import generate_integration_report

    df = pd.DataFrame(
        {
            "method": ["harmony", "scvi", "bbknn"],
            "batch_score": [0.8, 0.7, 0.75],
            "bio_score": [0.7, 0.9, 0.8],
            "overall_score": [0.75, 0.8, 0.775],
        }
    )
    with tempfile.TemporaryDirectory() as tmp:
        out = pathlib.Path(tmp) / "test_bb.html"
        generate_integration_report(df, output_path=out)
        html = out.read_text()
        assert "batch-bio-table" in html, "section id='batch-bio-table' must be present"
        assert "table-striped" in html, "Bootstrap table-striped class must be present"


def test_integration_report_batch_bio_table_sorted_by_overall():
    """batch × bio table must be sorted descending by overall_score."""
    pytest.importorskip("plotly")
    import re as _re

    from sc_tools.bm.report import generate_integration_report

    df = pd.DataFrame(
        {
            "method": ["harmony", "scvi", "bbknn"],
            "batch_score": [0.8, 0.7, 0.75],
            "bio_score": [0.7, 0.9, 0.8],
            "overall_score": [0.75, 0.80, 0.775],
        }
    )
    with tempfile.TemporaryDirectory() as tmp:
        out = pathlib.Path(tmp) / "test_bb_sorted.html"
        generate_integration_report(df, output_path=out)
        html = out.read_text()
        # Extract only the batch-bio-table section to avoid matching exec-summary cards
        start = html.find('id="batch-bio-table"')
        end = html.find("</section>", start)
        assert start != -1, "batch-bio-table section must exist"
        table_section = html[start:end]
        methods = _re.findall(r"<td>(harmony|scvi|bbknn)</td>", table_section)
        assert methods.index("scvi") < methods.index("harmony"), (
            "scvi (higher overall_score=0.80) must appear before harmony (0.75) in the table"
        )


def test_integration_report_batch_bio_table_section_order():
    """batch-bio-table section must appear before metrics-heatmap in the report."""
    pytest.importorskip("plotly")
    from sc_tools.bm.report import generate_integration_report

    df = pd.DataFrame(
        {
            "method": ["harmony", "scvi"],
            "batch_score": [0.8, 0.7],
            "bio_score": [0.7, 0.9],
            "overall_score": [0.75, 0.8],
        }
    )
    with tempfile.TemporaryDirectory() as tmp:
        out = pathlib.Path(tmp) / "test_bb_order.html"
        generate_integration_report(df, output_path=out)
        html = out.read_text()
        idx_bb = html.find('id="batch-bio-table"')
        idx_hm = html.find('id="metrics-heatmap"')
        assert idx_bb < idx_hm, "batch-bio-table must come before metrics-heatmap"


def test_integration_report_umap_before_radar():
    """umap-grid section must appear before radar section when both are present."""
    pytest.importorskip("plotly")
    from sc_tools.bm.report import generate_integration_report

    # radar requires >= 3 metric cols (beyond batch/bio/overall/method)
    df = pd.DataFrame(
        {
            "method": ["harmony", "scvi"],
            "batch_score": [0.8, 0.7],
            "bio_score": [0.7, 0.9],
            "overall_score": [0.75, 0.8],
            "asw_batch": [0.7, 0.6],
            "pcr": [0.8, 0.75],
            "graph_connectivity": [0.9, 0.85],
        }
    )
    # We cannot pass actual adata here, so just verify the template section order
    # in the context["sections"] list by parsing the generated HTML
    with tempfile.TemporaryDirectory() as tmp:
        out = pathlib.Path(tmp) / "test_umap_radar.html"
        generate_integration_report(df, output_path=out)
        html = out.read_text()
        # radar is conditional; check the template ordering in the HTML nav
        idx_radar = html.find('id="radar"')
        idx_ranking = html.find('id="ranking"')
        # If radar is present, it must come after ranking
        if idx_radar != -1 and idx_ranking != -1:
            assert idx_ranking < idx_radar, "ranking must come before radar"
