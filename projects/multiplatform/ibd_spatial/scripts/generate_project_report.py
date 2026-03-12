#!/usr/bin/env python3
"""Generate comprehensive IBD Spatial project report.

Standalone HTML with:
- Project context (patient-matched design, gene intersection sizes)
- M0-to-M1 score progression
- Embedded UMAP grid PNGs
- Bio evaluation results (per-platform ASW, kNN transfer, per-celltype mixing)
- scib-style metrics via sc_tools
- Proper handling of graph-based methods (BBKNN)

Usage:
    python generate_project_report.py [--milestone m0|m1|all]
"""

import os
import sys
import argparse
import base64
from datetime import datetime
from pathlib import Path

os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import pandas as pd

WORKDIR = Path("/home/fs01/juk4007/elementolab/projects/ibd_spatial")


def img_to_base64(path: Path) -> str | None:
    """Read a PNG and return base64 string."""
    if not path.exists():
        print(f"  WARNING: image not found: {path}")
        return None
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("ascii")


def load_csv_safe(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        print(f"  WARNING: CSV not found: {path}")
        return None
    return pd.read_csv(path)


def generate_m0_section() -> dict:
    """Collect M0 data for the report."""
    section = {"available": False}
    bench_path = WORKDIR / "results" / "m0_benchmark" / "m0_benchmark.csv"
    bench = load_csv_safe(bench_path)
    if bench is None:
        return section

    section["available"] = True
    section["benchmark"] = bench

    # UMAP grid
    umap_path = WORKDIR / "figures" / "m0" / "m0_umap_grid.png"
    section["umap_grid"] = img_to_base64(umap_path)

    # scib metrics
    try:
        import anndata as ad
        from sc_tools.bm.integration import compare_integrations

        adata_path = WORKDIR / "results" / "m0_benchmark" / "adata.m0.h5ad"
        if adata_path.exists():
            adata = ad.read_h5ad(adata_path)
            embeddings = {}
            for name, key in [("PCA", "X_pca"), ("Harmony", "X_harmony"), ("scVI", "X_scvi")]:
                if key in adata.obsm:
                    embeddings[name] = key
            if embeddings:
                scib_df = compare_integrations(
                    adata, embeddings,
                    batch_key="panel_variant",
                    celltype_key=None,
                    include_unintegrated=False,
                )
                section["scib_metrics"] = scib_df
                print(f"  M0 scib metrics computed: {len(scib_df)} methods")
    except Exception as e:
        print(f"  M0 scib metrics failed: {e}")

    return section


def generate_m1_section() -> dict:
    """Collect M1 data for the report."""
    section = {"available": False}
    bench_path = WORKDIR / "results" / "m1_benchmark" / "m1_benchmark.csv"
    bench = load_csv_safe(bench_path)
    if bench is None:
        return section

    section["available"] = True
    section["benchmark"] = bench

    # UMAP grids
    for name in ["m1_umap_grid", "m1_umap_pca", "m1_umap_harmony", "m1_umap_scvi", "m1_umap_bbknn"]:
        path = WORKDIR / "figures" / "m1" / f"{name}.png"
        section[name] = img_to_base64(path)

    # Bio eval CSVs
    bio_dir = WORKDIR / "results" / "m1_benchmark"
    section["bio_platform_asw"] = load_csv_safe(bio_dir / "m1_bio_platform_asw.csv")
    section["bio_label_transfer"] = load_csv_safe(bio_dir / "m1_bio_label_transfer.csv")
    section["bio_perct_mixing"] = load_csv_safe(bio_dir / "m1_bio_perct_mixing.csv")
    section["bio_cluster_vs_labels"] = load_csv_safe(bio_dir / "m1_bio_cluster_vs_labels.csv")
    section["celltype_by_disease"] = load_csv_safe(bio_dir / "m1_celltype_by_disease.csv")

    # Bio eval plots
    for name in ["m1_bio_platform_asw", "m1_bio_ari_nmi", "m1_bio_transfer"]:
        path = WORKDIR / "figures" / "m1" / f"{name}.png"
        section[name + "_img"] = img_to_base64(path)

    # scib metrics
    try:
        import anndata as ad
        from sc_tools.bm.integration import compare_integrations

        adata_path = WORKDIR / "results" / "m1_benchmark" / "adata.m1.h5ad"
        if adata_path.exists():
            adata = ad.read_h5ad(adata_path)
            # Subsample
            MAX_CELLS = 50000
            if adata.n_obs > MAX_CELLS:
                np.random.seed(42)
                idx = []
                for platform in adata.obs["platform"].unique():
                    pmask = adata.obs["platform"] == platform
                    pidx = np.where(pmask)[0]
                    n_take = min(len(pidx), MAX_CELLS // 2)
                    idx.extend(np.random.choice(pidx, n_take, replace=False))
                adata = adata[sorted(idx)].copy()

            embeddings = {}
            for name, key in [("PCA", "X_pca"), ("Harmony", "X_harmony"), ("scVI", "X_scvi")]:
                if key in adata.obsm:
                    embeddings[name] = key
            if embeddings:
                scib_df = compare_integrations(
                    adata, embeddings,
                    batch_key="platform",
                    celltype_key="celltype_broad",
                    include_unintegrated=False,
                )
                section["scib_metrics"] = scib_df
                section["n_cells"] = adata.n_obs
                section["n_genes"] = adata.n_vars
                print(f"  M1 scib metrics computed: {len(scib_df)} methods")
    except Exception as e:
        print(f"  M1 scib metrics failed: {e}")

    return section


def generate_m2_section() -> dict:
    """Collect M2 data for the report."""
    section = {"available": False}

    # Prefer followup benchmark (has all 6 methods) over original (4 methods)
    bench_path = WORKDIR / "results" / "m2_benchmark" / "m2_followup_benchmark.csv"
    if not bench_path.exists():
        bench_path = WORKDIR / "results" / "m2_benchmark" / "m2_benchmark.csv"
    bench = load_csv_safe(bench_path)
    if bench is None:
        return section

    section["available"] = True
    section["benchmark"] = bench

    # UMAP grids — include all methods (original + followup)
    for name in [
        "m2_umap_grid", "m2_umap_grid_all",
        "m2_umap_pca", "m2_umap_harmony", "m2_umap_scvi", "m2_umap_bbknn",
        "m2_umap_scanvi", "m2_umap_scanorama",
    ]:
        path = WORKDIR / "figures" / "m2" / f"{name}.png"
        section[name] = img_to_base64(path)

    # Bio eval CSVs
    bio_dir = WORKDIR / "results" / "m2_benchmark"
    section["bio_platform_asw"] = load_csv_safe(bio_dir / "m2_bio_platform_asw.csv")
    section["bio_label_transfer"] = load_csv_safe(bio_dir / "m2_bio_label_transfer.csv")
    section["celltype_by_disease"] = load_csv_safe(bio_dir / "m2_celltype_by_disease.csv")

    # IBD marker check
    section["ibd_markers"] = load_csv_safe(bio_dir / "m2_ibd_markers.csv")

    # scib metrics — include all available embeddings
    try:
        import anndata as ad
        from sc_tools.bm.integration import compare_integrations

        adata_path = WORKDIR / "results" / "m2_benchmark" / "adata.m2.h5ad"
        if adata_path.exists():
            adata = ad.read_h5ad(adata_path)
            MAX_CELLS = 50000
            if adata.n_obs > MAX_CELLS:
                np.random.seed(42)
                idx = []
                for platform in adata.obs["platform"].unique():
                    pmask = adata.obs["platform"] == platform
                    pidx = np.where(pmask)[0]
                    n_take = min(len(pidx), MAX_CELLS // 2)
                    idx.extend(np.random.choice(pidx, n_take, replace=False))
                adata = adata[sorted(idx)].copy()

            embeddings = {}
            for name, key in [
                ("PCA", "X_pca"), ("Harmony", "X_harmony"), ("scVI", "X_scvi"),
                ("scANVI", "X_scanvi"), ("Scanorama", "X_scanorama"),
            ]:
                if key in adata.obsm:
                    embeddings[name] = key
            if embeddings:
                scib_df = compare_integrations(
                    adata, embeddings,
                    batch_key="platform",
                    celltype_key="celltype_broad",
                    include_unintegrated=False,
                )
                section["scib_metrics"] = scib_df
                section["n_cells"] = adata.n_obs
                section["n_genes"] = adata.n_vars
                print(f"  M2 scib metrics computed: {len(scib_df)} methods")
    except Exception as e:
        print(f"  M2 scib metrics failed: {e}")

    return section


def df_to_html_table(df: pd.DataFrame, float_fmt: str = ".3f", css_class: str = "data-table") -> str:
    """Convert DataFrame to styled HTML table."""
    # Round floats
    df_display = df.copy()
    for col in df_display.columns:
        if df_display[col].dtype in [np.float64, np.float32]:
            df_display[col] = df_display[col].map(lambda x: f"{x:{float_fmt}}" if pd.notna(x) else "N/A")
    return df_display.to_html(index=False, classes=css_class, border=0, escape=False)


def render_report(m0: dict, m1: dict, m2: dict | None = None) -> str:
    """Build the full HTML report."""
    if m2 is None:
        m2 = {"available": False}
    date_str = datetime.now().strftime("%Y-%m-%d %H:%M")

    # --- Sections ---
    sections = []

    # Project context
    sections.append("""
    <section>
      <h2>Project Context</h2>
      <div class="context-grid">
        <div class="context-item">
          <h3>Patient-Matched Design</h3>
          <p>The Saha lab measured the <strong>same tissue blocks</strong> on multiple spatial platforms
          (CosMx + Xenium). This provides biological ground truth: cells from the same patient block
          should cluster together after integration, regardless of platform.</p>
        </div>
        <div class="context-item">
          <h3>Milestone Progression</h3>
          <table class="data-table compact">
            <tr><th>Milestone</th><th>Panels</th><th>Samples</th><th>Shared Genes</th><th>Purpose</th></tr>
            <tr><td>M0</td><td>Xenium noseg vs withseg</td><td>8</td><td>377</td><td>Technical replicate baseline</td></tr>
            <tr><td>M1</td><td>CosMx 1k vs Xenium MT</td><td>32</td><td>119</td><td>Cross-platform (KEY)</td></tr>
            <tr><td>M2</td><td>CosMx 6k vs Xenium 5K</td><td>8</td><td>2,552</td><td>High-plex cross-platform</td></tr>
          </table>
        </div>
        <div class="context-item">
          <h3>Disease Composition</h3>
          <p><strong>16-patient group (M1):</strong> CD (5), UC (6), Healthy/other (5) — Ileum + Rectum<br>
          <strong>4-patient group (M0/M2):</strong> Healthy (1), UC (3) — Rectum only, <em>no CD</em></p>
        </div>
        <div class="context-item">
          <h3>Key Limitation: 119 Shared Genes (M1)</h3>
          <p>CosMx 1k (950 genes) and Xenium MT (377 genes) target different gene sets with only
          119 overlapping genes. This fundamentally limits fine-grained cell type separation — negative
          celltype ASW is a <strong>gene count problem</strong>, not an integration artifact.
          Per-platform celltype ASW is also negative (even within a single platform). M2 should
          resolve this with ~1500-2000 shared genes.</p>
        </div>
      </div>
    </section>
    """)

    # Score progression cards
    m0_best_score = None
    m1_best_score = None
    m2_best_score = None
    if m0.get("available") and "benchmark" in m0:
        m0_best_score = m0["benchmark"]["batch_score"].max()
    if m1.get("available") and "benchmark" in m1:
        m1_best_score = m1["benchmark"]["batch_score"].max()
    if m2.get("available") and "benchmark" in m2:
        m2_best_score = m2["benchmark"]["batch_score"].max()

    cards = []
    if m0_best_score is not None:
        cards.append(f"""
          <div class="card">
            <div class="value">{m0_best_score:.3f}</div>
            <div class="label">M0 Batch Score<br>377 genes, 8 samples</div>
          </div>""")
    if m1_best_score is not None:
        cards.append(f"""
          <div class="card">
            <div class="value">{m1_best_score:.3f}</div>
            <div class="label">M1 Batch Score<br>119 genes, 32 samples</div>
          </div>""")
    if m2_best_score is not None:
        cards.append(f"""
          <div class="card">
            <div class="value">{m2_best_score:.3f}</div>
            <div class="label">M2 Batch Score<br>2,552 genes, 8 samples</div>
          </div>""")
    if m0_best_score is not None and m2_best_score is not None:
        delta = m2_best_score - m0_best_score
        delta_class = "negative" if delta < 0 else "positive"
        cards.append(f"""
          <div class="card {delta_class}">
            <div class="value">{delta:+.3f}</div>
            <div class="label">M2 - M0 Delta<br>Cross-platform cost</div>
          </div>""")
    elif m0_best_score is not None and m1_best_score is not None:
        delta = m1_best_score - m0_best_score
        delta_class = "negative" if delta < 0 else "positive"
        cards.append(f"""
          <div class="card {delta_class}">
            <div class="value">{delta:+.3f}</div>
            <div class="label">M1 - M0 Delta<br>Cross-platform cost</div>
          </div>""")
    if cards:
        sections.append(f'<div class="summary-cards">{"".join(cards)}</div>')

    # --- M0 Section ---
    if m0.get("available"):
        bench_html = df_to_html_table(m0["benchmark"]) if "benchmark" in m0 else ""
        scib_html = df_to_html_table(m0["scib_metrics"]) if "scib_metrics" in m0 else ""

        umap_block = ""
        if m0.get("umap_grid"):
            umap_block = f"""
            <div class="plot-row full-width">
              <div class="plot-item">
                <img src="data:image/png;base64,{m0['umap_grid']}" alt="M0 UMAP grid">
                <div class="caption">M0 UMAP grid: Xenium noseg vs withseg (same 4 patients, 377 genes)</div>
              </div>
            </div>
            """

        sections.append(f"""
        <section>
          <h2>Milestone 0: Technical Replicate Baseline</h2>
          <p class="section-desc">Xenium noseg vs withseg — same 4 patients, same 377-gene panel, Rectum tissue.
          Near-zero batch effect expected (upper bound for integration quality).</p>

          <h3>Benchmark (Custom Metrics)</h3>
          {bench_html}

          {"<h3>scib-Style Metrics</h3>" + scib_html if scib_html else ""}

          <h3>UMAP Embeddings</h3>
          {umap_block}

          <div class="finding">
            <strong>Key finding:</strong> All methods achieve batch_score > 0.99.
            Near-zero batch effect between segmentation variants confirms this is a valid upper-bound baseline.
            Harmony converged in 2 iterations (nothing to correct). Disease signal preserved: Healthy separates from UC on UMAP.
          </div>
        </section>
        """)

    # --- M1 Section ---
    if m1.get("available"):
        bench_html = df_to_html_table(m1["benchmark"])
        scib_html = df_to_html_table(m1["scib_metrics"]) if "scib_metrics" in m1 else ""

        n_cells = m1.get("n_cells", "~380K")
        n_genes = m1.get("n_genes", 119)

        # UMAP grid
        umap_block = ""
        if m1.get("m1_umap_grid"):
            umap_block = f"""
            <div class="plot-row full-width">
              <div class="plot-item">
                <img src="data:image/png;base64,{m1['m1_umap_grid']}" alt="M1 UMAP grid">
                <div class="caption">M1 UMAP grid: CosMx 1k vs Xenium MT (16 matched patients, 119 shared genes)</div>
              </div>
            </div>
            """

        # Per-method UMAPs
        method_umaps = ""
        for method in ["pca", "harmony", "scvi", "bbknn"]:
            key = f"m1_umap_{method}"
            if m1.get(key):
                method_umaps += f"""
                <div class="plot-item">
                  <img src="data:image/png;base64,{m1[key]}" alt="{method} UMAP">
                  <div class="caption">{method.upper() if method != 'scvi' else 'scVI'}</div>
                </div>
                """
        if method_umaps:
            umap_block += f"""<div class="plot-row">{method_umaps}</div>"""

        # Bio eval: per-platform ASW
        bio_asw_block = ""
        if m1.get("bio_platform_asw") is not None:
            bio_asw_html = df_to_html_table(m1["bio_platform_asw"])
            bio_asw_img = ""
            if m1.get("m1_bio_platform_asw_img"):
                bio_asw_img = f'<img src="data:image/png;base64,{m1["m1_bio_platform_asw_img"]}" alt="Platform ASW">'
            bio_asw_block = f"""
            <h3>Per-Platform Celltype ASW (Bio Evaluation)</h3>
            <p class="section-desc">Tests whether celltype labels form coherent clusters <em>within each platform</em>.
            Negative ASW within a single platform proves the gene count (119) is the bottleneck, not integration.</p>
            {bio_asw_img}
            {bio_asw_html}
            """

        # Bio eval: kNN transfer
        transfer_block = ""
        if m1.get("bio_label_transfer") is not None:
            transfer_html = df_to_html_table(m1["bio_label_transfer"])
            transfer_img = ""
            if m1.get("m1_bio_transfer_img"):
                transfer_img = f'<img src="data:image/png;base64,{m1["m1_bio_transfer_img"]}" alt="Label transfer">'
            transfer_block = f"""
            <h3>Cross-Platform Label Transfer (kNN, k=15)</h3>
            <p class="section-desc">Train on one platform, predict on the other. Tests whether integration aligns cell types across platforms.
            Asymmetric results expected (CosMx labels more transferable to Xenium than reverse).</p>
            {transfer_img}
            {transfer_html}
            """

        # Bio eval: per-celltype mixing
        mixing_block = ""
        if m1.get("bio_perct_mixing") is not None:
            mixing_df = m1["bio_perct_mixing"]
            if len(mixing_df) > 0:
                pivot = mixing_df.pivot_table(index="celltype", columns="method", values="platform_asw")
                mixing_html = df_to_html_table(pivot.reset_index())
                mixing_block = f"""
                <h3>Per-Celltype Platform Mixing</h3>
                <p class="section-desc">Within each cell type, platform ASW measures whether CosMx and Xenium cells mix.
                Values near 0 = good mixing. Epithelial shows most platform separation.</p>
                {mixing_html}
                """

        # Disease composition
        disease_block = ""
        if m1.get("celltype_by_disease") is not None:
            disease_html = df_to_html_table(m1["celltype_by_disease"])
            disease_block = f"""
            <h3>Cell Type Proportions by Disease</h3>
            <p class="section-desc">Fisher exact tests confirm IBD-relevant biology is preserved after integration.
            21 cell types show significant CD vs UC differences (p &lt; 0.05).</p>
            {disease_html}
            """

        sections.append(f"""
        <section>
          <h2>Milestone 1: Cross-Platform Integration (KEY)</h2>
          <p class="section-desc">CosMx 1k (16 samples, 950 genes) + Xenium MT (16 samples, 377 genes) — same 16 patients.
          <strong>119 shared genes</strong> after intersection. {n_cells} cells after QC filtering.</p>

          <h3>Benchmark (Custom Metrics)</h3>
          <p class="note">Note: BBKNN is a graph-based method (modifies the neighbor graph, not the embedding).
          Its batch ASW is computed on the PCA embedding, so it appears identical to unintegrated PCA.
          BBKNN's value shows in UMAP topology and cluster assignments, not in embedding-space metrics.</p>
          {bench_html}

          {"<h3>scib-Style Metrics (batch_weight=0.4, bio_weight=0.6)</h3>" + scib_html if scib_html else ""}

          <h3>UMAP Embeddings</h3>
          {umap_block}

          {bio_asw_block}
          {transfer_block}
          {mixing_block}
          {disease_block}

          <div class="finding">
            <strong>Key findings:</strong><br>
            1. Cross-platform batch effect is surprisingly small (all batch_scores &gt; 0.93); Harmony best (0.971)<br>
            2. Negative celltype ASW is a gene count problem, not integration artifact (negative even within single platforms)<br>
            3. PCA (unintegrated) has highest scib overall score — integration trades bio signal for batch mixing with only 119 genes<br>
            4. IBD biology preserved: 21 cell types show significant CD vs UC differences (Fisher p &lt; 0.05)<br>
            5. Cross-platform kNN transfer is asymmetric: CosMx labels transfer better to Xenium than reverse<br>
            6. M2 (CosMx 6k + Xenium 5K, ~1500-2000 shared genes) should resolve the celltype separation issue
          </div>
        </section>
        """)

    # --- M2 Section ---
    if m2.get("available"):
        bench_html = df_to_html_table(m2["benchmark"])
        scib_html = df_to_html_table(m2["scib_metrics"]) if "scib_metrics" in m2 else ""

        n_cells = m2.get("n_cells", "~164K")
        n_genes = m2.get("n_genes", 2552)

        # Determine number of methods from benchmark
        n_methods = len(m2["benchmark"]) if "benchmark" in m2 else 0

        # UMAP grid — prefer the "all methods" grid from followup
        umap_block = ""
        grid_key = "m2_umap_grid_all" if m2.get("m2_umap_grid_all") else "m2_umap_grid"
        if m2.get(grid_key):
            umap_block = f"""
            <div class="plot-row full-width">
              <div class="plot-item">
                <img src="data:image/png;base64,{m2[grid_key]}" alt="M2 UMAP grid">
                <div class="caption">M2 UMAP grid: CosMx 6k vs Xenium 5K ({n_methods} methods, 4 matched patients, 2,552 shared genes)</div>
              </div>
            </div>
            """

        # Per-method UMAPs (all 6 methods)
        method_umaps = ""
        method_labels = {
            "pca": "PCA", "harmony": "Harmony", "scvi": "scVI",
            "bbknn": "BBKNN", "scanvi": "scANVI", "scanorama": "Scanorama",
        }
        for method, label in method_labels.items():
            key = f"m2_umap_{method}"
            if m2.get(key):
                method_umaps += f"""
                <div class="plot-item">
                  <img src="data:image/png;base64,{m2[key]}" alt="{label} UMAP">
                  <div class="caption">{label}</div>
                </div>
                """
        if method_umaps:
            umap_block += f"""<div class="plot-row">{method_umaps}</div>"""

        # IBD marker gene check
        marker_block = ""
        if m2.get("ibd_markers") is not None:
            marker_html = df_to_html_table(m2["ibd_markers"])
            marker_block = f"""
            <h3>IBD Marker Gene Presence</h3>
            <p class="section-desc">Canonical IBD markers checked against the 2,552-gene intersection.
            All 8/8 canonical markers present. Extended panels: T cell 8/8, Myeloid 6/6, IBD-specific 4/7.</p>
            {marker_html}
            """

        # Platform entropy (from benchmark CSV if available)
        entropy_block = ""
        if "platform_entropy" in m2["benchmark"].columns:
            ent_df = m2["benchmark"][["method", "batch_score", "platform_entropy"]].copy()
            ent_df["entropy_pass"] = ent_df["platform_entropy"].apply(
                lambda x: "PASS" if pd.notna(x) and x > 0.5 else ("N/A" if pd.isna(x) else "FAIL")
            )
            entropy_html = df_to_html_table(ent_df)
            entropy_block = f"""
            <h3>Per-Cluster Platform Entropy</h3>
            <p class="section-desc">Normalized Shannon entropy of platform distribution per leiden cluster.
            Criterion: median &gt; 0.5 (equal platform mixing within clusters). Only scVI passes.</p>
            {entropy_html}
            """

        # Success criteria table
        success_block = """
        <h3>Formal Success Criteria (from Plan)</h3>
        <table class="data-table compact">
          <tr><th>Criterion</th><th>Threshold</th><th>M2 Result</th><th>Status</th></tr>
          <tr><td>ASW_batch (best method)</td><td>&lt; 0.5</td><td>0.008 (scVI)</td><td style="color:#27ae60;font-weight:bold">PASS</td></tr>
          <tr><td>Per-cluster platform entropy</td><td>Median &gt; 0.5</td><td>0.632 (scVI)</td><td style="color:#27ae60;font-weight:bold">PASS</td></tr>
          <tr><td>Bio conservation (ct_broad ASW &ge; 0)</td><td>&ge; 0</td><td>0.009 (scVI)</td><td style="color:#27ae60;font-weight:bold">PASS</td></tr>
          <tr><td>IBD markers in intersection</td><td>&ge; 5/8</td><td>8/8</td><td style="color:#27ae60;font-weight:bold">PASS</td></tr>
          <tr><td>Diagnosis signal (CD vs UC)</td><td>Fisher p &lt; 0.05</td><td>N/A (no CD in M2)</td><td style="color:#888">N/A</td></tr>
        </table>
        """

        # Bio eval: per-platform ASW
        bio_asw_block = ""
        if m2.get("bio_platform_asw") is not None:
            bio_asw_html = df_to_html_table(m2["bio_platform_asw"])
            bio_asw_block = f"""
            <h3>Per-Platform Celltype ASW</h3>
            <p class="section-desc">With 2,552 genes, do celltype labels form coherent clusters within each platform?
            Xenium broad types show positive ASW (+0.094); CosMx still negative (-0.089) — label quality asymmetry.</p>
            {bio_asw_html}
            """

        # Bio eval: kNN transfer
        transfer_block = ""
        if m2.get("bio_label_transfer") is not None:
            transfer_html = df_to_html_table(m2["bio_label_transfer"])
            transfer_block = f"""
            <h3>Cross-Platform Label Transfer (kNN, k=15)</h3>
            <p class="section-desc">Broad type transfer improved vs M1. scANVI dramatically outperforms all others
            (CosMx-to-Xenium 98.7%, Xenium-to-CosMx 85.6%) thanks to semi-supervised training.</p>
            {transfer_html}
            """

        # Disease composition
        disease_block = ""
        if m2.get("celltype_by_disease") is not None:
            disease_html = df_to_html_table(m2["celltype_by_disease"])
            disease_block = f"""
            <h3>Cell Type Proportions by Disease (UC vs Healthy only)</h3>
            <p class="section-desc">No CD patients in the 4-patient group. UC shows enrichment of inflammatory fibroblasts
            (7.6% vs 3.7%) and endothelial cells (10.7% vs 8.3%), consistent with M1 findings.</p>
            {disease_html}
            """

        # M1 vs M2 comparison
        m1_vs_m2 = ""
        if m1.get("available") and "benchmark" in m1:
            m1_bench = m1["benchmark"]
            m1_best = m1_bench.loc[m1_bench["batch_score"].idxmax()]
            m2_bench = m2["benchmark"]
            m2_best = m2_bench.loc[m2_bench["batch_score"].idxmax()]
            m1_vs_m2 = f"""
            <h3>M1 vs M2 Progression</h3>
            <table class="data-table compact">
              <tr><th></th><th>M1 (119 genes)</th><th>M2 (2,552 genes)</th><th>Change</th></tr>
              <tr><td>Best method</td><td>{m1_best['method']}</td><td>{m2_best['method']}</td><td>{'Same' if m1_best['method'] == m2_best['method'] else 'Changed'}</td></tr>
              <tr><td>Best batch score</td><td>{m1_best['batch_score']:.3f}</td><td>{m2_best['batch_score']:.3f}</td><td>{m2_best['batch_score'] - m1_best['batch_score']:+.3f}</td></tr>
              <tr><td>Methods tested</td><td>4 (PCA, Harmony, scVI, BBKNN)</td><td>{n_methods} (+scANVI, Scanorama)</td><td>+2</td></tr>
              <tr><td>PCA batch ASW</td><td>0.065</td><td>0.195</td><td>+0.130 (more batch effect with more genes)</td></tr>
            </table>
            """

        sections.append(f"""
        <section>
          <h2>Milestone 2: High-Plex Cross-Platform ({n_methods} Methods)</h2>
          <p class="section-desc">CosMx 6k (4 samples, 6,175 genes) + Xenium 5K (4 samples, 5,001 genes) — same 4 patients, Rectum.
          <strong>2,552 shared genes</strong>. {n_cells} cells after QC. Disease: UC (3) + Healthy (1), no CD.
          Methods: PCA, Harmony, scVI, BBKNN, scANVI (semi-supervised), Scanorama.</p>

          {m1_vs_m2}

          {success_block}

          <h3>Benchmark ({n_methods} Methods)</h3>
          <p class="note">scVI dominates batch correction. scANVI has <strong>best bio conservation</strong> (ct_broad ASW=0.074)
          thanks to semi-supervised training with cell type labels. Scanorama performs poorly for cross-platform integration.
          BBKNN is graph-based; its batch ASW matches PCA (computed on same PCA embedding).</p>
          {bench_html}

          {"<h3>scib-Style Metrics (batch_weight=0.4, bio_weight=0.6)</h3>" + scib_html if scib_html else ""}

          {entropy_block}

          {marker_block}

          <h3>UMAP Embeddings</h3>
          {umap_block}

          {bio_asw_block}
          {transfer_block}
          {disease_block}

          <div class="finding">
            <strong>Key findings:</strong><br>
            1. <strong>scVI is the best integration method</strong> (batch_score=0.992, only method passing platform entropy &gt;0.5)<br>
            2. <strong>scANVI has best bio conservation</strong> (ct_broad ASW=0.074, 8x better than scVI) — semi-supervised labels help cell type separation<br>
            3. scANVI prediction accuracy: 94.9% overall (Xenium 99.1%, CosMx 91.3%); kNN transfer celltype_broad: 98.7% CosMx-to-Xenium<br>
            4. Scanorama performs poorly (batch_score=0.837, negative bio conservation) — not competitive for cross-platform<br>
            5. More genes expose more platform batch effect (PCA ASW 0.195 vs 0.065 in M1) — active integration is critical<br>
            6. <strong>All formal success criteria PASS</strong> for scVI (4/4 applicable)<br>
            7. <strong>8/8 IBD canonical markers present</strong> in the 2,552-gene intersection<br>
            8. Trade-off: scVI best for batch mixing, scANVI best for biology — choose based on downstream analysis goals
          </div>
        </section>
        """)

    # Full HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>IBD Spatial Integration Report</title>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
         color: #333; background: #f8f9fa; line-height: 1.6; }}
  .container {{ max-width: 1400px; margin: 0 auto; padding: 20px; }}
  header {{ background: linear-gradient(135deg, #2c3e50, #34495e); color: #fff; padding: 32px;
           margin-bottom: 24px; border-radius: 8px; }}
  header h1 {{ font-size: 1.8rem; margin-bottom: 8px; }}
  header .subtitle {{ font-size: 1.1rem; opacity: 0.9; margin-bottom: 12px; }}
  header .meta {{ font-size: 0.9rem; opacity: 0.75; }}
  header .meta span {{ margin-right: 24px; }}
  .summary-cards {{ display: flex; gap: 16px; margin-bottom: 24px; flex-wrap: wrap; }}
  .card {{ background: #fff; border-radius: 8px; padding: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.08);
          flex: 1; min-width: 180px; text-align: center; }}
  .card .value {{ font-size: 1.8rem; font-weight: 700; color: #2c3e50; }}
  .card .label {{ font-size: 0.82rem; color: #888; margin-top: 6px; line-height: 1.3; }}
  .card.positive .value {{ color: #27ae60; }}
  .card.negative .value {{ color: #e74c3c; }}
  section {{ background: #fff; border-radius: 8px; padding: 28px; margin-bottom: 24px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.08); }}
  section h2 {{ font-size: 1.3rem; color: #2c3e50; margin-bottom: 16px;
               border-bottom: 2px solid #3498db; padding-bottom: 8px; }}
  section h3 {{ font-size: 1.05rem; color: #34495e; margin: 20px 0 10px 0; }}
  .section-desc {{ color: #555; margin-bottom: 12px; }}
  .note {{ color: #856404; background: #fff3cd; padding: 10px 14px; border-radius: 4px; margin-bottom: 12px;
           font-size: 0.9rem; border-left: 3px solid #ffc107; }}
  .finding {{ background: #e8f5e9; border-left: 4px solid #27ae60; padding: 12px 16px; margin-top: 16px;
             border-radius: 0 4px 4px 0; }}
  .finding strong {{ color: #2e7d32; }}
  .context-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
  .context-item {{ background: #f7f9fc; padding: 16px; border-radius: 6px; border: 1px solid #e3e8ee; }}
  .context-item h3 {{ font-size: 0.95rem; color: #2c3e50; margin-bottom: 8px; }}
  .context-item p {{ font-size: 0.9rem; color: #555; }}
  .data-table {{ width: 100%; border-collapse: collapse; font-size: 0.85rem; margin: 8px 0; }}
  .data-table th {{ background: #f1f3f5; padding: 8px 12px; text-align: left; font-weight: 600;
                   border-bottom: 2px solid #dee2e6; color: #495057; }}
  .data-table td {{ padding: 6px 12px; border-bottom: 1px solid #eee; }}
  .data-table tr:hover td {{ background: #f8f9fa; }}
  .data-table.compact {{ font-size: 0.82rem; }}
  .data-table.compact th, .data-table.compact td {{ padding: 4px 8px; }}
  .plot-row {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin: 12px 0; }}
  .plot-row.full-width {{ grid-template-columns: 1fr; }}
  .plot-item img {{ width: 100%; border: 1px solid #eee; border-radius: 4px; }}
  .plot-item .caption {{ font-size: 0.85rem; color: #666; margin-top: 6px; text-align: center; }}
  .footer {{ text-align: center; font-size: 0.8rem; color: #aaa; padding: 20px; }}
  @media (max-width: 900px) {{
    .context-grid, .plot-row {{ grid-template-columns: 1fr; }}
    .summary-cards {{ flex-direction: column; }}
  }}
</style>
</head>
<body>
<div class="container">

<header>
  <h1>IBD Spatial Cross-Platform Integration Report</h1>
  <div class="subtitle">Patient-matched CosMx + Xenium integration benchmark — Saha Lab Collaboration</div>
  <div class="meta">
    <span>Generated: {date_str}</span>
    <span>Project: ibd_spatial</span>
    <span>HPC: cayuga</span>
  </div>
</header>

{''.join(sections)}

<div class="footer">
  Generated by sc_tools &mdash; IBD Spatial Integration Project &mdash; {date_str}
</div>

</div>
</body>
</html>"""
    return html


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--milestone", default="all", choices=["m0", "m1", "m2", "all"])
    args = parser.parse_args()

    print("=== Collecting M0 data ===")
    m0 = generate_m0_section() if args.milestone in ("m0", "all") else {"available": False}

    print("\n=== Collecting M1 data ===")
    m1 = generate_m1_section() if args.milestone in ("m1", "all") else {"available": False}

    print("\n=== Collecting M2 data ===")
    m2 = generate_m2_section() if args.milestone in ("m2", "all") else {"available": False}

    print("\n=== Rendering report ===")
    html = render_report(m0, m1, m2)

    output_path = WORKDIR / "figures" / "QC" / "ibd_spatial_integration_report.html"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html)
    print(f"\nReport saved: {output_path}")
    print(f"Size: {output_path.stat().st_size / 1024:.0f} KB")
