import sthd

# Read data
adata_raw = sthd.io.read_multi_sample("data/")
sthd.qc.calculate_qc_metrics(adata_raw)
sthd.pl.qc_overview(adata_raw, save_name="qc_pre_filter.pdf")
sthd.io.save_h5ad(adata_raw, "output/raw_data.h5ad")

# Filtering & normalization
adata_filtered = sthd.pp.filter_cells_and_genes(adata_raw)
sthd.pp.normalize_and_log(adata_filtered)
sthd.qc.calculate_qc_metrics(adata_filtered)
sthd.pl.qc_overview(adata_filtered, save_name="qc_post_filter.pdf")
sthd.pp.remove_qc_genes(adata_filtered)
sthd.pp.filter_low_total_counts_per_sample(adata_filtered, fraction = 0.1)
sthd.io.save_h5ad(adata_filtered, "output/adata_filtered.h5ad")

# HVG & SVG
# adata_filtered = sthd.pp.highly_variable_genes(adata_filtered, batch_key="sample")
# sthd.pp.spatially_variable_genes(adata_filtered)
adata_filtered = sthd.pp.filter_by_consensus_hvg_svg(adata_filtered)

method = "scvi"

adata_corrected = sthd.pp.batch_integrate(adata_filtered.copy(), batch_key="sample", method=method)
sthd.pl.umap(adata_corrected, color="sample", save_name=f"umap_{method}.pdf")
sthd.io.save_h5ad(adata_corrected, f"output/adata_corrected_{method}.h5ad")


