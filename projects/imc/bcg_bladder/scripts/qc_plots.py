import os
from glob import glob
from pathlib import Path
import tifffile
from tqdm import tqdm
import yaml

import pandas as pd
import numpy as np

import scanpy as sc
import anndata

from skimage.exposure import adjust_gamma
from skimage import filters
import scipy.ndimage as ndi
import scipy

import seaborn as sns
import matplotlib.pyplot as plt

'''
USER DEFINED INPUTS
'''

sc.settings.set_figure_params(dpi=150, dpi_save=300, fontsize=12)

metadata_filename = 'metadata/bcg_config.yml'

with open(metadata_filename, "r") as stream:
    try:
        metadata = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# var_names = metadata['var_celltype_groups']

IMAGE_FILE_PATTERN = metadata['IMAGE_FILE_PATTERN']
MASK_FILE_PATTERN = metadata['MASK_FILE_PATTERN']
NUCMASK_FILE_PATTERN = metadata['NUCMASK_FILE_PATTERN']
CSV_FILE_PATERN = metadata['CSV_FILE_PATERN']

SINGLE_CELL_FILE = metadata['adata_labeled_file_name']

ROI_AREA_FILE = metadata['ROI_AREA_FILE']
ROI_MEAN_FILE = metadata['ROI_MEAN_FILE']
ROI_VAR_FILE = metadata['ROI_VAR_FILE']

ROI_NUC_FILE = metadata['ROI_NUC_FILE']
ROI_CYTO_FILE = metadata['ROI_CYTO_FILE']

ROI_SPILLOVER_FILE = metadata['ROI_SPILLOVER_FILE']

UMAP_FRACTION = metadata['UMAP_FRACTION']

fig_dict = {
    'nrow': [5, 6, 6, 8],
    'ncol': [8, 8, 10, 10],
    'figsize': [(25,15), (25,15), (25,15), (30,20)],
}


'''
FUNCTIONS
'''

def get_img_size(cell_mask_file_name):
    cell_mask = tifffile.imread(cell_mask_file_name)
    cell_area = cell_mask > 0
    
    return cell_area.sum() / 1e6

def get_marker_mean(image_file_name):
    img = tifffile.imread(image_file_name)
    return img.mean(axis = (1,2))

def get_marker_var(image_file_name):
    img = tifffile.imread(image_file_name)
    return img.var(axis = (1,2))

def get_TMA_size(cell_mask_file_name):
    cell_mask = tifffile.imread(cell_mask_file_name)
    smooth = filters.gaussian(cell_mask > 0, sigma=5)

    thresh_value = filters.threshold_triangle(cell_mask > 0)
    thresh = smooth > thresh_value

    fill = ndi.binary_fill_holes(thresh)
    
    plt.imshow(fill)
    plt.axis('off')
    plt.show()
    print(fill.shape)
    
    return fill.sum() / 1e6

def get_image_marker_mean_nuc_cyto(image_file_name, cell_mask_file, nuc_mask_file):
    image = tifffile.imread(image_file_name)
    nuc_mask = tifffile.imread(nuc_mask_file) > 0
    cyto_mask = np.bitwise_xor(tifffile.imread(cell_mask_file) > 0, nuc_mask)

    nuc_image = image * nuc_mask[np.newaxis, :, :]
    cyto_image = image * cyto_mask[np.newaxis, :, :]
    
    mean_nuc = nuc_image.mean(axis = (1,2))
    mean_cyto = cyto_image.mean(axis = (1,2))
    
    return mean_nuc, mean_cyto

def get_plot_subfig_param(
    num_var: int,
    fig_dict:dict = fig_dict
):

    for i in range(len(fig_dict['nrow'])):
        if num_var <= fig_dict['nrow'][i] * fig_dict['ncol'][i]:
            break

    return fig_dict['nrow'][i], fig_dict['ncol'][i], fig_dict['figsize'][i]

def roi_cellcount_plot(
    adata: anndata.AnnData,
    area_df: pd.DataFrame,
    title: str = 'ROI Cell Count Information',
    outdir: Path = 'figures',
    filename: Path = 'QC_count_plot.pdf'
):
    df = pd.DataFrame(adata.obs.groupby('roi').count()['sample'])
    df['cell count'] = df['sample']
    df = df.join(area_df.set_index('roi'))
    df['area (mm$^2$)'] = df['ROI_area']
    df['cell density (mm$^{-2}$)'] = df['cell count'] / df['area (mm$^2$)']

    del df['sample']
    del df['ROI_area']

    box_df = df.melt()

    g = sns.FacetGrid(
        box_df.sort_values('variable'),
        col = 'variable',
        sharey = False,
        aspect = 0.7,
        margin_titles = False)

    #g.map(sns.boxplot, 'variable', 'value', fliersize = 0)
    #g.map(sns.swarmplot, 'variable', 'value', color = 'black', size = 2)
    
    g.map(sns.violinplot, 'variable', 'value', color = 'grey', inner = 'box')
    g.set_axis_labels("", "")
    #g.set_axis_ticks("")
    g.set_titles(col_template="{col_name}", row_template="")
    g.set_xticklabels([])
    plt.suptitle(title)
    plt.savefig(Path(outdir) / Path(filename), bbox_inches = 'tight')
    plt.close()
    return box_df

def mv_ratio_plot(
    mean_df: pd.DataFrame,
    var_df: pd.DataFrame,
    title: str = 'log Var (Y) / log Mean (X) Ratio',
    outdir: Path = 'figures',
    filename: Path = 'QC_MA_plot.pdf'
) -> None:
    print('Plotting Mean Variance ratio per ROI')
    nrow, ncol, figsize = get_plot_subfig_param(num_var = mean_df.shape[1])
    fig, axs = plt.subplots(nrow, ncol, figsize = figsize, dpi = 300)

    for i, ax in enumerate(fig.axes):
        if i < len(var_df.columns):
            col = var_df.columns[i]
            sns.scatterplot(x = mean_df[col], y = var_df[col], ax = ax, s = 10, color = 'black', alpha = 0.5)
            ax.set_title(col)
            
            ax.set_xlabel('mean')
            ax.set_ylabel('variance')
            ax.set_yscale('log')
            ax.set_xscale('log')
            
            lims = np.array([
                np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
                np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
            ])
            
            ax.plot(lims, lims, 'k--', alpha=0.3, zorder=0)
            ax.plot(lims, lims*10, 'k--', alpha=0.3, zorder=0)
            ax.plot(lims, lims*100, 'k--', alpha=0.3, zorder=0)
            
        else:
            ax.axis('off')

    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(Path(outdir) / Path(filename))
    plt.close()

def epitope_evaluation_plot(
    mean_df: pd.DataFrame,
    var_df: pd.DataFrame,
    title: str = 'Epitope QC',
    outdir: Path = 'figures',
    filename: Path = 'QC_marker_specificity.pdf',
    title_fontsize: int = 14,
    text_fontsize: int = 12
):
    mean = mean_df.melt(ignore_index = False).reset_index().rename(columns = {'value':'mean'})
    var = var_df.melt(ignore_index = False).reset_index().rename(columns = {'value':'var'})
    mv = mean.merge(var, on = ['index', 'variable'])
    mv['log mean'] = np.log(mv['mean']) + 1
    mv['log var'] = np.log(mv['var']) + 1

    from scipy.stats import linregress

    grouped = mv.groupby('variable')
    result_df = pd.DataFrame(
        grouped.apply(lambda x: linregress(x['log mean'], x['log var']))
    ).reset_index()
    result_df[['slope', 'intercept', 'r_value', 'p_value', 'std_err']] = result_df[0].apply(pd.Series)
    del result_df[0]

    # result_df['Region Specificity Score'] = np.log(
    #     np.e**-1 + (result_df['slope'] - result_df['slope'].min()) / (result_df['slope'].max() - result_df['slope'].min())
    # )
    # result_df['General Abundance'] = np.log(
    #     np.e**-1 + (result_df['intercept'] - result_df['intercept'].min()) / (result_df['intercept'].max() - result_df['intercept'].min())
    # )

    result_df['Region Specificity Score'] = result_df['slope']
    # result_df['General Abundance'] = result_df['intercept']
    log_mean = mv.groupby('variable').mean()['log mean']
    log_mean = (log_mean - log_mean.min() + 1).tolist()
    result_df['General Abundance'] = log_mean

    plt.figure(figsize = (11,8), dpi = 300)
    sns.scatterplot(data = result_df, x = 'Region Specificity Score', y = 'General Abundance', color = 'red')
    for m in result_df.itertuples():
        plt.text(x = m[7], y = m[8], s = m[1], ha = 'left', va = 'baseline', fontsize = text_fontsize)

    #plt.yscale('log')
    plt.title(title, fontsize = title_fontsize)
    plt.savefig(Path(outdir) / Path(filename), bbox_inches = 'tight')
    plt.close()

    return mv, result_df

def mean_intensity_plot(
    mean_df: pd.DataFrame,
    log: bool = True,
    title: str = 'Marker Mean Intensity across ROI',
    outdir: Path = 'figures',
    filename: Path = 'QC_image_wide_marker_mean.pdf'
):
    if log:
        mean_df = pd.DataFrame(np.log(mean_df))
    # Mean Intensity Plot
    sns.clustermap(
        mean_df,
        cmap = 'Spectral_r',
        xticklabels = True,
        figsize = (15,10),
        vmax = 5,
        vmin = -5,
        cbar_kws = {"shrink": 0.5}
    )
    plt.title(title)
    plt.savefig(f'{outdir}/{filename}')
    plt.close()


def calculate_spillover(
    adata: anndata.AnnData,
    outdir: Path = 'metadata',
    filename: Path = 'ROI_spillover.csv'
) -> pd.DataFrame:
    ROI_SPILLOVER_PATH = Path('metadata') / Path(filename)
    print('Calculating Spillover Score...')
    if len(adata) > 100000:
        data = sc.pp.subsample(adata, n_obs = 100000, copy = True)
    else:
        data = adata

    if data.raw:
        feat_ = data.raw.var.index
        data = data.raw.X
    else:
        feat_ = data.var.index
        data = data.X
    
    spillover_df = scipy.stats.spearmanr(data)[0]
    #spillover_df = np.corrcoef(data)
    spillover_df = pd.DataFrame(spillover_df, index = feat_, columns = feat_)

    feat_ = spillover_df.columns.str.contains('(', regex = False) & spillover_df.columns.str.contains(')', regex = False)
    spillover_df = spillover_df.loc[feat_, feat_]
    spillover_df.to_csv(ROI_SPILLOVER_PATH)

    return spillover_df

def spillover_plot(
    spillover_df: pd.DataFrame,
    outdir: Path = 'figures',
    matrixplot_title: str = 'Spillover Score Matrix',
    matrixplot_filename: Path = 'QC_spillover_matrix.pdf',
    graph_title: str = 'QC Spillover Graph',
    graph_filename: Path = 'QC_spillover_graph.pdf',
):
    fig, ax = plt.subplots(figsize = (15, 13), dpi = 300)
    mask = np.triu(np.ones_like(spillover_df, ))
    sns.heatmap(
        spillover_df,
        mask = mask,
        cmap = 'Spectral_r',
        vmin = 0.2,
        vmax = 0.8,
        xticklabels = True,
        yticklabels = True,
        square = True,
        linewidths = 0.5,
        cbar_kws = {"shrink": 0.5}
    )

    import matplotlib.patches as patches
    rect = pd.DataFrame(
        np.eye(spillover_df.shape[0], dtype = bool),
        index = spillover_df.index,
        columns = spillover_df.columns
    )

    for i in range(1, spillover_df.shape[0]):
        for j in reversed(range(i)):
            if rect.iloc[i-1,j] and rect.iloc[i,j+1] and spillover_df.iloc[i,j] > 0.5:
                rect.iloc[i,j] = True
                ax.add_patch(
                    patches.Rectangle(
                        (j, i), 1.0, 1.0,
                        edgecolor='black',
                        fill=False,
                        lw=2
                    )
                )

    plt.title(matrixplot_title)
    plt.tight_layout()
    plt.savefig(f'{outdir}/{matrixplot_filename}')
    plt.close()


    fig, ax = plt.subplots(figsize = (9, 10), dpi = 300)
    import networkx as nx
    np.fill_diagonal(rect.values, False)
    adj = spillover_df * rect
    adj[adj.abs()<0.5] = 0
    spillover_graph = nx.from_pandas_adjacency(
        adj.fillna(0)
    )
    spillover_graph.remove_nodes_from(list(nx.isolates(spillover_graph)))
    pos = nx.nx_pydot.graphviz_layout(spillover_graph)
    edges,weights = zip(*nx.get_edge_attributes(spillover_graph,'weight').items())
    nx.draw(
        spillover_graph,
        with_labels=True,
        node_color = 'lightgray',
        #edge_color = 'lightblue',
        edgelist=edges,
        edge_color=weights,
        edge_cmap= plt.cm.Spectral_r,
        edge_vmin = 0.2,
        edge_vmax = 0.8,
        pos = pos,
        ax = ax)
    plt.title(graph_title)
    plt.tight_layout()
    plt.savefig(f'{outdir}/{graph_filename}')
    plt.close()


def nuc_cyto_ratio_plot(
    nuc_df: pd.DataFrame,
    cyto_df: pd.DataFrame,
    outdir: Path = 'figures',
    filename1: Path = 'QC_NUC_CYTO_ratio.pdf',
    filename2: Path = 'QC_NUC_CYTO_ratio_barplot.pdf'
):

    # Nucleus Cytoplasm ratio violinplot
    nuc_df_melt = nuc_df.melt(ignore_index = False)
    cyto_df_melt = cyto_df.melt(ignore_index = False)
    nuc_df_melt['compartment'] = 'nucleus'
    cyto_df_melt['compartment'] = 'cytoplasm'

    nuc_cyto_df = pd.concat([nuc_df_melt, cyto_df_melt])
    nuc_cyto_df['log expression'] = np.log(nuc_cyto_df['value'] + 1e-5)
    #print(nuc_cyto_df)
    fig, ax = plt.subplots(1,1,figsize=(15,5), dpi = 100)
    sns.violinplot(
        x="variable",
        y="log expression",
        hue = 'compartment',
        data=nuc_cyto_df,
        palette=['#0000FF','#00FF00'],
        ax = ax,
        split=True
    ) 
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
    plt.savefig(f'{outdir}/{filename1}', bbox_inches = 'tight')
    plt.close()

    
    # Nucleus Cytoplasm Ratio Plot
    nc_ratio = np.log(nuc_df) - np.log(cyto_df)
    nc_ratio = nc_ratio.loc[:,nc_ratio.columns.str.contains('DNA') | nc_ratio.columns.str.contains('HH')]
    features = nc_ratio.columns
    nc_ratio['order'] = nc_ratio.sum(axis = 1)
    nc_ratio = nc_ratio.reset_index()
    nc_ratio['index'] = pd.Categorical(nc_ratio['index'], categories = nc_ratio.sort_values('order')['index'])


    fig, ax = plt.subplots(
        len(features),
        1,
        figsize = (len(nc_ratio) / 5, len(features) * 3),
        dpi = 300
    )

    for i, f in enumerate(features):
        sns.barplot(
            data = nc_ratio,
            x = 'index',
            y = f,
            palette = sns.diverging_palette(220, 20, n = len(nc_ratio)),
            ax = ax[i]
        )
        ax[i].set_xlabel('')

        if i < len(features) - 1:
            ax[i].set_xticklabels([])
        else:
            ax[i].set_xticklabels(ax[i].get_xticklabels(), rotation = 90)
    plt.suptitle('Cell Segmentation Score Plot')
    plt.tight_layout()
    plt.savefig(f'{outdir}/{filename2}', bbox_inches = 'tight')
    plt.close()
    

'''
BASIC IO
'''
image_files = sorted(glob(IMAGE_FILE_PATTERN))
mask_files = sorted(glob(MASK_FILE_PATTERN))
nucmask_files = sorted(glob(NUCMASK_FILE_PATTERN))
csv_files = sorted(glob(CSV_FILE_PATERN))

image_files_series = pd.Series(
    image_files,
    index = [s.split('/')[-1].split('_full')[0] for s in image_files],
    name = 'image_dir'
)

mask_files_series = pd.Series(
    mask_files,
    index = [s.split('/')[-1].split('_full')[0] for s in mask_files],
    name = 'mask_dir'
)

nucmask_files_series = pd.Series(
    nucmask_files,
    index = [s.split('/')[-1].split('_full')[0] for s in nucmask_files],
    name = 'nucmask_dir'
)

csv_files_series = pd.Series(
    csv_files,
    index = [s.split('/')[-1].split('_full')[0] for s in csv_files],
    name = 'csv_dir'
)

files = pd.concat(
    [
        image_files_series,
        mask_files_series,
        nucmask_files_series,
        csv_files_series
    ],
    axis = 1
)


filtered = files[files.isna().any(axis=1)].index.tolist()
print(f'Filtering following ROIs due to missing data:\n{filtered}')
files = files.dropna()


adata = sc.read(SINGLE_CELL_FILE)

ROI_AREA_PATH = Path(ROI_AREA_FILE)
if ROI_AREA_FILE.replace('metadata/','') in os.listdir('metadata'):
    area_df = pd.read_csv(ROI_AREA_PATH, index_col = 0)
else:
    print('Calculating ROI Area...')
    TMA_area = dict()
    for file in tqdm(files['mask_dir'].tolist()):
        area = get_img_size(file)
        key = file.split('/')[-1].replace('_full_mask.tiff','')
        TMA_area[key] = area

    area_df = pd.DataFrame.from_records([TMA_area]).T.reset_index()
    area_df.columns = ['roi', 'ROI_area']
    area_df.to_csv(ROI_AREA_PATH)


ROI_MEAN_PATH = Path(ROI_MEAN_FILE)
if ROI_MEAN_FILE.replace('metadata/','') in os.listdir('metadata'):
    mean_df = pd.read_csv(ROI_MEAN_PATH, index_col = 0)
else:
    print('Calculating ROI mean...')
    TMA_mean = dict()
    for file in tqdm(files['image_dir'].tolist()):
        m = get_marker_mean(file)
        key = file.split('/')[-1].replace('_full.tiff','')
        TMA_mean[key] = m
    
    matrix = np.array(list(TMA_mean.values()))
    index = list(TMA_mean.keys())
    channels = pd.read_csv(csv_files[0], index_col = 0)['channel'].tolist()
    mean_df = pd.DataFrame(matrix, index = index, columns = channels)
    mean_df.to_csv(ROI_MEAN_PATH)


ROI_VAR_PATH = Path(ROI_VAR_FILE)
if ROI_VAR_FILE.replace('metadata/','') in os.listdir('metadata'):
    var_df = pd.read_csv(ROI_VAR_PATH, index_col = 0)
else:
    print('Calculating ROI variance...')
    TMA_var = dict()
    for file in tqdm(files['image_dir'].tolist()):
        v = get_marker_var(file)
        key = file.split('/')[-1].replace('_full.tiff','')
        TMA_var[key] = v
    
    matrix = np.array(list(TMA_var.values()))
    index = list(TMA_var.keys())
    channels = pd.read_csv(csv_files[0], index_col = 0)['channel'].tolist()
    var_df = pd.DataFrame(matrix, index = index, columns = channels)
    var_df.to_csv(ROI_VAR_PATH)


ROI_NUC_PATH = Path(ROI_NUC_FILE)
ROI_CYTO_PATH = Path(ROI_CYTO_FILE)
if ROI_NUC_FILE.replace('metadata/','') and ROI_CYTO_FILE.replace('metadata/','') in os.listdir('metadata'):
    nuc_df = pd.read_csv(ROI_NUC_PATH, index_col = 0)
    cyto_df = pd.read_csv(ROI_CYTO_PATH, index_col = 0)
else:
    print('Calculating Nucleus Cytoplasm Ratio')
    TMA_nuc_mean = dict()
    TMA_cyto_mean = dict()

    for file in tqdm(files.itertuples(), total = files.shape[0]):
        image_file, mask_file, nucmask_file, csv_file = file[1], file[2], file[3], file[4]
        nuc_mu, cyto_mu = get_image_marker_mean_nuc_cyto(image_file, mask_file, nucmask_file)
        
        key = image_file.split('/')[-1].replace('_full.tiff','')
        TMA_nuc_mean[key] = nuc_mu
        TMA_cyto_mean[key] = cyto_mu

    df = pd.read_csv(csv_file, index_col = 0)

    nuc_df = pd.DataFrame.from_records(TMA_nuc_mean, index = df['channel']).T
    cyto_df = pd.DataFrame.from_records(TMA_cyto_mean, index = df['channel']).T
    nuc_df.to_csv(ROI_NUC_PATH)
    cyto_df.to_csv(ROI_CYTO_PATH)


adata.var.index = adata.var.index.str.replace('[0-9]{2,3}[A-Z][a-z]', '', regex = True)
adata.raw.var.index = adata.raw.var.index.str.replace('[0-9]{2,3}[A-Z][a-z]', '', regex = True)
mean_df.columns = mean_df.columns.str.replace('[0-9]{2,3}[A-Z][a-z]', '', regex = True)
var_df.columns = var_df.columns.str.replace('[0-9]{2,3}[A-Z][a-z]', '', regex = True)
nuc_df.columns = nuc_df.columns.str.replace('[0-9]{2,3}[A-Z][a-z]', '', regex = True)
cyto_df.columns = cyto_df.columns.str.replace('[0-9]{2,3}[A-Z][a-z]', '', regex = True)

# Facet plot settings


'''
QC Plotting
'''

# plot basic information
roi_cellcount_plot(adata, area_df)

# Mean Intensity Plot
mean_intensity_plot(mean_df)

# Mean Variance Ratio plot
mv_ratio_plot(mean_df, var_df)

# Epitopte Eval Plot
epitope_evaluation_plot(mean_df, var_df)

# Quantify Spillover
ROI_SPILLOVER_PATH = Path(ROI_SPILLOVER_FILE)
if ROI_SPILLOVER_FILE in os.listdir('metadata'):
    spillover_df = pd.read_csv(ROI_SPILLOVER_PATH, index_col = 0)
else:
    spillover_df = calculate_spillover(adata)

spillover_plot(spillover_df)

nuc_cyto_ratio_plot(nuc_df, cyto_df)

# Plot UMAP Sample

# shuffle
umap_data = sc.pp.subsample(adata, fraction = UMAP_FRACTION, copy = True)
sc.pl.umap(
    umap_data,
    color = 'sample',
    size = 1,
    frameon = False,
    show = False,
    title = f'UMAP by samples\nPlotting {len(umap_data)} cells, {UMAP_FRACTION * 100}% of all cells'
)
filename = SINGLE_CELL_FILE.split('/')[-1].replace('.h5ad','')
plt.savefig(
    f'figures/QC_UMAP_sample_{filename}.pdf',
    bbox_inches = 'tight'
)
plt.close()

# Plot UMAP Marker
print('Plotting marker UMAP raw value')

nrow, ncol, figsize = get_plot_subfig_param(num_var = umap_data.shape[1])
fig, axs = plt.subplots(nrow, ncol, figsize = figsize, dpi = 300)

for i, ax in tqdm(enumerate(fig.axes)):
    if i < len(umap_data.var.index):
        var = umap_data.var.index[i]
        sc.pl.umap(
            umap_data,
            color = var,
            use_raw = True,
            size = 1,
            frameon = False,
            ax = ax,
            show = False,
            colorbar_loc = None,
        )
        
    else:
        ax.axis('off')
        
plt.suptitle(f'Cell Marker Raw Expression UMAP\nPlotting {len(umap_data)} cells, {UMAP_FRACTION * 100}% of all cells')
plt.tight_layout()
filename = SINGLE_CELL_FILE.split('/')[-1].replace('.h5ad','')
plt.savefig(
    f'figures/QC_umap_marker_raw_{filename}.pdf',
    bbox_inches = 'tight'
)
plt.close()


# Plot UMAP Marker
print('Plotting marker UMAP')
fig, axs = plt.subplots(nrow, ncol, figsize = figsize, dpi = 300)

for i, ax in tqdm(enumerate(fig.axes)):
    if i < len(umap_data.var.index):
        var = umap_data.var.index[i]
        sc.pl.umap(
            umap_data,
            color = var,
            use_raw = False,
            size = 1,
            frameon = False,
            ax = ax,
            show = False,
            colorbar_loc = None,
            vmin = 0,
            vmax = 3
        )
        
    else:
        ax.axis('off')
        
plt.suptitle(f'Cell Marker Scaled Expression UMAP\nPlotting {len(umap_data)} cells, {UMAP_FRACTION * 100}% of all cells')
plt.tight_layout()
filename = SINGLE_CELL_FILE.split('/')[-1].replace('.h5ad','')
plt.savefig(
    f'figures/QC_umap_marker_{filename}.pdf',
    bbox_inches = 'tight'
)
plt.close()

# Plot Raw Matrixplot
for c in adata.obs.columns[adata.obs.columns.str.contains('cluster_')]:
    mp = sc.pl.matrixplot(
        adata,
        groupby = c,
        var_names = var_names,
        log = True,
        cmap = 'Spectral_r',
        vmax = 1.5,
        return_fig = True
    )
    mp.add_totals().style(cmap = 'Spectral_r', edge_color='black')
    mp.savefig(f'figures/QC_matrixplot_raw_{c}.pdf')
    plt.close()


# Plot Normalized Matrixplot
for c in adata.obs.columns[adata.obs.columns.str.contains('cluster_')]:
    mp = sc.pl.matrixplot(
        adata,
        groupby = c,
        var_names = var_names,
        log = True,
        cmap = 'Spectral_r',
        standard_scale = 'var',
        return_fig = True
    )
    mp.add_totals().style(cmap = 'Spectral_r', edge_color='black')
    mp.savefig(f'figures/QC_matrixplot_normalized_{c}.pdf')
    plt.close()
