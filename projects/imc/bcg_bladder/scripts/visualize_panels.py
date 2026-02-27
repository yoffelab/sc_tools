from skimage.exposure import adjust_gamma
from tqdm import tqdm
from glob import glob
import tifffile
import matplotlib.pyplot as plt
import pandas as pd

def add_scale_box_to_fig(
    img,
    ax,
    box_width: int = 100,
    box_height: float = 3,
    color: str = 'white'
):    
    import matplotlib.patches as patches
    x = img.shape[1]
    y = img.shape[0]
    
    # Create a Rectangle patch
    rect = patches.Rectangle((x - box_width, y * (1-box_height/100)), box_width, y * (box_height/100), linewidth=0.1, edgecolor='black', facecolor=color)
    
    
    # Add the patch to the Axes
    ax.add_patch(rect)
    return ax

def visualize(
    img_name: str,
    csv_name: str,
    out_dir: str = 'images'
) -> None:
    fig_dict = {
        'nrow': [5, 6, 6, 8],
        'ncol': [8, 8, 10, 10],
        'figsize': [(20,15), (20,15), (20,15), (25,20)],
    }
    
    img = tifffile.imread(img_name)
    roi = img_name.split('/')[-1].replace('_full.tiff', '')
    df = pd.read_csv(csv_name, index_col = 0)
    df['channel'] = df['channel'].str.replace('[0-9]{2,3}[A-Z][a-z]', '', regex = True)
    n_feature = df.shape[0]

    for i in range(4):
        if n_feature < fig_dict['nrow'][i] * fig_dict['ncol'][i]:
            break

    fig, axs = plt.subplots(
        fig_dict['nrow'][i],
        fig_dict['ncol'][i],
        figsize = fig_dict['figsize'][i],
        dpi = 300
    )
    
    for i, ax in enumerate(fig.axes):

        if i < len(df):
            rescaled_img = adjust_gamma(img[i], gamma = 0.2, gain = 1)

            ax.imshow(rescaled_img, cmap = 'viridis')

            ax.set_title(df.iloc[i]['channel'], fontsize = 10)
            add_scale_box_to_fig(rescaled_img, ax, box_width = 200)

            ax.annotate(
                text = 'min: {:.2f}\nmax: {:.2f}'.format(img[i].min(), img[i].max()),
                xy = (rescaled_img.shape[1], 0),
                ha = 'right',
                va = 'bottom',
                fontsize = 6
            )
            
            ax.annotate(
                text = '200μm',
                xy = (rescaled_img.shape[1], rescaled_img.shape[0]),
                ha = 'right',
                va = 'top',
                fontsize = 6
            )
        ax.axis('off')

    plt.suptitle(roi)
    plt.tight_layout()
    plt.savefig(f'{out_dir}/{roi}.pdf', bbox_inches = 'tight')
    plt.close()

tiff_files = glob('processed/*/tiffs/*_full.tiff')
csv_files = glob('processed/*/tiffs/*_full.csv')

for img_name, csv_name in tqdm(zip(tiff_files, csv_files), total = len(tiff_files)):
    visualize(img_name, csv_name)
