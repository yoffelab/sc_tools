from glob import glob
import pandas as pd
import shutil
import os
import tifffile

mask_types = ['FullStack', 'Vessel', 'Epithleial', 'Nucleus', 'Cytoplasm', 'CellMask']
d = dict()
for mask_type in mask_types:
    files = sorted(glob(f'Images/*_{mask_type}*'))
    d[mask_type] = files

df = pd.DataFrame.from_records(d)
df = df.sample(frac = 1, random_state = 42)

train = df.iloc[:600,:]
test = df.iloc[600:,:]

def copy_record(rec, target_dir):
    for ind in rec.index:
        file = rec.loc[ind]
        file_target = rec.loc[ind].replace('Images/', '').split(ind)[0][:-1]

        if ind == 'FullStack':
            dest_fpath = f'{target_dir}/Input/{file_target}.tiff'
            os.makedirs(os.path.dirname(dest_fpath), exist_ok=True)
            #shutil.copy(file, dest_fpath)
            img = tifffile.imread(file)
            tifffile.imwrite(dest_fpath, img.transpose(1,2,0))
        else:
            dest_fpath = f'{target_dir}/Target/{ind}/{file_target}.tiff'
            os.makedirs(os.path.dirname(dest_fpath), exist_ok=True)
            #shutil.copy(file, dest_fpath)
            img = tifffile.imread(file)
            tifffile.imwrite(dest_fpath, img > 0)


train.apply(copy_record, target_dir = 'BRCA-danenberg', axis = 1)
test.apply(copy_record, target_dir = 'BRCA-danenberg/Test', axis = 1)
