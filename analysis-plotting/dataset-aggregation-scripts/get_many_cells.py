import h5py
import numpy as np
import os
from glob import glob
from tqdm import tqdm

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--path',dest='path',action='store',type=str)

args = parser.parse_args()

all_hdf5 = glob(os.path.join(args.path,'*.h5'))

model_id = args.path.strip(os.sep).split(os.sep)[-1]

try: os.makedirs('datasets/many-cells')
except: pass

wbc_output = open('datasets/many-cells/wbc-{}.csv'.format(model_id),'w')
rbc_output = open('datasets/many-cells/rbc-{}.csv'.format(model_id),'w')

for hdf5_file in tqdm(all_hdf5):
    slide_name = os.path.split(hdf5_file)[-1].split('.')[0]
    cell_type = slide_name[:3]
    slide_name = slide_name[4:]
    model_name = hdf5_file.split(os.sep)[-2]
    with h5py.File(hdf5_file,'r') as F:
        cells = F['cells']
        cell_keys = list(cells.keys())
        key_subset = np.random.choice(cell_keys,50,replace=False)
        for k in key_subset:
            ct = np.argmax(cells[k]['cell_type'])
            features = ','.join([str(x) for x in cells[k]['features'][()]])
            O = '{},{},{},{},{}\n'.format(model_name,slide_name,cell_type,ct,features)
            if cell_type == 'wbc':
                wbc_output.write(O)
            else:
                rbc_output.write(O)
