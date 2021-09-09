import h5py
import sys
import numpy as np
from tqdm import tqdm

H5_FILE = sys.argv[1]
OUTLIERS_FILE = sys.argv[2]

with open(OUTLIERS_FILE) as o:
    lines = [x.strip() for x in o.readlines()]

q1,q3,m,M = [np.array([float(y) for y in x.split(',')])
             for x in lines]
iqr = q3 - q1

with h5py.File(H5_FILE,'r') as F:
    all_keys = list(F.keys())
    with tqdm(all_keys,postfix='') as pbar:
        for k in all_keys:
            pbar.postfix = k
            pbar.update()
            all_cells = []
            for cell_idx in F[k]['cells']:
                all_cells.append(F[k]['cells'][cell_idx][::])
            all_cells = np.concatenate(all_cells,axis=0)
            all_cells_ = all_cells[np.all(all_cells >= (m),axis=1),]
            all_cells_ = all_cells_[np.all(all_cells_ <= (M),axis=1),]
            features = {
                'mean':np.mean(all_cells_,axis=0),
                'variance':np.var(all_cells_,axis=0),
                'q05':np.quantile(all_cells,0.05,axis=0),
                'q95':np.quantile(all_cells,0.95,axis=0),
                'median':np.quantile(all_cells,0.5,axis=0)
            }
            for n in features:
                line = features[n]
                print(k+','+','.join([str(x) for x in line])+','+n)
