import h5py
import sys
import numpy as np
from tqdm import tqdm

# defines outliers for MILe-ViCe datasets

H5_FILE = sys.argv[1]

all_cells = []
with h5py.File(H5_FILE,'r') as F:
    all_keys = list(F.keys())
    with tqdm(all_keys,postfix='') as pbar:
        for k in all_keys:
            pbar.postfix = k
            pbar.update()
            for cell_idx in F[k]['cells']:
                all_cells.append(F[k]['cells'][cell_idx][::])

all_cells = np.concatenate(all_cells,axis=0)
q1 = np.quantile(all_cells,0.25,axis=0)
q3 = np.quantile(all_cells,0.75,axis=0)
m = np.min(all_cells,axis=0)
M = np.max(all_cells,axis=0)
print(','.join([str(x) for x in q1]))
print(','.join([str(x) for x in q3]))
print(','.join([str(x) for x in m]))
print(','.join([str(x) for x in M]))
