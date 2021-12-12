import h5py
import sys
import numpy as np
from tqdm import tqdm

H5_FILE = sys.argv[1]
OUTPUT_FILE = sys.argv[2]
SUBSET = int(sys.argv[3])

O = open(OUTPUT_FILE,'w')
all_cell_sizes = []
with h5py.File(H5_FILE,'r') as F:
    keys = F.keys()
    for k in tqdm(keys):
        total_possibilities = []
        for cell_group in F[k]['cells']:
            for row in range(F[k]['cells'][cell_group].shape[0]):
                total_possibilities.append([cell_group,row])
        if len(total_possibilities) < SUBSET:
            s = [x for x in range(len(total_possibilities))]
        else:
            s = np.random.choice(len(total_possibilities),
                                 replace=False,size=SUBSET)

        for i in s:
            cell_group,row = total_possibilities[i]
            line = F[k]['cells'][cell_group][row,:]
            O.write(k + ',' + ','.join([str(x) for x in line]) + '\n')
