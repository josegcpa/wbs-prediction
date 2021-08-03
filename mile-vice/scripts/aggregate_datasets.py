from glob import glob
import os
import h5py
import numpy as np

try:
    os.makedirs('datasets')
except:
    pass

dataset_folders = [
    '/hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/_aggregates_wbc',
    '/hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/_aggregates_rbc'
]

dataset_output = [
    'datasets/wbc.h5',
    'datasets/rbc.h5'
]

for dataset_folder,dataset in zip(dataset_folders,dataset_output):
    with h5py.File(dataset,'w') as F:
        for sub_dataset in glob(os.path.join(dataset_folder,'*h5')):
            name = os.path.split(sub_dataset)[-1].split('.')[0]
            print(dataset,name)
            x = h5py.File(sub_dataset,'r')
            cells = x['cells']
            g = F.create_group(name)
            g_2 = g.create_group('cells')
            for C in cells:
                y = x['cells'][C][::]
                # avoid numerical imprecision issues on the long run
                y = y[~(np.sum(np.abs(y) > 1e5,axis=1)>0),:]
                #
                g_2[C] = y
            g['mean'] = x['mean'][::]
            g['variance'] = x['variance'][::]
