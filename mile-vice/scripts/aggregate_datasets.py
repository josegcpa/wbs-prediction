from glob import glob
import os
import h5py
import numpy as np
from tqdm import tqdm 

try:
    os.makedirs('datasets')
except:
    pass

with open("scripts/wbc_feature_subset","r") as o:
    try:
        wbc_feature_subset = [
            int(x.strip())-1 for x in o.read().strip().split(',')]
    except:
        wbc_feature_subset = []
    if len(wbc_feature_subset) == 0:
        wbc_feature_subset = None

with open("scripts/rbc_feature_subset","r") as o:
    try:
        rbc_feature_subset = [
            int(x.strip())-1 for x in o.read().strip().split(',')]
    except:
        rbc_feature_subset = []
    if len(wbc_feature_subset) == 0:
        rbc_feature_subset = None

dataset_folders = [
    '/hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/_aggregates_wbc',
    '/hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/_aggregates_rbc',
    '/hps/nobackup/research/gerstung/josegcpa/data/ADDEN_NDPI/_aggregates_wbc',
    '/hps/nobackup/research/gerstung/josegcpa/data/ADDEN_NDPI/_aggregates_rbc',
    '/hps/nobackup/research/gerstung/josegcpa/data/ADDEN_SVS_results/_aggregates_wbc',
    '/hps/nobackup/research/gerstung/josegcpa/data/ADDEN_SVS_results/_aggregates_rbc']

dataset_output_root = [
    'datasets/wbc',
    'datasets/rbc',
    'datasets/wbc_adden_1',
    'datasets/rbc_adden_1',
    'datasets/wbc_adden_2',
    'datasets/rbc_adden_2']

dataset_output = [x + '.h5' for x in dataset_output_root]

for i,(dataset_folder,dataset) in enumerate(zip(dataset_folders,dataset_output)):
    print(dataset)
    with h5py.File(dataset,'w') as F:
        for sub_dataset in tqdm(glob(os.path.join(dataset_folder,'*h5'))):
            name = os.path.split(sub_dataset)[-1].split('.')[0]
            try:
                x = h5py.File(sub_dataset,'r')
                cells = x['cells']
                g = F.create_group(name)
                g_2 = g.create_group('cells')
                for C in cells:
                    y = x['cells'][C][::]
                    # avoid numerical imprecision issues on the long run
                    y = y[~(np.sum(np.abs(y) > 1e5,axis=1)>0),:]
                    if 'wbc' in dataset and wbc_feature_subset:
                        y = y[:,wbc_feature_subset]
                    g_2[C] = y
                g['mean'] = x['mean'][::]
                g['variance'] = x['variance'][::]
                x.close()
            except:
                pass

dataset_output = [x + '_subset.h5' for x in dataset_output_root]

for i,(dataset_folder,dataset) in enumerate(zip(dataset_folders,dataset_output)):
    print(dataset)
    with h5py.File(dataset,'w') as F:
        for sub_dataset in tqdm(glob(os.path.join(dataset_folder,'*h5'))):
            name = os.path.split(sub_dataset)[-1].split('.')[0]
            try:
                x = h5py.File(sub_dataset,'r')
                cells = x['cells']
                g = F.create_group(name)
                g_2 = g.create_group('cells')
                for C in cells:
                    y = x['cells'][C][::]
                    # avoid numerical imprecision issues on the long run
                    y = y[~(np.sum(np.abs(y) > 1e5,axis=1)>0),:]
                    if 'wbc' in dataset and wbc_feature_subset:
                        y = y[:,wbc_feature_subset]
                    if 'rbc' in dataset and wbc_feature_subset:
                        y = y[:,rbc_feature_subset]
                    g_2[C] = y
                g['mean'] = x['mean'][::]
                g['variance'] = x['variance'][::]
                x.close()
            except: 
                pass