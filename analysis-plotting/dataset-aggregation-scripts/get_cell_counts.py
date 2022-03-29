from glob import glob
import os
import h5py
import numpy as np
from tqdm import tqdm

try:
    os.makedirs('datasets')
except:
    pass

dataset_folders = [
    '/nfs/research/gerstung/josegcpa/data/SLIDES/MLL_TIFF/_aggregates_wbc',
    '/nfs/research/gerstung/josegcpa/data/SLIDES/MLL_TIFF/_aggregates_rbc',
    '/nfs/research/gerstung/josegcpa/data/SLIDES/ADDEN_NDPI/_aggregates_wbc',
    '/nfs/research/gerstung/josegcpa/data/SLIDES/ADDEN_NDPI/_aggregates_rbc',
    '/nfs/research/gerstung/josegcpa/data/SLIDES/ADDEN_SVS_results/_aggregates_wbc',
    '/nfs/research/gerstung/josegcpa/data/SLIDES/ADDEN_SVS_results/_aggregates_rbc']

dataset_output_root = [
    'datasets/wbc',
    'datasets/rbc',
    'datasets/wbc_adden_1',
    'datasets/rbc_adden_1',
    'datasets/wbc_adden_2',
    'datasets/rbc_adden_2']

def generate_dataset(input_path,output_path,
                    m=None,M=None,N=100):
    with open(output_path,'w') as o:
        for sub_dataset in tqdm(glob(os.path.join(input_path,'*h5'))):
            name = os.path.split(sub_dataset)[-1].split('.')[0]
            x = h5py.File(sub_dataset,'r')
            cells = x['cells']
            n_cells = 0
            for C in cells:
                y = x['cells'][C][::]
                # avoid numerical imprecision issues on the long run
                y = y[~(np.sum(np.abs(y) > 1e5,axis=1)>0),:]

                if m is not None:
                    y = y[(np.all(y >= m,axis=1)==True),:]
                if M is not None:
                    y = y[(np.all(y <= M,axis=1)==True),:]

                n_cells += y.shape[0]

            o.write('{},{}\n'.format(name,n_cells))

dataset_output = [x + '_counts.csv' for x in dataset_output_root]

generate_dataset(dataset_folders[0],dataset_output[0])
generate_dataset(dataset_folders[1],dataset_output[1])
generate_dataset(dataset_folders[2],dataset_output[2])
generate_dataset(dataset_folders[3],dataset_output[3])
generate_dataset(dataset_folders[4],dataset_output[4])
generate_dataset(dataset_folders[5],dataset_output[5])