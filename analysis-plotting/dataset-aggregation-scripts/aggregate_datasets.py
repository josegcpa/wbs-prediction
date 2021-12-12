from glob import glob
import os
import h5py
import numpy as np
from tqdm import tqdm

try:
    os.makedirs('datasets')
except:
    pass

def generate_min_max(input_path):
    print("Calculating min and max for {}".format(input_path))
    m = []
    M = []
    for sub_dataset in tqdm(glob(os.path.join(input_path,'*h5'))):
        try:
            x = h5py.File(sub_dataset,'r')
            cells = x['cells']
            for C in cells:
                y = x['cells'][C][::]
                m.append(y.min(axis=0))
                M.append(y.max(axis=0))
        except:
            pass

    m,M = np.stack(m).min(axis=0),np.stack(M).max(axis=0)

    return m[np.newaxis,:],M[np.newaxis,:]

def generate_dataset(input_path,output_path,
                    m=None,M=None):
    with open(output_path,'w') as o:
        for sub_dataset in tqdm(glob(os.path.join(input_path,'*h5'))):
            name = os.path.split(sub_dataset)[-1].split('.')[0]
            x = h5py.File(sub_dataset,'r')
            cells = x['cells']
            all_cells = []
            for C in cells:
                y = x['cells'][C][::]
                # avoid numerical imprecision issues on the long run
                y = y[~(np.sum(np.abs(y) > 1e5,axis=1)>0),:]

                if m is not None:
                    y = y[(np.all(y >= m,axis=1)==True),:]
                if M is not None:
                    y = y[(np.all(y <= M,axis=1)==True),:]

                all_cells.append(y)
            all_cells = np.concatenate(all_cells,axis=0)
            mean = np.mean(all_cells,axis=0)
            variance = np.var(all_cells,axis=0)
            quantiles = np.quantile(
                all_cells,axis=0,q=[0.05,0.50,0.95])
            o.write('{},{},{}'.format(
                name,','.join([str(float(x)) for x in mean]),'mean\n'))
            o.write('{},{},{}'.format(
                name,','.join([str(float(x)) for x in variance]),'variance\n'))
            o.write('{},{},{}'.format(
                name,','.join([str(float(x)) for x in quantiles[0,:]]),'q05\n'))
            o.write('{},{},{}'.format(
                name,','.join([str(float(x)) for x in quantiles[1,:]]),'q95\n'))
            o.write('{},{},{}'.format(
                name,','.join([str(float(x)) for x in quantiles[2,:]]),'median\n'))
            x.close()

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

m_wbc,M_wbc = generate_min_max(dataset_folders[0])
m_rbc,M_rbc = generate_min_max(dataset_folders[1])

dataset_output = [x + '_summaries.csv' for x in dataset_output_root]

generate_dataset(dataset_folders[0],dataset_output[0],m_wbc,M_wbc)
generate_dataset(dataset_folders[1],dataset_output[1],m_rbc,M_rbc)
generate_dataset(dataset_folders[2],dataset_output[2],m_wbc,M_wbc)
generate_dataset(dataset_folders[3],dataset_output[3],m_rbc,M_rbc)
generate_dataset(dataset_folders[4],dataset_output[4],m_wbc,M_wbc)
generate_dataset(dataset_folders[5],dataset_output[5],m_rbc,M_rbc)