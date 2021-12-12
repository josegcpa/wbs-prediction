from glob import glob
import os
import h5py
import numpy as np
from tqdm import tqdm

try:
    os.makedirs('datasets')
except:
    pass

def get_feature_subset(file_path):
    try:
        with open(file_path,"r") as o:
            feature_subset = [
                int(x.strip())-1 for x in o.read().strip().split(',')]
    except:
            feature_subset = []
    if len(feature_subset) == 0:
        feature_subset = None
    return feature_subset

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
                # avoid numerical imprecision issues on the long run
                m.append(y.min(axis=0))
                M.append(y.max(axis=0))
        except:
            pass

    m,M = np.stack(m).min(axis=0),np.stack(M).max(axis=0)

    return m[np.newaxis,:],M[np.newaxis,:]

def generate_dataset(input_path,output_path,
                     m=None,M=None,
                     feature_subset=None):
    print("Generating dataset for {} in {}\n(feature subset: {})".format(
        input_path,output_path,feature_subset
    ))
    with h5py.File(output_path,'w') as F:
        for sub_dataset in tqdm(glob(os.path.join(input_path,'*h5'))):
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

                    if m is not None:
                        y = y[(np.all(y >= m,axis=1)==True),:]
                    if M is not None:
                        y = y[(np.all(y <= M,axis=1)==True),:]

                    if feature_subset:
                        y = y[:,feature_subset]
                    g_2[C] = y
                mean = x['mean'][::]
                variance = x['variance'][::]
                if feature_subset:
                    mean = mean[feature_subset]
                    variance = variance[feature_subset]
                g['mean'] = mean
                g['variance'] = variance
                x.close()
            except:
                pass

def generate_dataset(input_paths,output_path,
                     m=None,M=None,
                     feature_subset=None):
    print("Generating dataset for {} in {}\n(feature subset: {})".format(
        input_paths,output_path,feature_subset
    ))
    with h5py.File(output_path,'w') as F:
        for input_path in input_paths:
            for sub_dataset in tqdm(glob(os.path.join(input_path,'*h5'))):
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

                        if m is not None:
                            y = y[(np.all(y >= m,axis=1)==True),:]
                        if M is not None:
                            y = y[(np.all(y <= M,axis=1)==True),:]
                        if feature_subset:
                            y = y[:,feature_subset]
                        g_2[C] = y
                    mean = x['mean'][::]
                    variance = x['variance'][::]
                    if feature_subset:
                        mean = mean[feature_subset]
                        variance = variance[feature_subset]
                    g['mean'] = mean
                    g['variance'] = variance
                    x.close()
                except:
                    pass

wbc_feature_subset = get_feature_subset(
    "scripts/wbc_feature_subset")
rbc_feature_subset = get_feature_subset(
    "scripts/rbc_feature_subset")
wbc_feature_subset_unbiased = get_feature_subset(
    "scripts/wbc_feature_subset_unbiased")
rbc_feature_subset_unbiased = get_feature_subset(
    "scripts/rbc_feature_subset_unbiased")

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

dataset_output = [x + '.h5' for x in dataset_output_root]

generate_dataset([dataset_folders[0]],dataset_output[0],m_wbc,M_wbc)
generate_dataset([dataset_folders[1]],dataset_output[1],m_rbc,M_rbc)
generate_dataset([dataset_folders[2]],dataset_output[2],m_wbc,M_wbc)
generate_dataset([dataset_folders[3]],dataset_output[3],m_rbc,M_rbc)
generate_dataset([dataset_folders[4]],dataset_output[4],m_wbc,M_wbc)
generate_dataset([dataset_folders[5]],dataset_output[5],m_rbc,M_rbc)

dataset_output = [x + '_subset.h5' for x in dataset_output_root]

generate_dataset([dataset_folders[0]],dataset_output[0],
                 m_wbc,M_wbc,wbc_feature_subset)
generate_dataset([dataset_folders[1]],dataset_output[1],
                 m_rbc,M_rbc,rbc_feature_subset)
generate_dataset([dataset_folders[2]],dataset_output[2],
                 m_wbc,M_wbc,wbc_feature_subset)
generate_dataset([dataset_folders[3]],dataset_output[3],
                 m_rbc,M_rbc,rbc_feature_subset)
generate_dataset([dataset_folders[4]],dataset_output[4],
                 m_wbc,M_wbc,wbc_feature_subset)
generate_dataset([dataset_folders[5]],dataset_output[5],
                 m_rbc,M_rbc,rbc_feature_subset)

generate_dataset([dataset_folders[0],dataset_folders[2]],'datasets/wbc_mll_adden_1_subset.h5',
                 m_wbc,M_wbc,wbc_feature_subset)
generate_dataset([dataset_folders[1],dataset_folders[3]],'datasets/rbc_mll_adden_1_subset.h5',
                 m_rbc,M_rbc,rbc_feature_subset)

dataset_output = [x + '_subset_unbiased.h5' for x in dataset_output_root]

generate_dataset([dataset_folders[0]],dataset_output[0],
                 m_wbc,M_wbc,wbc_feature_subset_unbiased)
generate_dataset([dataset_folders[1]],dataset_output[1],
                 m_rbc,M_rbc,rbc_feature_subset_unbiased)
generate_dataset([dataset_folders[2]],dataset_output[2],
                 m_wbc,M_wbc,wbc_feature_subset_unbiased)
generate_dataset([dataset_folders[3]],dataset_output[3],
                 m_rbc,M_rbc,rbc_feature_subset_unbiased)
generate_dataset([dataset_folders[4]],dataset_output[4],
                 m_wbc,M_wbc,wbc_feature_subset_unbiased)
generate_dataset([dataset_folders[5]],dataset_output[5],
                 m_rbc,M_rbc,rbc_feature_subset_unbiased)
