from glob import glob
import os
import numpy as np
from tqdm import tqdm

try: os.makedirs('datasets')
except: pass
try: os.makedirs('datasets/histograms/')
except: pass

results_folder = '/nfs/research/gerstung/josegcpa/data/SLIDES/results'

dataset_folders = {
    'MLL':'/nfs/research/gerstung/josegcpa/data/SLIDES/MLL_TIFF/',
    'ADDEN1':'/nfs/research/gerstung/josegcpa/data/SLIDES/ADDEN_NDPI/',
    'ADDEN2':'/nfs/research/gerstung/josegcpa/data/SLIDES/ADDEN_SVS_results/'}

for dataset_folder in glob('{}/histograms/*'.format(results_folder)):
    dataset_key = dataset_folder.split(os.sep)[-1]
    for histogram_file in tqdm(glob('{}/*csv'.format(dataset_folder))):
        slide_id = histogram_file.split(os.sep)[-1][:-4]
        output = open('datasets/histograms/{}.csv'.format(slide_id),'w')
        with open(histogram_file,'r') as o:
            for l in o:
                if 'HIST' in l[:4]:
                    l = '{},{},{}'.format(slide_id,dataset_key,l)
                    output.write(l)

output = open('datasets/blur-maps.csv'.format(slide_id),'w')
for dataset_folder in glob('{}/blur-maps/*'.format(results_folder)):
    dataset_key = dataset_folder.split(os.sep)[-1]
    for blur_file in tqdm(glob('{}/*csv'.format(dataset_folder))):
        slide_id = blur_file.split(os.sep)[-1][:-4]
        with open(blur_file,'r') as o:
            for l in o:
                if 'BLUR' in l[:4]:
                    l = '{},{},{}'.format(slide_id,dataset_key,l)
                    output.write(l)

output_general = open('datasets/qc_summary.csv','w')
for dataset_key in dataset_folders:
    dataset_folder = dataset_folders[dataset_key]
    output = open('datasets/{}_qc.csv'.format(dataset_key),'w')
    for qc_file in tqdm(glob('{}/_quality_control/*'.format(dataset_folder))):
        slide_id = qc_file.split(os.sep)[-1]
        total = 0
        good = 0
        with open(qc_file,'r') as o:
            for l in o:
                if 'OUT' in l[:4]:
                    l = '{},{}'.format(slide_id,l)
                    output.write(l)
                    total += 1
                    if l.split(',')[-2] == '1':
                        good += 1
        output_general.write('{},{},{},{}\n'.format(slide_id,dataset_key,good,total))
