import h5py
import os
import numpy as np
import cv2
from tqdm import tqdm

try: os.makedirs('datasets')
except: pass

results_folder = '/nfs/research/gerstung/josegcpa/projects/01IMAGE/quality-net'

output = open('datasets/blur-quality-net.csv','w')
with h5py.File('{}/dataset_hdf5/train.h5'.format(results_folder)) as F:
    for k in tqdm(F):
        image = F[k]['image'][()]
        label = F[k]['class'][()]
        blur = cv2.Laplacian(np.uint8(image.mean(axis=-1)),
                             cv2.CV_64F).var()
        hist,_ = np.histogram(image,bins=256,range=[0,255])
        hist = ','.join([str(x) for x in hist])
        output.write('{},{},{},{}\n'.format(k,label,blur,hist))

with h5py.File('{}/dataset_hdf5/test.h5'.format(results_folder)) as F:
    for k in tqdm(F):
        image = F[k]['image'][()]
        label = F[k]['class'][()]
        blur = cv2.Laplacian(np.uint8(image.mean(axis=-1)),
                             cv2.CV_64F).var()
        hist,_ = np.histogram(image,bins=256,range=[0,255])
        hist = ','.join([str(x) for x in hist])
        output.write('{},{},{},{}\n'.format(k,label,blur,hist))
