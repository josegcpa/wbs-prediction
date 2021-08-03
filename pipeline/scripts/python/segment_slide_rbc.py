import sys
import argparse
import os
import numpy as np
import openslide
import time
from PIL import Image
from multiprocessing import Pool
import h5py
import xgboost
import pickle

from mask_rbc import wraper as mask_rbc
from image_generator import ImageGeneratorWithQueue
import MIA

from skimage import io

MIA_FEATURES = MIA.MIA_FEATURES

def fn(image):
    image,coords = image
    _,mask = mask_rbc(image)
    features = MIA.wrapper(image,mask)
    return features,coords

parser = argparse.ArgumentParser(
    prog = 'segmented_slide_rbc.py',
    description = 'segments all red blood cells in slide'
)

parser.add_argument('--csv_path',dest='csv_path',
                    action='store',
                    default=None,
                    help='Path to CSV file with the quality control.')
parser.add_argument('--slide_path',dest='slide_path',
                    action='store',
                    default=None,
                    help='Path to slide.')
parser.add_argument('--output_path',dest='output_path',
                    action='store',
                    default=False,
                    help='Output path for the hdf5 file.')
parser.add_argument('--n_processes_analysis',dest='n_processes_analysis',
                    action='store',
                    type=int,
                    default=1,
                    help='Number of processes for analysis.')
parser.add_argument('--n_processes_data',dest='n_processes_data',
                    action='store',
                    type=int,
                    default=1,
                    help='Number of processes for data loading.')

args = parser.parse_args()

times = []
image_list = []
pool = Pool(args.n_processes_analysis)
F = h5py.File(args.output_path,mode='w')
N = 0
i = 0

igwq = ImageGeneratorWithQueue(args.slide_path,args.csv_path,0,
                               maxsize=args.n_processes_data)
igwq.start()

osss = openslide.OpenSlide(args.slide_path)

for image in igwq.generate():
    a = time.time()
    image_list.append(image)
    i += 1
    if len(image_list) == args.n_processes_analysis:
        output = pool.map(fn,image_list)
        for im_dict_list in output:
            im_dict_list,coords = im_dict_list
            for obj in im_dict_list:
                y = obj['x'] + coords[1]
                x = obj['y'] + coords[0]
                features = []
                for k in MIA_FEATURES:
                    features.append(obj[k])
                if np.any([
                        np.any(x == 0),
                        np.any(x == (512+256)),
                        np.any(y == 0),
                        np.any(y == (512+256))
                ]):
                    pass
                else:
                    C = str([np.mean(x),np.mean(y)])
                    if C not in F:
                        N += 1
                        g = F.create_group(str(C))
                        g.create_dataset("X",x.shape,
                                         dtype=np.int32,data=x)
                        g.create_dataset("Y",y.shape,
                                         dtype=np.int32,data=y)
                        g.create_dataset("features",[len(features)],
                                         dtype=np.float64,data=features)
        times.append((time.time() - a)/args.n_processes_analysis)
        image_list = []
    if (i-1) % 100 == 0:
        print('Image {}. {} {}'.format(i-1,np.mean(times),N))

F.attrs['NCells'] = N
F.attrs['NImages'] = i
F.close()
