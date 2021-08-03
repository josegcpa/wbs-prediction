from glob import glob
import os
import h5py
import argparse
import openslide
import numpy as np
from tqdm import tqdm
from PIL import Image

if __name__ == "__main__":
    print("Reading cmd line arguments...")
    parser = argparse.ArgumentParser(description='Assemble dataset to train VAE.')

    parser.add_argument('--segmentations_folder',dest='segmentations_folder',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--slides_folder',dest='slides_folder',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--output_path',dest='output_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--size',dest='size',
                        action='store',
                        type=int,
                        default=64)
    parser.add_argument('--subset',dest='subset',
                        action='store',
                        type=int,
                        default=250)
    args = parser.parse_args()

    all_segmentations = glob(os.path.join(args.segmentations_folder,'*'))
    all_slides = glob(os.path.join(args.slides_folder,'*'))

    all_combinations = {}
    for s in all_segmentations:
        root = '.'.join(os.path.split(s)[-1].split('.')[:-1])
        all_combinations[root] = {'segmentation':s}
    for s in all_slides:
        root = '.'.join(os.path.split(s)[-1].split('.')[:-1])
        if root in all_combinations:
            all_combinations[root]['slide'] = s

    output = h5py.File(args.output_path,'w')
    n_images = 0
    for s in all_combinations:
        print(s)
        oss = openslide.OpenSlide(all_combinations[s]['slide'])
        f = h5py.File(all_combinations[s]['segmentation'],'r')
        all_keys = [x for x in f.keys()]
        if len(all_keys) > args.subset:
            K = np.random.choice(all_keys,args.subset,replace=False)
        else:
            K = all_keys
        im_list = []
        for k in K:
            try:
                x,y = [int(float(j)-args.size/2) for j in k[1:-1].split(',')]
                im = np.array(oss.read_region([x,y],0,[args.size,args.size]))
                im_list.append(im)
                n_images += 1
            except:
                pass

        im_array = np.stack(im_list,axis=0)
        _ = output.create_dataset(s,data=im_array)

    output.attrs['size'] = n_images
