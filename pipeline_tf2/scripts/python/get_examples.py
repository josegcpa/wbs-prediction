import argparse
import sys
import psutil
import pickle
import time
import numpy as np
import h5py
import openslide
import re
from tqdm import tqdm

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get cell examples.')

    parser.add_argument('--slide_path',dest='slide_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--aggregates_path',dest='aggregates_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--segmented_path',dest='segmented_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--output_path',dest='output_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--subset',dest='subset',
                        action='store',
                        type=int,
                        default=None)
    parser.add_argument('--flip',dest='flip',
                        action='store_true')
    parser.add_argument('--output_size',dest='output_size',
                        action='store',type=int,default=None)
    parser.add_argument('--random',dest='random',action='store_true')
    parser.add_argument('--box',dest='box',action='store_true')

    args = parser.parse_args()

    slide = openslide.OpenSlide(args.slide_path)
    aggregates = h5py.File(args.aggregates_path,'r')
    segmentations = h5py.File(args.segmented_path,'r')
    output = h5py.File(args.output_path,'w')

    subset_size = args.subset
    all_sizes = [
        aggregates['cells'][x].shape[0]
        for x in aggregates['cells']]
    n_cells = sum(all_sizes)
    if subset_size > n_cells:
        subset_size = n_cells
    subset_of_cells = np.random.choice(
        n_cells,
        size=subset_size,replace=False)
    subset_of_cells = np.sort(subset_of_cells)
    corresponding_subsets = np.floor(subset_of_cells / all_sizes[0])
    subset_dict = {}
    for c,s in zip(subset_of_cells,corresponding_subsets):
        s = int(s)
        c = c % all_sizes[0]
        if s not in subset_dict:
            subset_dict[s] = [c]
        else:
            subset_dict[s].append(c)
    for subset in subset_dict:
        a = all_sizes[0] * subset
        for cell_idx in tqdm(subset_dict[s]):
            cell = aggregates['cells'][str(subset)][cell_idx,:]
            i = cell_idx + a
            cell_center = aggregates['cell_centers'][i]
            cell_center_int = cell_center.decode('utf-8')[1:-1].split(',')
            cell_center_int = [float(cell_center_int[0]),float(cell_center_int[1])]
            cell_coordinates = segmentations[cell_center]
            if args.flip == True:
                x = cell_coordinates['Y'][::]
                y = cell_coordinates['X'][::]
            else:    
                x = cell_coordinates['X'][::]
                y = cell_coordinates['Y'][::]
            x_min,x_max = x.min(),x.max()
            y_min,y_max = y.min(),y.max()
            rx,ry = (x_max-x_min),(y_max-y_min)
            if args.output_size is None:
                location = [x.min()-3,y.min()-3]
                size = [x.max()-location[0] + 6,y.max()-location[1] + 6]
            else:
                if args.random == True:
                    shift_x,shift_y = [
                        np.random.randint(args.output_size)-rx,
                        np.random.randint(args.output_size)-ry]
                else:
                    shift_x = int(args.output_size//2 - rx//2)
                    shift_y = int(args.output_size//2 - ry//2)
                location = [x.min() - shift_x,y.min() - shift_y]
                size = [args.output_size,args.output_size]
            image = np.array(slide.read_region(location,0,size))[:,:,:3]
            if args.box == True: 
                A,B,C,D = shift_x,shift_y,shift_x+rx,shift_y+ry
                color = [160, 32, 240]
                image[A:C,B,:] = color
                image[A:C,D,:] = color
                image[A,B:D,:] = color
                image[C,B:D,:] = color
            segmented_image = np.zeros_like(image[:,:,0])
            segmented_image[(y-location[1],x-location[0])] = 255
            i += 1
            g = output.create_group(cell_center)
            g.create_dataset('image',data=image,dtype=np.uint8)
            g.create_dataset('isolated_image',data=segmented_image,dtype=np.uint8)
            g.create_dataset('cell_center_y',data=cell_center_int[0])
            g.create_dataset('cell_center_x',data=cell_center_int[1])
            g.create_dataset('features',data=cell)

