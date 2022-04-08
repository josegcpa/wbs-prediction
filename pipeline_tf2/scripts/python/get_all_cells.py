"""
Generate examples of cells in an HDF5 format.

Usage:
    python3 get_all_cells.py --help
"""

import argparse
import numpy as np
import h5py
import openslide
from tqdm import tqdm

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Get all cells.')

    parser.add_argument('--slide_path',dest='slide_path',
                        action='store',type=str,default=None,
                        help="Path to slide.")
    parser.add_argument('--aggregates_path',dest='aggregates_path',
                        action='store',type=str,default=None,
                        help="Path to cell characteristics HDF5")
    parser.add_argument('--output_path',dest='output_path',
                        action='store',type=str,default=None,
                        help="Path to output")
    parser.add_argument('--size',dest='size',
                        action='store',type=int,default=64,
                        help="Size of the crop")

    args = parser.parse_args()

    OS = openslide.OpenSlide(args.slide_path)
    aggregates = h5py.File(args.aggregates_path,'r')
    output = h5py.File(args.output_path,'w')

    m = args.size // 2
    all_cell_center_idx = aggregates['cell_centers']
    for center in tqdm(all_cell_center_idx):
        center = center.decode()
        center_int = center[1:-1].split(',')
        center_int = [float(center_int[0]),float(center_int[1])]
        loc = [int(center_int[0])-m,int(center_int[1])-m]
        image = OS.read_region(loc,0,[args.size,args.size])
        g = output.create_group(center)
        g.create_dataset('image',data=image,dtype=np.uint8)