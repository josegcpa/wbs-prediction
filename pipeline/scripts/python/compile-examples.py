import argparse
import sys
import time
import numpy as np
import h5py
from tqdm import tqdm
from glob import glob

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compile examples.')

    parser.add_argument('--subset',dest='subset',
                        action='store',
                        type=int,
                        default=50)
    parser.add_argument('--pattern',dest='pattern',
                        action='store',
                        type=str,
                        default="wbc")
    parser.add_argument('--output_path',dest='output_path',
                        action='store',
                        type=str,
                        default=None)
    args = parser.parse_args()

    out = h5py.File(args.output_path,"w")

    for examples_path in tqdm(glob("examples/*{}*h5".format(args.pattern))):
        slide_id = '_'.join(examples_path.split('_')[:-1])
        with h5py.File(examples_path,"r") as F:
            all_keys = [x for x in F.keys()]
            if args.subset > len(all_keys):
                S = len(all_keys)
            else:
                S = args.subset 

            ss = np.random.choice(all_keys,S,replace=False)

            for k in ss:
                cell = F[k]
                out_cell = out.create_group(k)
                for x in ['cell_center_x', 'cell_center_y', 'features', 'image']:
                    out_cell[x] = cell[x][()]
                out_cell['slide_id'] = slide_id
    
    out.close()