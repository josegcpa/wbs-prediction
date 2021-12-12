import h5py
from glob import glob
import os
import argparse
import numpy as np
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compile cell collections.')

    parser.add_argument('--collection_path',dest='collection_path',
                        action='store',
                        type=str,
                        default=None)

    args = parser.parse_args()

    all_h5 = glob(os.path.join(args.collection_path,"*h5"))

    for h5 in all_h5:
        F = h5py.File(h5,'r')
        cell_proportions = F['cell_proportions'][()]
        slide_id = os.path.split(h5)[-1].split('.')[0]
        cp_str = ','.join([str(x) for x in cell_proportions.tolist()])
        if 'rbc-' in slide_id:
            slide_id = slide_id[4:]
            cell_type = "cell_type_1"
        elif 'wbc-' in slide_id:
            slide_id = slide_id[4:]
            cell_type = "cell_type_0"
        print(slide_id + ',' + cell_type + ',' + cp_str)
