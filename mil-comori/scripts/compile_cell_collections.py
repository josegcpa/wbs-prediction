import h5py
from glob import glob
import os
import argparse
import numpy as np
import re
from tqdm import tqdm

if __name__ == "__main__":
    print("Reading cmd line arguments...")
    parser = argparse.ArgumentParser(description='Compile cell collections.')

    parser.add_argument('--collection_path',dest='collection_path',
                        action='store',type=str,default=None,
                        help="Path to folder containing collections")
    parser.add_argument('--pattern',dest='pattern',
                        action='store',type=str,default='wbc*h5',
                        help="Pattern to be used when globing collections")
    parser.add_argument('--output_path',dest='output_path',
                        action='store',type=str,default=None,
                        help="Output path for compile cell collection")
    parser.add_argument('--no_features',dest='no_features',
                        action='store_true',
                        help="If no features are to be stored")
    parser.add_argument('--simplify_cell_type',dest='simplify_cell_type',
                        action='store_true',
                        help="Virtual cell types are quantized")
    parser.add_argument('--n',dest='n',
                        action='store',type=int,default=50,
                        help="Number of cells from each cell collection")

    args = parser.parse_args()

    all_h5 = glob(os.path.join(args.collection_path,args.pattern))

    output = h5py.File(args.output_path,'w')
    i = 0
    for h5 in tqdm(all_h5):
        F_ = h5py.File(h5,'r')
        k = list(F_.keys())
        if 'cells' in k:
            F = F_['cells']
        else:
            F = F_
        all_keys = list(F.keys())
        S = np.minimum(len(all_keys),args.n)
        k_subset = np.random.choice(all_keys,S,replace=False)
        for k in k_subset:
            cell_center_int = str(k)[1:-1].split(',')
            cell_center_int = [float(cell_center_int[0]),float(cell_center_int[1])]
            g = output.create_group(str(i))
            g['slide_path'] = h5
            g.create_dataset('image',data=F[k]['image'][::],
                             dtype=np.uint8)
            if args.simplify_cell_type == True:
                cell_type = F[k]['cell_type'][::] > 0.4
                cell_type = cell_type.astype(np.int8)
                g.create_dataset('cell_type',data=cell_type,
                                 dtype=np.int8)
            else:
                g.create_dataset('cell_type',data=F[k]['cell_type'][::],
                                 dtype=np.float32)
            if args.no_features == False:
                g.create_dataset('features',data=F[k]['features'][::],
                                 dtype=np.float64)

            g.create_dataset('cell_center_x',data=cell_center_int[1])
            g.create_dataset('cell_center_y',data=cell_center_int[0])
            i += 1
