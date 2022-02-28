import argparse
import h5py
import numpy as np
from PIL import Image
import os
import re
from glob import glob
from tqdm import tqdm

def pad_to_size(image,size):
    h,w,_ = image.shape
    h_a = (size[0] - h) // 2
    w_a = (size[1] - w) // 2
    h_b = (size[0] - h_a - h)
    w_b = (size[1] - w_a - w)
    image = np.pad(image,((h_a, h_b), (w_a, w_b), (0, 0)),
                   mode='constant',constant_values=255)
    return image

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get virtual cell type examples.')

    parser.add_argument('--N_examples',dest='N_examples',
                        action='store',type=int,default=10)
    parser.add_argument('--size_wbc',dest='size_wbc',
                        action='store',type=int,default=128)
    parser.add_argument('--size_rbc',dest='size_rbc',
                        action='store',type=int,default=64)
    parser.add_argument('--dataset_id',dest='dataset_id',
                        action='store',type=str,default=None)
    parser.add_argument('--pattern',dest='pattern',
                        action='store',type=str,default="*")
    args = parser.parse_args()

    try: os.makedirs("virtual_cell_examples")
    except: pass
    try: os.makedirs(
        os.path.join("virtual_cell_examples",args.dataset_id))
    except: pass

    all_cell_collection_paths = glob(
        os.path.join('cell-collections',args.dataset_id,args.pattern+"h5"))

    for cell_collection_path in all_cell_collection_paths:
        print("Current cell collection: {}".format(cell_collection_path))
        root = cell_collection_path.split(os.sep)[-1][:-3]
        root_ = root.split('.')[0]
        root_dir = 'virtual_cell_examples/{}/{}'.format(
            args.dataset_id,root_)
        try: os.makedirs(root_dir)
        except: pass
        output_dir = "virtual_cell_examples/{}/{}/{}".format(
            args.dataset_id,root_,root)
        try: os.makedirs(output_dir)
        except: pass

        if 'wbc' in os.path.split(cell_collection_path)[-1][:10]:
            size = args.size_wbc
        else:
            size = args.size_rbc

        print("\tRetrieving all cells...")
        cell_type_mapping = {}
        cell_type_count = {}
        with h5py.File(cell_collection_path,'r') as F:
            all_keys = [k for k in F.keys()]
            for k in tqdm(all_keys):
                ct = F[k]['cell_type'][()]
                cell_type = np.argmax(ct) + 1
                if cell_type not in cell_type_mapping:
                    cell_type_mapping[cell_type] = [k]
                    cell_type_count[cell_type] = 1
                else:
                    cell_type_mapping[cell_type].append(k)
                    cell_type_count[cell_type] += 1

            print("\tCreating image...")
            M = max(cell_type_count.keys())
            for i in cell_type_mapping:
                t = args.N_examples ** 2
                if len(cell_type_mapping[i]) < t:
                    S = len(cell_type_mapping[i])
                    ss = cell_type_mapping[i]
                else:
                    S = t
                    options = [x for x in cell_type_mapping[i]]
                    ss = np.random.choice(options,S,replace=False)
                np.random.shuffle(ss)
                image_subset = [F[k]['image'] for k in ss]
                output_image = np.zeros((args.N_examples*size,args.N_examples*size,4),
                                        dtype=np.uint8)
                output_image[:,:,:] = 255
                for idx,image in enumerate(image_subset):
                    x,y = idx // args.N_examples,idx % args.N_examples
                    x,y = x*size,y*size
                    output_image[x:(x+size),y:(y+size),:3] = pad_to_size(image,(size,size))
                Image.fromarray(output_image).save(output_dir + '/' + str(i) + '.png')
