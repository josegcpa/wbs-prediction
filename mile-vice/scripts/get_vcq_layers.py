import argparse
import sys
import torch
import psutil
import pickle
import time

from networks import *
from data_generator import *
from metrics import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Test virtual cell classifier.')

    parser.add_argument('--model_path',dest='model_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--fold',dest='fold',
                        action='store',
                        type=int,
                        default=0)

    args = parser.parse_args()

    state_dict = torch.load(args.model_path,map_location='cpu')
    network_state_dict = state_dict[args.fold]['network']
    ks = network_state_dict.keys()
    ks = [x for x in ks if 'weight' in x and '.f.' in x]
    for dataset_idx,k in enumerate(ks):
        layer = np.squeeze(network_state_dict[k]).T
        for i,x in enumerate(layer):
            x = x.to('cpu').numpy()
            base_str = '{},{},{},{},'.format(
                args.model_path,str(args.fold),
                str(i),'dataset_'+str(dataset_idx))
            for f_no,b in enumerate(x):
                print(base_str+str(f_no)+','+str(b))
