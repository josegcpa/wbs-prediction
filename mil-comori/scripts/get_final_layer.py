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
    parser = argparse.ArgumentParser(
        description='Gets last layer of MIL-CoMorI.')

    parser.add_argument('--model_path',dest='model_path',
                        action='store',type=str,default=None,
                        help="Path to MIL-CoMorI model")
    parser.add_argument('--fold',dest='fold',
                        action='store',type=int,default=0,
                        help="Model fold")

    args = parser.parse_args()

    state_dict = torch.load(args.model_path,map_location='cpu')
    network_state_dict = state_dict[args.fold]['network']
    ks = network_state_dict.keys()
    ks = [x for x in ks if 'weight' in x and 'final_layers' in x]
    for task_idx,k in enumerate(ks):
        for i,x in enumerate(network_state_dict[k]):
            x = x.to('cpu').numpy()
            base_str = '{},{},{},{},'.format(
                args.model_path,str(args.fold),str(i),task_idx
            )
            for vc,b in enumerate(x):
                print(base_str+str(vc)+','+str(b))
