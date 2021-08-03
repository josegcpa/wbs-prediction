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

    print(network_state_dict['final_layers.0.0.weight'].shape)
    for i,x in enumerate(network_state_dict['final_layers.0.0.weight']):
        x = x.to('cpu').numpy()
        print(','.join([str(b) for b in x]))

    x = network_state_dict['final_layers.0.0.weight'].to('cpu').numpy()
    x = x[1,:] - x[0,:]
    print(','.join([str(b) for b in x]))
