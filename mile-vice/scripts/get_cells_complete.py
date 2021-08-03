print('-'*20)
import argparse
import sys
import torch
import psutil
import pickle
import time
import h5py
import openslide
import re

from networks import *
from data_generator import *
from metrics import *

if __name__ == "__main__":
    print("Reading cmd line arguments...")
    parser = argparse.ArgumentParser(description='Test virtual cell classifier.')

    parser.add_argument('--network_idx',dest='network_idx',
                        action='store',
                        type=int,
                        default=0)
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
    parser.add_argument('--model_path',dest='model_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--output_path',dest='output_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--fold',dest='fold',
                        action='store',
                        type=int,
                        default=0)
    parser.add_argument('--subset',dest='subset',
                        action='store',
                        type=int,
                        default=None)
    parser.add_argument('--flip',dest='flip',
                        default=False,
                        action='store_true')

    args = parser.parse_args()

    if torch.cuda.is_available():
          dev = "cuda:0"
    else:
          dev = "cpu"

    ni = args.network_idx
    state_dict = torch.load(args.model_path,map_location='cpu')
    network_state_dict = state_dict[args.fold]['network']
    n_features = [
        network_state_dict[k].shape[-1]
        for k in network_state_dict
        if len(network_state_dict[k].shape) == 4]
    n_virtual_cells = [
        int(network_state_dict[k].shape[0])
        for k in network_state_dict
        if len(network_state_dict[k].shape) == 4]
    final_layer_names = [x
                         for x in network_state_dict
                         if re.search("final.*weight",x)]
    final_layers= [network_state_dict[x]
                   for x in final_layer_names]
    n_classes = [x.shape[0] for x in final_layers]
    od_features = [int(x.shape[-1] - np.sum(n_virtual_cells))
                   for x in final_layers]

    stacked_network = VirtualCellClassifierStack(
        n_input_features=n_features,
        n_virtual_cells=n_virtual_cells,
        other_datasets_sizes=od_features,
        n_classes=n_classes).to(dev)

    stacked_network.load_state_dict(network_state_dict)

    means = np.squeeze(state_dict[args.fold]['means'][ni])
    stds = np.squeeze(state_dict[args.fold]['stds'][ni])

    slide = openslide.OpenSlide(args.slide_path)
    aggregates = h5py.File(args.aggregates_path,'r')
    segmentations = h5py.File(args.segmented_path,'r')
    output = h5py.File(args.output_path,'w')

    vcq = stacked_network.vcq_list[ni]

    if args.subset is not None:
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
            for cell_idx in subset_dict[s]:
                cell = aggregates['cells'][str(subset)][cell_idx,:]
                i = cell_idx + a
                cell_center = aggregates['cell_centers'][i]
                cell_center_int = cell_center.decode('utf-8')[1:-1].split(',')
                cell_center_int = [float(cell_center_int[0]),float(cell_center_int[1])]
                cell_coordinates = segmentations[cell_center]
                x = cell_coordinates['X'][::]
                y = cell_coordinates['Y'][::]
                location = [x.min()-3,y.min()-3]
                size = [x.max()-location[0] + 6,y.max()-location[1] + 6]
                if args.flip == True:
                    location = [location[1],location[0]]
                    size = [size[1],size[0]]
                image = np.array(slide.read_region(location,0,size))[:,:,:3]
                cell_transformed = (cell - means)/stds
                cell_transformed = cell_transformed[
                    np.newaxis,np.newaxis,np.newaxis,:]
                cell_transformed = torch.Tensor(np.float32(cell_transformed))
                vc = vcq.get_virtual_cells(cell_transformed).to('cpu').detach().numpy()
                i += 1
                g = output.create_group(cell_center)
                g.create_dataset('image',data=image,dtype=np.uint8)
                g.create_dataset('cell_type',data=vc,dtype=np.float32)
                g.create_dataset('features',data=cell,dtype=np.float64)
                g.create_dataset('cell_center_y',data=cell_center_int[0])
                g.create_dataset('cell_center_x',data=cell_center_int[1])

    else:
        i = 0
        for subset in aggregates['cells']:
            for cell in aggregates['cells'][subset]:
                cell_center = aggregates['cell_centers'][i]
                cell_center_int = cell_center.decode('utf-8')[1:-1].split(',')
                cell_center_int = [float(cell_center_int[0]),float(cell_center_int[1])]
                cell_coordinates = segmentations[cell_center]
                x = cell_coordinates['X'][::]
                y = cell_coordinates['Y'][::]
                location = [x.min()-3,y.min()-3]
                size = [x.max()-location[0] + 6,y.max()-location[1] + 6]
                if args.flip == True:
                    location = [location[1],location[0]]
                    size = [size[1],size[0]]
                image = np.array(slide.read_region(location,0,size))[:,:,:3]
                cell_transformed = (cell - means)/stds
                cell_transformed = cell_transformed[
                    np.newaxis,np.newaxis,np.newaxis,:]
                cell_transformed = torch.Tensor(np.float32(cell_transformed))
                vc = vcq.get_virtual_cells(cell_transformed).to('cpu').detach().numpy()
                i += 1
                g = output.create_group(cell_center)
                g.create_dataset('image',data=image,dtype=np.uint8)
                g.create_dataset('cell_type',data=vc,dtype=np.float32)
                g.create_dataset('features',data=cell,dtype=np.float64)
                g.create_dataset('cell_center_y',data=cell_center_int[0])
                g.create_dataset('cell_center_x',data=cell_center_int[1])
