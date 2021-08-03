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

    parser.add_argument('--slide_path',dest='slide_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--aggregates_path',dest='aggregates_path',
                        action='append',
                        type=str,
                        default=None)
    parser.add_argument('--segmented_path',dest='segmented_path',
                        action='append',
                        type=str,
                        default=None)
    parser.add_argument('--other_datasets',dest='other_datasets',
                        action='append',
                        type=str,
                        default=[])
    parser.add_argument('--dataset_id',dest='dataset_id',
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
                        default=500)
    parser.add_argument('--complete_summary',dest='complete_summary',
                        action='store_true',
                        default=False)
    parser.add_argument('--flip',dest='flip',
                        default=False,
                        action='store_true')

    args = parser.parse_args()

    if torch.cuda.is_available():
          dev = "cuda:0"
    else:
          dev = "cpu"

    n_od = len(args.other_datasets)
    n_d = len(args.aggregates_path)
    # loading state dict and inferring network architecture from
    # state dict
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

    # assembling network and opening files, loading data moments
    stacked_network = VirtualCellClassifierStack(
        n_input_features=n_features,
        n_virtual_cells=n_virtual_cells,
        other_datasets_sizes=[od_features[0]],
        n_classes=n_classes).to(dev)

    stacked_network.load_state_dict(network_state_dict)

    means = [np.squeeze(x) for x in state_dict[args.fold]['means']]
    stds = [np.squeeze(x) for x in state_dict[args.fold]['stds']]
    means_od = [np.squeeze(x).astype(np.float32)
                for x in state_dict[args.fold]['means'][n_d:]]
    stds_od = [np.squeeze(x).astype(np.float32)
               for x in state_dict[args.fold]['stds'][n_d:]]

    slide = openslide.OpenSlide(args.slide_path)
    aggregates = [h5py.File(x,'r') for x in args.aggregates_path]
    segmentations = [h5py.File(x,'r') for x in args.segmented_path]
    output = h5py.File(args.output_path,'w')

    # get data from cell collections here
    all_cell_sets = []
    for dataset_idx,(aggregate,segmentation) in enumerate(zip(aggregates,segmentations)):
        subset_size = args.subset
        all_sizes = [
            aggregate['cells'][x].shape[0]
            for x in aggregate['cells']]
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
        cells_transformed = []
        for subset in subset_dict:
            a = all_sizes[0] * subset
            for cell_idx in subset_dict[s]:
                cell = aggregate['cells'][str(subset)][cell_idx,:]
                if args.complete_summary:
                    i = cell_idx + a
                    cell_center = aggregate['cell_centers'][i]
                    cell_coordinates = segmentation[cell_center]
                    x = cell_coordinates['X'][::]
                    y = cell_coordinates['Y'][::]
                    location = [x.min()-3,y.min()-3]
                    size = [x.max()-location[0] + 6,y.max()-location[1] + 6]
                    if args.flip == True:
                        location = [location[1],location[0]]
                        size = [size[1],size[0]]
                    image = np.array(slide.read_region(location,0,size))[:,:,:3]
                    g = output.create_group('')
                    g.create_dataset('image',data=image,dtype=np.uint8)
                    g.create_dataset('cell_type',data=vc,dtype=np.float32)
                    g.create_dataset('features',data=cell,dtype=np.float64)
                    g.create_dataset('cell_center_y',data=cell_center_int[0])
                    g.create_dataset('cell_center_x',data=cell_center_int[1])

                cell_transformed = (cell - means[dataset_idx])/stds[dataset_idx]
                cell_transformed = cell_transformed[
                    np.newaxis,np.newaxis,np.newaxis,:]
                cells_transformed.append(cell_transformed)

        cells_transformed = np.concatenate(cells_transformed,axis=2)
        cells_transformed = torch.Tensor(np.float32(cells_transformed)).to(dev)
        all_cell_sets.append(cells_transformed)

    # get data from other datasets
    other_datasets = [torch.Tensor((CSVDataset(x)[[args.dataset_id]])).to(dev) for x in args.other_datasets]
    other_datasets = [(x-means_od[i])/stds_od[i] for i,x in enumerate(other_datasets)]
    # predictions only happen here
    for i in range(len(final_layers)):
        prediction = stacked_network([all_cell_sets,other_datasets],i)
        l = state_dict['args'].labels_path[i]
        print(l + ',' + ','.join([str(x) for x in prediction.detach().cpu().numpy().squeeze()]))
