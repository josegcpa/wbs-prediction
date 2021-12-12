import argparse
import torch
import re

from networks import *
from data_generator import *
from metrics import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Predict virtual cell dataset.')

    parser.add_argument('--dataset_path',dest='dataset_path',action='append',
                        type=str,default=None)
    parser.add_argument('--model_path',dest='model_path',action='store',
                        type=str,default=None)
    parser.add_argument('--fold',dest='fold',action='store',
                        type=int,default=0)

    args = parser.parse_args()

    if torch.cuda.is_available():
          dev = "cuda:0"
    else:
          dev = "cpu"

    n_d = len(args.dataset_path)

    # loading state dict and inferring network architecture from it
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
    final_layer_names = [x for x in network_state_dict
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
    stacked_network.train(False)

    means = [np.squeeze(x) for x in state_dict[args.fold]['means']]
    stds = [np.squeeze(x) for x in state_dict[args.fold]['stds']]

    all_datasets = [GenerateFromDataset(x) for x in args.dataset_path]

    all_sets = [set(x.keys) for x in all_datasets]
    all_keys = list(set.intersection(*all_sets))
    all_keys = ['XIII_17']
    # get data from cell collections here
    for key in all_keys:
        all_cell_sets = []
        for di,dataset in enumerate(all_datasets):
            cells = dataset[key].return_all_cells()
            cells_transformed = (cells-means[di])/stds[di]
            cells_transformed = cells_transformed[np.newaxis,np.newaxis,:,:]
            cells_transformed = torch.Tensor(np.float32(cells_transformed)).to(dev)
            all_cell_sets.append(cells_transformed)
        
        # predict here
        prediction = stacked_network.get_cells(all_cell_sets)
        for j,pred in enumerate(prediction):
            print(key + ',' + 'cell_type_{}'.format(j) + ',' + ','.join(
                [str(x) for x in pred.detach().cpu().numpy().squeeze()]))
