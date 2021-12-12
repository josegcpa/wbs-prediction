print('-'*20)
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
    print("Reading cmd line arguments...")
    parser = argparse.ArgumentParser(description='Test virtual cell classifier.')

    parser.add_argument('--n_virtual_cells',dest='n_virtual_cells',
                        action='store',
                        type=int,
                        default=20)
    parser.add_argument('--n_classes',dest='n_classes',
                        action='store',
                        type=int,
                        default=2)
    parser.add_argument('--dataset_path',dest='dataset_path',
                        action='append',
                        type=str,
                        default=None)
    parser.add_argument('--other_dataset_path',dest='other_dataset_path',
                        action='append',
                        type=str,
                        default=None)
    parser.add_argument('--batch_size',dest='batch_size',
                        action='store',
                        type=int,
                        default=32)
    parser.add_argument('--number_of_cells',dest='number_of_cells',
                        action='store',
                        type=int,
                        default=500)
    parser.add_argument('--model_path',dest='model_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--mode',dest='mode',
                        action='store',
                        type=str,
                        choices=['proportions','cells'],
                        default="proportions")
    parser.add_argument('--fold',dest='fold',
                        action='store',
                        type=int,
                        default=0)

    args = parser.parse_args()

    if torch.cuda.is_available():
          dev = "cuda:0"
    else:
          dev = "cpu"

    all_datasets = [GenerateFromDataset(x)  for x in args.dataset_path]
    if args.other_dataset_path is not None:
        all_other_datasets = [CSVDataset(x,handle_nans='remove')
                              for x in args.other_dataset_path]
        all_other_datasets_keys = [set(x.keys) for x in all_other_datasets]
    else:
        all_other_datasets = None

    all_lists = {}

    all_networks = []

    all_lists['queued_datasets'] = []

    if all_other_datasets is not None:
        other_datasets_size = [x.n_features
                               for x in all_other_datasets]
    else:
        other_datasets_size = [0]
    all_input_features = [d.n_features for d in all_datasets]
    all_n_virtual_cells = [args.n_virtual_cells for d in all_datasets]
    stacked_network = VirtualCellClassifierStack(
        all_input_features,all_n_virtual_cells,
        other_datasets_size,[args.n_classes]).to(dev)

    state_dict = torch.load(args.model_path,map_location='cpu')
    network_state_dict = state_dict[args.fold]['network']
    stacked_network.load_state_dict(network_state_dict)

    all_keys = [*[x.keys for x in all_datasets],*all_other_datasets_keys]
    mutual_keys = list(set.intersection(*[set(x) for x in all_keys]))

    all_datasets[0].mean = state_dict[args.fold]['means'][0]
    all_datasets[1].mean = state_dict[args.fold]['means'][1]
    all_datasets[0].std = state_dict[args.fold]['stds'][0]
    all_datasets[1].std = state_dict[args.fold]['stds'][1]
    if all_other_datasets is not None:
        for i in range(len(all_other_datasets)):
            all_other_datasets[i].mean = state_dict[
                args.fold]['means'][i+len(all_datasets)]
            all_other_datasets[i].std = state_dict[
                args.fold]['stds'][i+len(all_datasets)]

    if args.mode == 'proportions':
        for key in mutual_keys:
            T = []
            N = []
            for ds in all_datasets:
                n = np.minimum(ds.all_datasets[key].n_cells,
                               args.number_of_cells)
                T.append(torch.Tensor(ds.generate_n_cells([key],n)[:,np.newaxis,:,:]))
                N.append(n)
            OT = []
            for od in all_other_datasets:
                OT.append(torch.Tensor(od[[key]]))
            probs = stacked_network([T,OT]).to('cpu').detach().numpy().flatten()
            C = []
            for i in range(len(stacked_network.vcq_list)):
                cells = stacked_network[0].get_virtual_cells(T[0])
                cells = cells.to('cpu').detach().numpy().flatten()
                C.append(cells)
            output = []
            for n,c in zip(N,C):
                output.extend([str(x) for x in c])
                output.append(str(n))
            output.extend([str(x) for x in probs])

            print(key + ',' + ','.join(output))

    if args.mode == 'class':
        for key in mutual_keys:
            n0 = np.minimum(all_datasets[0].all_datasets[key].n_cells,
                            args.number_of_cells)
            n1 = np.minimum(all_datasets[1].all_datasets[key].n_cells,
                            args.number_of_cells)
            T0 = torch.Tensor(all_datasets[0].generate_n_cells([key],n0)[:,np.newaxis,:,:])
            T1 = torch.Tensor(all_datasets[1].generate_n_cells([key],n1)[:,np.newaxis,:,:])
            OT0 = torch.Tensor(all_other_datasets[0][[key]])
            OT1 = torch.Tensor(all_other_datasets[1][[key]])
            probs = stacked_network([[T0,T1],[OT0,OT1]]).to('cpu').detach().numpy().flatten()

            cells1 = cells1.to('cpu').detach().numpy().flatten()
            cells2 = cells2.to('cpu').detach().numpy().flatten()
            output = [
                *[str(x) for x in cells1],str(n0),
                *[str(x) for x in cells2],str(n1),
                *[str(x) for x in probs]]

            print(key + ',' + ','.join(output))

    elif args.mode == 'cells':
        for key in mutual_keys:
            n0 = np.minimum(all_datasets[0].all_datasets[key].n_cells,
                            args.number_of_cells)
            n1 = np.minimum(all_datasets[1].all_datasets[key].n_cells,
                            args.number_of_cells)
            T0 = torch.Tensor(all_datasets[0].generate_n_cells([key],n0)[:,np.newaxis,:,:])
            T1 = torch.Tensor(all_datasets[1].generate_n_cells([key],n1)[:,np.newaxis,:,:])

            OT0 = torch.Tensor(all_other_datasets[0][[key]])
            OT1 = torch.Tensor(all_other_datasets[1][[key]])
            probs = stacked_network([[T0,T1],[OT0,OT1]]).to('cpu').detach().numpy().flatten()
            cells1 = stacked_network[0].f(T0)
            cells2 = stacked_network[1].f(T1)

            cells1 = cells1.to('cpu').detach().numpy()[0,:,:,0]
            cells2 = cells2.to('cpu').detach().numpy()[0,:,:,0]

            for p,cells in zip(args.dataset_path,[cells1,cells2]):
                for cell in range(cells.shape[1]):
                    output = [str(x) for x in cells[:,cell]]
                    output.append(p)
                    print(key + ',' + ','.join(output))
