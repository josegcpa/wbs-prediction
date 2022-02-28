import argparse
import torch
import re
from tqdm import tqdm
from sklearn.model_selection import StratifiedKFold
from scipy.spatial.distance import cdist

from networks import *
from data_generator import *
from metrics import *

def obtain_model(network_state_dict):
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
    return VirtualCellClassifierStack(
        n_input_features=n_features,
        n_virtual_cells=n_virtual_cells,
        other_datasets_sizes=[od_features[0]],
        n_classes=n_classes).to(dev)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Test virtual cell classifier.')

    parser = argparse.ArgumentParser(description='Predict virtual cell dataset.')

    parser.add_argument('--dataset_path',dest='dataset_path',action='append',
                        type=str,default=None)
    parser.add_argument('--model_path_1',dest='model_path_1',action='store',
                        type=str,default=None)
    parser.add_argument('--fold_1',dest='fold_1',action='store',
                        type=int,default=0)
    parser.add_argument('--model_path_2',dest='model_path_2',action='store',
                        type=str,default=None)
    parser.add_argument('--fold_2',dest='fold_2',action='store',
                        type=int,default=0)
    parser.add_argument('--ob',dest='ob',action='store',
                        type=int,default=0)
    parser.add_argument('--excluded_ids',dest='excluded_ids',
                        nargs='+',action='store',type=str,default=None)
    parser.add_argument('--threshold',dest='threshold',
                        action='store',type=float,default=0.25)
    parser.add_argument('--n_cells_consensus',dest='n_cells_consensus',
                        action='store',type=int,default=250)

    args = parser.parse_args()

    if torch.cuda.is_available():
          dev = "cuda:0"
    else:
          dev = "cpu"

    n_d = len(args.dataset_path)

    # print("loading labels")

    # print("loading state dict and inferring network architecture from it")
    state_dict_1 = torch.load(args.model_path_1,map_location='cpu')
    network_state_dict_1 = state_dict_1[args.fold_1]['network']
    state_dict_2 = torch.load(args.model_path_2,map_location='cpu')
    network_state_dict_2 = state_dict_2[args.fold_2]['network']

    # print("assembling network and opening files, loading data moments")
    stacked_network_1 = obtain_model(network_state_dict_1)
    stacked_network_1.load_state_dict(network_state_dict_1)
    stacked_network_1.train(False)
    stacked_network_2 = obtain_model(network_state_dict_2)
    stacked_network_2.load_state_dict(network_state_dict_2)
    stacked_network_2.train(False)

    means_1 = [np.squeeze(x) for x in state_dict_1[args.fold_1]['means']]
    stds_1 = [np.squeeze(x) for x in state_dict_1[args.fold_1]['stds']]
    means_2 = [np.squeeze(x) for x in state_dict_2[args.fold_2]['means']]
    stds_2 = [np.squeeze(x) for x in state_dict_2[args.fold_2]['stds']]

    # print("loading datasets")
    all_datasets = [GenerateFromDataset(x,auto=False) for x in args.dataset_path]

    all_sets = [set(x.keys) for x in all_datasets]
    all_keys = list(set.intersection(*all_sets))
    all_keys = [k for k in all_keys if k not in args.excluded_ids]

    # print("building consensus set of cells")
    thr = args.threshold
    masks = []
    cells = []
    for dataset_obj in all_datasets:
        c = dataset_obj.generate_n_cells_flat(
            all_keys,args.n_cells_consensus)
        cells.append(c)

    cells_norm_1 = [(cells[i]-means_1[i])/stds_1[i] for i in range(len(cells))]
    cells_norm_1 = [torch.Tensor(x) for x in cells_norm_1]
    cells_norm_1 = [torch.unsqueeze(x,0) for x in cells_norm_1]
    cells_norm_1 = [torch.unsqueeze(x,0) for x in cells_norm_1]
    vc_1 = stacked_network_1.get_cells_individual(cells_norm_1)

    cells_norm_2 = [(cells[i]-means_2[i])/stds_2[i] for i in range(len(cells))]
    cells_norm_2 = [torch.Tensor(x) for x in cells_norm_2]
    cells_norm_2 = [torch.unsqueeze(x,0) for x in cells_norm_2]
    cells_norm_2 = [torch.unsqueeze(x,0) for x in cells_norm_2]
    vc_2 = stacked_network_2.get_cells_individual(cells_norm_2)

    for idx in range(len(all_datasets)):
        vc_pred_1 = np.squeeze(vc_1[idx].detach().numpy())
        vc_pred_2 = np.squeeze(vc_2[idx].detach().numpy())
        N_1 = np.zeros(vc_pred_1.shape[0],dtype=np.int32)
        N_2 = np.zeros(vc_pred_2.shape[0],dtype=np.int32)
        a,b = np.unique(np.argmax(vc_pred_1,axis=0),return_counts=True)
        N_1[a] += b
        a,b = np.unique(np.argmax(vc_pred_2,axis=0),return_counts=True)
        N_2[a] += b
        D = cdist(vc_pred_1,vc_pred_2,metric='correlation')
        Y = np.argmin(D,axis=1)
        D_min = D[[i for i in range(Y.size)],Y]
        thr = np.quantile(D_min,args.threshold)
        X = np.where(D_min <= thr)[0]
        Ds = D[X,Y[X]]
        for x,y,d in zip (X,Y[X],Ds):
            if N_1[x] > 100 and N_2[y] > 100:
                print('CELL_{},{},{},{},{},{}'.format(idx,x,y,(1-d)**2,N_1[x],N_2[y]))
