import argparse
import torch
import re
from tqdm import tqdm
from sklearn.model_selection import StratifiedKFold
from scipy.spatial.distance import cdist

from networks import *
from data_generator import *
from metrics import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Tests MIL-CoMorI models with consensual VCTs.')

    parser.add_argument('--dataset_path',dest='dataset_path',action='append',
                        type=str,default=None,
                        help="Dataset used for validation")
    parser.add_argument('--other_dataset_path',dest='other_dataset_path',
                        action='append',type=str,default=[],
                        help="Path for tabular datasets")
    parser.add_argument('--model_path',dest='model_path',action='store',
                        type=str,default=None,
                        help="Path for MIL-CoMorI model")
    parser.add_argument('--ob',dest='ob',action='store',
                        type=int,default=0,
                        help="Objective index (if MIL-CoMorI is multi-objective)")
    parser.add_argument('--subset',dest='subset',action='store',
                        type=int,default=500,
                        help="Subset cells for prediction")
    parser.add_argument('--labels_path',dest='labels_path',
                        action='store',type=str,default=None,
                        help="Path for labels")
    parser.add_argument('--excluded_ids',dest='excluded_ids',
                        nargs='+',action='store',type=str,default=None,
                        help="Excluded dataset IDs (separated by spaces)")
    parser.add_argument('--threshold',dest='threshold',
                        action='store',type=float,default=0.25,
                        help="Threshold for VCT corr. (for consensus)")
    parser.add_argument('--n_cells_consensus',dest='n_cells_consensus',
                        action='store',type=int,default=250,
                        help="Number of cells to calculate consensus")

    args = parser.parse_args()

    if torch.cuda.is_available():
          dev = "cuda:0"
    else:
          dev = "cpu"

    n_d = len(args.dataset_path)

    # print("loading labels")
    label_dict = get_classes(args.labels_path)
    labels = Labels(label_dict)

    # print("loading state dict and inferring network architecture from it")
    state_dict = torch.load(args.model_path,map_location='cpu')
    network_state_dict = state_dict[0]['network']

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

    # print("assembling network and opening files, loading data moments")
    stacked_network = VirtualCellClassifierStack(
        n_input_features=n_features,
        n_virtual_cells=n_virtual_cells,
        other_datasets_sizes=[od_features[0]],
        n_classes=n_classes).to(dev)
    stacked_network.load_state_dict(network_state_dict)
    stacked_network.train(False)

    # print("loading datasets")
    all_datasets = [GenerateFromDataset(x,auto=False) for x in args.dataset_path]
    all_datasets = [GenerateFromDataset(x,auto=False) for x in args.dataset_path]
    if args.other_dataset_path is not None:
        all_other_datasets = []
        for x in args.other_dataset_path:
            x = CSVDataset(x,handle_nans='remove')
            x.normalize = False
            all_other_datasets.append(x)

        all_other_datasets_keys = [set(x.keys) for x in all_other_datasets]
    else:
        all_other_datasets = None

    all_sets = [set(x.keys) for x in all_datasets]
    all_sets.append(labels.keys)
    if args.other_dataset_path is not None:
        all_sets.extend(all_other_datasets_keys)
    all_keys = list(set.intersection(*all_sets))
    all_keys = [k for k in all_keys if k not in args.excluded_ids]

    # print("building consensus set of cells")
    thr = args.threshold
    all_folds = [k for k in state_dict if k != 'args']
    masks = []
    cells = []
    for dataset_obj in all_datasets:
        c = dataset_obj.generate_n_cells_flat(
            all_keys,args.n_cells_consensus)
        cells.append(c)

    vc_predictions = [[] for _ in all_datasets]
    for k in range(len(all_folds)):
        network_state_dict = state_dict[k]['network']
        stacked_network.load_state_dict(network_state_dict)
        stacked_network.train = False
        means = [np.squeeze(x) for x in state_dict[k]['means']]
        stds = [np.squeeze(x) for x in state_dict[k]['stds']]
        cells_norm = [(cells[i]-means[i])/stds[i] for i in range(len(cells))]
        cells_norm = [torch.Tensor(x) for x in cells_norm]
        cells_norm = [torch.unsqueeze(x,0) for x in cells_norm]
        cells_norm = [torch.unsqueeze(x,0) for x in cells_norm]
        vc = stacked_network.get_cells_individual(cells_norm)
        for i in range(len(vc)):
            x = vc[i].detach().numpy()
            vc_predictions[i].append(np.squeeze(x))

    all_masks = []
    for curr_fold in range(len(all_folds)):
        consensus = [[] for _ in all_datasets]
        for idx in range(len(vc_predictions)):
            vc_pred = vc_predictions[idx]
            all_D_min = []
            for k in range(len(all_folds)):
                if k != curr_fold:
                    D = cdist(vc_pred[curr_fold],vc_pred[k],'correlation')
                    D_min = np.min(D,axis=1)
                    all_D_min.append(D_min)
            all_D_min = np.stack(all_D_min,axis=1)
            thr = args.threshold
            x = np.where(all_D_min <= thr)[0]
            x = [i for i,c in zip(*np.unique(x,return_counts=True)) if c >= 3]
            consensus[idx].append(set(x))

        consensus = [list(set.intersection(*x)) for x in consensus]
        for x in consensus:
            mask = np.zeros(n_virtual_cells[0])
            mask[x] = 1.
            masks.append(torch.Tensor(mask))
        all_masks.append(masks)

    metrics_workhorse = MetricFunction(
        {'Accuracy':accuracy,'Confusion Matrix':confusion_matrix,
        'AUC':auc,'Deviance':deviance,})

    all_metrics = {}
    for curr_fold in range(len(all_folds)):
        metrics_workhorse.reset()
        masks = all_masks[curr_fold]
        network_state_dict = state_dict[curr_fold]['network']
        stacked_network.load_state_dict(network_state_dict)
        stacked_network.train = False
        means = [np.squeeze(x) for x in state_dict[curr_fold]['means']]
        stds = [np.squeeze(x) for x in state_dict[curr_fold]['stds']]
        means_od = [np.squeeze(x).astype(np.float32)
                    for x in state_dict[curr_fold]['means'][n_d:]]
        stds_od = [np.squeeze(x).astype(np.float32)
                for x in state_dict[curr_fold]['stds'][n_d:]]
        key_set = list(set.intersection(
            set(state_dict[curr_fold]['testing_set']),set(all_keys)))
        for key in tqdm(key_set):
            all_cell_sets = []
            for dataset_idx,dataset in enumerate(all_datasets):
                subset_size = args.subset
                n_cells = dataset[key].n_cells
                if subset_size > n_cells:
                    subset_size = n_cells
                subset_size = n_cells
                cells = dataset[key].return_all_cells()
                cells_transformed = (cells - means[dataset_idx])/stds[dataset_idx]
                cells_transformed = cells_transformed[np.newaxis,np.newaxis,:,:]
                cells_transformed = torch.Tensor(np.float32(cells_transformed)).to(dev)
                all_cell_sets.append(cells_transformed)

            truth = labels[[key]]

            # get data from other datasets
            other_datasets = [torch.Tensor(x[[key]]).to(dev) 
                            for x in all_other_datasets]
            other_datasets = [(x-means_od[i])/stds_od[i] 
                            for i,x in enumerate(other_datasets)]
            # predict here
            output_prob = stacked_network.forward_w_mask(
                [all_cell_sets,other_datasets],args.ob,masks)

            output_onehot = torch.zeros_like(output_prob)
            output_onehot[(torch.arange(0,output_prob.shape[0]).long(),
                        torch.argmax(output_prob,dim=-1))] = 1

            truth_onehot = torch.zeros_like(output_prob)
            truth_onehot[(torch.arange(0,output_prob.shape[0]).long(),
                        truth)] = 1
            if truth_onehot.shape[1] == 2:
                # coherces to a binary classification problem
                truth_onehot = truth_onehot[:,1:]
                output_onehot = output_onehot[:,1:]
                output_prob = output_prob[:,1:]
            metrics_workhorse.update(
                truth_onehot.to('cpu').detach(),
                output_onehot.to('cpu').detach(),
                output_prob.to('cpu').detach())

        M = metrics_workhorse.compute()
        for metric in M:
            for i,m in enumerate(M[metric].flatten()):
                n = '{}_{}'.format(metric,i)
                if n not in all_metrics:
                    all_metrics[n] = [float(m)]
                else:
                    all_metrics[n].append(float(m))

    for metric in all_metrics:
        m = all_metrics[metric]
        print('{},{},{},{}'.format(
            args.model_path,'nan',metric,float(np.mean(m))))