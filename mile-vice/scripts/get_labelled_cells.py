import argparse
import torch
import h5py
import openslide
import re
from scipy.spatial.distance import cdist
from tqdm import tqdm

from networks import *
from data_generator import *
from metrics import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Test virtual cell classifier.')

    parser.add_argument('--cell_type',dest='cell_type',
                        action='store',type=str,default='wbc')

    parser.add_argument('--centers_file',dest='centers_file',
                        action='store',type=str,default=None)
    parser.add_argument('--labels_file',dest='labels_file',
                        action='store',type=str,default=None)
    parser.add_argument('--model_path',dest='model_path',
                        action='store',type=str,default=None)
    parser.add_argument('--feature_subset',dest='feature_subset',
                        default=None,
                        action='store')

    parser.add_argument('--slides_path',dest='slide_path',
                        action='store',type=str,default=None)
    parser.add_argument('--aggregates_path',dest='aggregates_path',
                        action='store',type=str,default=None)
    parser.add_argument('--segmented_path',dest='segmented_path',
                        action='store',type=str,default=None)

    parser.add_argument('--fold',dest='fold',
                        action='store',type=int,default=0)

    args = parser.parse_args()

    if torch.cuda.is_available():
          dev = "cuda:0"
    else:
          dev = "cpu"

    try:
        with open(args.feature_subset,"r") as o:
            feature_subset = [
                int(x.strip())-1 for x in o.read().strip().split(',')]
    except:
        feature_subset = []
    if len(feature_subset) == 0:
        feature_subset = None

    with open(args.centers_file) as o:
        lines = [x.strip().split(',') for x in o.readlines()]
        centers = {}
        for i,s,x,y in lines:
            if s not in centers:
                centers[s] = [{'cell_idx':i,'center':(x,y)}]
            else:
                centers[s].append(
                    {'cell_idx':i,'center':(float(x),float(y))})

    with open(args.labels_file) as o:
        lines = [x.strip().split(',') for x in o.readlines()]
        labels = {}
        for x in lines:
            if x[1] not in labels:
                labels[x[1]] = {'labels':[x[2]],'user_idx':[x[3]]}
            else:
                labels[x[1]]['labels'].append(x[2])
                labels[x[1]]['user_idx'].append(x[3])

    ni = 0 if args.cell_type == 'wbc' else 1
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
    od_features = [int(final_layers[0].shape[-1] - np.sum(n_virtual_cells))]

    stacked_network = VirtualCellClassifierStack(
        n_input_features=n_features,
        n_virtual_cells=n_virtual_cells,
        other_datasets_sizes=od_features,
        n_classes=n_classes).to(dev)

    stacked_network.load_state_dict(network_state_dict)

    means = np.squeeze(state_dict[args.fold]['means'][ni])
    stds = np.squeeze(state_dict[args.fold]['stds'][ni])

    vcq = stacked_network.vcq_list[ni]
    TT = tqdm()
    for slide_id in centers:
        TT.update()
        slide_path = '{}/{}.tiff'.format(args.slide_path,slide_id)
        aggregates_path = '{}/{}.h5'.format(args.aggregates_path,slide_id)
        segmented_path = '{}/{}.h5'.format(args.segmented_path,slide_id)
        slide = openslide.OpenSlide(slide_path)
        aggregates = h5py.File(aggregates_path,'r')
        segmentations = h5py.File(segmented_path,'r')

        cell_centers = aggregates['cell_centers'][()].tolist()
        cell_centers = [x.decode() for x in cell_centers]
        cell_centers_num = [x[1:-1].split(',') for x in cell_centers]
        cell_centers_num = [
            [float(x[1]),float(x[0])] for x in cell_centers_num]
        centers_num = np.array(cell_centers_num)

        labelled_centers_num = np.array(
            [x['center'] for x in centers[slide_id]])
        dist_mat = cdist(centers_num,labelled_centers_num)
        dist_mat = np.where(dist_mat < 5,True,False)
        centers_original,centers_labelled = np.where(dist_mat)
        for c,a in zip(centers_original,centers_labelled):
            cells_idx = '0'
            if c > aggregates['cells'][cells_idx].shape[0]:
                c = c - aggregates['cells'][cells_idx].shape[0]
                cells_idx = '1'
            features = aggregates['cells'][cells_idx][c,:]
            if feature_subset is not None:
                features = features[feature_subset]
            features = (features - means)/stds
            features = features[np.newaxis,np.newaxis,np.newaxis,:]
            features = torch.Tensor(features)
            vc = vcq.get_virtual_cells(features).detach().cpu().numpy()
            vc_ = np.where(vc > 0.5,1,0)
            if np.any(vc_ > 0):
                vc_max = np.argmax(vc_)
                vc_prob = vc[0,vc_max]
                c_idx = centers[slide_id][a]['cell_idx']
                if c_idx in labels:
                    label = labels[c_idx]
                    for l,u in zip(label['labels'],label['user_idx']):
                        print('{},{},{},{},{},{},{},{}'.format(
                            c_idx,slide_id,
                            centers[slide_id][a]['center'][0],
                            centers[slide_id][a]['center'][1],
                            vc_max,vc_prob,l,u))
