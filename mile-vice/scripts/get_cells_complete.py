print('-'*20)
import argparse
import torch
import h5py
import openslide
import re
from tqdm import tqdm

from networks import *
from data_generator import *
from metrics import *

if __name__ == "__main__":
    print("Reading cmd line arguments...")
    parser = argparse.ArgumentParser(
        description='''Predicts VCT for files from pipeline_tf2 output,
        gathers image of each cell from slide and coordinates''')

    parser.add_argument('--network_idx',dest='network_idx',
                        action='store',type=int,default=0,
                        help="Cell type index (wbc=0,rbc=1)")
    parser.add_argument('--slide_path',dest='slide_path',
                        action='store',type=str,default=None,
                        help="Path for the slide")
    parser.add_argument('--aggregates_path',dest='aggregates_path',
                        action='store',type=str,default=None,
                        help="Path for hdf5 containing cell characteristics")
    parser.add_argument('--segmented_path',dest='segmented_path',
                        action='store',type=str,default=None,
                        help="Path for hdf5 containing cell segmentations")
    parser.add_argument('--model_path',dest='model_path',
                        action='store',type=str,default=None,
                        help="Path for MILe-ViCe model")
    parser.add_argument('--output_path',dest='output_path',
                        action='store',type=str,default=None,
                        help="HDF5 output path for VCT class + images")
    parser.add_argument('--fold',dest='fold',
                        action='store',type=int,default=0,
                        help="Model fold")
    parser.add_argument('--subset',dest='subset',
                        action='store',type=int,default=100,
                        help="How many cells to subset")
    parser.add_argument('--flip',dest='flip',
                        default=False,action='store_true',
                        help="Whether the coordinates should be flipped")
    parser.add_argument('--feature_subset',dest='feature_subset',
                        default=None,action='store',
                        help="File to feature subset")
    parser.add_argument('--size',dest='size',type=int,
                        default=None,action='store',
                        help="Size of the images to be stored")

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

    if args.size is not None:
        S = args.size
        HS = args.size//2

    ni = args.network_idx
    
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
    
    slide = openslide.OpenSlide(args.slide_path)
    aggregates = h5py.File(args.aggregates_path,'r')
    segmentations = h5py.File(args.segmented_path,'r')
    output = h5py.File(args.output_path,'w')

    vcq = stacked_network.vcq_list[ni]

    subset_size = args.subset
    
    all_cells = []
    all_cell_centers = aggregates['cell_centers'][()]
    for cell_group in aggregates['cells']:
        all_cells.append(aggregates['cells'][cell_group][()])
    all_cells = np.concatenate(all_cells,axis=0)
    if feature_subset is not None:
        all_cells = all_cells[:,feature_subset]
    cell_transformed = (all_cells-means[ni])/stds[ni]
    cell_transformed = cell_transformed[np.newaxis,np.newaxis,:,:]
    cell_transformed = torch.Tensor(np.float32(cell_transformed))
    all_vc = vcq.f(cell_transformed)
    cell_proportions = np.squeeze(vcq.g(all_vc).detach().numpy())
    all_vc = all_vc.to('cpu').detach().numpy()
    all_vc = all_vc[0,:,:,0].T

    output['cell_proportions'] = cell_proportions
    cell_idx_subset = np.random.choice(
        all_cells.shape[0],np.minimum(subset_size,all_cells.shape[0]),
        replace=False)

    cell_group_h5 = output.create_group('cells')
    for cell_idx in cell_idx_subset:
        cell_center = all_cell_centers[cell_idx]
        vc = all_vc[cell_idx,:]
        cell = all_cells[cell_idx,:]
        cell_center_int = cell_center.decode('utf-8')[1:-1].split(',')
        cell_center_int = [float(cell_center_int[0]),float(cell_center_int[1])]
        cell_coordinates = segmentations[cell_center]
        x = cell_coordinates['X'][::]
        y = cell_coordinates['Y'][::]
        if args.size is None:
            location = [x.min()-3,y.min()-3]
            size = [x.max()-location[0] + 6,y.max()-location[1] + 6]
        else:
            cx,cy = (x.max()+x.min())/2,(y.max()+y.min())/2
            cx,cy = [np.round(el) for el in [cx,cy]]
            location = int(cx-HS),int(cy-HS)
            size = S,S
        if args.flip == True:
            location = [location[1],location[0]]
            size = [size[1],size[0]]
        image = np.array(slide.read_region(location,0,size))[:,:,:3]
        g = cell_group_h5.create_group(cell_center)
        g.create_dataset('image',data=image,dtype=np.uint8)
        g.create_dataset('cell_type',data=vc,dtype=np.float32)
        g.create_dataset('features',data=cell,dtype=np.float64)
        g.create_dataset('cell_center_y',data=cell_center_int[0])
        g.create_dataset('cell_center_x',data=cell_center_int[1])
