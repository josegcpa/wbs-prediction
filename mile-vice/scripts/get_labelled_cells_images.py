import argparse
import h5py
import openslide
import re
import os
import numpy as np
from scipy.spatial.distance import cdist
from tqdm import tqdm

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Get images for the expert annotated cell collections.')

    parser.add_argument('--centers_file',dest='centers_file',
                        action='store',type=str,default=None)
    parser.add_argument('--labels_file',dest='labels_file',
                        action='store',type=str,default=None)
    parser.add_argument('--output_path',dest='output_path',
                        action='store',type=str,default=None)

    parser.add_argument('--slides_path',dest='slide_path',
                        action='store',type=str,default=None)
    parser.add_argument('--aggregates_path',dest='aggregates_path',
                        action='store',type=str,default=None)
    parser.add_argument('--segmented_path',dest='segmented_path',
                        action='store',type=str,default=None)

    args = parser.parse_args()

    try:
        with open(args.feature_subset,"r") as o:
            feature_subset = [
                int(x.strip())-1 for x in o.read().strip().split(',')]
    except:
        feature_subset = []
    try:
        os.makedirs(args.output_path)
    except:
        pass

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
            c = aggregates['cell_centers'][c].decode()
            x,y = segmentations[c]['X'][()],segmentations[c]['Y'][()]
            A,B = x.min()-8,y.min()-8
            C,D = x.max()+8,y.max()+8
            h,w = C-A,D-B
            c_idx = centers[slide_id][a]['cell_idx']
            im = slide.read_region([int(A),int(B)],0,[int(h),int(w)])
            output_path = os.path.join(args.output_path,'{}.png'.format(c_idx))
            im.save(output_path)
