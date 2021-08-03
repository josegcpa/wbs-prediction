import sys
import argparse
import os
import numpy as np
import time
import openslide
import h5py
import cv2

import MIA
wrapper_single_image_separate = MIA.wrapper_single_image_separate
MIA_FEATURES = MIA.MIA_FEATURES
MIA_FEATURES_NUCLEUS = MIA.MIA_FEATURES_NUCLEUS

def segment_nucleus(image,mask):
    idxs = np.where(mask>0)
    pixel_values = image[idxs]
    pixel_values_ = np.float32(pixel_values)
    _, labels, (centers) = cv2.kmeans(
        pixel_values_, 2, None,
        (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 100, 0.2),
        10, cv2.KMEANS_RANDOM_CENTERS)
    nucleus_center = np.argmin(np.mean(centers,axis=-1))

    label_mask = np.zeros_like(mask)
    label_mask[idxs] = labels[:,0] + 1
    label_mask = np.where(label_mask == (nucleus_center+1),1,0)
    label_mask = label_mask.astype(np.uint8)
    labels, stats = cv2.connectedComponentsWithStats(label_mask, 4)[1:3]
    largest_label = 1 + np.where(stats[1:, cv2.CC_STAT_AREA] > 100)[0]
    label_mask = np.where(np.isin(labels,largest_label),1,0).astype(np.uint8)
    cnt = cv2.findContours(label_mask,cv2.RETR_EXTERNAL,
                           cv2.CHAIN_APPROX_NONE)[0]
    return label_mask,cnt

parser = argparse.ArgumentParser(
    prog = 'aggregate_hdf5.py',
    description = 'creates dataset and summary statistics for hdf5 file'
)

parser.add_argument('--hdf5_path',dest='hdf5_path',
                    action='store',
                    default=None,
                    help='Path to hdf5 file with segmentation.')
parser.add_argument('--slide_path',dest='slide_path',
                    action='store',
                    default=None,
                    help='Path to slide.')
parser.add_argument('--output_path',dest='output_path',
                    action='store',
                    default=None,
                    help='Output path for the hdf5 file.')

args = parser.parse_args()

OS = openslide.OpenSlide(args.slide_path)
F = h5py.File(args.hdf5_path,mode='r')
F_out = h5py.File(args.output_path,mode='w')
F_out_cells = F_out.create_group('cells')
N_accumulator = 0
sum_accumulator = np.zeros([len(MIA_FEATURES)+len(MIA_FEATURES_NUCLEUS)],
                           dtype=np.float64)
squared_sum_accumulator = np.zeros([len(MIA_FEATURES)+len(MIA_FEATURES_NUCLEUS)],
                                   dtype=np.float64)
cells = []
nuclei = []
cell_centers = []
N_groups = 0
for item in F:
    cell = F[item]['features']
    ### process to segment the nucleus ###
    X = F[item]['X']
    Y = F[item]['Y']
    x,y = np.min(X)-5,np.min(Y)-5
    h,w = np.max(X)-x+10,np.max(Y)-y+10
    cell_image = np.array(OS.read_region((x,y),0,(h,w)))
    mask = np.zeros([w,h],dtype=np.uint8)
    mask[Y-y,X-x] = 1
    features = []
    for i,feature in enumerate(MIA_FEATURES):
        features.append(cell[i])
    if np.any(np.isnan(features)): pass
    elif np.any(np.isinf(features)): pass
    else:
        nucleus_mask,nucleus_cnt = segment_nucleus(
            np.array(cell_image),mask)
        if len(nucleus_cnt) > 0:
            nucleus = wrapper_single_image_separate(
                cell_image,nucleus_mask,nucleus_cnt)
            for i,feature in enumerate(MIA_FEATURES_NUCLEUS):
                features.append(nucleus[feature])
            features = np.array(features,dtype=np.float64)
            cell_centers.append(item)
            N_accumulator += 1
            sum_accumulator += features
            squared_sum_accumulator += np.square(features)
            cells.append(features)

    if len(cells) >= 100000:
        cells = np.array(cells)
        F_out_cells.create_dataset(
            str(N_groups),cells.shape,
            dtype=np.float64,data=cells)
        cells = []
        N_groups += 1

if len(cells) > 0:
    cells = np.array(cells)
    F_out_cells.create_dataset(
        str(N_groups),cells.shape,
        dtype=np.float64,data=cells)

means = sum_accumulator / N_accumulator
variances = squared_sum_accumulator / N_accumulator - means

F_out.create_dataset(
    'mean',means.shape,
    dtype=np.float64,data=means)
F_out.create_dataset(
    'variance',variances.shape,
    dtype=np.float64,data=variances)
F_out.create_dataset(
    'cell_centers',
    data=np.array(cell_centers,dtype="S"))

F.close()
F_out.close()
