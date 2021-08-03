import sys
import numpy as np
import argparse
import time
import h5py
import openslide
import cv2
from PIL import Image

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
      description = 'Characterise dataset of WBC nuclei.'
)

parser.add_argument('--slide_path',dest='slide_path',
                    action='store',
                    default=None)
parser.add_argument('--feature_file_path',dest='feature_file_path',
                    action='store',
                    default=None)
parser.add_argument('--output_file_path',dest='output_file_path',
                    action='store',
                    default=None)

args = parser.parse_args()

F = h5py.File(args.feature_file_path,'r')
O = h5py.File(args.output_file_path,'w')
OS = openslide.OpenSlide(args.slide_path)

for i,C in enumerate(F):
    if i % 100 == 0:
        print(i)
    if C not in O:
        X = F[C]['X']
        Y = F[C]['Y']
        x,y = np.min(X)-5,np.min(Y)-5
        h,w = np.max(X)-x+10,np.max(Y)-y+10

        cell_image = np.array(OS.read_region((x,y),0,(h,w)))
        mask = np.zeros([w,h],dtype=np.uint8)
        mask[Y-y,X-x] = 1
        nucleus_mask,nucleus_cnt = segment_nucleus(
            np.array(cell_image),mask)
        nucleus_features = wrapper_single_image_separate(
            cell_image,nucleus_mask,nucleus_cnt)
        feature_array = np.zeros([len(MIA_FEATURES_NUCLEUS)],
                                 dtype=np.float64)
        for i,f in enumerate(MIA_FEATURES_NUCLEUS):
            feature_array[i] = nucleus_features[f]
        g = O.create_group(C)
        g["features_nucleus"] = feature_array
