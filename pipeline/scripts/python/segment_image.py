import sys
import numpy as np
import argparse
import time
import h5py
import openslide
import cv2
from PIL import Image
from scipy.ndimage import convolve
from skimage.filters import apply_hysteresis_threshold
from skimage import io
from scipy.ndimage.morphology import binary_fill_holes
import tensorflow as tf

from image_generator import ImageGeneratorWithQueue
import MIA
MIA_FEATURES = MIA.MIA_FEATURES

sys.path.append('/nfs/research1/gerstung/josegcpa/projects/01IMAGE/tf_implementations')

import unet_utilities

def convolve_n(image,F,n=1,thr=0.6):
    for i in range(n):
        image = convolve(image.astype(np.float32),F) > thr
    return image.astype(np.bool)

def draw_hulls(image):
    contours,hierarchy = cv2.findContours(image.astype(np.uint8),2,1)
    cnt = contours[-1]
    hull = cv2.convexHull(cnt,returnPoints = True)
    defects = cv2.convexityDefects(cnt, cv2.convexHull(cnt,returnPoints = False))
    if defects is not None:
        defects = defects[defects[:,0,-1]>2000,:,:]
        if defects.size > 0:
            for defect in defects:
                a = tuple([x for x in cnt[defect[0,0],0,:]])
                b = tuple([x for x in cnt[defect[0,1],0,:]])
                output = cv2.line(image.astype(np.uint8),a,b,1)
                output = binary_fill_holes(output)
        else:
            output = image
    else:
        output = image
    return output,cnt

def refine_prediction(mask):
    F_size = 9
    Filter = np.ones([F_size,F_size]) / (F_size**2)
    mask = mask[0,:,:,0]
    mask_binary = apply_hysteresis_threshold(mask,0.45,0.5)
    mask_binary_holes = binary_fill_holes(mask_binary)

    num_labels, labels_im = cv2.connectedComponents(
        np.uint8(mask_binary_holes))
    for i in range(1,num_labels):
        try:
            mask_binary_holes = np.zeros_like(labels_im)
            mask_binary_holes[labels_im == i] = 1
            S = mask_binary_holes.sum()
            if (S > 1000) and (S < 10000):
                mask_convolved = convolve_n(
                    mask_binary_holes.astype(np.float32),
                    Filter,3,thr=0.5)
                mask_hulls,cnt = draw_hulls(mask_convolved)
                mask_hulls = binary_fill_holes(mask_hulls)
                yield np.where(mask_hulls > 0)
        except:
            pass

def image_to_patches(image,size=(512,512),stride=None):
    if stride is None:
        stride = size
    x,y,_ = image.shape
    for i in range(0,x,stride[0]):
        if i + size[0] > x:
            i = x - size[0]
        for j in range(0,y,stride[1]):
            if j + size[1] > y:
                j = y - size[1]
            a,b = i,i+size[0]
            c,d = j,j+size[1]
            yield image[a:b,c:d,:],((a,b),(c,d))

parser = argparse.ArgumentParser(
      prog = 'u-net.py',
      description = 'Multi-purpose U-Net implementation.')

parser.add_argument('--image_path',dest='image_path',
                    action='store',
                    default=None)
parser.add_argument('--checkpoint_path',dest='checkpoint_path',
                    action='store',
                    default=None,
                    help='Path to U-Net checkpoint.')
parser.add_argument('--depth',dest='depth',
                    action='store',
                    type=float,
                    default=1.0,
                    help='Depth of the U-Net.')
parser.add_argument('--squeeze_and_excite',dest='squeeze_and_excite',
                    action='store_true',
                    default=False,
                    help='Whether squeeze-and-excite layers should be used.')
parser.add_argument('--output_path',dest='output_path',
                    action='store',
                    default=False,
                    help='Output path for the mask.')
parser.add_argument('--n_processes_data',dest='n_processes_data',
                    action='store',
                    type=int,
                    default=1,
                    help='Number of processes for data loading.')

args = parser.parse_args()

inputs = tf.placeholder(tf.uint8,
                        [None,128,128,3],
                        name='Input_Tensor')
inputs = tf.image.convert_image_dtype(inputs,tf.float32)

flipped_inputs = tf.image.flip_left_right(inputs)
inputs = tf.concat(
    [inputs,
     tf.image.rot90(inputs,1),
     tf.image.rot90(inputs,2),
     tf.image.rot90(inputs,3)],
  axis=0
)

unet = unet_utilities.u_net(
  inputs=inputs,
  padding='SAME',
  is_training=False,
  depth_mult=args.depth,
  squeeze_and_excite=args.squeeze_and_excite
)[0]
prediction_network = tf.expand_dims(
  tf.nn.softmax(unet,axis=-1)[:,:,:,1],
  axis=-1)

pred_list = [
  prediction_network[0,:,:,:],
  tf.image.rot90(prediction_network[1,:,:,:],-1),
  tf.image.rot90(prediction_network[2,:,:,:],-2),
  tf.image.rot90(prediction_network[3,:,:,:],-3)
]

prediction_network = tf.stack(pred_list,axis=0)
prediction_network = tf.reduce_mean(prediction_network,
                                    axis=0,
                                    keepdims=True)

image_input = np.array(Image.open(args.image_path))[:,:,:]
output = np.zeros(image_input.shape[:2])
denominator = np.zeros(image_input.shape[:2])
size = [128,128]
stride = [64,64]
saver = tf.train.Saver()

with tf.Session() as sess:
    saver.restore(sess,args.checkpoint_path)
    for image,coords in image_to_patches(image_input,size=size,stride=stride):
        (a,b),(c,d) = coords
        A = time.time()
        segmented_image,image = sess.run(
            [prediction_network,inputs],
            feed_dict={'Input_Tensor:0':image[np.newaxis,:,:,:]})
        denominator[a:b,c:d] +=1
        output[a:b,c:d] += segmented_image[0,:,:,0]

output = output / denominator

output_ = np.zeros_like(output)
for x,y in refine_prediction(output[np.newaxis,:,:,np.newaxis]):
    if np.any([np.any(x==0),np.any(x==image.shape[0]),
               np.any(y==0),np.any(y==image.shape[1])]) == True:
        pass
    else:
        output_[x,y] = 1

io.imsave(args.output_path,output_)
