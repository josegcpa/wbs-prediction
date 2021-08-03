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

def refine_prediction_characterise(image,mask):
    F_size = 9
    Filter = np.ones([F_size,F_size]) / (F_size**2)
    image = image[0,:,:,:]
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
            if (S > 1000) and (S < 8000):
                mask_convolved = convolve_n(
                    mask_binary_holes.astype(np.float32),
                    Filter,3,thr=0.5)
                mask_hulls,cnt = draw_hulls(mask_convolved)
                mask_hulls = binary_fill_holes(mask_hulls)

                features = MIA.wrapper_single_image(image,mask_hulls,cnt)
                yield features
        except:
            pass

parser = argparse.ArgumentParser(
      prog = 'u-net.py',
      description = 'Multi-purpose U-Net implementation.'
)

parser.add_argument('--csv_path',dest='csv_path',
                    action='store',
                    default=None,
                    help='Path to CSV file with the quality file.')
parser.add_argument('--slide_path',dest='slide_path',
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
                    help='Output path for the hdf5 file.')
parser.add_argument('--n_processes_data',dest='n_processes_data',
                    action='store',
                    type=int,
                    default=1,
                    help='Number of processes for data loading.')

args = parser.parse_args()

inputs = tf.placeholder(tf.uint8,
                        [None,512+256,512+256,3],
                        name='Input_Tensor')
inputs = tf.image.convert_image_dtype(inputs,tf.float32)

flipped_inputs = tf.image.flip_left_right(inputs)
inputs = tf.concat(
    [inputs,
     tf.image.rot90(inputs,1),
     tf.image.rot90(inputs,2),
     tf.image.rot90(inputs,3),
     #flipped_inputs,
     #tf.image.rot90(flipped_inputs,1),
     #tf.image.rot90(flipped_inputs,2),
     #tf.image.rot90(flipped_inputs,3)
     ],
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

flipped_prediction = tf.image.flip_left_right(
  prediction_network[4:,:,:,:])
prediction_network = prediction_network[:4,:,:,:]

pred_list = [
  prediction_network[0,:,:,:],
  tf.image.rot90(prediction_network[1,:,:,:],-1),
  tf.image.rot90(prediction_network[2,:,:,:],-2),
  tf.image.rot90(prediction_network[3,:,:,:],-3),
  #flipped_prediction[0,:,:,:],
  #tf.image.rot90(flipped_prediction[1,:,:,:],1),
  #tf.image.rot90(flipped_prediction[2,:,:,:],-2),
  #tf.image.rot90(flipped_prediction[3,:,:,:],-1)
]

prediction_network = tf.stack(pred_list,axis=0)
prediction_network = tf.reduce_mean(prediction_network,
                                    axis=0,
                                    keepdims=True)

igwq = ImageGeneratorWithQueue(args.slide_path,args.csv_path,
                               maxsize=args.n_processes_data)
igwq.start()

saver = tf.train.Saver()
X = []
Y = []
F = h5py.File(args.output_path,mode='w')
N = 0
i = 0

with tf.Session() as sess:
    saver.restore(sess,args.checkpoint_path)
    for image,coords in igwq.generate():
        i += 1
        A = time.time()
        segmented_image,image = sess.run(
            [prediction_network,inputs],
            feed_dict={'Input_Tensor:0':image[np.newaxis,:,:,:]})
        for obj in refine_prediction_characterise(image,segmented_image):
            y,x = obj['x'],obj['y']
            if np.any([
                    np.any(x == 0),
                    np.any(x == (512+256)),
                    np.any(y == 0),
                    np.any(y == (512+256))
            ]):
                pass
            else:
                features = []
                for k in MIA_FEATURES:
                    features.append(obj[k])
                x += coords[0]
                y += coords[1]
                C = str([np.mean(x),np.mean(y)])
                if C not in F:
                    N += 1
                    g = F.create_group(str(C))
                    g.create_dataset("X",x.shape,dtype=np.int32,data=x)
                    g.create_dataset("Y",y.shape,dtype=np.int32,data=y)
                    g.create_dataset("features",[len(features)],
                                     dtype=np.float64,data=features)
        B = time.time()
        if i % 100 == 0:
            sys.stderr.write('Image {}. '.format(i))
            sys.stderr.write('{} {}'.format(B-A,N))
            sys.stderr.write('\n')

F.attrs['NCells'] = N
F.attrs['NImages'] = i
F.close()
