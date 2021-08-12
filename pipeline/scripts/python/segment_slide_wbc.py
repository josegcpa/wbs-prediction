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

import unet_utilities

F_size = 9
Filter = np.ones([F_size,F_size]) / (F_size**2)

def convolve_n(image,F,n=1,thr=0.6):
    for i in range(n):
        image = convolve(image.astype(np.float32),F) > thr
    return image.astype(np.bool)

def draw_hulls(image):
    contours,hierarchy = cv2.findContours(image.astype(np.uint8),2,1)
    cnt = contours[-1]
    hull = cv2.convexHull(cnt,returnPoints = True)
    defects = cv2.convexityDefects(
        cnt, cv2.convexHull(cnt,returnPoints = False))
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

def characterise_cell(image_labels_im_i):
    image,labels_im,i = image_labels_im_i
    m = 16
    sh = image.shape
    x,y = np.where(labels_im == i)
    R = x.min(),y.min(),x.max(),y.max()
    sx,sy = R[2]-R[0]+m*2,R[3]-R[1]+m*2
    cc = [R[0]-m,R[1]-m]
    cc.extend([cc[0]+sx,cc[1]+sy])
    features = None
    if np.any(
            [cc[0]<0,cc[1]<0,
             cc[2]>sh[0],cc[3]>sh[1],
             sx>128,sy>128]):
        pass
    else:
        x = x - cc[0]
        y = y - cc[1]
        S = len(x)
        if (S > 1000) and (S < 8000):
            sub_image = image[cc[0]:cc[2],cc[1]:cc[3],:]
            mask_binary_holes = np.zeros([sx,sy])
            mask_binary_holes[(x,y)] = 1
            mask_convolved = convolve_n(
                mask_binary_holes.astype(np.float32),
                Filter,3,thr=0.5)
            mask_hulls,cnt = draw_hulls(mask_convolved)
            mask_hulls = binary_fill_holes(mask_hulls)
            s = np.stack([mask_hulls for _ in range(3)],axis=2)
            features = MIA.wrapper_single_image(
                sub_image,mask_hulls,cnt)
    return features


def refine_prediction_characterise(image,mask):
    mask_binary = apply_hysteresis_threshold(mask,0.45,0.5)
    mask_binary_holes = binary_fill_holes(mask_binary)

    num_labels, labels_im = cv2.connectedComponents(
        np.uint8(mask_binary_holes))

    for i in range(1,num_labels):
        #try:
            features = characterise_cell([image,labels_im,i])
            if features is not None:
                yield features
        #except:
        #    pass

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

h,w,extra = 512,512,128
n_i = 1
inputs = tf.placeholder(tf.uint8,
                        [None,h+extra,w+extra,3],
                        name='Input_Tensor')
inputs = tf.image.convert_image_dtype(inputs,tf.float32)

#flipped_inputs = tf.image.flip_left_right(inputs)
inputs_ = []
for i in range(n_i):
    inputs_curr = inputs[i,:,:,:]
    inputs_.extend([
        inputs_curr,
        tf.image.rot90(inputs_curr,1),
        tf.image.rot90(inputs_curr,2),
        tf.image.rot90(inputs_curr,3)
    ])
inputs = tf.stack(inputs_,axis=0)

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

#flipped_prediction = tf.image.flip_left_right(
#  prediction_network[4:,:,:,:])
#prediction_network = prediction_network[:4,:,:,:]

pred_list = []
prediction_network_ = []
for i in range(n_i):
    pred_curr = prediction_network[(4*i):(4*i+4),:,:,:]
    pred_list.extend([
        pred_curr[0,:,:,:],
        tf.image.rot90(pred_curr[1,:,:,:],-1),
        tf.image.rot90(pred_curr[2,:,:,:],-2),
        tf.image.rot90(pred_curr[3,:,:,:],-3)
    ])
    pred_list = tf.stack(pred_list,axis=0)
    prediction_network_.append(
        tf.reduce_mean(pred_list,axis=0,keepdims=True)
    )
    pred_list = []

prediction_network = tf.concat(prediction_network_,axis=0)

igwq = ImageGeneratorWithQueue(args.slide_path,args.csv_path,
                               extra_padding=extra//2,
                               maxsize=args.n_processes_data)
igwq.start()

saver = tf.train.Saver()
X = []
Y = []
F = h5py.File(args.output_path,mode='w')
N = 0
i = 0
times = []
Centers = {}

with tf.Session() as sess:
    saver.restore(sess,args.checkpoint_path)
    images = []
    for image,coords in igwq.generate():
        i += 1
        images.append(image)
        if len(images) == n_i:
            images = np.stack(images,axis=0)
            A = time.time()
            segmented_images,images = sess.run(
                [prediction_network,inputs],
                feed_dict={'Input_Tensor:0':images})
            for image_idx in range(n_i):
                image = images[image_idx,:,:,:]
                segmented_image = segmented_images[image_idx,:,:,:]
                for obj in refine_prediction_characterise(
                        image,segmented_image):
                    y,x = obj['x'],obj['y']
                    if np.any([
                            x.min() == 0,
                            x.max() == (h+extra),
                            y.min() == 0,
                            y.max() == (w+extra)
                    ]):
                        pass
                    else:
                        features = []
                        for k in MIA_FEATURES:
                            features.append(obj[k])
                        x += coords[0]
                        y += coords[1]
                        C = str([np.mean(x),np.mean(y)])
                        if C not in Centers:
                            Centers[C] = 1
                            N += 1
                            g = F.create_group(str(C))
                            g.create_dataset(
                                "X",x.shape,dtype=np.int32,data=x)
                            g.create_dataset(
                                "Y",y.shape,dtype=np.int32,data=y)
                            g.create_dataset(
                                "features",[len(features)],
                                dtype=np.float64,data=features)
            B = time.time()
            times.append(B-A)
            images = []
        if i % 100 == 0:
            sys.stdout.write('Image {}. '.format(i))
            sys.stdout.write(
                'Av. time/image={:.3f}s '.format(
                    float(np.mean(times))/n_i))
            sys.stdout.write(
                '(for last 100 images={:.3f}s); '.format(
                    float(np.mean(times[-100:]))/n_i))
            sys.stdout.write(
                'No. detected cells={}'.format(N))
            sys.stdout.write('\n')

F.attrs['NCells'] = N
F.attrs['NImages'] = i
F.close()
