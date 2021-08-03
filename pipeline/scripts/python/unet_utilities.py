import os
import numpy as np
from math import inf
import cv2
import tifffile as tiff
import tensorflow as tf
import h5py
from scipy.spatial import distance
from PIL import Image

import tf_da
from data_generators import *

try:
    tf.logging.set_verbosity(tf.logging.ERROR)
except:
    pass

slim = tf.contrib.slim

"""
Deep learning/TF-related operations.
"""

def safe_log(tensor):
    """
    Prevents log(0)

    Arguments:
    * tensor - tensor
    """
    return tf.log(tf.clip_by_value(tensor,1e-32,tf.reduce_max(tensor)))

def variables(vs):
    """
    For logging purposes.

    Arguments:
    * vs - variables from tf.all_variables() or similar functions
    """
    return int(np.sum([np.prod(v.get_shape().as_list()) for v in vs]))

def u_net(inputs,
          final_endpoint = None,
          padding = 'VALID',
          factorization = False,
          beta = 0,
          residuals = False,
          n_classes = 2,
          depth_mult = 1,
          is_training = True,
          aux_node = False,
          squeeze_and_excite = False):

    """
    Implementation of a standard U-net with some tweaks, namely:
        * Possibility of choosing whether padding should be SAME or VALID (for
        SAME it is necessary that the height % 16 == 0 and width % 16 == 0)
        * Possibility of factorizing the convolutions (a 3x3 convolution
        becomes a 1x3 and a 3x1 convolutions, turning the number of
        operations/filter from 9 to 6)
        * beta is the beta parameter in l2 regularization
        * residuals activates residual blocks in the links

    Arguments [default]:
    * inputs - input tensor;
    * final_endpoint - final layer to return ['Final']
    * padding - whether VALID or SAME padding should be used ['VALID']
    * factorization - whether convolutions should be factorized (3x3conv become
    sequential 1x3 and 3x1 convolutions) [False]
    * beta - L2 regularization factor [0]
    * residuals - whether to use residual linkers in the shortcut links [False]
    * n_classes - number of classes in the output layer (only works with 2) [2]
    * depth_mult - factor to increase or decrease the depth in each layer
    """

    def pixel_normalization(inputs):
        """
        The pixel normalization presented in the ProGAN paper. This
        normalization is basically a local response normalization performed
        depth-wise. It offers the advantage of having no learnable
        parameters, which makes training faster, and allows the usage of
        the WGAN-GP loss function.
        This is used to normalize the layers from both the generator (as
        per the ProGAN paper) and the discriminator (as the WGAN-GP paper
        suggests).
        """
        norm_factor = tf.sqrt(
            tf.reduce_mean(tf.square(inputs),axis=-1) + 1e-16
        )
        norm_factor = tf.expand_dims(norm_factor,-1)
        return inputs / norm_factor

    def conv2d(x,depth,size,stride,scope,
               factorization = False,padding = padding):
        if factorization == False:
            scope = scope + str(size) + 'x' + str(size)
            x = slim.conv2d(x,depth,[size,size],stride = stride,
                            scope = scope,padding = padding)
        else:
            scope_1 = scope + str(size) + 'x' + '1'
            scope_2 = scope + '1' + 'x' + str(size)
            x = slim.conv2d(x,depth,[size,1],stride = 1,scope = scope_1,
                            padding = padding)
            x = slim.conv2d(x,depth,[1,size],stride = 1,scope = scope_2,
                            padding = padding)
            if stride > 1:
                scope = scope + 'maxpool2d_' + str(stride) + 'x' + str(stride)
                x = slim.max_pool2d(x, [stride,stride], stride = stride,
                                    padding = 'SAME',scope = scope)
        x = slim.dropout(x,keep_prob=0.8,is_training=is_training)
        return x

    def block(x,depth,size,stride,padding,factorization = False):
        x = conv2d(x,depth,size,stride,'conv2d_1_',factorization)
        x = conv2d(x,depth,size,stride,'conv2d_2_',factorization)
        return x

    def red_block_wrapper(net,depth,factorization,residuals,endpoints,
                          name):
        net = block(net,depth,3,1,factorization)
        if residuals == True:
            endpoints[name] = residual_block(
                net,depth,3,factorization
            )
        else:
            endpoints[name] = net

        net = slim.max_pool2d(net,[2,2],stride = 2,
                              scope = 'maxpool2d_2x2')
        return net,endpoints

    def rec_block_wrapper(net,depth,factorization,
                          residuals,endpoints,padding,
                          name,red_block_name):


        if padding == 'VALID':
            new_size = net.get_shape().as_list()
            new_size = [new_size[1] * 2,new_size[2] * 2]
            net = tf.image.resize_nearest_neighbor(net,new_size)
            net = slim.conv2d(net,depth,[2,2],scope = 'conv2d_0_',
                              padding = 'SAME')
            n = endpoints[red_block_name].get_shape().as_list()[1]
            n = (n - net.get_shape().as_list()[1])/2
            concat_net = crop(endpoints[red_block_name],n)
        else:
            net = slim.conv2d_transpose(net,depth,[2,2],stride=2)
            concat_net = endpoints[red_block_name]
        net = tf.concat([net,concat_net],
                        axis = 3)
        net = block(net,depth,3,1,padding,factorization)
        endpoints[name] = net

        return net,endpoints

    def crop(concat_net,n):
        n = int(n)
        slice_size = concat_net.get_shape().as_list()
        slice_size = [
            - 1,
            int(slice_size[1] - n * 2),
            int(slice_size[2] - n * 2),
            slice_size[3]
        ]

        concat_net = tf.slice(concat_net,
                              [0,n,n,0],
                              slice_size
                              )
        return concat_net

    def fc_layer_classification(net,depth,n_classes,
                                **kwargs):
        # This node is used to classify the presence/absence of objects
        # in the bottleneck. This forces the network to collect
        # relevant features along the whole network.
        with tf.variable_scope('Aux_Node',None,[net]):
            with slim.arg_scope([slim.fully_connected],
                                **kwargs):
                pre_class = slim.fully_connected(
                    net,
                    num_outputs=depth,
                    activation_fn=None,
                    scope='pre_class'
                )
                classification = slim.fully_connected(
                    pre_class,
                    num_outputs=int(n_classes),
                    activation_fn=None,
                    scope='classification'
                )
        return classification

    def residual_block(input_n,depth,size,factorization):
        r = conv2d(input_n,depth,size,1,'res_conv2d_0_',
                   factorization = factorization,padding = 'SAME')
        r = conv2d(r,depth,size,1,'res_conv2d_1_',
                   factorization = factorization,padding = 'SAME')
        return tf.add(input_n,r)

    def c_squeeze_and_excite(net,r=1):
        with tf.variable_scope('ChannelSqueezeAndExcite'):
            squeeze = tf.reduce_mean(net,axis=[1,2])
            squeeze = tf.layers.dense(squeeze,
                                      units=net.get_shape().as_list()[-1]//r,
                                      activation=tf.nn.relu,
                                      name='channel_squeeze_and_excite_1')
            excite = tf.layers.dense(squeeze,
                                     units=net.get_shape().as_list()[-1],
                                     activation=tf.sigmoid,
                                     name='channel_squeeze_and_excite_2')
            excite = tf.expand_dims(
                tf.expand_dims(excite,axis=1),
                axis=1)
        return net * excite

    def s_squeeze_and_excite(net):
        with tf.variable_scope('SpatialSqueezeAndExcite'):
            squeeze = slim.conv2d(net,
                                  num_outputs=1,
                                  kernel_size=[1,1],
                                  stride=1,
                                  activation_fn=tf.nn.sigmoid,
                                  scope = 'spatial_squeeze_and_excite',
                                  padding = 'SAME')
        return net * squeeze

    def sc_squeeze_and_excite(net,r=1):
        with tf.variable_scope('ConcurrentSCSqueezeAndExcite'):
            c = c_squeeze_and_excite(net,r=r)
            s = s_squeeze_and_excite(net)
        return c + s

    endpoints = {}

    if beta > 0:
        weights_regularizer = tf.contrib.layers.l2_regularizer(beta)
    else:
        weights_regularizer = None

    classifications = []

    with slim.arg_scope(
        [slim.conv2d],
        activation_fn=tf.nn.elu,
        padding=padding,
        weights_regularizer=weights_regularizer,
        #normalizer_fn=pixel_normalization
        normalizer_fn=slim.batch_norm,
        normalizer_params={"is_training":is_training,
                           "decay":0.9,
                           "zero_debias_moving_mean":True}
        ):

        with tf.variable_scope('U-net', None, [inputs]):
            with tf.variable_scope('Red_Operations',None,[inputs]):
                with tf.variable_scope('Red_Block_1',None):
                    net,endpoints = red_block_wrapper(
                        inputs,
                        depth=int(64 * depth_mult),
                        factorization=factorization,
                        residuals=residuals,
                        endpoints=endpoints,
                        name='Red_Block_1')
                    if squeeze_and_excite == True:
                        net = sc_squeeze_and_excite(net,r=2)

                with tf.variable_scope('Red_Block_2',None,[net]):
                    net,endpoints = red_block_wrapper(
                        net,
                        depth=int(128 * depth_mult),
                        factorization=factorization,
                        residuals=residuals,
                        endpoints=endpoints,
                        name='Red_Block_2')
                    if squeeze_and_excite == True:
                        net = sc_squeeze_and_excite(net,r=2)

                with tf.variable_scope('Red_Block_3',None,[net]):
                    net,endpoints = red_block_wrapper(
                        net,
                        depth=int(256 * depth_mult),
                        factorization=factorization,
                        residuals=residuals,
                        endpoints=endpoints,
                        name='Red_Block_3')
                    if squeeze_and_excite == True:
                        net = sc_squeeze_and_excite(net,r=2)

                with tf.variable_scope('Red_Block_4',None,[net]):
                    net,endpoints = red_block_wrapper(
                        net,
                        depth=int(512 * depth_mult),
                        factorization=factorization,
                        residuals=residuals,
                        endpoints=endpoints,
                        name='Red_Block_4')
                    if squeeze_and_excite == True:
                        net = sc_squeeze_and_excite(net,r=2)

                with tf.variable_scope('Red_Block_5',None,[net]):
                    d_current = int(1024 * depth_mult)
                    net = block(net,d_current,3,1,factorization)
                    endpoints['Red_Block_5'] = net
                    if final_endpoint == 'Red_Block_5':
                        return net,endpoints

                    if aux_node == True:
                        flat_bottleneck = tf.reduce_max(net,axis=[1,2])
                        classification = fc_layer_classification(
                            flat_bottleneck,
                            depth=256,
                            n_classes=1,
                            activation_fn=tf.nn.relu,
                            weights_regularizer=weights_regularizer,
                            normalizer_fn=slim.batch_norm,
                            normalizer_params={"is_training":is_training}
                        )
                if aux_node == True:
                    endpoints['Classification'] = classification
                    classifications.append(classification)
                    if final_endpoint == 'Classification':
                        return net,endpoints

            with tf.variable_scope('Rec_Operations',None,[net]):
                with tf.variable_scope('Rec_Block_1',None,[net]):
                    net,endpoints = rec_block_wrapper(
                        net,
                        depth=int(512 * depth_mult),
                        factorization=factorization,
                        residuals=residuals,
                        endpoints=endpoints,
                        padding=padding,
                        name='Rec_Block_1',
                        red_block_name='Red_Block_4')
                    if squeeze_and_excite == True:
                        net = sc_squeeze_and_excite(net,r=2)

                with tf.variable_scope('Rec_Block_2',None,[net]):
                    net,endpoints = rec_block_wrapper(
                        net,
                        depth=int(256 * depth_mult),
                        factorization=factorization,
                        residuals=residuals,
                        endpoints=endpoints,
                        padding=padding,
                        name='Rec_Block_2',
                        red_block_name='Red_Block_3')
                    if squeeze_and_excite == True:
                        net = sc_squeeze_and_excite(net,r=2)

                with tf.variable_scope('Rec_Block_3',None,[net]):
                    net,endpoints = rec_block_wrapper(
                        net,
                        depth=int(128 * depth_mult),
                        factorization=factorization,
                        residuals=residuals,
                        endpoints=endpoints,
                        padding=padding,
                        name='Rec_Block_3',
                        red_block_name='Red_Block_2')
                    if squeeze_and_excite == True:
                        net = sc_squeeze_and_excite(net,r=2)

                with tf.variable_scope('Rec_Block_4',None,[net]):
                    net,endpoints = rec_block_wrapper(
                        net,
                        depth=int(64 * depth_mult),
                        factorization=factorization,
                        residuals=residuals,
                        endpoints=endpoints,
                        padding=padding,
                        name='Rec_Block_4',
                        red_block_name='Red_Block_1')
                    if squeeze_and_excite == True:
                        net = sc_squeeze_and_excite(net,r=2)

                with tf.variable_scope('Final',None,[net]):
                    net = slim.conv2d(net, n_classes, [1, 1],
                                      normalizer_fn=None,
                                      activation_fn=None,
                                      scope='conv2d_0_sigmoid')

                    endpoints['Final'] = net
                    if final_endpoint == 'Final':
                        return net,endpoints
    for endpoint in endpoints:
        print(endpoint,endpoints[endpoint])
    return net,endpoints,classifications

def iglovikov_loss(truth,network):
    """
    Loss function suggested by Iglovikov to include the Intersection of Union
    metric, a non-differentiable segmentation metric, with softmax cross
    entropy.

    Arguments [default]:
    * truth - tensor with ground truth
    * network - tensor with the output from the U-Net
    """

    network_softmax = tf.nn.softmax(network,axis = -1)

    j = tf.divide(
        tf.multiply(truth,network_softmax),
        tf.subtract(truth + network_softmax,tf.multiply(truth,network_softmax))
        )

    h = tf.losses.softmax_cross_entropy(truth,network)

    loss = tf.subtract(h,tf.reduce_mean(j))

    return loss

def active_contour_loss(truth,network):
    """
    Loss function implemented in "Learning Active Contour Models for Medical
    Image Segmentation".

    Arguments [default]:
    * truth - tensor with ground truth
    * network - tensor with the output from the U-Net
    """

    # length term

    x = network[:,1:,:,:] - network[:,:-1,:,:]
    y = network[:,:,1:,:] - network[:,:,:-1,:]

    delta_x = x[:,1:,:-2,:]**2
    delta_y = y[:,:-2,1:,:]**2
    delta_u = tf.math.abs(delta_x + delta_y)

    length = tf.reduce_mean(tf.math.sqrt(delta_u+1e-8),axis=[1,2,3])

    # region term

    region_in = tf.math.abs(
        tf.reduce_mean(network*(1.-truth)**2,axis=[1,2,3]))
    region_out = tf.math.abs(
        tf.reduce_mean((1.-network)*(-truth)**2,axis=[1,2,3]))

    mu = 1.
    lambda_p = 1.

    return length + lambda_p * (region_in * mu + region_out)

def weighted_softmax_cross_entropy(prediction,truth,weights):
    """
    Function that calculates softmax cross-entropy and weighs it using a
    pixel-wise map.

    Arguments:
    * prediction - output from the network
    * truth - ground truth
    * weights - weight map
    """
    prediction = tf.nn.softmax(prediction,axis = -1)
    h = tf.add(
        tf.multiply(truth,safe_log(prediction)),
        tf.multiply(1 - truth,safe_log(1 - prediction))
        )

    weighted_cross_entropy = tf.reduce_mean(h * weights)

    return weighted_cross_entropy

"""
Image and image processing-related operations. Image generators.
"""

def masked_mean(arr,mask):
    masked_arr = arr * mask
    masked_mean = np.sum(masked_arr) / np.sum(mask)
    return masked_mean

def image_to_array(image_path):
    """Opens image in image_path and returns it as a numpy array.

    Arguments:
    * image_path - the path to an image
    """

    with Image.open(image_path) as o:
        arr = np.array(o)
        if len(arr.shape) == 3:
            if arr.shape[2] == 4:
                arr = arr[:,:,:-1]
        return arr

def normal_image_generator(image_list,*mask_list,
                           classification_list=None,
                           random=True,eternal=True):

    """
    A standard image generator, to be used for testing purposes only. Goes over
    a list of image paths and outputs an image and its corresponding segmented
    image.

    Arguments:
    * image_path_list - a list of image paths
    * *mask_list - a list of ground truth images or any other masks
    * classification_list - a list containing classifications
    * random - whether to randomize the list of images each time
    * eternal - whether to continue indefinitely
    """

    cont = True

    while cont == True:
        indexes = np.arange(0,len(image_list))
        if random == True:
            np.random.shuffle(indexes)

        for index in indexes:
            image = image_list[index]

            masks = [mask[index] for mask in mask_list]

            if classification_list != None:
                classification = classification_list[index]
            else:
                classification = None

            im_shape = image.shape

            if im_shape[0] == im_shape[1]:
                yield image, (*masks), classification

        if eternal == False:
            cont = False

def prediction_image_generator(image_path_list):
    for image_path in image_path_list:
        yield np.array(Image.open(image_path)),os.path.basename(image_path)

def generate_tiles(large_image,
                   input_height = 256,input_width = 256,
                   padding = 'VALID'):
    """Uses a large image to generate smaller tiles for prediction.

    Arguments [default]:
    * large_image - a numpy array containing a large image
    * input_height - input height [256]
    * input_width - input width [256]
    * padding - whether VALID or SAME padding should be used ['VALID']
    """

    if padding == 'VALID':
        extra = 92

    stride_height = input_height - extra * 2
    stride_width = input_width - extra * 2

    stride_height,stride_width = input_height,input_width

    height,width = large_image.shape[0:2]
    h_tiles,w_tiles = floor(height/stride_height),floor(width/stride_width)
    for i in range(h_tiles - 1):
        h = i * stride_height
        for j in range(w_tiles - 1):
            w = j * stride_width
            tile = large_image[h:h + input_height,w:w + input_width,:]
            yield tile,h + extra,w + extra
        w = (j + 1) * stride_width
        if w + input_width > width:
            w = width - input_width
        tile = large_image[h:h + input_height,w:w + input_width,:]
        yield tile,h + extra,w + extra
    h = (i + 1) * stride_height
    if h + input_height > height:
        h = stride_height - input_height
    for j in range(w_tiles - 1):
        w = j * stride_width
        tile = large_image[h:h + input_height,w:w + input_width,:]
        yield tile,h + extra,w + extra
    w = (j + 1) * stride_width
    tile = large_image[h:h + input_height,w:w + input_width,:]
    yield tile,h + extra,w + extra

def remap_tiles(mask,division_mask,h_1,w_1,tile):
    """
    Function used to remap the tiles to the original image size. Currently, it
    is working with one class prediction/"channel". The division mask will be
    used in the end to divide the final_prediction and use a sort of a
    "consensus" approach for overlapping regions.

    * mask - mask with the same size of a large image
    * division_mask - mask that sums the number of positive pixels in
    overlapping regions
    * h_1 - height for the tile insertion
    * w_1 - width for the tile insertion
    * tile - numpy array with the output from U-Net
    """

    x,y = tile.shape[0:2]
    mask[h_1:h_1 + x,w_1:w_1 + y,:] += tile
    division_mask[h_1:h_1 + x,w_1:w_1 + y,:] += np.ones(tile.shape)
    return mask,division_mask

def generate_images(image_path_list,truth_path,
                    input_height = 256,input_width = 256,
                    padding = 'VALID',n_classes = 2,
                    truth_only = False,
                    weight_maps = True,
                    mode = 'train'):
    """
    Multi-purpose image generator.

    Arguments [default]:
    * image_path_list - a list of image paths
    * truth_path - a list of ground truth image paths
    * net_x - output height for the network [None]
    * net_y - output width for the network [None]
    * input_height - input height for the network [256]
    * input_width - input width for the network [256]
    * padding - whether VALID or SAME padding should be used ['VALID']
    * n_classes - no. of classes [2]
    * truth_only - whether only positive images should be used [False]
    * mode - algorithm mode [train].
    """

    def check_size(image_path):
        size = Image.open(image_path).size
        return np.all([size[0] == input_height,size[1] == input_width])

    image_path_list = [x for x in image_path_list if check_size(x) == True]

    if mode == 'train':
        image_list = []
        truth_list = []
        weight_map_list = []

        for i,image_path in enumerate(image_path_list):
            class_array = np.array([-1 for i in range(n_classes)])
            #if i % 5 == 0: print(i)
            image_name = image_path.split(os.sep)[-1]
            truth_image_path = truth_path + os.sep + image_name
            truth_img = image_to_array(truth_image_path)

            #assigns a class for edges
            if len(truth_img.shape) == 3:
                truth_img = np.mean(truth_img,axis=2)
            if n_classes == 3:
                edges = cv2.Canny(truth_img.astype('uint8'),100,200)
                edges = cv2.dilate(edges,
                                   kernel=np.ones([3,3]),
                                   iterations=1)
                #assigns one channel per class
                classes = np.unique(truth_img)
                for index,label in enumerate(classes):
                    class_array[index] = label
                mask = np.zeros((truth_img.shape[0],
                                 truth_img.shape[1],
                                 n_classes))
                for j in range(len(classes)):
                    mask[:,:,j][truth_img == classes[j]] = 1
                mask[:,:,len(classes) - 1][[edges > 0]] = 0
                mask[:,:,len(classes)][edges > 0] = 1

            if n_classes == 2:
                #assigns one channel per class
                classes = np.unique(truth_img)
                for index,label in enumerate(classes):
                    class_array[index] = label
                classes = class_array
                mask = np.zeros((truth_img.shape[0],
                                 truth_img.shape[1],
                                 n_classes))
                for j in range(n_classes):
                    mask[:,:,j][truth_img == classes[j]] = 1
            truth_img = mask
            weight_map = get_weight_map(truth_img)
            dist_weight_map = get_near_weight_map(truth_img,w0=2,sigma=20)
            weight_map = weight_map + dist_weight_map
            image_list.append(image_to_array(image_path))
            truth_list.append(truth_img)
            weight_map_list.append(weight_map[:,:,np.newaxis])

        generator = normal_image_generator(
            image_list,truth_list,weight_map_list,
            classification_list=None,
            random=True,eternal=True
        )

    elif mode == 'test':
        image_list = []
        truth_list = []

        for i,image_path in enumerate(image_path_list):
            #if i % 5 == 0: print(i)
            image_name = image_path.split(os.sep)[-1]
            truth_image_path = truth_path + os.sep + image_name
            truth_img = image_to_array(truth_image_path)

            #assigns a class for edges
            if len(truth_img.shape) == 3:
                truth_img = np.mean(truth_img,axis=2)
            if n_classes == 3:
                edges = cv2.Canny(truth_img.astype('uint8'),100,200)
                edges = cv2.dilate(edges,
                                   kernel=np.ones([3,3]),
                                   iterations=1)
                #assigns one channel per class
                classes = np.unique(truth_img)
                mask = np.zeros((truth_img.shape[0],
                                 truth_img.shape[1],
                                 n_classes))
                for j in range(len(classes)):
                    mask[:,:,j][truth_img == classes[j]] = 1
                mask[:,:,len(classes) - 1][[edges > 0]] = 0
                mask[:,:,len(classes)][edges > 0] = 1

            if n_classes == 2:
                #assigns one channel per class
                classes = np.unique(truth_img)
                mask = np.zeros((truth_img.shape[0],
                                 truth_img.shape[1],
                                 n_classes))
                for j in range(len(classes)):
                    mask[:,:,j][truth_img == classes[j]] = 1
            truth_img = mask

            image_list.append(image_to_array(image_path))
            truth_list.append(truth_img)

        generator = normal_image_generator(
            image_list,truth_list,
            classification_list=None,
            random=False,eternal=False
        )

    elif mode == 'predict' or mode == 'tumble_predict':
        generator = prediction_image_generator(image_path_list)
        for image,image_path in generator:
            yield image,image_path

    elif mode == 'large_predict':
        batch_paths = []
        batch_coord = []
        batch_shape = []
        for large_image_path in image_path_list:
            large_image = np.array(
                Image.open(large_image_path)
                )
            generator = generate_tiles(
                large_image,
                input_height = input_height,
                input_width = input_width,
                padding = padding
                )

            for img,x,y in generator:
                img = (img - np.min(img)) / (np.max(img) - np.min(img))
                if len(batch) >= batch_size:
                    batch = []
                    batch_paths = []
                    batch_coord = []
                    batch_shape = []

                shape = img.shape
                yield image,large_image_path,(x,y),batch_shape

    if mode == 'train':
        while True:
            for element in generator:
                element = tuple(x for x in element[:-1])
                yield element
    else:
        for element in generator:
            element = tuple(x for x in element[:-1])
            yield element

def generate_images_affinity(
    image_path_list,truth_path,
    input_height = 256,input_width = 256,
    padding = 'VALID',n_classes = 2,
    truth_only = False,
    weight_maps = True,
    mode = 'train'):

    """
    Multi-purpose image generator for affinity calculation.

    Arguments [default]:
    * image_path_list - a list of image paths
    * truth_path - a list of ground truth image paths
    * net_x - output height for the network [None]
    * net_y - output width for the network [None]
    * input_height - input height for the network [256]
    * input_width - input width for the network [256]
    * padding - whether VALID or SAME padding should be used ['VALID']
    * n_classes - no. of classes [2]
    * truth_only - whether only positive images should be used [False]
    * mode - algorithm mode [train].
    """

    def image_to_nearby_pixels(image,directions=[1]):
        image_shape = image.shape
        output = np.zeros([image_shape[0],image_shape[1],len(directions) * 4])

        for i,direction in enumerate(directions):
            padded_mask = np.pad(image,pad_width=direction,mode='reflect')
            index_list = [
                [0,image_shape[0]],
                [2*direction,image_shape[0] + 2*direction],
                [0,image_shape[1]],
                [2*direction,image_shape[1] + 2*direction]
            ]
            original_image_index_x = [direction,image_shape[0]+direction]
            original_image_index_y = [direction,image_shape[1]+direction]
            output[:,:,0 * i] = padded_mask[
                index_list[0][0]:index_list[0][1],
                original_image_index_y[0]:original_image_index_y[1]]
            output[:,:,1 * i] = padded_mask[
                index_list[1][0]:index_list[1][1],
                original_image_index_y[0]:original_image_index_y[1]]
            output[:,:,2 * i] = padded_mask[
                original_image_index_x[0]:original_image_index_x[1],
                index_list[2][0]:index_list[2][1]]
            output[:,:,3 * i] = padded_mask[
                original_image_index_x[0]:original_image_index_x[1],
                index_list[3][0]:index_list[3][1]]

        return output

    def nearby_pixels_to_affinity(image,nearby_pixels_image):
        if len(image.shape) == 2:
            image=np.expand_dims(image,axis=2)

        affinity_mask = np.float32(nearby_pixels_image == image_shape)

        return affinity_mask

    def check_size(image_path):
        size = Image.open(image_path).size
        return np.all([size[0] == input_height,size[1] == input_width])

    image_path_list = [x for x in image_path_list if check_size(x) == True]

    if mode == 'train':
        image_list = []
        truth_list = []

        for i,image_path in enumerate(image_path_list):
            class_array = np.array([-1 for i in range(n_classes)])
            #if i % 5 == 0: print(i)
            image_name = image_path.split(os.sep)[-1]
            truth_image_path = truth_path + os.sep + image_name
            truth_img = image_to_array(truth_image_path)

            #assigns a class for edges
            if len(truth_img.shape) == 3:
                truth_img = np.mean(truth_img,axis=2)

            nearby_pixels_image = image_to_nearby_pixels(truth_img,
                                                         [1,3,5,7,9,11])
            truth_img = nearby_pixels_to_affinity(truth_img,
                                                  nearby_pixels_image)
            truth_img = np.concatenate(
                [truth_img, 1 - truth_img],
                axis=-1
            )

            image_list.append(image_to_array(image_path))
            truth_list.append(truth_img)

        generator = normal_image_generator(
            image_list,truth_list,
            classification_list=None,
            random=True,eternal=True
        )

    elif mode == 'test':
        image_list = []
        truth_list = []

        for i,image_path in enumerate(image_path_list):
            class_array = np.array([-1 for i in range(n_classes)])
            #if i % 5 == 0: print(i)
            image_name = image_path.split(os.sep)[-1]
            truth_image_path = truth_path + os.sep + image_name
            truth_img = image_to_array(truth_image_path)

            #assigns a class for edges
            if len(truth_img.shape) == 3:
                truth_img = np.mean(truth_img,axis=2)

            nearby_pixels_image = image_to_nearby_pixels(truth_img,
                                                         [1,3,5,7,9,11])
            truth_img = nearby_pixels_to_affinity(truth_img,
                                                  nearby_pixels_image)
            truth_img = np.concatenate(
                [truth_img, 1 - truth_img],
                axis=-1
            )

            image_list.append(image_to_array(image_path))
            truth_list.append(truth_img)

        generator = normal_image_generator(
            image_list,truth_list,
            classification_list=None,
            random=False,eternal=False
        )

    elif mode == 'predict' or mode == 'tumble_predict':
        generator = prediction_image_generator(image_path_list)
        for image,image_path in generator:
            yield image,image_path

    elif mode == 'large_predict':
        batch_paths = []
        batch_coord = []
        batch_shape = []
        for large_image_path in image_path_list:
            large_image = np.array(
                Image.open(large_image_path)
                )
            generator = generate_tiles(
                large_image,
                input_height = input_height,
                input_width = input_width,
                padding = padding
                )

            for img,x,y in generator:
                img = (img - np.min(img)) / (np.max(img) - np.min(img))
                if len(batch) >= batch_size:
                    batch = []
                    batch_paths = []
                    batch_coord = []
                    batch_shape = []

                shape = img.shape
                yield image,large_image_path,(x,y),batch_shape

    if mode == 'train':
        while True:
            for element in generator:
                element = tuple(x for x in element[:-1])
                yield element
    else:
        for element in generator:
            element = tuple(x for x in element[:-1])
            yield element

def generate_images_propagation(
    image_path_list,truth_path,propagation_path,
    batch_size=4,
    net_x = None,net_y = None,
    input_height = 256,input_width = 256,resize = False,
    padding = 'VALID',truth_only = False,weight_maps = True):
    """
    Multi-purpose image generator.

    Arguments [default]:
    * image_path_list - a list of image paths
    * truth_path - a list of ground truth image paths
    * batch_size - the size of the batch
    * input_height - input height for the network [256]
    * input_width - input width for the network [256]
    * padding - whether VALID or SAME padding should be used ['VALID']
    * n_classes - no. of classes [2]
    * truth_only - whether only positive images should be used [False]
    * mode - algorithm mode [train].
    """

    if mode == 'train' or mode == 'test':
        image_list = []
        truth_list = []
        weight_map_list = []

        for i,image_path in enumerate(image_path_list):
            class_array = np.array([-1 for i in range(n_classes)])
            #if i % 5 == 0: print(i)
            image_name = image_path.split(os.sep)[-1]
            truth_image_path = os.path.join(truth_path,image_name)
            propagation_image_path = os.path.join(propagation_path,image_name)
            truth_img = image_to_array(truth_image_path)
            propa_img = image_to_array(propagation_image_path)

            if len(truth_img.shape) == 3:
                truth_img = np.mean(truth_img,axis=2)

            if len(propa_img.shape) == 3:
                propa_img = np.mean(propa_img,axis=2)

            truth_img = np.where(truth_img > 0.5,1.,0.)
            truth_img = np.stack([1. - truth_img,truth_img],axis=2)

            propa_img = np.where(propa_img > 0,1.,0.)

            weight_map = get_weight_map(truth_img)
            dist_weight_map = get_near_weight_map(truth_img,w0=5,sigma=20)
            weight_map = (weight_map + dist_weight_map) * propa_img
            image_list.append(image_to_array(image_path))
            truth_list.append(truth_img)
            weight_map_list.append(weight_map[:,:,np.newaxis])

        generator = normal_image_generator(image_list,truth_list,
                                           weight_map_list)

        while True:
            for element in generator:
                img,truth_img,weight_map,_ = element
                yield img,truth_img,weight_map

    if mode == 'predict':
        a = True

        image_list = []
        weight_map_list = []

        for i,image_path in enumerate(image_path_list):
            class_array = np.array([-1 for i in range(n_classes)])
            #if i % 5 == 0: print(i)
            image_name = image_path.split(os.sep)[-1]
            propagation_image_path = os.path.join(propagation_path,image_name)
            propa_img = image_to_array(propagation_image_path)

            if len(propa_img.shape) == 3:
                propa_img = np.mean(propa_img,axis=2)

            propa_img = np.where(propa_img > 0,1.,0.)

            image_list.append(image_to_array(image_path))
            weight_map_list.append(propa_img)

        generator = normal_image_generator(image_list,weight_map_list)

        batch = []
        weight_batch = []
        while a == True:
            for element in generator:
                img,weight_map,_ = element
                yield img,weight_map

def generate_images_h5py_dataset(h5py_path,
                                input_height=256,
                                input_width=256,
                                key_list = None):
    segmentation_dataset = SegmentationDataset(
        hdf5_file=h5py_path,
        dimensions=(0,0,input_height,input_width),
        mode='segmentation',
        rel_keys=['image','mask','weight_map'],
        rotate_record=True,
        transform=None
    )
    size = len(segmentation_dataset)
    key_list = [x for x in key_list if x in segmentation_dataset.hf_keys]
    while True:
        if key_list is None:
            random_idx = np.random.randint(0,size)
            rr = segmentation_dataset[random_idx]
            mask = np.concatenate([1-rr['mask'],rr['mask']],axis=2)
            yield rr['image'],mask,rr['weight_map']
        else:
            random_key = np.random.choice(key_list)
            rr = segmentation_dataset[random_key]
            mask = np.concatenate([1-rr['mask'],rr['mask']],axis=2)
            yield rr['image'],mask,rr['weight_map']

def classification_generator(image_path_list,classification_list,
                             chances = [0,0,0],
                             input_height = 256,input_width = 256,
                             padding = 'VALID',n_classes = 2,
                             mode = 'train',batch_size = 4):
    """
    Realtime image generator for encoder pre-training.
    """

    batch = []
    classification_batch = []
    a = True

    image_list = [image_to_array(im_path) for im_path in image_path_list]
    tmp_images = []
    tmp_classifications = []
    for image,classification in zip(image_list,classification_list):
        if image.shape[0] == input_height and image.shape[1] == input_width:
            tmp_images.append(image)
            tmp_classifications.append(classification)

            #if len(tmp_images) % 5 == 0: print(len(tmp_images))
    image_list = tmp_images
    classification_list = tmp_classifications

    if mode == 'train':
        generator = realtime_image_augmentation(
            image_list=image_list,
            truth_list=None,
            weight_map_list=None,
            classification_list=classification_list,
            noise_chance=chances[0],
            blur_chance=chances[1],
            flip_chance=chances[2])

    else:
        generator = normal_image_generator(image_list,classification_list)

    while True:
        for element in generator:
            if mode == 'train':
                img,_,_,classification = element
            else:
                img,classification = element
            yield batch,classification

def tf_dataset_from_generator(generator,generator_params,
                              output_types,output_shapes,
                              is_training,buffer_size,batch_size):
    """
    Returns the next element
    """
    dataset = tf.data.Dataset.from_generator(
        generator=lambda: generator(**generator_params),
        output_types=output_types,
        output_shapes=output_shapes
    )

    if is_training == True:
        dataset = dataset.repeat()
        dataset = dataset.shuffle(buffer_size=buffer_size)

    dataset = dataset.batch(batch_size)
    iterator = dataset.make_one_shot_iterator()
    return iterator.get_next()

"""
Weight map calculation from segmentation maps operations.
"""

def get_poormans_weight_map(truth_image,w0=0.5,convolution_size=9):
    """
    Instead of defining a weight map using distances to nearest positive pixels
    (computationally heavy), use the product of a large convolution to weigh
    the pixels (pixels closer to a lot of white pixels will have higher weight
    values, while those further away will be lower).

    EDIT: turns out this is a not very successful idea unless it is done with
    very large convolutions and, ultimately, it does not lead to better results
    than the distance weight map. Keeping this here for some possible
    usefulness it might have.

    Arguments [default]:
    * truth_image - numpy array with ground truth image
    * w0 - lowest weight value possible
    * convolution_size - size of the convolutional filter to be applied
    """
    truth_image = truth_image.copy()
    truth_image = truth_image[:,:,1]
    if truth_image.max() > 0:
        w = cv2.GaussianBlur(
            truth_image,
            (convolution_size,convolution_size),
            0)
        w = (w - w.min())/(w.max() - w.min())
        w = w * (1 - w0) + w0
        w[truth_image > 0] = 1
    else:
        w = np.ones(truth_image.shape) * w0
    return w

def get_near_weight_map(truth_image,w0=10,sigma=5):
    """
    Calculates the weight map for pixels between nearby cells.

    Arguments [default]:
    * truth_image - numpy array with ground truth image
    * w0 - w0 value for the weight map [10]
    * sigma - sigma value for the weight map [5]
    """
    truth_image = truth_image.copy()
    truth_image = truth_image[:,:,1]
    sh = truth_image.shape
    truth_image[truth_image > 0] = 255
    kernel = np.ones((3,3))
    truth_image = cv2.morphologyEx(truth_image.astype('uint8'),
                                   cv2.MORPH_OPEN, kernel)
    if truth_image.max() > 0:
        n_cc,cc = cv2.connectedComponents(truth_image)
        zero_pixels = np.where(truth_image == 0)
        zero_pixels = np.array(zero_pixels).T
        edges = cv2.Canny(truth_image,100,200)
        edges = cv2.dilate(edges,kernel)
        edges[edges > 0] = 1
        accumulator = np.zeros((zero_pixels.shape[0],n_cc))
        for i in range(1,n_cc):
            arr = np.zeros(truth_image.shape)
            arr[cc == i] = 1
            arr = arr * edges
            one_pixels = np.where(arr > 0)
            one_pixels = np.array(one_pixels).T
            d = distance.cdist(zero_pixels,one_pixels)
            d = np.amin(d,axis=1)
            accumulator[:,i] = d
        accumulator = np.sort(accumulator,axis=1,kind='mergesort')
        if accumulator.shape[1] > 2:
            weight_values = np.exp(
                - np.square(
                    accumulator[:,1] + accumulator[:,2]
                ) / (2 * (sigma ** 2))
            )
            output = np.zeros(sh)
            output[zero_pixels[:,0],zero_pixels[:,1]] = w0 * weight_values
        else:
            output = np.zeros(sh)
    else:
        output = np.zeros(sh)
    return output

def get_weight_map(truth_image,w0=0.5,sigma=10):
    """
    Has to be one channel only.

    Arguments [default]:
    * truth_image - numpy array with ground truth image
    * w0 - w0 value for the weight map [10]
    * sigma - sigma value for the weight map [5]
    """
    stride = 128
    size = 256
    sh = truth_image.shape
    max_x = sh[0] // stride
    max_y = sh[1] // stride
    if max_x % stride != 0: max_x += 1
    if max_y % stride != 0: max_y += 1
    kernel = np.ones((3,3))
    truth_image = truth_image.copy()
    truth_image = truth_image[:,:,1]
    truth_image[truth_image > 0] = 255
    truth_image = cv2.morphologyEx(truth_image.astype('uint8'),
                                   cv2.MORPH_OPEN, kernel)
    truth_image_ = truth_image.copy()
    edges = cv2.Canny(truth_image,100,200)

    if 255 in truth_image:
        final_weight_mask = np.zeros((sh[0],sh[1]),dtype = 'float32')
        for i in range(max_x):
            for j in range(max_y):
                sub_image = truth_image[i*stride:(i+1)*size,
                                        j*stride:(j+1)*size]
                sub_edges = edges[i*stride:(i+1)*size,
                                  j*stride:(j+1)*size]
                tmp_weight_mask = final_weight_mask[i*stride:(i+1)*size,
                                                    j*stride:(j+1)*size]
                ssh = sub_image.shape
                pixel_coords = np.where(sub_image == 0)
                pixel_coords_t = np.transpose(pixel_coords)
                weight_mask = np.zeros((ssh[0],ssh[1]),dtype = 'float32')

                cp_coords = np.transpose(np.where(sub_edges > 0))
                distances = distance.cdist(pixel_coords_t,cp_coords)
                distances[distances == 0] = inf
                if distances.any():
                    mins = np.array(np.min(distances,axis = 1))
                    weight_map = ((mins * 2) ** 2)/(2 * sigma ** 2)
                    weight_map = w0 + ((1 - w0) * np.exp(-weight_map))

                    weight_mask[pixel_coords[0],
                                pixel_coords[1]] = weight_map

                else:
                    weight_mask = np.ones((ssh[0],ssh[1]),
                                          dtype = 'float32') * w0
                tmp_weight_mask = np.where(
                    weight_mask > tmp_weight_mask,
                    weight_mask,
                    tmp_weight_mask
                )
                final_weight_mask[i*stride:(i+1)*size,
                                  j*stride:(j+1)*size] = tmp_weight_mask

    else:
        final_weight_mask = np.ones(truth_image.shape) * w0

    final_weight_mask[truth_image > 0] = 1
    return final_weight_mask

"""
Miscellaneous operations.
"""

def log_write_print(log_file,to_write):
    """
    Convenience funtion to write something in a log file and also print it.

    Arguments:
    * log_file - path to the file
    * to_write - what should be written/printed
    """
    with open(log_file,'a') as o:
        o.write(to_write + '\n')
    print(to_write)
