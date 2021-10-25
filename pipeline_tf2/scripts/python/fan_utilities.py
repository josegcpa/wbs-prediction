import os
import numpy as np
import h5py
from PIL import Image
from glob import glob

import tensorflow as tf
from tensorflow import keras

output_layer_ids = ['block_1_expand','block_3_expand',
                    'block_5_expand','block_7_expand']

class InvertibleConv2D(keras.layers.Layer):
    def __init__(self,n_input_features,n_output_features,batch_size):
        super(InvertibleConv2D,self).__init__()
        self.n_input_features = n_input_features
        self.n_output_features = n_output_features
        self.batch_size = batch_size
        self.w_shape = (1,1,self.n_input_features,self.n_output_features)
        self.init = keras.initializers.GlorotUniform()
        self.W = tf.Variable(self.init(self.w_shape))

    def encode(self,x):
        sh = x.shape.as_list()
        return tf.nn.conv2d(x,self.W,[1,1],'SAME')
    
    def decode(self,x):
        sh = x.shape.as_list()
        return tf.nn.conv2d_transpose(
            x,self.W,[self.batch_size,sh[1],sh[2],self.n_input_features],
            strides=[1,1])
    
    def update_batch_size(self,batch_size):
        self.batch_size = batch_size

    def call(self,x):
        if x.shape[-1] == self.n_input_features:
            return self.encode(x)
        else:
            return self.decode(x)

class FANLayer(keras.layers.Layer):
    def __init__(self,n_features,ratio,epsilon=1e-8,
                 upscaling='standard'):
        super(FANLayer,self).__init__()
        self.n_features = n_features
        self.ratio = ratio
        self.epsilon = epsilon
        self.upscaling = upscaling
        self.setup_layers()
    
    def setup_layers(self):
        self.sigmoid = keras.Sequential()
        self.sigmoid.add(keras.layers.Conv2D(self.n_features,1))
        self.sigmoid.add(keras.layers.Activation('sigmoid'))
        
        self.relu = keras.Sequential()
        self.relu.add(keras.layers.Conv2D(self.n_features,1))
        self.relu.add(keras.layers.Activation('relu'))

        if self.upscaling == 'standard':
            self.sigmoid.add(keras.layers.UpSampling2D(self.ratio))
            self.relu.add(keras.layers.UpSampling2D(self.ratio))
        elif self.upscaling == 'transpose':
            self.sigmoid.add(
                keras.layers.Conv2DTranspose(
                    self.n_features,self.ratio,padding='same',
                    strides=self.ratio))
            self.relu.add(
                keras.layers.Conv2DTranspose(
                    self.n_features,self.ratio,padding='same',
                    strides=self.ratio))
    
    def call(self,x,z):
        mean_x = tf.reduce_mean(x,axis=[1,2],keepdims=True)
        std_x =  tf.math.reduce_std(x,axis=[1,2],keepdims=True)
        x = (x - mean_x) / (std_x**2 + self.epsilon)
        sigmoid_output = self.sigmoid(z)
        relu_output = self.relu(z)
        return x * sigmoid_output + relu_output

class FAN(keras.Model):
    def __init__(self,sub_model,n_features,
                 input_height,input_width,epsilon=1e-8,
                 batch_size=1,upscaling='standard'):
        super(FAN,self).__init__()
        self.n_features = n_features
        self.epsilon = epsilon
        self.sub_model = sub_model
        self.input_height = input_height
        self.input_width = input_width
        self.batch_size = batch_size
        self.upscaling = upscaling
        self.n_z = len(self.sub_model.output)
        self.infer_ratios()
        self.setup_layers()
    
    def setup_layers(self):
        self.encoder_decoder = InvertibleConv2D(
            3,self.n_features,self.batch_size)
        self.linear_layers = {
            'fan_'+str(i): FANLayer(
                self.n_features,ratio,upscaling=self.upscaling) 
            for i,ratio in zip(range(self.n_z),self.ratios[-1::-1])}
        
    def infer_ratios(self):
        tmp = tf.ones([1,self.input_height,self.input_width,3])
        outs = self.sub_model(tmp)
        self.ratios = [self.input_height//x.shape[1] for x in outs]

    def call(self,x):
        feature_space = self.encoder_decoder.encode(x)
        zs = self.sub_model(x)
        tr = feature_space
        for lin,z in zip(self.linear_layers,zs[-1::-1]):
            tr = self.linear_layers[lin](tr,z)
        return self.encoder_decoder.decode(tr)

class ColourAugmentation(keras.layers.Layer):
    def __init__(self,
                 brightness_delta,
                 contrast_lower,contrast_upper,
                 hue_delta,
                 saturation_lower,saturation_upper,
                 probability=0.1):
        super(ColourAugmentation,self).__init__()
        self.probability = probability
        self.brightness_delta = brightness_delta
        self.contrast_lower = contrast_lower
        self.contrast_upper = contrast_upper
        self.hue_delta = hue_delta
        self.saturation_lower = saturation_lower
        self.saturation_upper = saturation_upper
        
    def brightness(self,x):
        return tf.image.random_brightness(
            x,self.brightness_delta)
    
    def contrast(self,x):
        return tf.image.random_contrast(
            x,self.contrast_lower,self.contrast_upper)
    
    def hue(self,x):
        return tf.image.random_hue(
            x,self.hue_delta)
    
    def saturation(self,x):
        return tf.image.random_saturation(
            x,self.saturation_lower,self.saturation_upper)
    
    def call(self,x):
        fn_list = [self.brightness,self.contrast,
                   self.hue,self.saturation]
        np.random.shuffle(fn_list)
        for fn in fn_list:
            if np.random.uniform() < self.probability:
                x = fn(x)
        x = tf.clip_by_value(x,0,1)
        return x
    
class Flipper(keras.layers.Layer):
    def __init__(self,probability=0.1):
        super(Flipper,self).__init__()
        self.probability = probability
            
    def call(self,x):
        if np.random.uniform() < self.probability:
            x = tf.image.flip_left_right(x)
        if np.random.uniform() < self.probability:
            x = tf.image.flip_up_down(x)
        return x

class ImageCallBack(keras.callbacks.Callback):
    def __init__(self,save_every_n,tf_dataset,log_dir):
        super(ImageCallBack, self).__init__()
        self.save_every_n = save_every_n
        self.tf_dataset = iter(tf_dataset)
        self.log_dir = log_dir
        self.writer = tf.summary.create_file_writer(self.log_dir)
        self.count = 0

    def on_train_batch_end(self, batch, logs=None):
        if self.count % self.save_every_n == 0:
            batch = next(self.tf_dataset)
            y_augmented,y_true = batch
            prediction = self.model.predict(y_augmented)
            with self.writer.as_default():
                tf.summary.image("0:InputImage",y_augmented,self.count)
                tf.summary.image("1:GroundTruth",y_true,self.count)
                tf.summary.image("2:Prediction",prediction,self.count)
                tf.summary.scalar("Loss",logs['loss'],self.count)
                tf.summary.scalar("MAE",logs['mean_absolute_error'],self.count)
        self.count += 1

class DataGenerator:
    def __init__(self,image_folder_path,shuffle=True,transform=None):
        self.image_folder_path = image_folder_path
        self.image_paths = glob('{}/*'.format(self.image_folder_path))
        self.shuffle = shuffle
        self.transform = transform
        self.n_images = len(self.image_paths)
    
    def generate(self,with_path=False):
        image_idx = [x for x in range(len(self.image_paths))]
        if self.shuffle == True:
            np.random.shuffle(image_idx)
        for idx in image_idx:
            P = self.image_paths[idx]
            x = np.array(Image.open(P))[:,:,:3]
            x = tf.convert_to_tensor(x) / 255
            if self.transform is not None:
                x = self.transform(x)
            if with_path == True:
                yield x,P
            else:
                yield x

class DataGeneratorHDF5:
    def __init__(self,hdf5_path,shuffle=True,transform=None):
        self.hdf5_path = hdf5_path
        self.h5 = h5py.File(self.hdf5_path,'r')
        self.shuffle = shuffle
        self.transform = transform
        self.all_keys = list(self.h5.keys())
        self.n_images = len(self.all_keys)
    
    def generate(self,with_path=False):
        image_idx = [x for x in range(self.n_images)]
        if self.shuffle == True:
            np.random.shuffle(image_idx)
        for idx in image_idx:
            P = self.all_keys[idx]
            x = self.h5[P]['image']
            x = tf.convert_to_tensor(x) / 255
            if self.transform is not None:
                x = self.transform(x)
            if with_path == True:
                yield x,P
            else:
                yield x

class LargeImage:
    def __init__(self,image,tile_size=[512,512],
                 output_channels=3,offset=0):
        """
        Class facilitating the prediction for large images by 
        performing all the necessary operations - tiling and 
        reconstructing the output.
        """
        self.image = image
        self.tile_size = tile_size
        self.output_channels = output_channels
        self.offset = offset
        self.h = self.tile_size[0]
        self.w = self.tile_size[1]
        self.sh = self.image.shape[:2]
        self.output = np.zeros([self.sh[0],self.sh[1],self.output_channels])
        self.denominator = np.zeros([self.sh[0],self.sh[1],1])

    def tile_image(self):
        for x in range(0,self.sh[0],self.h):
            x = x - self.offset
            if x + self.tile_size[0] > self.sh[0]:
                x = self.sh[0] - self.tile_size[0]
            for y in range(0,self.sh[1],self.w):
                y = y - self.offset
                if y + self.tile_size[1] > self.sh[1]:
                    y = self.sh[1] - self.tile_size[1]
                x_1,x_2 = x, x+self.h
                y_1,y_2 = y, y+self.w
                yield self.image[x_1:x_2,y_1:y_2,:],((x_1,x_2),(y_1,y_2))

    def update_output(self,image,coords):
        (x_1,x_2),(y_1,y_2) = coords
        self.output[x_1:x_2,y_1:y_2,:] += image
        self.denominator[x_1:x_2,y_1:y_2,:] += 1

    def return_output(self):
        return self.output/self.denominator
