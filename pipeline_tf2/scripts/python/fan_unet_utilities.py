import os
import numpy as np
import h5py
from PIL import Image
from glob import glob

import tensorflow as tf
from tensorflow import keras

from fan_utilities import *
from unet_utilities import *

class FANUNet(keras.Model):
    def __init__(self,fan_model,unet_model):
        super(FANUNet,self).__init__()
        self.fan_model = fan_model
        self.unet_model = unet_model

    def call(self,x):
        prediction = self.unet_model(self.fan_model(x))

        return prediction

    def train_step(self, data):
        x, y, w = data
        with tf.GradientTape() as tape:
            y_pred = self(x, training=True)
            loss = self.loss_fn(
                y,y_pred,w,regularization_losses=self.losses)

        self.compiled_loss(
            y,y_pred,sample_weight=None,regularization_losses=self.losses) 

        trainable_vars = self.trainable_variables
        gradients = tape.gradient(loss, trainable_vars)
        self.optimizer.apply_gradients(zip(gradients, trainable_vars))
        self.compiled_metrics.update_state(y, y_pred)
        return {m.name: m.result() for m in self.metrics}