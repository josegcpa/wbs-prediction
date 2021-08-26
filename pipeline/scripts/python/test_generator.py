"""
Implementation of a DNN using the attention modules described in [1].
Generalized for a folder containing images and a csv containing its
classes.

[1] http://openaccess.thecvf.com/content_ECCV_2018/papers/Pau_Rodriguez_Lopez_Attend_and_Rectify_ECCV_2018_paper.pdf
"""

import sys
import os
import time
from glob import glob
from math import ceil
import numpy as np
import cv2
from PIL import Image
import openslide

from image_generator import ImageGeneratorWithQueue

height = 1024
width = 1024
n_channels = 3

csv_path = sys.argv[1]

image_path_list = sorted(glob(csv_path + '/*'))
igwq = ImageGeneratorWithQueue(csv_path,None,maxsize=16,extra_padding=32,height=height,width=width)
igwq.start()
for image in igwq.generate():
    print(image[0].shape,image[1])

