import openslide
import os
import sys
from math import floor
from random import shuffle
import numpy as np
from PIL import Image

class OpenslideSplitter:
    """General class for splitting a large OpenSlide image into smaller images,
    doing some operations to remove images which are too blurry or contain no
    apparent cells.
    """
    def __init__(self,openslide_path,big_frame_size,random = False):
        self.openslide_path = openslide_path
        self.big_frame_size = big_frame_size
        self.random = random
        self.openslide_handle = openslide.OpenSlide(self.openslide_path)
        self.identifier = self.get_identifier()
        self.coord_list = self.get_coords()
        self.images = self.image_generator()

    def __repr__(self):
        return '<OpenslideSplitter for {0}>'.format(self.openslide_path)

    def get_identifier(self):
        identifier = self.openslide_path.split(os.sep)[-1]
        return ''.join(identifier.split('.')[:-1])

    def get_coords(self):

        x,y = self.openslide_handle.dimensions
        no_shift_x = floor(x/self.big_frame_size)
        no_shift_y = floor(y/self.big_frame_size)

        coord_list = []

        for i in range(0,no_shift_x):
            for j in range(0,no_shift_y):
                x1,y1 = (self.big_frame_size * i,
                         self.big_frame_size * j)

                coord_list.append(((x1,y1),(self.big_frame_size,
                                            self.big_frame_size),
                                   (i,j)))

            x1,y1 = (self.big_frame_size * i,self.
                     big_frame_size * no_shift_y)

            coord_list.append(((x1,y1),
                               (self.big_frame_size,
                                y - self.big_frame_size * no_shift_y),
                               (i,j + 1)))

        for j in range(0,no_shift_y):
            x1,y1 = (self.big_frame_size * no_shift_x,
                     self.big_frame_size * j)

            coord_list.append(((x1,y1),
                               (x-self.big_frame_size*no_shift_x,
                                self.big_frame_size),
                               (i + 1,j)))

        x1,y1 = (self.big_frame_size * no_shift_x,
                 self.big_frame_size * no_shift_y)

        coord_list.append(((x1,y1),(x - self.big_frame_size * no_shift_x,
                                    y - self.big_frame_size * no_shift_y),
                          (i + 1,j + 1)))

        self.image_no = len(coord_list)
        if self.random == True:
            shuffle(coord_list)
        return coord_list

    def image_generator(self):
        filtered_coord_list = []
        i = 0
        for coord in self.coord_list:
            coords,frame,mapping = coord
            try:
                curr_frame = np.array(self.openslide_handle.read_region(
                    coords,
                    0,
                    frame
                    ))[:,:,:3]
                yield coords,curr_frame,mapping
            except:
                yield None,None,None

def save_image(im,name):
    """Saves an image im as name in folder."""
    Image.fromarray(im).save(name)

slide_path = os.path.abspath(sys.argv[1])
output_path = os.path.abspath(sys.argv[2])

print(slide_path,output_path)

slide_path_split = slide_path.split(os.sep)
slide_name = slide_path_split[-1][:-5]
batch_name = slide_path_split[-2]
oss = OpenslideSplitter(slide_path,512)
count = 0
folder_name = '0'
try:
    os.makedirs('{}/{}'.format(output_path,folder_name))
except:
   pass
for coords,curr_frame,mapping in oss.images:
    try:
        name = '{0}/{1}/{2}_{3}_{4}_{5}_{6}_{7}.png'.format(
            output_path,
            folder_name,
            batch_name,
            slide_name,
            coords[0],
            coords[1],
            mapping[0],
            mapping[1]
        )
        print(name)
        if coords != None:
            save_image(curr_frame,name)
        if count % 1000 == 0:
            folder_name = str(count // 1000)
            try:
                os.makedirs('{}/{}'.format(output_path,folder_name))
            except:
                pass
        count += 1
    except Exception as e:
      pass
