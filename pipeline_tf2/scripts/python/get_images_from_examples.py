import h5py
import numpy as np
from PIL import Image
import os
import re

from glob import glob

try:
    os.makedirs("images_from_examples")
except:
    pass

all_examples = glob("examples/*wbc.h5")

mll_examples = [x for x in all_examples if re.search('.*[A-Z]+_[0-9]+_wbc',x)][:10]
adden_examples = [x for x in all_examples if not re.search('.*[A-Z]+_[0-9]+_wbc',x)]

N = 20
H,W = 210,105

all_isolated_images_mll = []

for file_path in mll_examples:
    
    with h5py.File(file_path,"r") as F:
        root = file_path.split(os.sep)[-1]

        print(file_path,len(F))

        patch_work = np.zeros([H*N,W*N,3],dtype=np.uint8)

        for i,k in enumerate(F):
            image = F[k]['image'][::]
            seg_image = F[k]['isolated_image'][::]
            isolated_image = np.where(
                seg_image[:,:,np.newaxis] > 0,image,np.uint8(image * 0.5)
            )
            x,y = np.where(seg_image > 0)

            isolated_image = np.concatenate([image,isolated_image])
            sh = isolated_image.shape

            all_isolated_images_mll.append((isolated_image,F,isolated_image.shape))

            if i < (N**2):
                place_x,place_y = int(i // N),(i % N)
                place_x_,place_y_ = place_x*H,place_y*W
                place_x = int(place_x_ + int(H/2 - sh[0] / 2))
                place_y = int(place_y_ + int(W/2 - sh[1] / 2))
                try:
                    patch_work[place_x:(place_x+sh[0]),place_y:(place_y+sh[1])] = isolated_image
                except:
                    pass

    output_image_path = "images_from_examples/MLL_wbc_{}.png".format(root)
    if os.path.exists(output_image_path) == False:
        patch_work = patch_work[:(place_x_ + H),:]
        Image.fromarray(patch_work).save(output_image_path)

N = 25

np.random.seed(4422)

isolated_image_idx = np.random.choice(len(all_isolated_images_mll),N*N,replace=False)

isolated_images = [all_isolated_images_mll[i] 
                   for i in isolated_image_idx]

H,W = [np.max([x[2][0] for x in isolated_images])+2,
       np.max([x[2][1] for x in isolated_images])+2]

patch_work = np.zeros([H*N,W*N,3],dtype=np.uint8)

for i,k in enumerate(isolated_images):
    isolated_image = k[0]
    sh = isolated_image.shape
    place_x,place_y = int(i // N),(i % N)
    place_x_,place_y_ = place_x*H,place_y*W
    place_x = int(place_x_ + int(H/2 - sh[0] / 2))
    place_y = int(place_y_ + int(W/2 - sh[1] / 2))

    patch_work[place_x:(place_x+sh[0]),place_y:(place_y+sh[1])] = isolated_image

Image.fromarray(patch_work).save('images_from_examples/MLL_wbc.png')

N = 20
H,W = 210,105

all_isolated_images_adden = []

for file_path in adden_examples:
    with h5py.File(file_path,"r") as F:
        root = file_path.split(os.sep)[-1]
        print(file_path,len(F))

        patch_work = np.zeros([H*N,W*N,3],dtype=np.uint8)

        for i,k in enumerate(F):
            image = F[k]['image'][::]
            seg_image = F[k]['isolated_image'][::]
            isolated_image = np.where(
                seg_image[:,:,np.newaxis] > 0,image,np.uint8(image * 0.5)
            )
            x,y = np.where(seg_image > 0)

            isolated_image = np.concatenate([image,isolated_image])
            sh = isolated_image.shape

            all_isolated_images_adden.append((isolated_image,F,isolated_image.shape))

            if i < N**2:
                place_x,place_y = int(i // N),(i % N)
                place_x_,place_y_ = place_x*H,place_y*W
                place_x = int(place_x_ + int(H/2 - sh[0] / 2))
                place_y = int(place_y_ + int(W/2 - sh[1] / 2))

                patch_work[place_x:(place_x+sh[0]),place_y:(place_y+sh[1])] = isolated_image

    output_image_path = "images_from_examples/AC2_wbc_{}.png".format(root)
    if os.path.exists(output_image_path) == False:
        patch_work = patch_work[:(place_x_ + H),:]
        Image.fromarray(patch_work).save(output_image_path)

N = 30

np.random.seed(4242)

isolated_image_idx = np.random.choice(len(all_isolated_images_adden),N*N,replace=False)

isolated_images = [all_isolated_images_adden[i] 
                   for i in isolated_image_idx]

H,W = [np.max([x[2][0] for x in isolated_images])+2,
       np.max([x[2][1] for x in isolated_images])+2]

patch_work = np.zeros([H*N,W*N,3],dtype=np.uint8)

for i,k in enumerate(isolated_images):
    isolated_image = k[0]
    sh = isolated_image.shape
    place_x,place_y = int(i // N),(i % N)        
    place_x_,place_y_ = place_x*H,place_y*W
    place_x = int(place_x_ + int(H/2 - sh[0] / 2))
    place_y = int(place_y_ + int(W/2 - sh[1] / 2))

    patch_work[place_x:(place_x+sh[0]),place_y:(place_y+sh[1])] = isolated_image

Image.fromarray(patch_work).save('images_from_examples/AC2_wbc.png')
