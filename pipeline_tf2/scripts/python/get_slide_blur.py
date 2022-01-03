import argparse
import numpy as np
import cv2

from image_generator import image_generator_slide

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get cell examples.')

    parser.add_argument('--slide_path',dest='slide_path',
                        action='store',
                        type=str,
                        default=None)
    args = parser.parse_args()

    for tile,coords in image_generator_slide(args.slide_path,512,512):
        blur_var = cv2.Laplacian(
            np.uint8(np.mean(tile,axis=-1)),cv2.CV_64F).var()

        print('BLUR,{},{}'.format(coords,blur_var))
