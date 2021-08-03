import h5py
import numpy as np
import torch
from PIL import Image

class DataGenerator:
    def __init__(self,path,transform=None):
        self.path = path
        self.transform = transform
        self.dataset = h5py.File(path,'r')
        self.keys = [x for x in self.dataset.keys()]
        self.get_match()

    def get_match(self):
        self.match = []
        for k in self.keys:
            n = [[k,i] for i in range(self.dataset[k].shape[0])]
            self.match.extend(n)

    def iterate(self):
        K = self.keys.copy()
        size = len(K)
        while True:
            np.random.shuffle(K)
            for key in K:
                d = self.dataset[key]
                n_images = d.shape[0]
                x = np.random.choice(n_images)
                o = d[x,:,:,:]
                yield o

    def __getitem__(self,i):
        k = self.match[i][0]
        j = self.match[i][1]
        o = self.dataset[k][j,:,:,:3]
        if self.transform:
            o = Image.fromarray(o)
            o = self.transform(o)
        return o,torch.Tensor([1.])

    def __len__(self):
        return len(self.match)
