import h5py
import numpy as np

def gen(i):
    with h5py.File('test/{}.h5'.format(i),'w') as f:
        f['mean'] = np.ones([1])
        g = f.create_group('cells')
        g['0'] = np.ones([1000,1]) * i

for x in range(10):
    gen(x)
