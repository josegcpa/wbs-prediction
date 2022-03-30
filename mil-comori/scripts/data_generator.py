import h5py
from glob import glob
import os
import numpy as np
from multiprocessing import Queue,Process

# Contains scripts for data generation

class HDF5Dataset:
    def __init__(self,F,feature_subset=None):
        self.F = F
        self.get_n_features()
        self.get_n_cells()

    def get_n_features(self):
        self.n_features = self.F['cells']['0'].shape[1]

    def get_n_cells(self):
        self.n_subsets = len(self.F['cells'])
        self.subset_size = self.F['cells']['0'].shape[0]
        self.n_cells = sum([self.F['cells'][x].shape[0] for x in self.F['cells']])

    def return_n_cells(self,n):
        output = np.zeros([n,self.n_features])
        S = np.random.choice(self.n_cells,size=n,
                             replace=n > self.n_cells)
        S = np.sort(S)
        S_subset = S // self.subset_size
        S_index = S % self.subset_size
        acc = 0
        for _s in np.unique(S_subset):
            idx = S_index[S_subset == _s]
            if len(idx) > 0:
                tmp = self.F['cells'][str(_s)][()]
                output[acc:(acc+len(idx)),:] = tmp[idx,:]
                acc += len(idx)
        return output

    def return_all_cells(self):
        all_cells = []
        for s in self.F['cells']:
            all_cells.append(self.F['cells'][str(s)][()])
        all_cells = np.concatenate(all_cells,axis=0)
        return all_cells

class BatchGenerator:
    def __init__(self,datasets,labels,id_subset,
                 other_datasets=None,
                 n_cells=1000,batch_size=25,
                 normalize=True,normalize_range=False):
        self.datasets = datasets
        self.labels = labels
        self.id_subset = id_subset
        self.other_datasets = other_datasets
        self.n_cells = n_cells
        self.batch_size = batch_size
        self.normalize = normalize
        self.normalize_range = normalize_range

    def fetch_batch(self,ids=None):
        if ids is None:
            ids = np.random.choice(
                self.id_subset,size=self.batch_size,
                replace=self.batch_size>len(self.datasets[0].keys))
        x = [d.generate_n_cells(
             ids,self.n_cells,
             norm=self.normalize,normalize_range=self.normalize_range)
             for d in self.datasets]
        n_cells = [d.get_n_cells(ids) for d in self.datasets]
        if self.other_datasets is not None:
            other_x = [d[ids] for d in self.other_datasets]
        else:
            other_x = None
        return {'datasets':x,
                'other_datasets':other_x,
                'labels':self.labels[ids],
                'slide_ids':ids,
                'n_cells':n_cells}

class QueueGenerator:
    def __init__(self,datasets,labels,id_subset,
                 n_cells,batch_size,
                 maxsize=1,normalize=True):
        self.datasets = datasets
        self.labels = labels
        self.id_subset = id_subset
        self.n_cells = n_cells
        self.batch_size = batch_size
        self.maxsize = maxsize
        self.normalize = normalize
        self.go = True

        self.q = Queue(self.maxsize)
        self.p = Process(target=self.gen,args=((self.q),))
        self.p.daemon = True
        self.p.start()

    def gen(self,q):
        while self.go == True:
            random_subset = np.random.choice(
                self.id_subset,size=self.batch_size,
                replace=self.batch_size>len(self.datasets[0].keys))
            x = [d.generate_n_cells(random_subset,self.n_cells,self.normalize)
                 for d in self.datasets]
            q.put({'datasets':x,
                   'labels':self.labels[random_subset],
                   'slide_ids':random_subset})
        self.q.put(None)

    def fetch_from_queue(self):
        x = self.q.get()
        return x

class GenerateFromDataset:
    def __init__(self,path,maxsize=1,min_cells=None,feature_subset=None,auto=True):
        self.maxsize = maxsize
        self.min_cells = min_cells
        self.mean = None
        self.std = None
        self.path = path
        self.dataset = h5py.File(self.path,'r')
        self.keys = [x for x in self.dataset.keys()]
        self.read_all_datasets()
        self.n_datasets = len(self.keys)
        self.n_features = self.all_datasets[self.keys[0]].n_features
        if auto == True:
            self.get_moments()

    def __getitem__(self,idx):
        return self.all_datasets[idx]

    def read_all_datasets(self):
        self.all_datasets = {}
        for k in self.keys:
            try:
                d = HDF5Dataset(self.dataset[k])
                if self.min_cells:
                    if d.n_cells >= self.min_cells:
                        self.all_datasets[k] = d
                else:
                    self.all_datasets[k] = d
            except:
                pass
        self.keys = [x for x in self.all_datasets.keys()]

    def generate_n_cells(self,S,n_cells,norm=True,normalize_range=False):
        n_datasets = len(S)
        output = np.zeros([n_datasets,n_cells,self.n_features])
        for d in range(n_datasets):
            x = self.all_datasets[S[d]].return_n_cells(n_cells)
            output[d,:,:] = x
        if norm == True:
            output = output - self.mean
            output /= self.std
        if normalize_range == True:
            output = (output-self.minimum)/(self.maximum-self.minimum)
        return output

    def generate_n_cells_flat(self,S,n_cells):
        n_datasets = len(S)
        output = []
        for d in range(n_datasets):
            x = self.all_datasets[S[d]].return_n_cells(n_cells)
            output.append(x)
        output = np.concatenate(output,axis=0)
        return output

    def get_n_cells(self,S):
        """
        Returns the number of cells for a given sample.
        """
        n_datasets = len(S)
        output = np.zeros([n_datasets])
        for d in range(n_datasets):
            x = self.all_datasets[S[d]].n_cells
            output[d] = x
        return output

    def get_moments(self):
        self.minimum = np.ones(self.n_features) * np.inf
        self.maximum = np.ones(self.n_features) * -np.inf
        self.s = np.zeros(self.n_features,dtype=np.float64)
        self.ss = np.zeros(self.n_features,dtype=np.float64)
        n = 0
        for key in self.keys:
            x = self.all_datasets[key].return_n_cells(
                self.all_datasets[key].n_cells)
            n += x.shape[0]
            m = x.min(axis=0)
            M = x.max(axis=0)
            s = x.sum(axis=0)
            ss = np.sum(x**2,axis=0)
            # update
            self.minimum = np.where(m < self.minimum,m,self.minimum)
            self.maximum = np.where(M > self.maximum,M,self.maximum)
            self.s += s
            self.ss += ss
        self.mean = self.s/n
        self.var = (self.ss - (self.s**2)/n)/n
        self.std = np.sqrt(self.var)
        self.minimum = (self.minimum - self.mean) / self.std
        self.maximum = (self.maximum - self.mean) / self.std

        self.mean = self.mean[np.newaxis,np.newaxis,:]
        self.var = self.var[np.newaxis,np.newaxis,:]
        self.std = self.std[np.newaxis,np.newaxis,:]
        self.minimum = self.minimum[np.newaxis,np.newaxis,:]
        self.maximum = self.maximum[np.newaxis,np.newaxis,:]

    def get_moments_old(self):
        big_mat = []
        for key in self.keys:
            big_mat.append(self.all_datasets[key].return_n_cells(
                self.all_datasets[key].n_cells))
        big_mat = np.concatenate(big_mat,axis=0)
        self.mean = np.mean(big_mat,axis=0)[np.newaxis,np.newaxis,:]
        self.var = np.var(big_mat,axis=0)[np.newaxis,np.newaxis,:]
        self.std = np.sqrt(self.var)
        big_mat = (big_mat - self.mean) / self.std
        big_mat = big_mat[0,:,:]
        self.maximum = np.max(big_mat,axis=0)
        self.maximum = self.maximum[np.newaxis,np.newaxis,:]
        self.minimum = np.min(big_mat,axis=0)
        self.minimum = self.minimum[np.newaxis,np.newaxis,:]

    def get_min_max(self,keys=None):
        big_mat = []
        if keys is None:
            keys = self.keys
        for key in keys:
            big_mat.append(self.all_datasets[key].return_n_cells(
                self.all_datasets[key].n_cells))
        big_mat = np.concatenate(big_mat,axis=0)
        big_mat = (big_mat - self.mean) / self.std
        big_mat = big_mat[0,:,:]
        self.maximum = np.max(big_mat,axis=0)
        self.maximum = self.maximum[np.newaxis,np.newaxis,:]
        self.minimum = np.min(big_mat,axis=0)
        self.minimum = self.minimum[np.newaxis,np.newaxis,:]

class Labels:
    def __init__(self,class_dict):
        self.class_dict = class_dict
        self.labels = [x[1] for x in self.class_dict.items()]
        self.keys = [x for x in self.class_dict.keys()]
        self.keys_correspondence = {x:i for i,x in enumerate(self.keys)}
        self.make_one_hot()

    def make_one_hot(self):
        unique_labels = np.sort(np.unique(self.labels))
        self.one_hot = np.zeros([len(self.keys)])
        for i,label in enumerate(self.labels):
            self.one_hot[i] = np.where(unique_labels == label)[0]
        self.unique_labels = unique_labels
        self.n_classes = len(unique_labels)

    def __getitem__(self,ids):
        ids_idx = [self.keys_correspondence[x] for x in ids]
        return self.one_hot[ids_idx]

class CSVDataset:
    def __init__(self,
                 csv_path,
                 prop_num=0.5,
                 handle_nans=None):
        # data has to be in format id,feature1,feature2,...
        self.csv_path = csv_path
        self.pn = prop_num
        self.get_features()
        self.normalize = False
        self.handle_nans = handle_nans
        if self.handle_nans == 'remove':
            self.remove_nans()
            self.get_moments(self.keys)
        elif self.handle_nans == 'median_impute':
            self.medians = np.ones([1,self.final_array.shape[1]])*np.nan
            self.get_moments(self.keys)
        else:
            raise Exception(
                "{} not implemented for handle_nans in CSVDataset!".format(
                    handle_nans))

    def is_numeric(self,s):
        n_dot = 0
        if s == '':
            return False
        for i,c in enumerate(s):
            if c == '-':
                if i != 0:
                    return False
            if c.isdigit() == True:
                pass
            elif c == '.':
                n_dot += 1
            else:
                return False
            if n_dot > 1:
                return False
        return True

    def make_one_hot(self,ids):
        unique_ids = [x for x in set(ids) if x != 'nan']
        unique_ids = np.sort(unique_ids)
        one_hot = np.zeros([len(ids),len(unique_ids)])
        for i,label in enumerate(ids):
            if label == 'nan':
                one_hot[i,:] = np.nan
            else:
                one_hot[i,np.where(unique_ids == label)] = 1
        return one_hot

    def get_features(self):
        with open(self.csv_path) as o:
            lines = [x.strip().split(',') for x in o.readlines()]

        self.columns_to_normalize = []
        self.n_cols = len(lines[0])
        self.n_feature_cols = self.n_cols - 1
        self.ids = []

        all_columns = [[] for x in range(self.n_feature_cols)]

        for line in lines:
            self.ids.append(line[0])
            for i,f in enumerate(line[1:]):
                all_columns[i].append(f)

        self.ids = {x:i for i,x in enumerate(self.ids)}
        self.keys = self.ids.keys()

        self.num_cols_ext = [list(map(self.is_numeric,x))
                             for x in all_columns]
        self.num_cols = [sum(x)>(len(self.ids)*self.pn)
                         for x in self.num_cols_ext]

        final_columns = []
        a = 0
        for i,s in enumerate(self.num_cols):
            if self.num_cols[i] == True:
                col = map(lambda x: x[0] if x[1] == True else 'nan',
                          zip(all_columns[i],self.num_cols_ext[i]))
                arr = np.array(list(col),dtype=np.float32)[:,np.newaxis]
                self.columns_to_normalize.append(a)
                final_columns.append(arr)
                a += 1
            else:
                oh = self.make_one_hot(all_columns[i])
                final_columns.append(oh)
                a += oh.shape[1]

        self.final_array = np.concatenate(final_columns,1)
        self.n_features = self.final_array.shape[1]

    def remove_nans(self):
        nan_idx = np.any(np.isnan(self.final_array),axis=1)
        nan_idx = np.where(np.logical_not(nan_idx))[0]
        nan_idx = np.unique(nan_idx)
        self.final_array = self.final_array[nan_idx,:]
        self.ids = [x for i,x in enumerate(self.ids) if i in nan_idx]
        self.ids = {x:i for i,x in enumerate(self.ids)}
        self.keys = self.ids.keys()

    def get_moments(self,id_list):
        self.normalize = False
        arr = self.__getitem__(id_list)
        if self.handle_nans == "median_impute":
            self.medians = np.zeros([1,arr.shape[1]])
            for i in range(arr.shape[1]):
                self.medians[0,i] = np.nanmedian(arr[:,i])
        self.mean = np.zeros([1,arr.shape[1]])
        self.std = np.ones([1,arr.shape[1]])
        for i in self.columns_to_normalize:
            self.mean[0,i] = np.nanmean(arr[:,i])
            self.std[0,i] = np.nanstd(arr[:,i])
        self.normalize = True

    def __getitem__(self,id_list):
        id_list = [self.ids[x] for x in id_list]
        o = self.final_array[id_list,:]
        if self.handle_nans == "median_impute":
            for i in range(o.shape[1]):
                o[np.isnan(o[:,i]),i] = self.medians[0,i]
        if self.normalize == True:
            o = (o - self.mean)/self.std
        return o

def split_dataset(dataset,p=0.7):
    train_set = np.random.choice(len(dataset),size=[int(p*len(dataset))],
                                 replace=False)
    test_set = [x for x in range(len(dataset)) if x not in train_set]
    return train_set,test_set

def split_dataset_stratified(dataset,labels,p=0.7):
    train_set = []
    test_set = []
    for l in np.unique(labels):
        dataset_label = [dataset[i]
                         for i in range(len(dataset)) if labels[i] == l]
        tr,te = split_dataset(dataset_label,p)
        train_set.extend(tr)
        test_set.extend(te)
    return train_set,test_set

def get_classes(path):
    with open(path,'r') as o:
        tmp = [x.strip().split(',') for x in o.readlines()[1:]]
    return {x[0]:x[1] for x in tmp}
