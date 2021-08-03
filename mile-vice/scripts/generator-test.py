import argparse
import sys
import torch
import psutil

from networks import VirtualCellClassifier
from data_generator import *
from metrics import *

def split_dataset(dataset,p=0.7):
    train_set = np.random.choice(len(dataset),size=[int(p*len(dataset))],
                                 replace=False)
    test_set = [x for x in range(len(dataset)) if x not in train_set]
    return train_set,test_set

def get_classes(path):
    with open(path,'r') as o:
        tmp = [x.strip().split(',') for x in o.readlines()[1:]]
    return {x[0]:x[1] for x in tmp}

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

if __name__ == "__main__":
    print("Training")

    parser = argparse.ArgumentParser(description='Test virtual cell classifier.')

    parser.add_argument('--dataset_path',dest='dataset_path',
                        action='append',
                        type=str,
                        default=None)
    parser.add_argument('--labels_path',dest='labels_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--number_of_steps',dest='number_of_steps',
                        action='store',
                        type=int,
                        default=10000)
    parser.add_argument('--batch_size',dest='batch_size',
                        action='store',
                        type=int,
                        default=100)
    parser.add_argument('--n_folds',dest='n_folds',
                        action='store',
                        type=int,
                        default=10)

    args = parser.parse_args()

    if torch.cuda.is_available():
          dev = "cuda:0"
    else:
          dev = "cpu"
    print(dev)

    all_datasets = [GenerateFromDataset(x) for x in args.dataset_path]

    label_dict = get_classes(args.labels_path)
    labels = Labels(label_dict)
    all_sets = [set(x.keys) for x in all_datasets]
    all_sets.append(set(labels.keys))
    all_keys = list(set.intersection(*all_sets))
    folds = [split_dataset(all_keys) for _ in range(args.n_folds)]

    for fold in range(args.n_folds):
        print('Fold: {}'.format(fold))

        training_set = [all_keys[i] for i in folds[fold][0]]
        for D in all_datasets:
            D.keys = training_set
            D.get_moments()
            D.keys = all_keys

        for dataset in all_datasets:
            for key in training_set:
                D = dataset.all_datasets[key]
                for _ in range(10):
                    batch = D.return_n_cells(10000)
                    d = np.isfinite(batch) == False
                    post_batch = (batch - dataset.minimum)/(dataset.maximum - dataset.minimum)
                    post_d = np.isfinite(post_batch) == False
                    L = labels[[key]]
                    x = post_d.sum(0).sum(0)
                    if np.all(np.isfinite(d)) == False or np.all(np.isfinite(post_d)) == False:
                        print(np.isnan(d).sum(),np.isnan(post_d).sum())

                        print(dataset.path,batch.min() - batch.max(),
                              post_batch.min() - post_batch.max())
