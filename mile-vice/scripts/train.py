print('-'*20)
import argparse
import sys
import torch
import psutil
import pickle

from networks import VirtualCellClassifier
from data_generator import *
from metrics import *

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
    print("Reading cmd line arguments...")
    parser = argparse.ArgumentParser(description='Test virtual cell classifier.')

    parser.add_argument('--n_virtual_cells',dest='n_virtual_cells',
                        action='store',
                        type=int,
                        default=20)
    parser.add_argument('--n_classes',dest='n_classes',
                        action='store',
                        type=int,
                        default=2)

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
                        default=32)
    parser.add_argument('--number_of_cells',dest='number_of_cells',
                        action='store',
                        type=int,
                        default=500)
    parser.add_argument('--learning_rate',dest='learning_rate',
                        action='store',
                        type=float,
                        default=0.0001)
    parser.add_argument('--n_folds',dest='n_folds',
                        action='store',
                        type=int,
                        default=1)
    parser.add_argument('--model_path',dest='model_path',
                        action='store',
                        type=str,
                        default=None)

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

    folds = [split_dataset_stratified(all_keys,labels[all_keys]) for _ in range(args.n_folds)]

    for fold in range(args.n_folds):
        print("\tFold: {}".format(fold),file=sys.stderr)
        metrics_workhorse = MetricFunction(
            {
                'Precision':precision,
                'Recall':recall,
                'Accuracy':accuracy,
                'F1-score':f1_score,
                'Confusion Matrix':confusion_matrix,
                'AUC':auc
            }
        )
        training_set = [all_keys[i] for i in folds[fold][0]]
        testing_set = [all_keys[i] for i in folds[fold][1]]
        all_networks = []

        for D in all_datasets:
            D.keys = training_set
            D.get_moments()
            D.keys = all_keys

        queued_dataset = BatchGenerator(
            all_datasets,labels,training_set,
            n_cells=args.number_of_cells,batch_size=args.batch_size)

        for dataset in all_datasets:
            networks = VirtualCellClassifier(
                n_input_features=dataset.n_features,
                n_virtual_cells=args.n_virtual_cells,
                n_classes=labels.n_classes).to(dev)
            all_networks.append(networks)

        final_layer = torch.nn.Sequential(
            torch.nn.Linear(args.n_virtual_cells*len(all_networks),
                            labels.n_classes),
            torch.nn.Softmax(dim=-1)).to(dev)

        def param_generator():
            for network in all_networks:
                for p in network.parameters():
                    yield p
            for p in final_layer.parameters():
                yield p

        weights = np.zeros([len(training_set),labels.n_classes])
        weights[(np.arange(len(training_set)).astype(np.int32),
                 labels[training_set].astype(np.int32))] = 1
        C = torch.Tensor(1.-np.sum(weights,axis=0)/np.sum(weights))

        criterion = torch.nn.CrossEntropyLoss(weight=C).to(dev)
        optimizer = torch.optim.AdamW(
            param_generator(),
            lr=args.learning_rate,
            weight_decay=5e-2)

        schedule = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer,factor=0.8,min_lr=0.000001)

        S = 0
        print("Training now...")
        while True:
            if S > args.number_of_steps:
                break

            outputs = []
            B = queued_dataset.fetch_batch()
            sampled_datasets = B['datasets']
            truth = torch.LongTensor(B['labels']).to(dev)
            slide_ids = B['slide_ids']
            for batch,network in zip(sampled_datasets,all_networks):
                batch = batch[:,np.newaxis,:,:]
                batch = torch.Tensor(batch).to(dev)
                outputs.append(network.get_virtual_cells(batch))
            output = torch.cat(tuple(outputs),dim=-1)
            output_prob = final_layer(output)
            loss = criterion(output_prob,truth)
            loss.backward()
            optimizer.step()
            if torch.isnan(loss):
                print('!')
                for i in range(args.batch_size):
                    if torch.any(torch.isnan(output_prob[i,:])):
                        print(output_prob[i,:],truth[i],slide_ids[i])
                        nf = [(1-np.isfinite(x)).astype(np.bool) for x in B['datasets']]
                        print(B['datasets'][0][nf[0]],B['datasets'][1][nf[1]])
                        w = [np.where(x.sum(0).sum(0)) for x in nf]
                        print(w)
                        print([x[n] for x,n in zip(B['datasets'],nf)])
                print(telomeroror)

            if S % 50 == 0:
                network.train(False)
                loss_acc = []
                for _ in range(10):
                    random_subset = np.random.choice(training_set,size=args.batch_size)
                    outputs = []
                    truth = torch.LongTensor(labels[random_subset]).to(dev)
                    for dataset,network in zip(all_datasets,all_networks):
                        batch = dataset.generate_n_cells(random_subset,
                                                         args.number_of_cells)
                        batch = batch[:,np.newaxis,:,:]
                        batch = torch.Tensor(batch).to(dev)
                        outputs.append(network.get_virtual_cells(batch))
                    output = torch.cat(tuple(outputs),dim=-1)
                    output_prob = final_layer(output)
                    loss = criterion(output_prob,truth)
                    loss_acc.append(loss.to('cpu').detach().numpy())
                L = np.mean(loss_acc)
                schedule.step(loss)
                print(S,L,"(current learning rate={})".format(schedule._last_lr))
                network.train(True)
            S += 1

        truth = torch.LongTensor(labels[training_set])
        outputs = []
        for dataset,network in zip(all_datasets,all_networks):
            batch = dataset.generate_n_cells(training_set,
                                             args.number_of_cells)
            batch = batch[:,np.newaxis,:,:]
            batch = torch.Tensor(batch).to(dev)
            outputs.append(network.get_virtual_cells(batch))
        output = torch.cat(tuple(outputs),dim=-1)
        output_prob = final_layer(output).to('cpu').detach()
        output_onehot = torch.zeros_like(output_prob)
        output_onehot[(torch.arange(0,output_prob.shape[0]).long(),
                       torch.argmax(output_prob,dim=-1))] = 1

        truth_onehot = torch.zeros_like(output_prob)
        truth_onehot[(torch.arange(0,output_prob.shape[0]).long(),
                      truth)] = 1
        metrics_workhorse.update(truth_onehot.to('cpu').detach(),
                                 output_onehot.to('cpu').detach())
        M = metrics_workhorse.compute()
        for metric in M:
            for i,m in enumerate(M[metric].flatten()):
                print('TRAIN,{},{}_{},{},{},{}'.format(
                    fold,metric,i,float(m),args.n_virtual_cells,
                    args.number_of_cells))
        metrics_workhorse.reset()

        truth = torch.LongTensor(labels[testing_set])
        outputs = []
        for dataset,network in zip(all_datasets,all_networks):
            batch = dataset.generate_n_cells(testing_set,
                                             args.number_of_cells)
            batch = batch[:,np.newaxis,:,:]
            batch = torch.Tensor(batch).to(dev)
            outputs.append(network.get_virtual_cells(batch))
        output = torch.cat(tuple(outputs),dim=-1)
        output_prob = final_layer(output).to('cpu').detach()
        output_onehot = torch.zeros_like(output_prob)
        output_onehot[(torch.arange(0,output_prob.shape[0]).long(),
                       torch.argmax(output_prob,dim=-1))] = 1

        truth_onehot = torch.zeros_like(output_prob)
        truth_onehot[(torch.arange(0,output_prob.shape[0]).long(),
                      truth)] = 1
        metrics_workhorse.update(truth_onehot.to('cpu').detach(),
                                 output_onehot.to('cpu').detach())
        M = metrics_workhorse.compute()
        for metric in M:
            for i,m in enumerate(M[metric].flatten()):
                print('TEST,{},{}_{},{},{},{}'.format(
                    fold,metric,i,float(m),args.n_virtual_cells,
                    args.number_of_cells))
        metrics_workhorse.reset()

        if args.model_path is not None:
            model_path_fold = '{}_{}'.format(args.model_path,fold)
            state_dict = {}
            for i,network in enumerate(all_networks):
                state_dict['network_{}'.format(i)] = network.state_dict()
            state_dict['final_layer'] = final_layer.state_dict()
            with open(model_path_fold,'wb') as o:
                pickle.dump(state_dict,o)
