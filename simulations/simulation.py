import numpy as np
import argparse
import os
from multiprocessing import Pool
import time
import pickle
from scipy import stats
from sklearn.datasets import make_spd_matrix
from sklearn import metrics
import ignite

from networks import *
from metrics import *
from data_simulation import *

parser = argparse.ArgumentParser(description='Test virtual cell classifier.')

parser.add_argument('--n_features',dest='n_features',
                    action='store',
                    type=int,
                    default=100)
parser.add_argument('--n_virtual_cells',dest='n_virtual_cells',
                    action='store',
                    type=int,
                    default=20)
parser.add_argument('--n_virtual_cells_model',dest='n_virtual_cells_model',
                    action='store',
                    type=int,
                    default=None)
parser.add_argument('--n_classes',dest='n_classes',
                    action='store',
                    type=int,
                    default=2)
parser.add_argument('--standard_deviation_centers',dest='standard_deviation_centers',
                    action='store',
                    type=float,
                    default=1.0)
parser.add_argument('--no_shared_variance',dest='no_shared_variance',
                    action='store_false')
parser.add_argument('--multivariate_normal',dest='multivariate_normal',
                    action='store_true')

parser.add_argument('--number_of_steps',dest='number_of_steps',
                    action='store',
                    type=int,
                    default=10000)
parser.add_argument('--batch_size',dest='batch_size',
                    action='store',
                    type=int,
                    default=100)
parser.add_argument('--number_of_cells',dest='number_of_cells',
                    action='store',
                    type=int,
                    default=500)
parser.add_argument('--learning_rate',dest='learning_rate',
                    action='store',
                    type=float,
                    default=0.0001)
parser.add_argument('--n_cpus',dest='n_cpus',
                    action='store',
                    type=int,
                    default=1)
parser.add_argument('--n_replicates',dest='n_replicates',
                    action='store',
                    type=int,
                    default=1)

args = parser.parse_args()

if not args.n_virtual_cells_model:
    args.n_virtual_cells_model = args.n_virtual_cells

N_FEATURES = args.n_features
N_VIRTUAL_CELLS = args.n_virtual_cells
N_CLASSES = args.n_classes

try:
    os.makedirs('simulation_outputs')
except:
    pass

metrics_workhorse = MetricFunction(
    {
        'Precision':precision,
        'Recall':recall,
        'Accuracy':accuracy,
        'F1-score':f1_score,
        'Confusion Matrix':confusion_matrix
    }
)

data_generator = DataGenerator(
    n_features=N_FEATURES,
    n_virtual_cells=N_VIRTUAL_CELLS,
    n_classes=N_CLASSES,
    n_cells=args.number_of_cells,
    batch_size=args.batch_size,
    standard_deviation_centers=args.standard_deviation_centers,
    shared_variance=args.no_shared_variance,
    multivariate_normal=args.multivariate_normal)

data_loader = torch.utils.data.DataLoader(
    data_generator,
    batch_size=args.batch_size,
    shuffle=True,
    num_workers=args.n_cpus)

output_dict = {}

for replicate in range(args.n_replicates):
    network = VirtualCellClassifier(
        n_input_features=N_FEATURES,
        n_virtual_cells=args.n_virtual_cells_model,
        n_classes=N_CLASSES)

    criterion = torch.nn.CrossEntropyLoss()
    optimizer = torch.optim.AdamW(network.parameters(),
                                  lr=args.learning_rate,
                                  weight_decay=1e-4)

    times = []
    loss_store = []
    print('Training replicate {}...'.format(replicate))
    for step,(batch,truth) in enumerate(data_loader):
        a = time.time()
        truth = torch.LongTensor(truth)
        batch = np.expand_dims(batch,1)
        batch = torch.Tensor(batch)
        output = network(batch)
        loss = criterion(output,truth)
        loss_store.append(float(loss))
        loss.backward()
        optimizer.step()
        times.append(time.time()-a)
        if args.n_classes <= 2:
            binary_output = torch.argmax(output,axis=-1)
        else:
            binary_output = torch.round(output)
            truth = torch.nn.functional.one_hot(truth,args.n_classes)
        metrics_workhorse.update(binary_output,truth)

        if step % 100 == 0:
            print('')
            print(''.join(['#' for _ in range(50)]))
            print(step,loss)
            print(''.join(['#' for _ in range(50)]))

            computed_metrics = metrics_workhorse.compute()
            for x in computed_metrics:
                m = computed_metrics[x]
                try:
                    print(x,'\n',m)
                except:
                    print(x,'\n',np.nan)
            metrics_workhorse.reset()
            print('Time (training_step):',np.mean(times))
        if step >= args.number_of_steps:
            break

    output_dict[replicate] = {
        'model':network.state_dict(),
        'metrics':metrics_workhorse,
        'data_generator':data_generator,
        'loss':loss_store
    }

output_name = 'simulation_outputs/NF{}_NVC{}_NVCM{}_NC{}_MVN{}_SV{}.pkl'.format(
    args.n_features,
    args.n_virtual_cells,
    args.n_virtual_cells_model,
    args.n_classes,
    args.multivariate_normal,
    args.no_shared_variance
)

with open(output_name,'wb') as o:
    pickle.dump(output_dict,o)
