print('-'*20)
import argparse
import sys
import torch
import time
from sklearn.model_selection import StratifiedKFold

from networks import *
from data_generator import *
from metrics import *

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
    parser.add_argument('--other_dataset_path',dest='other_dataset_path',
                        action='append',
                        type=str,
                        default=None)
    parser.add_argument('--labels_path',dest='labels_path',
                        action='append',
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
    parser.add_argument('--weight_decay',dest='weight_decay',
                        action='store',
                        type=float,
                        default=0.005)
    parser.add_argument('--dropout_rate',dest='dropout_rate',
                        action='store',
                        type=float,
                        default=0.)
    parser.add_argument('--n_folds',dest='n_folds',
                        action='store',
                        type=int,
                        default=1)
    parser.add_argument('--model_path',dest='model_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--random_state',dest='random_state',
                        action='store',
                        type=int,
                        default=1)
    parser.add_argument('--excluded_ids',dest='excluded_ids',
                        action='store',nargs='+',
                        type=str,
                        default=[])
    parser.add_argument('--median_impute',dest='median_impute',
                        action='store_true',
                        default=False)
    parser.add_argument('--range',dest='range',
                        action='store_true',
                        default=False)
    parser.add_argument('--min_cells',dest='min_cells',
                        action='store',type=int,
                        default=None)

    args = parser.parse_args()

    if torch.cuda.is_available():
          dev = "cuda:0"
    else:
          dev = "cpu"

    if args.median_impute == True:
        handle_nans = 'median_impute'
    else:
        handle_nans = 'remove'

    all_datasets = [GenerateFromDataset(x) for x in args.dataset_path]
    if args.other_dataset_path is not None:
        all_other_datasets = [CSVDataset(x,handle_nans=handle_nans)
                              for x in args.other_dataset_path]
        all_other_datasets_keys = [set(x.keys) for x in all_other_datasets]
    else:
        all_other_datasets = None

    loss_history = []
    all_lists = {}
    all_lists['labels'] = []
    all_lists['keys'] = []
    for label_path in args.labels_path:
        label_dict = get_classes(label_path)
        labels = Labels(label_dict)
        all_lists['labels'].append(labels)
        all_lists['keys'].extend(labels.keys)
    all_lists['keys'] = set(all_lists['keys'])
    all_lists['multiclass'] = {}
    for k in all_lists['keys']:
        mc = []
        for label in all_lists['labels']:
            if k in label.keys:
                mc.append(str(int(label[[k]])))
            else:
                mc.append('X')
        mc = '_'.join(mc)
        all_lists['multiclass'][k] = mc

    all_sets = [set(x.keys) for x in all_datasets]
    if args.other_dataset_path is not None:
        all_sets.extend(all_other_datasets_keys)
    all_sets.append(all_lists['keys'])
    all_keys = list(set.intersection(*all_sets))
    all_keys = [k for k in all_keys if k not in args.excluded_ids]
    all_mc = [all_lists['multiclass'][k] for k in all_keys]
    SF = StratifiedKFold(args.n_folds,shuffle=True,
                         random_state=args.random_state)
    sf = SF.split(all_keys,all_mc)
    folds = [f for f in sf]

    state_dict = {}
    state_dict['args'] = args

    # used to calculate CV AUC
    all_probs = [[] for l in args.labels_path]
    all_classes = [[] for l in args.labels_path]
    for fold in range(args.n_folds):
        print("\tFold: {}".format(fold),file=sys.stderr)
        metrics_workhorse = MetricFunction(
            {
                'Accuracy':accuracy,
                'Confusion Matrix':confusion_matrix,
                'AUC':auc,
                'Deviance':deviance,
            }
        )
        training_set = [all_keys[i] for i in folds[fold][0]]
        testing_set = [all_keys[i] for i in folds[fold][1]]
        all_networks = []

        for D in all_datasets:
            D.keys = training_set
            D.get_moments()
            D.keys = all_keys
        if all_other_datasets is not None:
            for D in all_other_datasets:
                D.get_moments(training_set)

        all_lists['queued_datasets'] = []
        for labels in all_lists['labels']:
            ts = [x for x in training_set if x in labels.keys]
            if args.range == False:
                Q = BatchGenerator(
                    all_datasets,labels,ts,
                    all_other_datasets,
                    n_cells=args.number_of_cells,
                    batch_size=args.batch_size)
            else:
                Q = BatchGenerator(
                    all_datasets,labels,ts,
                    all_other_datasets,
                    n_cells=args.number_of_cells,
                    batch_size=args.batch_size,
                    normalize=True,normalize_range=True)
            all_lists['queued_datasets'].append(Q)

        all_n_classes = [l.n_classes for l in all_lists['labels']]
        if all_other_datasets is not None:
            other_datasets_size = [x.n_features
                                   for x in all_other_datasets]
        else:
            other_datasets_size = [0]
        all_input_features = [d.n_features for d in all_datasets]
        all_n_virtual_cells = [args.n_virtual_cells for _ in all_datasets]
        stacked_network = VirtualCellClassifierStack(
            all_input_features,all_n_virtual_cells,
            other_datasets_size,all_n_classes,
            dropout_rate=args.dropout_rate).to(dev)
        def param_generator():
            for p in stacked_network.parameters():
                yield p

        all_lists['criterions'] = []
        for labels in all_lists['labels']:
            ts = [x for x in training_set if x in labels.keys]
            weights = np.zeros([len(ts),labels.n_classes])
            weights[(np.arange(len(ts)).astype(np.int32),
                    labels[ts].astype(np.int32))] = 1
            C = torch.Tensor(1.-np.sum(weights,axis=0)/np.sum(weights))
            criterion = torch.nn.CrossEntropyLoss(
                reduction='none',weight=C).to(dev)
            all_lists['criterions'].append(criterion)
        optimizer = torch.optim.AdamW(
            param_generator(),
            lr=args.learning_rate,
            weight_decay=args.weight_decay)

        schedule = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer,patience=100,
            factor=0.5,min_lr=args.learning_rate/10000)

        S = 0
        min_lr_counter = 0
        print("Training now...")
        times = []
        while True:
            if S > args.number_of_steps:
                break
            time_a = time.time()
            loss = torch.zeros([len(all_lists['criterions'])]).to(dev)
            for i in range(len(all_lists['criterions'])):
                criterion = all_lists['criterions'][i]
                queued_dataset = all_lists['queued_datasets'][i]
                B = queued_dataset.fetch_batch()
                sampled_datasets = [x[:,np.newaxis,:,:]
                                    for x in B['datasets']]
                truth = torch.LongTensor(B['labels']).to(dev)
                slide_ids = B['slide_ids']
                d = [torch.Tensor(x).to(dev) for x in sampled_datasets]
                if B['other_datasets'] is not None:
                    od = [torch.Tensor(x).to(dev) for x in B['other_datasets']]
                else:
                    od = None
                output_prob = stacked_network([d,od],i,args.dropout_rate>0.)
                L = criterion(output_prob,truth)
                # loss weight based on amount of evidence (number of cells)
                n_cells = np.array(B['n_cells'])
                W = np.minimum(args.number_of_cells,n_cells)
                W = np.sum(W/(1+args.number_of_cells*2),axis=0)
                W = torch.Tensor(1 + W)
                loss[i] = torch.mean(L*W)
            if len(all_lists['criterions']) > 1:
                loss = torch.Tensor(np.random.dirichlet(
                    [10. for _ in range(len(all_lists['criterions']))])) * loss
            loss = loss.sum()
            loss.backward()
            optimizer.step()
            times.append(time.time() - time_a)
            lll = loss.to('cpu').detach().numpy()
            loss_history.append(lll)
            schedule.step(lll)
            if schedule._last_lr[0] == args.learning_rate/10000:
                min_lr_counter += 1
            if S % 100 == 0:
                stacked_network.train(False)
                loss_acc = [torch.zeros([1])
                            for _ in range(len(all_lists['criterions']))]
                for _ in range(5):
                    loss = torch.zeros([1]).to(dev)
                    for i in range(len(all_lists['criterions'])):
                        criterion = all_lists['criterions'][i]
                        queued_dataset = all_lists['queued_datasets'][i]
                        B = queued_dataset.fetch_batch()
                        sampled_datasets = [x[:,np.newaxis,:,:]
                                            for x in B['datasets']]
                        truth = torch.LongTensor(B['labels']).to(dev)
                        slide_ids = B['slide_ids']
                        d = [torch.Tensor(x).to(dev) for x in sampled_datasets]
                        if B['other_datasets'] is not None:
                            od = [torch.Tensor(x).to(dev) for x in B['other_datasets']]
                        else:
                            od = None
                        output_prob = stacked_network([d,od],i)
                        l = criterion(output_prob,truth).to('cpu').detach().numpy()
                        l = np.mean(l)
                        loss_acc[i] += l / 10

                L = [x for x in loss_acc]
                print(S,L,"(curr. lr={}; av.time={:.2f}s)".format(
                    schedule._last_lr[0],float(np.mean(times))))
                stacked_network.train(True)
                times = []
                if min_lr_counter > 250:
                    break
            S += 1

        stacked_network.train(False)

        # calculating training set metrics
        for ob in range(len(all_lists['criterions'])):
            queued_dataset = all_lists['queued_datasets'][ob]
            labels = all_lists['labels'][ob]
            ts = [x for x in training_set if x in labels.keys]
            B = queued_dataset.fetch_batch(ts)
            sampled_datasets = [x[:,np.newaxis,:,:]
                                for x in B['datasets']]
            truth = torch.LongTensor(B['labels']).to(dev)
            slide_ids = B['slide_ids']
            d = [torch.Tensor(x).to(dev) for x in sampled_datasets]
            if B['other_datasets'] is not None:
                od = [torch.Tensor(x).to(dev) for x in B['other_datasets']]
            else:
                od = None
            output_prob = stacked_network([d,od],ob)

            output_onehot = torch.zeros_like(output_prob)
            output_onehot[(torch.arange(0,output_prob.shape[0]).long(),
                           torch.argmax(output_prob,dim=-1))] = 1

            truth_onehot = torch.zeros_like(output_prob)
            truth_onehot[(torch.arange(0,output_prob.shape[0]).long(),
                         truth)] = 1
            if truth_onehot.shape[1] == 2:
                # coherces to a binary classification problem
                truth_onehot = truth_onehot[:,1:]
                output_onehot = output_onehot[:,1:]
                output_prob = output_prob[:,1:]
            metrics_workhorse.update(
                truth_onehot.to('cpu').detach(),
                output_onehot.to('cpu').detach(),
                output_prob.to('cpu').detach())
            M = metrics_workhorse.compute()
            for metric in M:
                for i,m in enumerate(M[metric].flatten()):
                    print('TRAIN,{},{}_{},{},{},{},{}'.format(
                        fold,metric,i,float(m),args.n_virtual_cells,
                        args.number_of_cells,ob))
            metrics_workhorse.reset()

        # calculating testing set metrics
        for ob in range(len(all_lists['criterions'])):
            queued_dataset = all_lists['queued_datasets'][ob]
            labels = all_lists['labels'][ob]
            ts = [x for x in testing_set if x in labels.keys]
            if args.range == True:
                for dataset in queued_dataset.datasets:
                    dataset.get_min_max()
            B = queued_dataset.fetch_batch(ts)
            sampled_datasets = [x[:,np.newaxis,:,:]
                                for x in B['datasets']]
            truth = torch.LongTensor(B['labels']).to(dev)
            slide_ids = B['slide_ids']
            d = [torch.Tensor(x).to(dev) for x in sampled_datasets]
            if B['other_datasets'] is not None:
                od = [torch.Tensor(x).to(dev) for x in B['other_datasets']]
            else:
                od = None
            output_prob = stacked_network([d,od],ob)

            output_onehot = torch.zeros_like(output_prob)
            output_onehot[(torch.arange(0,output_prob.shape[0]).long(),
                           torch.argmax(output_prob,dim=-1))] = 1

            truth_onehot = torch.zeros_like(output_prob)
            truth_onehot[(torch.arange(0,output_prob.shape[0]).long(),
                         truth)] = 1
            if truth_onehot.shape[1] == 2:
                # coherces to a binary classification problem
                truth_onehot = truth_onehot[:,1:]
                output_onehot = output_onehot[:,1:]
                output_prob = output_prob[:,1:]

                all_probs[ob].append(
                    [float(x) for x in output_prob.cpu().detach().numpy()])
                all_classes[ob].append(
                    [float(x) for x in truth_onehot.cpu().detach().numpy()])
            else:
                arr_list = output_prob.cpu().detach().numpy().tolist()
                arr_list = [
                    ','.join([str(y) for y in x]) for x in arr_list
                ]
                all_probs[ob].append(arr_list)
                all_classes[ob].append(
                    [float(np.argmax(x)) 
                     for x in truth_onehot.cpu().detach().numpy()])

            metrics_workhorse.update(
                truth_onehot.to('cpu').detach(),
                output_onehot.to('cpu').detach(),
                output_prob.to('cpu').detach())
            M = metrics_workhorse.compute()
            for metric in M:
                for i,m in enumerate(M[metric].flatten()):
                    print('TEST,{},{}_{},{},{},{},{}'.format(
                        fold,metric,i,float(m),args.n_virtual_cells,
                        args.number_of_cells,ob))
            metrics_workhorse.reset()

        if args.model_path is not None:
            state_dict[fold] = {}
            state_dict[fold]['network'.format(i)] = stacked_network.to(
                'cpu').state_dict()

            state_dict[fold]['means'] = []
            state_dict[fold]['stds'] = []
            state_dict[fold]['minimum'] = []
            state_dict[fold]['maximum'] = []
            state_dict[fold]['training_set'] = training_set
            state_dict[fold]['testing_set'] = testing_set
            for D in all_datasets:
                state_dict[fold]['means'].append(D.mean)
                state_dict[fold]['stds'].append(D.std)
                state_dict[fold]['minimum'].append(D.minimum)
                state_dict[fold]['maximum'].append(D.maximum)
            if all_other_datasets is not None:
                for D in all_other_datasets:
                    state_dict[fold]['means'].append(D.mean)
                    state_dict[fold]['stds'].append(D.std)

    for ob in range(len(all_probs)):
        for fold_idx,(p,c) in enumerate(zip(all_probs[ob],all_classes[ob])):
            for p_,c_ in zip(p,c):
                print('PROBS_CLASS,{},{},{},{},{}'.format(
                    fold_idx,ob,p_,c_,labels.unique_labels[int(c_)]))
    if args.model_path is not None:
        torch.save(state_dict,args.model_path)
