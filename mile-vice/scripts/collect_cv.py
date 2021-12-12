import os
from glob import glob
import argparse

parser = argparse.ArgumentParser(description='Collect metrics.')
parser.add_argument('--log_files',dest='log_files',
                    nargs='+')
parser.add_argument('--subset',dest='subset',
                    action='store_true',default='full')
args = parser.parse_args()

all_csv_files = args.log_files

all_metrics = {}

all_possible_ids = [
    "anemia","binary","class","cv","disease","mds","multi","objective"]

for F in all_csv_files:
    d = os.path.split(F)[-1].split('.')
    datasets = ['cells']
    name = d[0]

    if args.subset == True:
        subset = "subset_features"
    else:
        subset = "all_features"
    name = '_'.join([x for x in name.split('_') if x in all_possible_ids])

    if 'bc' in d[-1]:
        datasets.append('bc')
    if 'dem' in d[-1]:
        datasets.append('dem')
    datasets = '_'.join(datasets)
    all_metrics[F] = {'TRAIN':{},'TEST':{}}
    with open(F) as o:
        for line in o:
            line = line.strip()
            if 'TRAIN' in line or 'TEST' in line:
                tt,fold,metric = line.split(',')[0:3]
                task = line.split(',')[-1]
                line = '{},{},{},{}'.format(line,name,datasets,subset)
                if task not in all_metrics[F][tt]:
                    all_metrics[F][tt][task] = {}
                if fold not in all_metrics[F][tt][task]:
                    all_metrics[F][tt][task][fold] = {metric:line}
                else:
                    all_metrics[F][tt][task][fold][metric] = line

for F in all_metrics:
    for tt in ['TRAIN','TEST']:
        for task in all_metrics[F][tt]:
            for fold in all_metrics[F][tt][task]:
                for metric in all_metrics[F][tt][task][fold]:
                    print(all_metrics[F][tt][task][fold][metric])
