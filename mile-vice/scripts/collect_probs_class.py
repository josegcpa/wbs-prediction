import os
from glob import glob

all_csv_files = glob('logs/cv*o')

all_metrics = {}

for F in all_csv_files:
    d = os.path.split(F)[-1].split('.')
    datasets = ['cells']
    name = d[0]
    if 'bc' in d[-2]:
        datasets.append('bc')
    if 'dem' in d[-2]:
        datasets.append('dem')
    datasets = '_'.join(datasets)
    all_metrics[F] = {'PROBS_CLASS':{}}
    with open(F) as o:
        for line in o:
            line = line.strip()
            if 'PROBS_CLASS' in line:
                tt,fold,metric = line.split(',')[0:3]
                task = line.split(',')[-1]
                line = '{},{},{}'.format(line,name,datasets)
                if task not in all_metrics[F][tt]:
                    all_metrics[F][tt][task] = {}
                if fold not in all_metrics[F][tt][task]:
                    all_metrics[F][tt][task][fold] = {metric:line}
                else:
                    all_metrics[F][tt][task][fold][metric] = line

for F in all_metrics:
    for tt in ['PROBS_CLASS']:
        for task in all_metrics[F][tt]:
            for fold in all_metrics[F][tt][task]:
                for metric in all_metrics[F][tt][task][fold]:
                    pass #print(all_metrics[F][tt][task][fold][metric])
