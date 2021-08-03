import os
from glob import glob

all_csv_files = glob('logs/cv*o')

all_metrics = {}

for F in all_csv_files:
    d = os.path.split(F)[-1].split('_')
    datasets = ['cells']
    try:
        float(d[2])
        name = d[1]
    except:
        name = '{}_{}'.format(d[1],d[2])
    if 'bc' in F:
        datasets.append('bc')
    if 'dem' in F:
        datasets.append('dem')
    datasets = '_'.join(datasets)
    all_metrics[F] = {'TRAIN':{},'TEST':{}}
    with open(F) as o:
        for line in o:
            line = line.strip()
            if 'TRAIN' in line or 'TEST' in line:
                tt,fold,metric = line.split(',')[0:3]
                line = '{},{},{}'.format(line,name,datasets)
                if fold not in all_metrics[F][tt]:
                    all_metrics[F][tt][fold] = {metric:line}
                else:
                    all_metrics[F][tt][fold][metric] = line

for F in all_metrics:
    for tt in ['TRAIN','TEST']:
        for fold in all_metrics[F][tt]:
            for metric in all_metrics[F][tt][fold]:
                print(all_metrics[F][tt][fold][metric])
