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
    try:
        nvc = int(d[-2].split('_')[0])
    except:
        nvc = int(d[-3])
    with open(F) as o:
        for line in o:
            line = line.strip()
            if 'PROBS_CLASS' in line:
                tt,fold,metric = line.split(',')[0:3]
                task = line.split(',')[-1]
                line = '{},{},{},{}'.format(line,name,nvc,datasets)
                print(line)
