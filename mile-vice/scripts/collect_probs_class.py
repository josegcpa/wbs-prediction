import os
import argparse
from glob import glob

parser = argparse.ArgumentParser(description='Collect metrics.')
parser.add_argument('--log_files',dest='log_files',
                    nargs='+')
parser.add_argument('--subset',dest='subset',
                    action='store_true',default='full')
parser.add_argument('--multiclass',dest='multiclass',
                    action='store_true')
args = parser.parse_args()

all_csv_files = args.log_files

if args.multiclass == True:
    all_csv_files = [x for x in all_csv_files if 'multi_class' in x]
else:
    all_csv_files = [x for x in all_csv_files if 'multi_class' not in x]

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
    try:
        nvc = int(d[-1].split('_')[0])
    except:
        nvc = int(d[-2])

    with open(F) as o:
        for line in o:
            line = line.strip()
            if 'PROBS_CLASS' in line:
                tt,fold,metric = line.split(',')[0:3]
                task = line.split(',')[-1]
                line = '{},{},{},{},{}'.format(line,name,nvc,datasets,subset)
                print(line)
