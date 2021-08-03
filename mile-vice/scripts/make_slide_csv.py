import h5py
import numpy as np

all_ids = {}
all_datasets = {}
all_datasets["wbc"] = {}
all_ids['wbc'] = []
with h5py.File('datasets/wbc.h5','r') as x:
    for ID in x:
        if ID not in all_ids['wbc']:
            all_ids['wbc'].append(ID)
        all_datasets['wbc'][ID] = np.concatenate([x[ID]['mean'],x[ID]['variance']])

all_datasets["rbc"] = {}
all_ids['rbc'] = []
with h5py.File('datasets/rbc.h5','r') as x:
    for ID in x:
        if ID not in all_ids['rbc']:
            all_ids['rbc'].append(ID)
        all_datasets['rbc'][ID] = np.concatenate([x[ID]['mean'],x[ID]['variance']])

all_datasets["blood_counts"] = {}
all_ids['blood_counts'] = []
with open("datasets/blood_counts_no_header.csv",'r') as o:
    for line in o:
        d = line.strip().split(',')
        if d[0] not in all_ids['blood_counts']:
            all_ids['blood_counts'].append(d[0])
        all_datasets['blood_counts'][d[0]] = d[1:]

all_datasets["demographics"] = {}
all_ids['demographics'] = []
with open("datasets/demographics_no_header.csv",'r') as o:
    for line in o:
        d = line.strip().split(',')
        if d[0] not in all_ids['demographics']:
            all_ids['demographics'].append(d[0])
        all_datasets['demographics'][d[0]] = d[1:]

all_datasets["class"] = {}
all_ids['class'] = []
with open("labels/all_classes.csv") as o:
    for line in o:
        d = line.strip().split(',')
        if d[0] not in all_ids['class']:
            all_ids['class'].append(d[0])
        all_datasets["class"][d[0]] = d[1:]

common_ids = list(set.intersection(*[set(all_ids[x]) for x in all_ids]))

for ID in common_ids:
    data = ','.join([','.join([str(x) for x in all_datasets['wbc'][ID]]),
                     ','.join([str(x) for x in all_datasets['rbc'][ID]]),
                     ','.join([str(x) for x in all_datasets['blood_counts'][ID]]),
                     ','.join([str(x) for x in all_datasets['demographics'][ID]]),
                     ','.join([str(x) for x in all_datasets['class'][ID]])])
    data = ID + ',' + data
    print(data)


