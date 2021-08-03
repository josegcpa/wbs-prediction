print('-'*20)
import numpy as np

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

C = [0,1,2,3,4,5]
P = np.array([0.5,0.2,0.2,0.1,0.3,0.05])
P = P / P.sum()
ids = np.arange(300)
classes = np.random.choice(C,size=500,replace=True,
                           p=P)

train,test = split_dataset_stratified(ids,classes)
train_set = [classes[i] for i in train]
test_set = [classes[i] for i in test]

Tr,Tr_c = np.unique(train_set,return_counts = True)
Te,Te_c = np.unique(test_set,return_counts = True)
full,full_c = np.unique(classes,return_counts = True)

for i in C:
    if i in Tr:
        a = Tr_c[i]
    else:
        a = 0
    if i in Te:
        b = Te_c[i]
    else:
        b = 0
    print('proportion for {}: {} (train); {} (test); {} (true)'.format(
        i,a/np.sum(Tr_c),b/np.sum(Te_c),full_c[i] / np.sum(full_c)))
