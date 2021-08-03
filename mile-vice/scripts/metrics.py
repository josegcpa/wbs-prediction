import torch
from sklearn import metrics

def precision(metric_dict):
    tp = metric_dict['tp'].sum(axis=0).float()
    fp = metric_dict['fp'].sum(axis=0).float()
    return tp / (tp+fp)

def deviance(metric_dict):
    pred = metric_dict['probs']
    truth = metric_dict['y']
    return 2*metrics.log_loss(truth,pred,normalize=False)

def auc(metric_dict):
    pred = metric_dict['probs']
    truth = metric_dict['y']
    return metrics.roc_auc_score(truth,pred,average='weighted')

def recall(metric_dict):
    tp = metric_dict['tp'].sum(axis=0).float()
    fn = metric_dict['fn'].sum(axis=0).float()
    return tp / (tp+fn)

def confusion_matrix(metric_dict):
    if metric_dict['y_hat'].shape[1] > 1:
        pred = torch.argmax(metric_dict['y_hat'],dim=1).numpy()
        truth = torch.argmax(metric_dict['y'],dim=1).numpy()
    else:
        pred = metric_dict['y_hat'].numpy()
        truth = metric_dict['y'].numpy()
    return metrics.confusion_matrix(truth,pred)

def dice_coefficient(metric_dict):
    tp = metric_dict['tp'].sum(axis=0).float()
    fp = metric_dict['fp'].sum(axis=0).float()
    fn = metric_dict['fn'].sum(axis=0).float()
    return tp / (tp+fp+fn)

def accuracy(metric_dict):
    cond_1 = len(metric_dict['tp'].shape) == 1

    try:
        cond_2 = metric_dict['tp'].shape[1] == 1
    except:
        cond_2 = False

    if cond_1 or cond_2:
        tp = metric_dict['tp'].sum().float()
        tn = metric_dict['tn'].sum().float()
        fp = metric_dict['fp'].sum().float()
        fn = metric_dict['fn'].sum().float()
        return (tp+tn)/(tp+tn+fp+fn)
    else:
        tp = metric_dict['tp'].sum().float()
        fp = metric_dict['fp'].sum().float()
        fn = metric_dict['fn'].sum().float()

        return (tp)/(tp+fp+fn)

def f1_score(metric_dict):
    prec = precision(metric_dict).float()
    rec = recall(metric_dict).float()
    return 2 / (1/prec + 1/rec)

class MetricFunction:
    def __init__(self,metric_fn_dict):
        self.tp = []
        self.tn = []
        self.fp = []
        self.fn = []
        self.y_hat = []
        self.y = []
        self.probs = []
        self.metric_fn_dict = metric_fn_dict
        self.metric_store = {
            x:[] for x in self.metric_fn_dict
        }

    def update(self,y,y_hat,probs=None):
        tp = torch.logical_and(y==1,y==y_hat)
        tn = torch.logical_and(y==0,y==y_hat)
        fp = torch.logical_and(y==0,y!=y_hat)
        fn = torch.logical_and(y==1,y!=y_hat)
        self.tp.append(tp)
        self.tn.append(tn)
        self.fp.append(fp)
        self.fn.append(fn)
        if probs is not None:
            self.probs.append(probs)
        else:
            self.probs.append(y_hat)
        self.y_hat.append(y_hat)
        self.y.append(y)

    def compute(self):
        output_dict = {}
        tp = torch.cat(self.tp,axis=0)
        tn = torch.cat(self.tn,axis=0)
        fp = torch.cat(self.fp,axis=0)
        fn = torch.cat(self.fn,axis=0)
        y = torch.cat(self.y,axis=0)
        y_hat = torch.cat(self.y_hat,axis=0)
        probs = torch.cat(self.probs,axis=0)
        for metric_fn_key in self.metric_fn_dict:
            metric_fn = self.metric_fn_dict[metric_fn_key]
            output_dict[metric_fn_key] = metric_fn({
                'tp':tp,
                'tn':tn,
                'fp':fp,
                'fn':fn,
                'y':y,
                'y_hat':y_hat,
                'probs':probs
            })
            self.metric_store[metric_fn_key].append(output_dict[metric_fn_key])

        return output_dict

    def reset(self):
        self.tp = []
        self.tn = []
        self.fp = []
        self.fn = []
        self.y = []
        self.y_hat = []
        self.probs = []
