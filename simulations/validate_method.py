from glob import glob
import pickle
import os
import gc

from networks import *
from data_simulation import *
from metrics import *

exec_path = ' '.join([
    "bsub -n {CPU} -M 8000 -o /dev/null -e /dev/null",
    "-J EX_NET_{NF}_{NVC}_{NC}_{NVCM}",
    "python3 simulation.py --n_features {NF}",
    "--n_virtual_cells {NVC} --n_classes {NC} --n_cpus {CPU}",
    "--learning_rate 0.002 --n_replicates 10",
    "--number_of_steps 2500 --n_virtual_cells_model {NVCM}"
])

N_VAL_ITERATIONS = 100
BATCH_SIZE = 100
N_ITERATIONS = 20

def read_pickle(pickle_path):
    return pickle.load(open(pickle_path,'rb'))

def get_info(simulation_path):
        nf,nvc,nvcm,nc,mvn,sv = simulation_path.split(os.sep)[-1][:-4].split('_')
        nf = int(nf[2:])
        nvc = int(nvc[3:])
        nvcm = int(nvcm[4:])
        nc = int(nc[2:])
        mvn = bool(mvn[3:])
        sv = bool(sv[2:])
        return nf,nvc,nvcm,nc,mvn,sv

def print_simulation_validation(simulation_path,
                                n_iterations,
                                batch_size):
    metric_machine = MetricFunction(
        {
            'Precision':precision,
            'Recall':recall,
            'Accuracy':accuracy,
            'F1-score':f1_score,
            'Confusion Matrix':confusion_matrix
        }
    )
    simulation_info = get_info(simulation_path)
    nf,nvc,nvcm,nc,mvn,sv = simulation_info
    try:
        simulation_data = read_pickle(simulation_path)
        model = VirtualCellClassifier(nf,nvcm,nc)
        n_replicates = len(simulation_data)
        #n_replicates = 1
        for replicate in range(n_replicates):
            gc.collect()
            metric_machine.reset()
            CENTROIDS_COS_DIS = {i:[] for i in range(nc)}
            N_COS_DIS = {i:0 for i in range(nc)}
            model.load_state_dict(simulation_data[replicate]['model'])
            generator = simulation_data[replicate]['data_generator']
            for iteration in range(n_iterations):
                batch,truth = generator.generate_batch(batch_size)
                batch = np.expand_dims(batch,1)
                batch = torch.Tensor(batch)
                truth = torch.Tensor(truth).to(torch.int64)
                output = model(batch)
                if nvcm == nvc:
                    proportions = model.get_virtual_cells(batch).detach().numpy()
                    for i,C in enumerate(generator.class_proportions):
                        tmp = proportions[truth == i,:].sum(axis=0)
                        CENTROIDS_COS_DIS[i].append(tmp)
                        N_COS_DIS[i] += proportions[truth==i,:].shape[0]
                if nc <= 2:
                    binary_output = torch.argmax(output,axis=-1)
                else:
                    binary_output = torch.round(output)
                    truth = torch.nn.functional.one_hot(truth,nc)

                metric_machine.update(binary_output,truth)
            if nvcm == nvc:
                for i,C in enumerate(generator.class_proportions):
                    centroid = np.divide(
                        np.concatenate(CENTROIDS_COS_DIS[i],axis=0).sum(axis=0),
                        N_COS_DIS[i]
                    )
                    num = np.sum(centroid*C)
                    den = np.multiply(
                        np.sum(np.sqrt(centroid**2)),
                        np.sum(np.sqrt(C**2))
                    )
                    cosine_distance = num / den
                    output_string = ",".join([str(xx) for xx in [
                        nf,nvc,nvcm,nc,mvn,sv,
                        replicate,'CosineDistance',i,cosine_distance
                    ]])
                    print(output_string)

            M = metric_machine.compute()
            for m in ['Precision','Recall','F1-score']:
                for i in range(len(M[m])):
                    output_string = ",".join([str(xx) for xx in [
                        nf,nvc,nvcm,nc,mvn,sv,
                        replicate,m,i,float(M[m][i])
                        ]])
                    print(output_string)
            output_string = ",".join([str(xx) for xx in [
                nf,nvc,nvcm,nc,mvn,sv,
                replicate,'Accuracy',"ALL",float(M['Accuracy'])
                ]])
            print(output_string)
    except Exception:
        tmp = exec_path.format(NF=nf,NVC=nvc,NVCM=nvcm,NC=nc,CPU=8)
        if mvn == True:
            tmp += ' --multivariate_normal'
        if sv == False:
            tmp += ' --no_shared_variance'
        print(tmp)

def get_all_simulations(folder_path,
                        n_iterations,
                        batch_size):
    output_string = ",".join([str(xx) for xx in [
        "NF","NVC","NVCM","NC","MVN","SV",
        "REPLICATE","METRIC","CLASS","METRIC_VALUE"
    ]])
    print(output_string)
    all_simulations = glob(folder_path + '/*.pkl')
    for simulation_path in all_simulations:
        gc.collect()
        print_simulation_validation(simulation_path)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Test virtual cell classifier.')

    parser.add_argument('--simulation_path',dest='simulation_path',
                        action='store',
                        type=str,
                        default=None)
    parser.add_argument('--batch_size',dest='batch_size',
                        action='store',
                        type=int,
                        default=100)
    parser.add_argument('--n_iterations',dest='n_iterations',
                        action='store',
                        type=int,
                        default=20)
    args = parser.parse_args()

    if args.simulation_path:
        print_simulation_validation(args.simulation_path,
                                    args.n_iterations,
                                    args.batch_size)
    else:
        get_all_simulations('simulation_outputs',
                            args.n_iterations,
                            args.batch_size)
