import sys
import argparse
import os
import numpy as np
import time
import openslide
import h5py
import pickle
import xgboost
from sklearn.preprocessing import StandardScaler

import MIA

MIA_FEATURES = MIA.MIA_FEATURES

parser = argparse.ArgumentParser(
    prog = 'aggregate_hdf5.py',
    description = 'creates dataset and summary statistics for hdf5 file'
)

parser.add_argument('--hdf5_path',dest='hdf5_path',
                    action='store',
                    default=None,
                    help='Path to hdf5 file.')
parser.add_argument('--output_path',dest='output_path',
                    action='store',
                    default=None,
                    help='Output path for the hdf5 file.')

args = parser.parse_args()

standardized_scaler = StandardScaler()
fitted_model = xgboost.XGBClassifier(n_estimators=50,eval_metric="logloss")
with open("scripts/python/rbc_scaler_params","rb") as o:
    sc_params = pickle.load(o)


standardized_scaler.n_features_in_ = len(sc_params['mean'])
standardized_scaler.mean_ = sc_params['mean']
standardized_scaler.scale_ = sc_params['std']
standardized_scaler.var_ = sc_params['std']**2
fitted_model.load_model("scripts/python/rbc_xgb_model")

F = h5py.File(args.hdf5_path,mode='r')
F_out = h5py.File(args.output_path,mode='w')
F_out_cells = F_out.create_group('cells')
N_accumulator = 0
sum_accumulator = np.zeros([len(MIA_FEATURES)],dtype=np.float64)
squared_sum_accumulator = np.zeros([len(MIA_FEATURES)],dtype=np.float64)
cells = []
cell_centers = []
pre_cells = []
pre_cell_centers = []
N_groups = 0
times = []
for item in F:
    cell = F[item]['features']
    features = []
    for i,feature in enumerate(MIA_FEATURES):
        features.append(cell[i])
    features = np.array(features)
    if np.any(np.isnan(features)):
        pass
    elif np.any(np.isinf(features)):
        pass
    else:
        pre_cells.append(features)
        pre_cell_centers.append(item)

    if len(pre_cells) % 1000 == 0:
        N_accumulator += len(pre_cells)
        pre_cells = np.array(pre_cells)
        pre_cells_ = standardized_scaler.transform(pre_cells)
        predictions = fitted_model.predict(pre_cells_)
        pre_cells = pre_cells[predictions,:]
        pre_cell_centers = [
            pre_cell_centers[idx]
            for idx in range(len(pre_cell_centers))
            if predictions[idx] == True]
        sum_accumulator += pre_cells.sum(axis=0)
        squared_sum_accumulator += np.square(pre_cells).sum(axis=0)
        for idx in range(pre_cells.shape[0]):
            cells.append(pre_cells[idx,:])
            cell_centers.append(pre_cell_centers[idx])
        pre_cells = []
        pre_cell_centers = []

    if len(cells) >= 100000:
        cells = np.stack(cells,axis=0)
        F_out_cells.create_dataset(
            str(N_groups),cells.shape,
            dtype=np.float64,data=cells)
        cells = []
        N_groups += 1

if len(pre_cells) > 0:
    N_accumulator += len(pre_cells)
    pre_cells = np.array(pre_cells)
    pre_cells_ = standardized_scaler.transform(pre_cells)
    predictions = fitted_model.predict(pre_cells_)
    pre_cells = pre_cells[predictions,:]
    pre_cell_centers = [
        pre_cell_centers[idx]
        for idx in range(len(pre_cell_centers))
        if predictions[idx] == True]
    sum_accumulator += pre_cells.sum(axis=0)
    squared_sum_accumulator += np.square(pre_cells).sum(axis=0)
    for idx in range(pre_cells.shape[0]):
        cells.append(pre_cells[idx,:])
        cell_centers.append(pre_cell_centers[idx])
    pre_cells = []
    pre_cell_centers = []

print(len(cells))

if len(cells) > 0:
    cells = np.array(cells)
    F_out_cells.create_dataset(
        str(N_groups),cells.shape,
        dtype=np.float64,data=cells)
    cells = []
    N_groups += 1

means = sum_accumulator / N_accumulator
variances = squared_sum_accumulator / N_accumulator - means

F_out.create_dataset(
    'mean',means.shape,
    dtype=np.float64,data=means)
F_out.create_dataset(
    'variance',variances.shape,
    dtype=np.float64,data=variances)
F_out.create_dataset(
    'cell_centers',
    data=np.array(cell_centers,dtype="S"))

F.close()
F_out.close()
