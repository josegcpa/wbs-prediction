# Weakly supervised learning for virtual cell quantification

Uses a multi-instance learning approach to deduce virtual cells (or cellular archetypes) that can be relevant for specific classifications. The script responsible for most of the training is `scripts/train_stack.py`, which uses the models, metrics and data generators specified in `scripts/networks.py`, `scripts/metrics.py` and `scripts/data_generators.py`, respectively. 

We use the hierarchical format HDF5 for each one of the datasets. Particularly, we structure these files such that a file is composed of multiple groups, each of which corresponds to an array with $n$ rows (the number of detected cells) and $f$ columns (the features characterising each cell). Given that we are interested in multiple types of cells (RBC and WBC), we have two separate hdf5 files specified in such a way. The names of the groups in the HDF5 file (slide IDs) are important as this is what is used to match groups of cells between different files.

The script MILe-ViCe (through `train_stack.py`) allows multiple cell datasets to be provided (`dataset_path`), as well as additionaly tabular datasets (`other_dataset_path`), provided they are CSV files where the first row corresponds to the slide ID. If a variable is categorical (i.e. most elements are non-numerical) this will be converted to a one-hot encoded variable. MILe-ViCe also allows for multiple objectives through the specificaiton of multiple `labels_path`. More parameters can be accessed through `python3 train_stack.py --help`. An example is provided below:

```
    python3 scripts/train_stack.py \
        --n_virtual_cells 25 \
        --n_classes 2 \
        --dataset_path datasets/wbc.h5 \
        --dataset_path datasets/rbc.h5 \
        --labels_path labels/A.csv \
        --labels_path labels/B.csv \
        --number_of_steps 10000 \
        --batch_size 32 \
        --number_of_cells 500 \
        --learning_rate 0.01 \
        --n_folds 5 \
        --weight_decay 0.2 \
        --other_dataset_path datasets/other_data.csv \
        --model_path models/model \
        --excluded_ids A B C DE F \
        --median_impute \
        --min_cells 50 
```