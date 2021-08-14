# A complete computational assessment of the cytomorphological determinants of myelodyplastic syndromes

## Motivation

This is the repository for [PLACEHOLDER](). In this work, we use the whole blood slides of >300 individuals with myelodyplastic syndromes and anaemias and use them to develop a method that is capable of predicting a disease and retrieving examples of cells which are relevant for each classification.

## Code map

### Requirements

#### Software

* `python`
* `snakemake`
* `R` (analysis and plotting)

#### Required python packages

`opencv-python`, `tensorflow==1.12`, `scikit-image`, `h5py`, `albumentations`, `psutil`, `pytorch`, `tifffile`

### Project enumeration

1. `pipeline` - contains the pipeline for WBC and RBC detection and characterisation from WBS
2. `simulations` - contains simulations validating MILe-ViCe
3. `mile-vice` - contains the code to train and run MILe-ViCe on the output from `pipeline`
4. `rbc-segmentation` - contains the code to train a predictor that filters poorly predictions for detected RBC
5. (STILL TESTING) `vae-characterisation` - characterisation of blood cells using a beta-variational autoencoder
