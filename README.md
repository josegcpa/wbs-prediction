# Code repository for "Computational cytomorphology reveals morphotypes in whole blood slides"

## Motivation

This is the repository for "Computational cytomorphology reveals morphotypes in whole blood slides" ([available as a preprint on XXXXX]()). In this work, we use the whole blood slides of >300 individuals with myelodyplastic syndromes and anaemias and use them to develop a method that is capable of automatically detecting cells in whole blood slides, predicting a disease and retrieving examples of cells which are relevant for each classification.

## Code map

### Requirements

#### Software

* `Python`
* `Snakemake`
* `R` (analysis and plotting)
* `conda` (recommended for virtual environment management)

#### Required python packages

Each folder requires different packages and a standard environment manager is advised (be it `conda` or `virtualenv`). Within each folder, a specification for the required packages is provided.

### Project enumeration

1. `pipeline_tf2` - contains Haemorasis, the pipeline for WBC and RBC detection and characterisation from WBS. Also available as a Docker container in [Docker hub](https://hub.docker.com/repository/docker/josegcpa/blood-cell-detection)
2. `mile-vice` - contains the code to train and run MILe-ViCe on the output from [Haemorasis](https://github.com/josegcpa/haemorasis)
3. `analysis-plotting` - contains the code to analyse and plot the results from the previous processes

## Running only the `analysis-plotting` scripts

This requires only a few steps, particularly:

1. Downloading and unzipping the `mile-vice` analysis output inside of the `mile-vice` folder ([download link]())
2. Downloading and unzipping the `data` and `datasets` required for this analysis inside the `analysis-plotting` folder ([download-link]())

Inside the `analysis-plotting` folder, a `Snakefile` is available that will execute all of the scripts and handle dependencies (provided all necessary packages are installed)