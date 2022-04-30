# Code repository for "Computational analysis of peripheral blood smears detects disease-associated cytomorphologies"

## Motivation

This is the repository for "Computational analysis of peripheral blood smears detects disease-associated cytomorphologies" ([available as a preprint on MedRxiv](https://www.medrxiv.org/content/10.1101/2022.04.19.22273757v1)). In this work, we use the whole blood slides of >300 individuals with myelodyplastic syndromes and anaemias and use them to develop a method that is capable of automatically detecting cells in whole blood slides, predicting a disease and retrieving examples of cells which are relevant for each classification.

## Code map

### Requirements

#### Software

* `Python`
* `Snakemake`
* `R` (analysis and plotting)
* `conda` (recommended for virtual environment management)

#### Required python packages

Each folder requires different packages and a standard environment manager is advised (be it `conda` or `virtualenv`). Within each folder, a specification for the required packages is provided.

#### Required R packages

In each folder a short description is provided delineating the required packages to run each analysis.

### Folder enumeration

1. [`pipeline_tf2`](https://github.com/josegcpa/wbs-prediction/tree/main/pipeline_tf2) - contains Haemorasis, the pipeline for WBC and RBC detection and characterisation from WBS. Also available as a Docker container in [Docker hub](https://hub.docker.com/repository/docker/josegcpa/blood-cell-detection)
2. [`mil-comori`](https://github.com/josegcpa/wbs-prediction/tree/main/mil-comori) - contains the code to train and run Morphotype analysis on the output from [Haemorasis](https://github.com/josegcpa/haemorasis)
3. [`analysis-plotting`](https://github.com/josegcpa/wbs-prediction/tree/main/analysis-plotting) - contains the code to analyse and plot the results from the previous processes

## Running the analyses 

It should be noted that, due to the high data volume required for this work, we have chose to make the necessary inputs and outputs for `mil-comori`, as well as the necessary inputs for `analysis-plotting`, available through a [Figshare project](https://figshare.com/projects/Computational_analysis_of_peripheral_blood_smears_detects_disease-associated_cytomorphologies/132443). Instructions detailing how and what to download are provided in README files available in each folder.

## Running the entire analysis (morphotype analysis + analysis-plotting)

1. 

## Running only the `analysis-plotting` scripts

This requires only a few steps, particularly:

1. Downloading and unzipping the Morphotype analysis output inside of the `mil-comori` folder ([download link](https://figshare.com/account/projects/132443/articles/19369391))
2. Downloading and unzipping the `data` and `datasets` required for this analysis inside the `analysis-plotting` folder ([download-link](https://figshare.com/account/projects/132443/articles/19371008))

Inside the `analysis-plotting` folder, a `Snakefile` is available that will execute all of the scripts and handle dependencies (provided all necessary packages are installed).