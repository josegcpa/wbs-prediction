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

It should be noted that, due to the high data volume required for this work, we have chose to make the necessary inputs and outputs for `mil-comori`, as well as the necessary inputs for `analysis-plotting`, available through a [Figshare project](https://figshare.com/projects/Computational_analysis_of_peripheral_blood_smears_detects_disease-associated_cytomorphologies/132443). We provide simplified instructions on how to reproduce the results from the Morhpotype analysis and the downstream statistical analysis below. In any case, README files explaining how to run each analysis are contained in each folder.

This was developed and tested using Python 3.6.8 and on 8GB RAM on CentOS Linux 8 (kernel: Linux 4.18.0-240.22.1.el8_3.x86_64), and for Haemorasis, the blood cell detection pipeline (`pipeline_tf2`), we used NVidia Quadro M6000 GPUs. Approximately 3GB of free disk space are required to run both the Morphotype analysis and the downstream statistical analysis of the results (if you want to run only the latter approximately only 200MB of disk space are required). We also recommend having at least 8GB of RAM and a CentOS-based machine, and using `conda` to manage packages as this greatly simplifies dependencies.

### Simplified instructions to run Morphotype analysis and analysis-plotting

1. Clone this Github repository to your local machine 
    * `git clone https://github.com/josegcpa/wbs-prediction.git`
2. Download and unzip the necessary data from Figshare using the download script
    * `sh download-and-unzip-data-full.sh`
3. Enter the Morphotype analysis directory (`cd mil-comori`) and follow the instructions presented there in the README.md file
4. After having performed the analysis in step 3., enter the statistical analysis and figure generation folder (`cd analysis-plotting`) and follow the instructions presented there

### Simplified instructions to run only the statistical analysis and figure generation scripts

1. Clone this Github repository to your local machine 
    * `git clone https://github.com/josegcpa/wbs-prediction.git`
2. Download and unzip the necessary data from Figshare using the download script containing the necessary analysis data and the Morphotype analysis output
    * `sh download-and-unzip-analysis-plotting-only.sh`
3. Enter the statistical analysis and figure generation folder (`cd analysis-plotting`) and follow the instructions presented there