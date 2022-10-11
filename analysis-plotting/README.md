# Statistical analysis for "Computational analysis of peripheral blood smears detects disease-associated cytomorphologies"

This folder contains the `R` scripts necessary to conduct the analysis of the results from "Computational cytomorphology reveals morphotypes in whole blood slides".

* Scripts starting with `A-` correspond to a more general analysis of the cohort and of the blood cell detection pipeline metrics.
* Scripts starting with `B-` correspond to the analysis of the discriminatory power of morphometric features in standard classification tasks
* Scripts starting with `C-` correspond to the analysis of single-cell morphometry and MILe-ViCe results

## Software and package details

The R version used is `3.6.1` and the installed packages, as well as their respective versions, are provided in the `packages.R` file and below:

```
P <- c("caret==6.0-83","cowplot==1.1.1","dendextend==1.15.2","dunn.test==1.3.5","ggplot2==3.3.5","ggpubr==0.4.0","ggrepel==0.9.1","ggsci==2.9","glmnet==2.0-16","MASS==7.3-51.3","MLmetrics==1.1.1","pROC==1.18.0","progress==1.2.2","RANN==2.6.1","RRF==1.9.1","tidyverse==1.3.1","umap==0.2.7.0","WRS2==1.1-3")
```

## Running these analyses

### Necessary downloads (in case these haven't been done as specified in the [wbs-prediction repository](https://github.com/josegcpa/wbs-prediction))

#### If you have run Morphotype analysis

1. Download and unzip the data available in https://doi.org/10.6084/m9.figshare.19371008.v1.

#### If you have not run Morphotype analysis

1. Download the data available in https://doi.org/10.6084/m9.figshare.19369391.v2 and unzip to the folder called `mil-comori` in the main project folder (`../mil-comori`)
2. Download and unzip the data available in https://doi.org/10.6084/m9.figshare.19371008.v1.

### Statistical analysis and figure generation

To help simplify matters, we have setup a `Snakefile` pipeline that coordinates all the necessary steps to generate figures. To run it, please run the command `snakemake -c 1`, where the flag `-c` specifies the number of available cores (and concurrent processes). Running this step with a single available core takes approximately 50 minutes.
