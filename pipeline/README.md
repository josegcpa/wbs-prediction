# Blood cell segmentation pipeline

Implemented in a combination of `python` and `bash`, with workflow management in `snakemake`.

## Setup

Create a virtual environment using Python3.6 and install the packages available in `requirements.txt` (`pip3 install requirements.txt`). CUDA 9 (and a CUDA 9-supporting GPU) should be available. If you wish to run this exclusively on the CPU please install the packages in `requirements.txt` and uninstall `tensorflow-gpu` (`pip3 uninstall tensorflow-gpu`) and install `tensorflow` (`pip3 install tensorflow==1.12`)
