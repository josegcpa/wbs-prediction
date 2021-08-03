bsub\
    -M 16000\
    -P gpu\
    -gpu "num=1:j_exclusive=yes"\
    -q research-rh74\
    -m "gpu-009 gpu-011"\
    -o logs/train-mse-wbc.o\
    -e logs/train-mse-wbc.e\
    -J TRAIN_VAE_WBC\
    python3 run.py --config configs/wbc.yaml

bsub\
    -M 16000\
    -P gpu\
    -gpu "num=1:j_exclusive=yes"\
    -q research-rh74\
    -m "gpu-009 gpu-011"\
    -o logs/train-mse-wbc-96px.o\
    -e logs/train-mse-wbc-96px.e\
    -J TRAIN_VAE_WBC_96\
    python3 run.py --config configs/wbc-96px.yaml

bsub\
    -M 16000\
    -P gpu\
    -gpu "num=1:j_exclusive=yes"\
    -q research-rh74\
    -m "gpu-009 gpu-011"\
    -o logs/train-mse-rbc.o\
    -e logs/train-mse-rbc.e\
    -J TRAIN_VAE_RBC\
    python3 run.py --config configs/rbc.yaml

