source /hps/nobackup/research/gerstung/josegcpa/virtual_envs/tf_gpu_cuda9/bin/activate

SLIDE_PATH=$1
OUTPUT_PATH=$2
UNET_CKPT=$3
QC_CKPT=$4

mkdir -p logs

snakemake \
    --config slide_path=$SLIDE_PATH output_path=$OUTPUT_PATH unet_ckpt=$UNET_CKPT qc_ckpt=$QC_CKPT \
    --latency-wait 3000 \
    --jobs 2 \
    -k \
    -p \
     --cluster 'bsub -e logs/{params.log_id}.e -o logs/{params.log_id}.o -M 16000 -n 16' \
    --rerun-incomplete $5
