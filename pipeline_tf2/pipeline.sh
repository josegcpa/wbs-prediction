source /hps/research/gerstung/josegcpa/projects/01IMAGE/tf2/bin/activate

SLIDE_PATH=$1
OUTPUT_PATH=$2
QC_CKPT=$3
FAN_CKPT=$4
UNET_CKPT=$5
DATASET=$6
RF=$7

mkdir -p logs

snakemake \
    --config slide_path=$SLIDE_PATH output_path=$OUTPUT_PATH fan_ckpt=$FAN_CKPT unet_ckpt=$UNET_CKPT qc_ckpt=$QC_CKPT dataset=$DATASET rescale_factor=$RF\
    -s Snakefile \
    --latency-wait 3000 \
    --jobs 3 \
    -k \
    -p \
    --rerun-incomplete $8
     #--cluster 'bsub -q research-rh74 -e logs/{params.log_id}.e -o logs/{params.log_id}.o -M {params.n_cores} -n {params.mem}' \
