source /hps/nobackup/research/gerstung/josegcpa/virtual_envs/tf_gpu_cuda9/bin/activate

SLIDE_PATH=$1

mkdir -p logs

snakemake \
    --config slide_path=$1 output_path=/hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF \
    --latency-wait 3000 \
    --jobs 2 \
    -k \
     --cluster 'bsub -e logs/{params.log_id}.e -o logs/{params.log_id}.o -M 16000 -n 16' \
    --rerun-incomplete
