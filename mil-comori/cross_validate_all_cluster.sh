CONFIG_FILE=$1
EXTRA_ARG=$2

mkdir -p logs
mkdir -p models
mkdir -p best_models
mkdir -p metrics
mkdir -p logs_cluster

snakemake --latency-wait 30000 \
    -k \
    -p \
    --rerun-incomplete \
    --configfile $CONFIG_FILE -c 1 \
    --cluster 'bsub -J MILE_VICE_{params.log_id} -e logs_cluster/{params.log_id}.e -o logs_cluster/{params.log_id}.o -M 4000 -n 16' \
    --jobs 10000 $EXTRA_ARG
