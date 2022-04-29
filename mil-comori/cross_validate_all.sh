CONFIG_FILE=$1

mkdir -p logs
mkdir -p models
mkdir -p best_models
mkdir -p metrics
mkdir -p logs_cluster

snakemake --latency-wait 30000 \
    -k \
    -p \
    --rerun-incomplete \
    --configfile $CONFIG_FILE -c 1
