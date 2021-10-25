source ./config

N=500

mkdir -p examples_large

for file in $MLL_DIR/*$MLL_FMT
do
    root=$(basename $file | cut -d '.' -f 1)
    echo "Retrieving examples for $root..."
    segmented_file_rbc=$MLL_DIR/_segmented_rbc/$root.h5
    aggregate_file_rbc=$MLL_DIR/_aggregates_rbc/$root.h5
    bsub -n 1 -M 1000 -o /dev/null -e /dev/null \
        -J EXAMPLES_WBC_$root \
        python3 scripts/python/get_examples.py\
        --slide_path $file\
        --aggregates_path $aggregate_file_rbc\
        --segmented_path $segmented_file_rbc\
        --output_path examples_large/"$root"_rbc.h5\
        --subset $N\
        --output_size 256\
        --box
done