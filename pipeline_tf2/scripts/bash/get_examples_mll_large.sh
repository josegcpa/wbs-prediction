source ./config

N=500

output_dir=examples_large/mll

mkdir -p $output_dir

for file in $MLL_DIR/*$MLL_FMT
do
    root=$(basename $file | cut -d '.' -f 1)
    echo "Retrieving examples for $root..."
    segmented_file_wbc=$MLL_DIR/_segmented_wbc/$root.h5
    segmented_file_rbc=$MLL_DIR/_segmented_rbc/$root.h5
    aggregate_file_wbc=$MLL_DIR/_aggregates_wbc/$root.h5
    aggregate_file_rbc=$MLL_DIR/_aggregates_rbc/$root.h5
    bsub -n 1 -M 1000 -o /dev/null -e /dev/null \
        -J EXAMPLES_WBC_$root \
        python3 scripts/python/get_examples.py\
        --slide_path $file\
        --aggregates_path $aggregate_file_wbc\
        --segmented_path $segmented_file_wbc\
        --output_path $output_dir/"$root"_wbc.h5\
        --subset $N\
        --output_size 512\
        --random
done