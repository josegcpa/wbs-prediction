source ./config

N=1000

mkdir -p examples

for file in $ADDEN_2_DIR/*$ADDEN_2_FMT
do
    root=$(basename $file | cut -d '.' -f 1)
    n_checkpoints=$(find $ADDEN_2_DIR_OUT/_checkpoints -name "$root*" | wc -l)
    if [[ $n_checkpoints -eq 5 ]]
    then
        echo "Retrieving examples for $root..."
        segmented_file_wbc=$ADDEN_2_DIR_OUT/_segmented_wbc/$root.h5
        segmented_file_rbc=$ADDEN_2_DIR_OUT/_segmented_rbc/$root.h5
        aggregate_file_wbc=$ADDEN_2_DIR_OUT/_aggregates_wbc/$root.h5
        aggregate_file_rbc=$ADDEN_2_DIR_OUT/_aggregates_rbc/$root.h5
        if [[ ! -f examples/"$root"_wbc.h5 ]]
        then
            bsub -n 2 -M 2000 -o /dev/null -e /dev/null \
                python3 scripts/python/get_examples.py\
                --slide_path $file\
                --aggregates_path $aggregate_file_wbc\
                --segmented_path $segmented_file_wbc\
                --output_path examples/"$root"_wbc.h5\
                --subset $N
        fi

        if [[ ! -f examples/"$root"_rbc.h5 ]]
        then
            bsub -n 2 -M 2000 -o /dev/null -e /dev/null \
                python3 scripts/python/get_examples.py\
                --slide_path $file\
                --aggregates_path $aggregate_file_rbc\
                --segmented_path $segmented_file_rbc\
                --output_path examples/"$root"_rbc.h5\
                --subset $N
        fi
    fi
done

for file in $MLL_DIR/*$MLL_FMT
do
    root=$(basename $file | cut -d '.' -f 1)
    n_checkpoints=$(find $MLL_DIR/_checkpoints -name "$root*" | wc -l)
    if [[ $n_checkpoints -eq 5 ]]
    then
        echo "Retrieving examples for $root..."
        segmented_file_wbc=$MLL_DIR/_segmented_wbc/$root.h5
        segmented_file_rbc=$MLL_DIR/_segmented_rbc/$root.h5
        aggregate_file_wbc=$MLL_DIR/_aggregates_wbc/$root.h5
        aggregate_file_rbc=$MLL_DIR/_aggregates_rbc/$root.h5
        if [[ ! -f examples/"$root"_wbc.h5 ]]
        then
            bsub -n 2 -M 2000 -o /dev/null -e /dev/null \
                python3 scripts/python/get_examples.py\
                --slide_path $file\
                --aggregates_path $aggregate_file_wbc\
                --segmented_path $segmented_file_wbc\
                --output_path examples/"$root"_wbc.h5\
                --subset $N
        fi

        if [[ ! -f examples/"$root"_rbc.h5 ]]
        then
            bsub -n 2 -M 2000 -o /dev/null -e /dev/null \
                python3 scripts/python/get_examples.py\
                --slide_path $file\
                --aggregates_path $aggregate_file_rbc\
                --segmented_path $segmented_file_rbc\
                --output_path examples/"$root"_rbc.h5\
                --subset $N
        fi
    fi
done
