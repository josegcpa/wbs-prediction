DATASET=${1:-MLL}
GPU=${2:-yes}

source ./config

if [[ $DATASET == 'MLL' ]]
then
    DIR=$MLL_DIR
    FMT=$MLL_FMT
    OUT=$MLL_OUTPUT_DIR
    RS=$RESCALE_FACTOR_MLL
elif [[ $DATASET == 'Adden1' ]]
then
    DIR=$ADDEN1_DIR
    FMT=$ADDEN1_FMT
    OUT=$ADDEN1_OUTPUT_DIR
    RS=$RESCALE_FACTOR_ADDEN1
elif [[ $DATASET == 'Adden2' ]]
then
    DIR=$ADDEN2_DIR
    FMT=$ADDEN2_FMT
    OUT=$ADDEN2_OUTPUT_DIR
    RS=$RESCALE_FACTOR_ADDEN2
else
    echo "error: first argument $DATASET should be one of MLL, Adden1 or Adden2"
    exit 1
fi

PIPELINE_SCRIPT=pipeline.sh

for slide_path in $DIR/*$FMT
do
    EXEC_STR="$PIPELINE_SCRIPT -i $slide_path -o $OUT -q $QC_CKPT -x $XGB_PATH -u $UNET_CKPT -r $RS"
    R=$(basename $slide_path)
    if [[ $GPU == "yes" ]]
    then
        bsub \
            -o logs/$R.o -e logs/$R.e \
            -M 8000 -q gpu \
            -gpu "num=1:j_exclusive=yes" \
            -E 'if [[ $(nvidia-smi | grep "No runnin" | wc -l) == 1 ]]; then X=0; else X=1; fi; echo $X' \
            -J SLIDE_$R \
            -g /SLIDE/GPU \
            -W 24:00 \
            sh $EXEC_STR "${@:3}"
    elif [[ $GPU == "smart" ]]
    then
        slide_id=$(basename $slide_path | cut -d '.' -f 1)
        if [[ -f $OUT/_checkpoints/"$slide_id"_seg ]]
        then
            bsub \
                -o logs/$R.o -e logs/$R.e \
                -M 8000 -n 8 \
                -J SLIDE_$R \
                -g /SLIDE/CPU \
                -W 24:00 \
                sh $EXEC_STR "${@:3}"
        else
            bsub \
                -o logs/$R.o -e logs/$R.e \
                -M 8000 -q gpu \
                -gpu "num=1:j_exclusive=yes" \
                -E 'if [[ $(nvidia-smi | grep "No runnin" | wc -l) == 1 ]]; then X=0; else X=1; fi; echo $X' \
                -J SLIDE_$R \
                -g /SLIDE/GPU \
                -W 24:00 \
                sh $EXEC_STR "${@:3}"
        fi
    elif [[ $GPU == "no" ]]
    then
        bsub \
            -o logs/$R.o -e logs/$R.e \
            -M 8000 -n 8 \
            -J SLIDE_$R \
            -g /SLIDE/CPU \
            -W 24:00 \
            sh $EXEC_STR "${@:3}"
    fi
done
