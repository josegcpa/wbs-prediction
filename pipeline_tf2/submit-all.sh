DATASET=${1:-MLL}
FAN=${2:-yes}
GPU=${3:-yes}

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

if [[ $FAN == "yes" ]]
then
    PIPELINE_SCRIPT=pipeline-fan.sh
    UNET_CKPT=$UNET_CKPT_FAN
elif [[ $FAN == "no" ]]
then
    PIPELINE_SCRIPT=pipeline.sh
else
    echo "second argument $FAN should be either yes or no"
fi

for slide_path in $DIR/*$FMT
do
    R=$(basename $slide_path)
    if [[ $GPU == "yes" ]]
    then
        echo bsub \
            -o logs/$R.o -e logs/$R.e \
            -M 16000 -P gpu \
            -gpu "num=1:j_exclusive=yes" \
            -E 'if [[ $(nvidia-smi | grep "No runnin" | wc -l) == 1 ]]; then X=0; else X=1; fi; echo $X' \
            -J SLIDE_$R \
            -W 24:00 \
            sh $PIPELINE_SCRIPT $slide_path $OUT $QC_CKPT $FAN_CKPT $UNET_CKPT $DATASET $RS "${@:4}"
    else
        echo bsub \
            -o logs/$R.o -e logs/$R.e \
            -M 16000 -n 8 \
            -J SLIDE_$R \
            -W 24:00 \
            sh $PIPELINE_SCRIPT $slide_path $OUT $QC_CKPT $FAN_CKPT $UNET_CKPT $DATASET $RS "${@:4}"
    fi
done
