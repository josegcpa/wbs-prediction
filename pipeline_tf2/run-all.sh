DATASET=${1:-MLL}
FAN=${2:-yes}

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
elif [[ $FAN == "no" ]]
then
    PIPELINE_SCRIPT=pipeline.sh
else
    echo "second argument $FAN should be either yes or no"
fi

for slide_path in $DIR/*$FMT
do
    R=$(basename $slide_path)
    sh $PIPELINE_SCRIPT $slide_path $OUT $QC_CKPT $FAN_CKPT $UNET_CKPT $DATASET $RS "${@:3}"
done
