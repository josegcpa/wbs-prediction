source ./config

out_path=$ADDEN_2_DIR

for file in $(find $ADDEN_2_DIR -name *$ADDEN_2_FMT)
do
    slide_name=$(basename $file)
    log_file=logs/salim_$slide_name.o
    log_file_err=logs/salim_$slide_name.e
    bsub\
        -P gpu\
        -M 16000\
        -gpu "num=1:j_exclusive=yes"\
        -g /salim_gpu\
        -J SALIM_$(basename $file | cut -d '.' -f 1)\
        -o logs/salim_$(basename $slide_name).o\
        -e logs/salim_$(basename $slide_name).e\
        -W 24:00\
        sh pipeline.sh $file $out_path $UNET_CKPT_ADDEN_2 $QC_CKPT
done
