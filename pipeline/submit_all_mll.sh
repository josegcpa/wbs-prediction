source ./config

out_path=$MLL_DIR

for file in $(find $MLL_DIR -name *$MLL_FMT)
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
    sh pipeline.sh $file $out_path $UNET_CKPT $QC_CKPT
done
