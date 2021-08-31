source ./config

out_path=$ADDEN_1_DIR

for file in $(find $ADDEN_1_DIR -name *$ADDEN_1_FMT)
do
    slide_name=$(basename $file)
    log_file=logs/salim_$slide_name.o
    log_file_err=logs/salim_$slide_name.e
    job_name=SALIM_$(basename $file | rev | cut -d '.' -f 2- | rev)
    bsub\
        -P gpu\
        -M 16000\
        -gpu "num=1:j_exclusive=yes"\
        -g /salim_gpu\
        -J $job_name\
        -o logs/salim_$(basename $slide_name).o\
        -e logs/salim_$(basename $slide_name).e\
        -W 24:00\
        sh pipeline.sh $file $out_path $UNET_CKPT $QC_CKPT
done
