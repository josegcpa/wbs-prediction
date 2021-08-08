out_path=/hps/nobackup/research/gerstung/josegcpa/data/ADDEN_SVS/

for file in $(find /hps/nobackup/research/gerstung/josegcpa/data/ADDEN_SVS/ -name *svs)
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
        -W 48:00\
        sh pipeline.sh $file $out_path
done
