N_CPUS=$1
for N_FEATURES in 20 50 100 200 500
do
    for N_VIRTUAL_CELLS in 2 5 7 10 15 20 25
    do
        for N_CLASSES in 2 4 5 10 20
        do
            if [ $N_CLASSES -lt $N_VIRTUAL_CELLS ]
            then
                for N_VIRTUAL_CELLS_MODEL in 5 10 20 50 100
                do
                    bsub\
                        -n $N_CPUS\
                        -M 8000\
                        -J EX_NET_"$N_FEATURES"_"$N_VIRTUAL_CELLS"_"$N_CLASSES"_"$N_VIRTUAL_CELLS_MODEL"\
                        -o /dev/null\
                        -e /dev/null\
                        python3 simulation.py\
                        --n_features $N_FEATURES\
                        --n_virtual_cells $N_VIRTUAL_CELLS\
                        --n_classes $N_CLASSES\
                        --n_cpus $N_CPUS\
                        --learning_rate 0.002\
                        --n_replicates 10\
                        --number_of_steps 2500\
                        --n_virtual_cells_model $N_VIRTUAL_CELLS_MODEL $2 $3 $4
                    sleep 0.01
                done
            fi
        done
    done
done
