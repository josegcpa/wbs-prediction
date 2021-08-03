for simulation in simulation_outputs/*pkl
do
    bsub -n 8 -M 16000\
        -o logs/VALIDATE_$(basename ${simulation:1:$((${#simulation}-5))}).o\
        -e logs/VALIDATE_$(basename ${simulation:1:$((${#simulation}-5))}).e\
        python3 validate_method.py\
        --simulation_path $simulation\
        --n_iterations 50
done
