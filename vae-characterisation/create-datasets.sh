# training
echo bsub -n 4 -M 8000 python3 scripts/assemble-dataset.py --segmentations_folder /hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/_segmented_rbc/ --slides_folder /hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/ --output_path datasets/rbc.h5 --size 64
echo bsub -n 4 -M 8000 python3 scripts/assemble-dataset.py --segmentations_folder /hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/_segmented_wbc/ --slides_folder /hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/ --output_path datasets/wbc.h5 --size 128
bsub -n 4 -M 8000 python3 scripts/assemble-dataset.py --segmentations_folder /hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/_segmented_wbc/ --slides_folder /hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/ --output_path datasets/wbc_96.h5 --size 96

# validation
echo bsub -n 4 -M 8000 python3 scripts/assemble-dataset.py --segmentations_folder /hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/_segmented_rbc/ --slides_folder /hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/ --output_path datasets/rbc_val.h5 --size 64 --subset 100
echo bsub -n 4 -M 8000 python3 scripts/assemble-dataset.py --segmentations_folder /hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/_segmented_wbc/ --slides_folder /hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/ --output_path datasets/wbc_val.h5 --size 128 --subset 100
bsub -n 4 -M 8000 python3 scripts/assemble-dataset.py --segmentations_folder /hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/_segmented_wbc/ --slides_folder /hps/nobackup/research/gerstung/josegcpa/data/MLL_TIFF/ --output_path datasets/wbc_96_val.h5 --size 96 --subset 100


