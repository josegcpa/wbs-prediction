for cell_type in wbc rbc
do
    if [[ $cell_type == wbc ]]
    then
        size=128
    else
        size=64
    fi

    F=/nfs/research/gerstung/josegcpa/data/SLIDES/ADDEN_NDPI
    mkdir -p /nfs/research/gerstung/josegcpa/data/all_cells/ADDEN_NDPI
    for slide in $F/*ndpi
    do
        root=$(basename $slide | cut -d '.' -f 1)
        seg_path=/nfs/research/gerstung/josegcpa/data/SLIDES/ADDEN_NDPI/_segmented_$cell_type/$root.h5
        agg_path=/nfs/research/gerstung/josegcpa/data/SLIDES/ADDEN_NDPI/_aggregates_$cell_type/$root.h5
        bsub -n 1 -M 2000 -o /dev/null -e python3 scripts/python/get_all_cells.py \
            --slide_path $slide \
            --segmented_path $seg_path \
            --aggregates_path $agg_path \
            --size $size \
            --output_path /nfs/research/gerstung/josegcpa/data/all_cells/ADDEN_NDPI/$cell_type-$root.h5
    done

    F=/nfs/research/gerstung/josegcpa/data/SLIDES/MLL_TIFF
    mkdir -p /nfs/research/gerstung/josegcpa/data/all_cells/MLL_TIFF
    for slide in $F/*tiff
    do
        root=$(basename $slide | cut -d '.' -f 1)
        seg_path=/nfs/research/gerstung/josegcpa/data/SLIDES/MLL_TIFF/_segmented_$cell_type/$root.h5
        agg_path=/nfs/research/gerstung/josegcpa/data/SLIDES/MLL_TIFF/_aggregates_$cell_type/$root.h5
        bsub -n 1 -M 2000 -o /dev/null -e python3 scripts/python/get_all_cells.py \
            --slide_path $slide \
            --segmented_path $seg_path \
            --aggregates_path $agg_path \
            --size $size \
            --output_path /nfs/research/gerstung/josegcpa/data/all_cells/MLL_TIFF/$cell_type-$root.h5
    done

    F=/nfs/research/gerstung/josegcpa/data/SLIDES/ADDEN_SVS
    mkdir -p /nfs/research/gerstung/josegcpa/data/all_cells/ADDEN_SVS
    for slide in $F/*svs
    do
        root=$(basename $slide | cut -d '.' -f 1)
        seg_path=/nfs/research/gerstung/josegcpa/data/SLIDES/ADDEN_SVS_results/_segmented_$cell_type/$root.h5
        agg_path=/nfs/research/gerstung/josegcpa/data/SLIDES/ADDEN_SVS_results/_aggregates_$cell_type/$root.h5
        bsub -n 1 -M 2000 -o /dev/null -e /dev/null python3 scripts/python/get_all_cells.py \
            --slide_path $slide \
            --segmented_path $seg_path \
            --aggregates_path $agg_path \
            --size $size \
            --output_path /nfs/research/gerstung/josegcpa/data/all_cells/ADDEN_SVS/$cell_type-$root.h5
    done

done
