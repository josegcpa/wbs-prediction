mkdir -p checkpoints

for S in A-0-cohort-demographics-blood-counts.R A-1-describe-cell-populations.R A-2-qc-analysis.R A-3-u-net-metrics.R B-0-feature-viz.R B-1-glmnet.R C-0-u-map.R C-3-annotated-cell-prediction-analysis.R
do
    ckpt=checkpoints/$S.check
    if [[ ! -f $ckpt ]]
    then
        echo running $S
        Rscript $S >/dev/null 2>/dev/null && touch $ckpt && echo     success $S
    fi
done

for S in C-1-analyse-metrics.R C-2-virtual-cell-analysis-consensus.R C-2-virtual-cell-analysis.R C-3-annotated-cell-prediction-analysis.R C-4-analyse-metrics-validation-consensus.R C-4-analyse-metrics-validation.R
do 
    for model_id in full full_dropout subset subset_dropout
    do
        ckpt=checkpoints/$S-$model_id.check
        if [[ ! -f $ckpt ]]
        then
            echo running $S for model $model_id
            Rscript $S $model_id  >/dev/null 2>/dev/null && touch $ckpt && echo     success $S
        fi
    done
done
