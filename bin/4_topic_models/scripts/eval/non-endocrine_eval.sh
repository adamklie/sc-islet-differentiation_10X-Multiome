#!/bin/bash
# Evaluate cisTopic models
# Submit: /cellar/users/aklie/opt/SLURM/cpu.sh -s <this>.sh -j cistopic_eval -m 32G -t 04:00:00 -c 1

date
echo -e "Job ID: $SLURM_JOB_ID\n"

conda deactivate 2>/dev/null
source /cellar/users/aklie/opt/deeptopic/.venv/bin/activate

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/pycistopic/evaluate_models.py \
    --cistopic_obj /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/non-endocrine/objects/sc-islet-differentiation_10X-Multiome_non-endocrine.pkl \
    --models /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/non-endocrine/models_mallet/merged \
    --output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/non-endocrine/evaluation \
    --groupby cell_type"
echo -e "Running:\n$cmd\n"
eval $cmd

date
