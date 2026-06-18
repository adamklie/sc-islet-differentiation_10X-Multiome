#!/bin/bash
# Select best model, binarize topics, and export BED files for CREsted
# Submit: /cellar/users/aklie/opt/SLURM/cpu.sh -s <this>.sh -j binarize -m 32G -t 01:00:00 -c 1

date
echo -e "Job ID: $SLURM_JOB_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv/bin/activate

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/pycistopic/select_and_binarize.py \
    --cistopic_obj /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/endocrine_replicate/objects/sc-islet-differentiation_10X-Multiome_endocrine_replicate.pkl \
    --models /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/endocrine_replicate/models_mallet/merged \
    --output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/endocrine_replicate/binarized \
    "
echo -e "Running:\n$cmd\n"
eval $cmd

date
