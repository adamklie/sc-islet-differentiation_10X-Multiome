#!/bin/bash
# Select 100-topic model, binarize topics, and export BED files for CREsted
# Submit: /cellar/users/aklie/opt/SLURM/cpu.sh -s <this>.sh -j binarize_100 -m 64G -t 02:00:00 -c 1

date
echo -e "Job ID: $SLURM_JOB_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv/bin/activate

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/pycistopic/select_and_binarize.py \
    --cistopic_obj /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/objects/sc-islet-differentiation_10X-Multiome_all.pkl \
    --models /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/models_mallet/merged \
    --output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all_100/binarized \
    --n_topics 100 \
    --smooth_topics"
echo -e "Running:\n$cmd\n"
eval $cmd

date
