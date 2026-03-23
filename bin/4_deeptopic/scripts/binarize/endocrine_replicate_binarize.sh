#!/bin/bash
# Select n_topics=50 model, binarize topics, and export BED files for CREsted
# Submit: /cellar/users/aklie/opt/SLURM/cpu.sh -s <this>.sh -j binarize_endocrine_replicate -m 32G -t 01:00:00 -c 1

date
echo -e "Job ID: $SLURM_JOB_ID\n"

source /cellar/users/aklie/projects/igvf/single_cell_utilities/pycistopic/.venv/bin/activate

cmd="python /cellar/users/aklie/projects/ML4GLand/deeptopic/scripts/pycistopic/select_and_binarize.py \
    --cistopic_obj /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_deeptopic/endocrine_replicate/objects/sc-islet-differentiation_10X-Multiome_endocrine_replicate.pkl \
    --models /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_deeptopic/endocrine_replicate/models_mallet/merged \
    --output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_deeptopic/endocrine_replicate/binarized \
    --n_topics 50 \
    --method otsu \
    --smooth_topics"
echo -e "Running:\n$cmd\n"
eval $cmd

date
