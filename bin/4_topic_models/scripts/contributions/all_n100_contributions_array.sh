#!/bin/bash
#####
# CREsted contributions for all/n100 — ARRAY JOB (1 topic per GPU)
# Each array task processes a single topic. Max 1,000 regions per topic
# (trimmed from 5k 2026-04-13 — MoDISco converges well below 5k and cuts GPU time).
#
# Submit: gpu_array.sh -s <this>.sh -j contrib_n100 -m 64G -t 14-00:00:00 -c 4 -g a30:1 -n 100 -x 4
#####

date
echo -e "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID\n"

conda deactivate 2>/dev/null
source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
export KERAS_BACKEND=torch

TOPIC_NUM=$SLURM_ARRAY_TASK_ID
TOPIC="Topic${TOPIC_NUM}"
echo "Processing: $TOPIC"

BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/n100

python /cellar/users/aklie/opt/deeptopic/scripts/crested/compute_contribution_scores.py \
    --training_dir ${BASE}/crested_model \
    --beds_dir ${BASE}/binarized/beds \
    --regions_file ${BASE}/binarized/consensus_regions.bed \
    --genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --output_dir ${BASE}/crested_model/contributions \
    --topics ${TOPIC} \
    --max_regions 1000 \
    --batch_size 64 \
    --method expected_integrated_grad \
    --backend torch

date
