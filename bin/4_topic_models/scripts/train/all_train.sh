#!/bin/bash
# Train DeepTopic CNN using CREsted
# Submit: /cellar/users/aklie/opt/SLURM/gpu.sh -s <this>.sh -j crested_train -g a30:1 -m 32G -t 08:00:00 -c 4

date
echo -e "Job ID: $SLURM_JOB_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv/bin/activate

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/crested/train_topic_model.py \
    --beds_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/binarized/beds \
    --regions_file /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/binarized/consensus_regions.bed \
    --genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/crested_model \
    --project_name sc-islet-differentiation_10X-Multiome_all \
    "
echo -e "Running:\n$cmd\n"
eval $cmd

date
