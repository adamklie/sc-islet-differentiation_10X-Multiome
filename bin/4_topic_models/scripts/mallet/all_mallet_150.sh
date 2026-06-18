#!/bin/bash
# Run cisTopic MALLET model for 150 topics
# Submit: /cellar/users/aklie/opt/SLURM/cpu.sh -s <this>.sh -j mallet_150 -m 128G -t 14-00:00:00 -c 8

date
echo -e "Job ID: $SLURM_JOB_ID\n"

conda deactivate 2>/dev/null
source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/pycistopic/run_models.py \
    --input /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/objects/sc-islet-differentiation_10X-Multiome_all.pkl \
    --output /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/models_mallet/topic_150 \
    --n_topics 150 \
    --method mallet \
    --n_cpu 8 \
    --n_iter 150 \
    --alpha 50.0 \
    --eta 0.1 \
    --seed 555 \
    --mallet_path /cellar/users/aklie/opt/Mallet-202108/bin/mallet \
    --mallet_memory 120G"
echo -e "Running:\n$cmd\n"
eval $cmd

date
