#!/bin/bash
# Run cisTopic models (CGS)
# Submit: /cellar/users/aklie/opt/SLURM/cpu.sh -s <this>.sh -j cistopic_cgs -m 64G -t 24:00:00 -c 8

date
echo -e "Job ID: $SLURM_JOB_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv/bin/activate

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/pycistopic/run_models.py \
    --input /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/endocrine/objects/sc-islet-differentiation_10X-Multiome_endocrine.pkl \
    --output /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/endocrine/models_cgs \
    --n_topics 2 5 10 12 14 16 18 20 22 24 26 28 30 35 40 45 50 \
    --method cgs \
    --n_cpu 8 \
    --n_iter 150 \
    --alpha 50.0 \
    --eta 0.1 \
    --seed 555"
echo -e "Running:\n$cmd\n"
eval $cmd

date
