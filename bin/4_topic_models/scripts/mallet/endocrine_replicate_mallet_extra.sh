#!/bin/bash
# Run cisTopic models (MALLET) — one topic per SLURM array task
# Submit: /cellar/users/aklie/opt/SLURM/cpu_array.sh -s <this>.sh -j cistopic_mallet -m 128G -t 24:00:00 -c 8 -n <num_topics> -x <num_topics>

date
echo -e "Job ID: $SLURM_JOB_ID (array task: $SLURM_ARRAY_TASK_ID)\n"

conda deactivate 2>/dev/null
source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate

# Select one n_topics value for this array task (1-indexed)
N_TOPICS=(55 60 75 100)
topic=${N_TOPICS[$SLURM_ARRAY_TASK_ID-1]}
echo -e "Running topic: $topic\n"

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/pycistopic/run_models.py \
    --input /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/endocrine_replicate/objects/sc-islet-differentiation_10X-Multiome_endocrine_replicate.pkl \
    --output /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/endocrine_replicate/models_mallet/topic_${topic} \
    --n_topics $topic \
    --method mallet \
    --n_cpu 8 \
    --n_iter 150 \
    --alpha 50.0 \
    --eta 0.1 \
    --seed 555 \
    --mallet_path /cellar/users/aklie/opt/Mallet-202108/bin/mallet \
    --mallet_memory 30G"
echo -e "Running:\n$cmd\n"
eval $cmd

date
