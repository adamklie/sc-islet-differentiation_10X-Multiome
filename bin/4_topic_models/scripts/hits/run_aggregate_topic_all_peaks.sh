#!/bin/bash
#####
# Aggregate per-topic finemo hits.tsv from the topic all-peaks run into
# topic x motif matrices. Topic-side analog of
# bin/4_multi_task_model/scripts/hits/run_aggregate_all_peaks.sh.
#
# Submit (CPU):
#   /cellar/users/aklie/opt/SLURM/cpu.sh -s <this>.sh -j topic_allpeak_agg \
#     -m 16G -t 04:00:00 -c 4 \
#     -o bin/slurm_logs/4_topic_models/hits_agg
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

source /cellar/users/aklie/opt/miniconda3/etc/profile.d/conda.sh
conda activate eugene_tools

BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome
SCRIPT=${BASE}/bin/4_topic_models/scripts/hits/aggregate_hits_topic_all_peaks.py
HITS_ROOT=${BASE}/results/4_topic_models/all/n100/crested_model/motifs/hits_all_peaks
OUT=${BASE}/results/4_topic_models/all/n100/crested_model/motifs/hits_all_peaks

python ${SCRIPT} \
    --hits_root ${HITS_ROOT} \
    --output_dir ${OUT}

date
