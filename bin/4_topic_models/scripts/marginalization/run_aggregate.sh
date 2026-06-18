#!/bin/bash
# Aggregate per-topic DeepTopic marginalization NPZs → motif × topic matrices.
# Reads from crested_model/marginalization/Topic{N}/marginalization_data.counts.npz
# Writes to   crested_model/motifs/analysis/marginalization_matrix{,_zscore}.csv
#
# Submit with afterok dependency on the marg array job:
#   sbatch --dependency=afterok:<MARG_ARRAY_ID> --time=00:30:00 \
#          --partition=carter-compute -c 2 --mem=8G \
#          --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_topic_models/marginalization/topic_marg_agg.%j.out \
#          /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_topic_models/scripts/marginalization/run_aggregate.sh
set -euo pipefail
PY=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python
BIN=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_topic_models/scripts/marginalization
BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/n100
$PY $BIN/aggregate_marginalization.py \
    --marginalization_root $BASE/crested_model/marginalization \
    --output_dir $BASE/crested_model/motifs/analysis
