#!/bin/bash
# Aggregate per-CT v3 marginalization NPZs → motif × CT matrices (v3).
# Reads from crested/marginalization_v3/<CT>/marginalization_data.counts.npz
# Writes to   crested/motifs_v3/analysis/marginalization_matrix{,_zscore}.csv
#
# Submit with afterok dependency on the marg array job:
#   sbatch --dependency=afterok:<MARG_V3_ARRAY_ID> --time=00:30:00 \
#          --partition=carter-compute -c 2 --mem=8G \
#          --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_multi_task_model/pr_marg_agg_v3.%j.out \
#          run_aggregate_v3.sh
set -euo pipefail
PY=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python
BASE=/cellar/users/aklie/projects/islet_organoid_differentiation/sc-islet-differentiation_10X-Multiome
$PY $BASE/bin/4_multi_task_model/scripts/marginalization/aggregate_marginalization.py \
    --marginalization_root $BASE/results/4_multi_task_model/crested/marginalization_v3 \
    --output_dir $BASE/results/4_multi_task_model/crested/motifs_v3/analysis
