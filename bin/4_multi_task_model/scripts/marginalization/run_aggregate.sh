#!/bin/bash
# Aggregate per-CT marginalization NPZs → motif × CT matrices.
# Submit with afterok dependency on the marg array job:
#   sbatch --dependency=afterok:<MARG_ARRAY_ID> --time=00:30:00 -p carter-compute \
#          -c 2 --mem=8G --output=slurm_logs/marg_agg.%j.out run_aggregate.sh
set -euo pipefail
PY=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python
BASE=/cellar/users/aklie/projects/islet_organoid_differentiation/sc-islet-differentiation_10X-Multiome
$PY $BASE/bin/4_multi_task_model/scripts/marginalization/aggregate_marginalization.py \
    --marginalization_root $BASE/results/4_multi_task_model/crested/marginalization \
    --output_dir $BASE/results/4_multi_task_model/crested/motifs/analysis
