#!/bin/bash
# Aggregate per-CT v3 finemo hits → motif × CT matrices (v3).
# Reads from motifs_v3/hits/<CT>/hits.tsv, writes to motifs_v3/hits/.
#
# Submit with afterok dependency on the hits array job:
#   sbatch --dependency=afterok:<HITS_V3_ARRAY_ID> --time=00:30:00 \
#          --partition=carter-compute -c 2 --mem=8G \
#          --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_multi_task_model/pr_hits_agg_v3.%j.out \
#          run_aggregate_v3.sh
set -euo pipefail
PY=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python
BASE=/cellar/users/aklie/projects/islet_organoid_differentiation/sc-islet-differentiation_10X-Multiome
$PY $BASE/bin/4_multi_task_model/scripts/hits/aggregate_hits.py \
    --hits_root $BASE/results/4_multi_task_model/crested/motifs_v3/hits \
    --output_dir $BASE/results/4_multi_task_model/crested/motifs_v3/hits
