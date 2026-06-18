#!/bin/bash
# Aggregate per-CT finemo hits → motif × CT matrices.
set -euo pipefail
PY=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python
BASE=/cellar/users/aklie/projects/islet_organoid_differentiation/sc-islet-differentiation_10X-Multiome
$PY $BASE/bin/4_multi_task_model/scripts/hits/aggregate_hits.py \
    --hits_root $BASE/results/4_multi_task_model/crested/motifs/hits \
    --output_dir $BASE/results/4_multi_task_model/crested/motifs/hits
