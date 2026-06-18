#!/bin/bash
#####
# Aggregate per-CT finemo hits.tsv from the all-peaks run into motif x CT matrices.
# Mirrors run_aggregate_v3.sh but points at hits_all_peaks/.
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Use eugene_tools python directly — conda activate is unreliable in batch jobs.
PY=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python

BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome
SCRIPT=${BASE}/bin/4_multi_task_model/scripts/hits/aggregate_hits.py
# Updated to the all-peak dual-track catalog hits path (was motifs_v3/hits_all_peaks/).
HITS_ROOT=${BASE}/results/4_multi_task_model/crested/motifs_all_peaks/hits
OUT=${BASE}/results/4_multi_task_model/crested/motifs_all_peaks/hits
META=${BASE}/config/cell_type_metadata.tsv

${PY} ${SCRIPT} \
    --hits_root ${HITS_ROOT} \
    --output_dir ${OUT} \
    --metadata ${META}

date
