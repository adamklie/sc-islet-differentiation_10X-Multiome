#!/bin/bash
#####
# Topic n100 v2 pilot — contributions + MoDISco at max_regions=2000 for 5
# confident topics (Topic1, Topic9, Topic16, Topic44, Topic80) to vet whether
# a full v2 rebuild (vs current v1 at max_regions=1000) brings meaningful
# catalog uplift before committing 24-36h of GPU compute.
#
# Outputs land in parallel _pilot_v2/ directories so v1 stays untouched.
#
# Submit:
#   sbatch --partition=carter-gpu --account=carter-gpu --gres=gpu:a30:1 \
#     --cpus-per-task=4 --mem=64G --time=6:00:00 \
#     --array=1,9,16,44,80%5 --job-name=topic_pilot_v2 \
#     --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_topic_models/pilot_v2/%x.%A.%a.out \
#     <this>.sh
#####
set -uo pipefail
date
echo "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID"

source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
export KERAS_BACKEND=torch

TOPIC="Topic${SLURM_ARRAY_TASK_ID}"
echo "Processing: $TOPIC"

BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/n100
CONTRIB_OUT=${BASE}/crested_model/contributions_pilot_v2
MODISCO_OUT=${BASE}/crested_model/modisco_pilot_v2
mkdir -p ${CONTRIB_OUT} ${MODISCO_OUT}/${TOPIC}

# --- Step 1: contributions at max_regions=2000 ---
echo "=== Step 1: contributions (max_regions=2000) for $TOPIC ==="
python /cellar/users/aklie/opt/deeptopic/scripts/crested/compute_contribution_scores.py \
    --training_dir ${BASE}/crested_model \
    --beds_dir ${BASE}/binarized/beds \
    --regions_file ${BASE}/binarized/consensus_regions.bed \
    --genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --output_dir ${CONTRIB_OUT} \
    --topics ${TOPIC} \
    --max_regions 2000 \
    --batch_size 64 \
    --method integrated_grad \
    --backend torch
date

# --- Step 2: convert npz -> modisco h5 ---
echo "=== Step 2: convert npz -> modisco h5 for $TOPIC ==="
PY=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python
${PY} -c "
import numpy as np, h5py, sys
topic = '${TOPIC}'
src = '${CONTRIB_OUT}/${TOPIC}'
out = '${MODISCO_OUT}/${TOPIC}/${TOPIC}_contribs.h5'
proj = np.load(f'{src}/projected_contribs.npz')['arr_0']
hypo = np.load(f'{src}/hypothetical_contribs.npz')['arr_0']
oh   = np.load(f'{src}/one_hot_seqs.npz')['arr_0']
print(f'  {topic}: shape={proj.shape}')
with h5py.File(out, 'w') as f:
    f.create_dataset('projected_shap/seq', data=proj.astype(np.float16))
    f.create_dataset('shap/seq',          data=hypo.astype(np.float16))
    f.create_dataset('raw/seq',           data=oh.astype(np.int8))
print(f'  wrote {out}')
"
date

# --- Step 3: MoDISco motifs (matching v1 params: w=400, n=50000, l=2) ---
echo "=== Step 3: modisco motifs for $TOPIC ==="
MODISCO_BIN=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/modisco
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK
${MODISCO_BIN} motifs \
    -i ${MODISCO_OUT}/${TOPIC}/${TOPIC}_contribs.h5 \
    -n 50000 \
    -o ${MODISCO_OUT}/${TOPIC}/${TOPIC}.modisco.h5 \
    -l 2 \
    -w 400 \
    -v
date

echo "=== DONE for $TOPIC ==="
