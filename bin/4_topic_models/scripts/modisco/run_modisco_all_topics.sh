#!/bin/bash

#####
# Run MoDISco on all completed CREsted contribution topics
# Step 1: Convert npz -> h5 for each topic
# Step 2: modisco motifs (de novo motif discovery)
# Step 3: modisco report (TomTom annotation + HTML)
# Step 4: Extract PFMs from modisco h5
#
# USAGE: sbatch \
#   --job-name=modisco_topics \
#   --partition=carter-compute \
#   --output=slurm_logs/4_topic_models/modisco/%x.%A.%a.out \
#   --array=1-50%8 \
#   --mem=128G \
#   --cpus-per-task=16 \
#   -t 14-00:00:00 \
#   run_modisco_all_topics.sh <subset>
#
# Example: sbatch ... run_modisco_all_topics.sh endocrine_replicate
#####

date
echo -e "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID\n"

# Config
PYTHON=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python
MODISCO_BIN=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/modisco
PFM_SCRIPT=/cellar/users/aklie/projects/ML4GLand/chrombpnet/scripts/modisco_to_pfm.py
MEME_DB=/cellar/users/aklie/data/ref/motifs/jvierstra/motif-clustering-v2.0beta/motifs.meme
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Subset from command line arg
SUBSET=${1:-endocrine_replicate}
BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/${SUBSET}
CONTRIBS=${BASE}/crested_model/contributions
MODISCO_OUT=${BASE}/crested_model/modisco

TOPIC="Topic${SLURM_ARRAY_TASK_ID}"
TOPIC_CONTRIBS=${CONTRIBS}/${TOPIC}
TOPIC_OUT=${MODISCO_OUT}/${TOPIC}

# Skip if contributions not ready
if [ ! -f "${TOPIC_CONTRIBS}/projected_contribs.npz" ]; then
    echo "Contributions not ready for ${TOPIC}, skipping."
    exit 0
fi

# Skip if already done
if [ -f "${TOPIC_OUT}/${TOPIC}.modisco.h5" ]; then
    echo "MoDISco already complete for ${TOPIC}, skipping."
    exit 0
fi

mkdir -p ${TOPIC_OUT}

# Step 1: Convert npz to h5 (channels-first format matching ChromBPNet)
echo "Step 1: Converting npz -> h5 for ${TOPIC}..."
$PYTHON -c "
import numpy as np, h5py, sys
topic, contribs_dir, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
proj = np.load(f'{contribs_dir}/projected_contribs.npz')['arr_0']
hypo = np.load(f'{contribs_dir}/hypothetical_contribs.npz')['arr_0']
one_hot = np.load(f'{contribs_dir}/one_hot_seqs.npz')['arr_0']
print(f'  {topic}: {proj.shape[0]} regions, shape={proj.shape}')
h5_path = f'{out_dir}/{topic}_contribs.h5'
with h5py.File(h5_path, 'w') as f:
    f.create_dataset('projected_shap/seq', data=proj.astype(np.float16))
    f.create_dataset('shap/seq', data=hypo.astype(np.float16))
    f.create_dataset('raw/seq', data=one_hot.astype(np.int8))
print(f'  Saved: {h5_path}')
" ${TOPIC} ${TOPIC_CONTRIBS} ${TOPIC_OUT}

# Step 2: MoDISco motif discovery
echo "Step 2: Running modisco motifs for ${TOPIC}..."
$MODISCO_BIN motifs \
    -i ${TOPIC_OUT}/${TOPIC}_contribs.h5 \
    -n 50000 \
    -o ${TOPIC_OUT}/${TOPIC}.modisco.h5 \
    -l 2 \
    -w 400 \
    -v

# Step 3: MoDISco report (TomTom annotation)
echo "Step 3: Running modisco report for ${TOPIC}..."
$MODISCO_BIN report \
    -i ${TOPIC_OUT}/${TOPIC}.modisco.h5 \
    -o ${TOPIC_OUT}/${TOPIC}.modisco.report/ \
    -s ${TOPIC_OUT}/${TOPIC}.modisco.report/ \
    -m ${MEME_DB}

# Step 4: Extract PFMs
echo "Step 4: Extracting PFMs for ${TOPIC}..."
mkdir -p ${TOPIC_OUT}/pfms
$PYTHON ${PFM_SCRIPT} \
    -m ${TOPIC_OUT}/${TOPIC}.modisco.h5 \
    -o ${TOPIC_OUT}/pfms \
    -op ${TOPIC} \
    -f 2

# Consolidate PFMs into single file per topic
temp_file="${TOPIC_OUT}/pfms/${TOPIC}.tmp.pfm"
final_file="${TOPIC_OUT}/pfms/${TOPIC}.pfm"
rm -f "$final_file" "$temp_file"
for x in ${TOPIC_OUT}/pfms/${TOPIC}*.pfm; do
    [ "$x" = "$final_file" ] && continue
    [ "$x" = "$temp_file" ] && continue
    echo ">$(basename "$x")" >> "$temp_file"
    cat "$x" >> "$temp_file"
done
[ -f "$temp_file" ] && mv "$temp_file" "$final_file"

echo "Done with ${TOPIC}!"
date
