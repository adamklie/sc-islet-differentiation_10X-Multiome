#!/bin/bash

#####
# Submit CREsted contribution jobs for all subsets WITH subsampling (5K regions max)
#
# This uses the existing compute_contribution_scores.py with --max_regions 5000
# to subsample the top 5K highest-confidence positive regions per topic.
# This is ~4x faster than computing on all positive regions (~20K).
#
# USAGE: bash submit_contributions_subsampled.sh [subset1 subset2 ...]
# Default: all 6 subsets except endocrine_replicate (already done unsubsampled)
#
# Example:
#   bash submit_contributions_subsampled.sh endocrine non-endocrine_replicate
#   bash submit_contributions_subsampled.sh  # all 5 remaining subsets
#####

SCRIPT=/cellar/users/aklie/opt/deeptopic/scripts/crested/compute_contribution_scores.py
PYTHON=/cellar/users/aklie/opt/deeptopic/.venv-keras/bin/python
GENOME=/cellar/users/aklie/data/ref/genomes/hg38/hg38.fa
BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models
LOGDIR=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_topic_models/contributions

MAX_REGIONS=5000
BATCH_SIZE=64
METHOD=expected_integrated_grad
BACKEND=torch

# Default subsets (skip endo_rep — already done)
if [ $# -eq 0 ]; then
    SUBSETS=(endocrine non-endocrine_replicate non-endocrine all_replicate all)
else
    SUBSETS=("$@")
fi

mkdir -p $LOGDIR

for SUBSET in "${SUBSETS[@]}"; do
    TRAINING_DIR=${BASE}/${SUBSET}/crested_model
    BEDS_DIR=${BASE}/${SUBSET}/binarized/beds
    REGIONS_FILE=${BASE}/${SUBSET}/binarized/consensus_regions.bed
    OUTPUT_DIR=${BASE}/${SUBSET}/crested_model/contributions

    # Check that training is done
    if [ ! -f "${TRAINING_DIR}/evaluation/best_checkpoint.txt" ]; then
        echo "SKIP $SUBSET: no best_checkpoint.txt found"
        continue
    fi

    # Check beds dir exists
    if [ ! -d "$BEDS_DIR" ]; then
        echo "SKIP $SUBSET: beds dir not found at $BEDS_DIR"
        continue
    fi

    # Check if regions file exists, make it optional
    REGIONS_ARG=""
    if [ -f "$REGIONS_FILE" ]; then
        REGIONS_ARG="--regions_file $REGIONS_FILE"
    fi

    echo "Submitting $SUBSET (max_regions=$MAX_REGIONS)..."

    sbatch \
        --job-name=contrib_${SUBSET}_5k \
        --partition=carter-gpu \
        --gres=gpu:a30:1 \
        --mem=64G \
        --cpus-per-task=4 \
        --time=14-00:00:00 \
        --output=${LOGDIR}/contrib_${SUBSET}_5k_%j.out \
        --wrap="
export KERAS_BACKEND=${BACKEND}
${PYTHON} ${SCRIPT} \
    --training_dir ${TRAINING_DIR} \
    --beds_dir ${BEDS_DIR} \
    ${REGIONS_ARG} \
    --genome ${GENOME} \
    --output_dir ${OUTPUT_DIR}_5k \
    --max_regions ${MAX_REGIONS} \
    --batch_size ${BATCH_SIZE} \
    --method ${METHOD} \
    --backend ${BACKEND}
"
done

echo ""
echo "Note: Subsampled outputs go to contributions_5k/ (separate from unsubsampled contributions/)"
echo "      This preserves the partial unsubsampled runs."
