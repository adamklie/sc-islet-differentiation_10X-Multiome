#!/bin/bash
set -euo pipefail
# Predict n100 topic probabilities and save AnnData with predictions + combined layer + QC plot.
# Run ONCE before per-topic contribution scoring.
#
# Submit: sbatch --partition=carter-gpu --account=carter-gpu --gres=gpu:a30:1 \
#   --cpus-per-task=4 --mem=64G --time=04:00:00 --job-name=n100_predict \
#   --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/%x.%j.out \
#   <this>.sh

date
echo -e "Job ID: $SLURM_JOB_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
export KERAS_BACKEND=torch

SCRIPT_DIR=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_topic_models/scripts/contributions
BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/n100

python ${SCRIPT_DIR}/predict_and_save_adata_topics.py \
    --checkpoint ${BASE}/crested_model/checkpoints/60.keras \
    --beds_dir ${BASE}/binarized/beds \
    --regions_file ${BASE}/binarized/consensus_regions.bed \
    --genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --chrom_sizes /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
    --output_h5ad ${BASE}/crested_model/contributions_specific/adata_with_predictions.h5ad \
    --qc_plot ${BASE}/crested_model/contributions_specific/qc_sort_and_filter_cutoff.pdf

date
