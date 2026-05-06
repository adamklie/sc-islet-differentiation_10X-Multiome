#!/bin/bash
set -euo pipefail
#####
# Per-topic contribution scoring for the n100 CREsted topic model (all subset).
# Follows CREsted enhancer_code_analysis recipe:
#   - inputs the predictions+combined-layer AnnData built by all_n100_predict_and_save_adata.sh
#   - sort_and_filter_regions_on_specificity(top_k=1000, model_name="combined")
#   - contribution_scores_specific (writes one NPZ per topic for tfmodisco-lite)
#
# Submit (parallelize 100 topics, 4 concurrent):
#   sbatch --partition=carter-gpu --account=carter-gpu --gres=gpu:a30:1 \
#     --cpus-per-task=4 --mem=64G --time=14-00:00:00 \
#     --array=1-100%4 --job-name=n100_contrib_spec \
#     --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/%x.%A_%a.out \
#     <this>.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
export KERAS_BACKEND=torch

TOPIC="Topic${SLURM_ARRAY_TASK_ID}"
echo "Topic: $TOPIC"

TOOLS=/cellar/users/aklie/opt/crested/scripts
BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/n100

python ${TOOLS}/compute_contributions_specific.py \
    --input_h5ad ${BASE}/crested_model/contributions_specific/adata_with_predictions.h5ad \
    --checkpoint ${BASE}/crested_model/checkpoints/60.keras \
    --genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --chrom_sizes /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
    --output_dir ${BASE}/crested_model/contributions_specific \
    --cell_types ${TOPIC} \
    --top_k_regions 1000 \
    --method integrated_grad \
    --batch_size 128 \
    --backend torch

date
