#!/bin/bash
#####
# Per-cell-type contribution scoring on finetuned_v2 (peak regression).
# Follows CREsted enhancer_code_analysis recipe:
#   - inputs the predictions+combined-layer AnnData built by predict_and_save_adata.sh
#   - sort_and_filter_regions_on_specificity(top_k=1000, model_name="combined")
#   - contribution_scores_specific (writes one NPZ per class for tfmodisco-lite)
#
# Submit (parallelize 22 cell types, 4 concurrent):
#   sbatch --partition=carter-gpu --account=carter-gpu --gres=gpu:a30:1 \
#     --cpus-per-task=4 --mem=64G --time=14-00:00:00 \
#     --array=1-22%4 --job-name=pr_contrib_spec_v2 \
#     --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/%x.%A_%a.out \
#     <this>.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
NVIDIA_LIBS=$VIRTUAL_ENV/lib/python3.11/site-packages/nvidia
export LD_LIBRARY_PATH=$NVIDIA_LIBS/cublas/lib:$NVIDIA_LIBS/cuda_runtime/lib:$NVIDIA_LIBS/cudnn/lib:$NVIDIA_LIBS/cufft/lib:$NVIDIA_LIBS/curand/lib:$NVIDIA_LIBS/cusolver/lib:$NVIDIA_LIBS/cusparse/lib:$NVIDIA_LIBS/nccl/lib:$NVIDIA_LIBS/nvjitlink/lib:$LD_LIBRARY_PATH
export KERAS_BACKEND=tensorflow

celltypes=(
    DE ENP_phase1 FB_FLT1 PFG1 PFG2 PGT1 PGT2 PGT3 PP1 PP2
    SC_delta_GHRL early_ENP early_SC_EC early_SC_alpha early_SC_beta
    exocrine late_ENP late_SC_EC late_SC_alpha late_SC_beta liver
    proliferating_endocrine
)
celltype=${celltypes[$SLURM_ARRAY_TASK_ID - 1]}
echo "Cell type: $celltype"

SCRIPT_DIR=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_multi_task_model/scripts/contributions
BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model

python ${SCRIPT_DIR}/compute_contributions_specific_v2.py \
    --input_h5ad ${BASE}/crested/contributions_specific_v2/adata_with_predictions.h5ad \
    --checkpoint ${BASE}/crested/finetuned_v2/checkpoints/28.keras \
    --genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --chrom_sizes /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
    --output_dir ${BASE}/crested/contributions_specific_v2 \
    --cell_types ${celltype} \
    --top_k_regions 1000 \
    --method integrated_grad \
    --batch_size 128

date
