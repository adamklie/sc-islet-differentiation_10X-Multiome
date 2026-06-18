#!/bin/bash
#####
# Multi-task (peak regression) motif marginalization, per cell type.
# Uses the unified 15-motif catalog (combined_mapped.meme) and the finetuned_v2
# CREsted model. Background sequences are drawn from each CT's narrowPeak file
# (same backgrounds the ChromBPNet pipeline uses).
#
# Submit as a 22-CT array (concurrency 4 to keep GPU busy without hogging):
#   /cellar/users/aklie/opt/SLURM/gpu_array.sh \
#     -s /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_multi_task_model/scripts/marginalization/sc-islet-differentiation_10X-Multiome_motif_marginalization.sh \
#     -j pr_marg -g a30:1 -c 4 -m 32G -t 14-00:00:00 -n 22 -x 4 \
#     -o /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/
#####

date
echo -e "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID\n"

# CREsted Keras venv (same env used for contributions_specific_v2)
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

SCRIPT=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_multi_task_model/scripts/marginalization/marginalize_multitask.py
BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model
PEAK_DIR=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all

n_seqs=100

python ${SCRIPT} \
    --model ${BASE}/crested/finetuned_v2/checkpoints/28.keras \
    --motif_file ${BASE}/crested/motifs/combined_mapped.meme \
    --genome_fasta /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --peaks ${PEAK_DIR}/${celltype}.narrowPeak \
    --adata ${BASE}/crested/contributions_specific_v2/adata_with_predictions.h5ad \
    --celltype ${celltype} \
    --output_dir ${BASE}/crested/marginalization \
    --num_seqs ${n_seqs} \
    --batch_size 64

date
