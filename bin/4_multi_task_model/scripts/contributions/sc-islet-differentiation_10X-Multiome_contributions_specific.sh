#!/bin/bash
#####
# Compute contributions on CELL-TYPE-SPECIFIC regions (Gini-filtered + z-scored)
# Uses modified script that selects top 5K by z-score instead of raw count.
#
# Submit: gpu_array.sh -s <this>.sh -j pr_contrib_spec -m 64G -t 14-00:00:00 -c 4 -g a30:1 -n 22 -x 4
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

python ${SCRIPT_DIR}/compute_contributions_specific.py \
    --checkpoint $(cat ${BASE}/crested/evaluation/best_checkpoint.txt) \
    --bigwigs_dir ${BASE}/bigwigs \
    --regions_file /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/consensus_peaks.bed \
    --genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --chrom_sizes /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
    --output_dir ${BASE}/crested/contributions_specific/${celltype} \
    --splits_file /cellar/users/aklie/data/ref/genomes/hg38/chrombpnet/splits/fold_0.json \
    --cell_types ${celltype} \
    --method expected_integrated_grad \
    --top_k_regions 5000 \
    --backend tensorflow \
    --batch_size 64

date
