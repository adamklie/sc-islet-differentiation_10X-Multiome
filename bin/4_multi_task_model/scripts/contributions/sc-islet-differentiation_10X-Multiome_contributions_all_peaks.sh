#!/bin/bash
#####
# Per-cell-type contribution scoring on finetuned_v2 (peak regression),
# computed over ALL peaks (no specificity filter). Writes per-CT NPZ pairs
# at contributions_all_peaks/<CT>_{oh,contrib}.npz — same schema as the
# specific_v2/v3 outputs so the existing finemo hit pipeline reuses unchanged.
#
# Submit via /cellar/users/aklie/opt/SLURM/gpu_array.sh
#   -s <this>.sh -j mt_allpeak_contrib -g a30:1 -c 4 -m 64G -t 14-00:00:00
#   -n 22 -x 4 -o bin/slurm_logs/4_multi_task_model/contributions_all_peaks
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

mkdir -p ${BASE}/crested/contributions_all_peaks

NARROWPEAK=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/${celltype}.narrowPeak

python ${SCRIPT_DIR}/compute_contributions_all_peaks.py \
    --input_h5ad ${BASE}/crested/contributions_specific_v2/adata_with_predictions.h5ad \
    --checkpoint ${BASE}/crested/finetuned_v2/checkpoints/28.keras \
    --genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --chrom_sizes /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
    --output_dir ${BASE}/crested/contributions_all_peaks \
    --cell_types ${celltype} \
    --peaks_bed ${NARROWPEAK} \
    --method integrated_grad \
    --batch_size 128 \
    --chunk_size 20000 \
    --skip_if_exists

date
