#! /bin/bash

#####
# Script to fine-tune a CREsted peak regression model on cell-type-specific regions
# USAGE: sbatch \
#--job-name=crested_finetune \
#--partition=carter-gpu \
#--account=carter-gpu \
#--gpus=a30:1 \
#--output slurm_logs/%x.%A.out \
#--mem=64G \
#-n 4 \
#-t 02-00:00:00 \
#crested_finetune.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Set-up env
source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
# TF needs pip-installed nvidia libs on LD_LIBRARY_PATH
NVIDIA_LIBS=$VIRTUAL_ENV/lib/python3.11/site-packages/nvidia
export LD_LIBRARY_PATH=$NVIDIA_LIBS/cublas/lib:$NVIDIA_LIBS/cuda_runtime/lib:$NVIDIA_LIBS/cudnn/lib:$NVIDIA_LIBS/cufft/lib:$NVIDIA_LIBS/curand/lib:$NVIDIA_LIBS/cusolver/lib:$NVIDIA_LIBS/cusparse/lib:$NVIDIA_LIBS/nccl/lib:$NVIDIA_LIBS/nvjitlink/lib:$LD_LIBRARY_PATH

# Run cmd
cmd="python /cellar/users/aklie/opt/crested/scripts/finetune_peak_regression.py \
--base_training_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/basemodel \
--bigwigs_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/bigwigs \
--regions_file /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/consensus_peaks.bed \
--genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
--chrom_sizes /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
--output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/finetuned \
--splits_file /cellar/users/aklie/data/ref/genomes/hg38/chrombpnet/splits/fold_0.json \
--project_name sc-islet-differentiation_10X-Multiome \
--run_name finetuned \
--backend tensorflow \
--gini_std_threshold 1.0 \
--learning_rate 0.0001 \
--epochs 60 \
--batch_size 64 \
--seq_len 2114 \
--top_k_percent 0.03 \
--seed 42 \
--mixed_precision"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Date
date
