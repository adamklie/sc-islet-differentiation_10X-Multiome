#! /bin/bash

#####
# Script to evaluate a CREsted peak regression model
# USAGE: sbatch \
#--job-name=crested_evaluate \
#--partition=carter-gpu \
#--account=carter-gpu \
#--gpus=a30:1 \
#--output slurm_logs/%x.%A.out \
#--mem=32G \
#-n 4 \
#-t 00-04:00:00 \
#crested_evaluate.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Set-up env
source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
# TF needs pip-installed nvidia libs on LD_LIBRARY_PATH
NVIDIA_LIBS=$VIRTUAL_ENV/lib/python3.11/site-packages/nvidia
export LD_LIBRARY_PATH=$NVIDIA_LIBS/cublas/lib:$NVIDIA_LIBS/cuda_runtime/lib:$NVIDIA_LIBS/cudnn/lib:$NVIDIA_LIBS/cufft/lib:$NVIDIA_LIBS/curand/lib:$NVIDIA_LIBS/cusolver/lib:$NVIDIA_LIBS/cusparse/lib:$NVIDIA_LIBS/nccl/lib:$NVIDIA_LIBS/nvjitlink/lib:$LD_LIBRARY_PATH

# Run cmd
cmd="python /cellar/users/aklie/opt/crested/scripts/evaluate_peak_regression.py \
--training_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/finetuned \
--bigwigs_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/bigwigs \
--regions_file /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/consensus_peaks.bed \
--genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
--chrom_sizes /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
--output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/evaluation \
--splits_file /cellar/users/aklie/data/ref/genomes/hg38/chrombpnet/splits/fold_0.json \
--backend tensorflow \
--batch_size 64"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Date
date
