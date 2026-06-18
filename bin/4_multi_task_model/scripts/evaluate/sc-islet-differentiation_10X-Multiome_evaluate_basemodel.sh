#!/bin/bash
# Evaluate the BASE peak regression model (for comparison with finetuned)
date
echo -e "Job ID: $SLURM_JOB_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
NVIDIA_LIBS=$VIRTUAL_ENV/lib/python3.11/site-packages/nvidia
export LD_LIBRARY_PATH=$NVIDIA_LIBS/cublas/lib:$NVIDIA_LIBS/cuda_runtime/lib:$NVIDIA_LIBS/cudnn/lib:$NVIDIA_LIBS/cufft/lib:$NVIDIA_LIBS/curand/lib:$NVIDIA_LIBS/cusolver/lib:$NVIDIA_LIBS/cusparse/lib:$NVIDIA_LIBS/nccl/lib:$NVIDIA_LIBS/nvjitlink/lib:$LD_LIBRARY_PATH

cmd="python /cellar/users/aklie/opt/crested/scripts/evaluate_peak_regression.py \
--training_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/basemodel \
--bigwigs_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/bigwigs \
--regions_file /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/consensus_peaks.bed \
--genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
--chrom_sizes /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
--output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/evaluation_basemodel \
--splits_file /cellar/users/aklie/data/ref/genomes/hg38/chrombpnet/splits/fold_0.json \
--backend tensorflow \
--batch_size 64"
echo -e "Running command:\n$cmd\n"
eval $cmd

date
