#!/bin/bash
# Fine-tune peak regression model v2 — lower LR (1e-5), tighter early stopping
# Based on CREsted tutorial: finetuning should converge quickly (2-5 epochs)
# Previous run with lr=1e-4 caused catastrophic forgetting (29 epochs, best at 22)

date
echo -e "Job ID: $SLURM_JOB_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
NVIDIA_LIBS=$VIRTUAL_ENV/lib/python3.11/site-packages/nvidia
export LD_LIBRARY_PATH=$NVIDIA_LIBS/cublas/lib:$NVIDIA_LIBS/cuda_runtime/lib:$NVIDIA_LIBS/cudnn/lib:$NVIDIA_LIBS/cufft/lib:$NVIDIA_LIBS/curand/lib:$NVIDIA_LIBS/cusolver/lib:$NVIDIA_LIBS/cusparse/lib:$NVIDIA_LIBS/nccl/lib:$NVIDIA_LIBS/nvjitlink/lib:$LD_LIBRARY_PATH

cmd="python /cellar/users/aklie/opt/crested/scripts/finetune_peak_regression.py \
--base_training_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/basemodel \
--bigwigs_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/bigwigs \
--regions_file /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/consensus_peaks.bed \
--genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
--chrom_sizes /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
--output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/finetuned_v2 \
--splits_file /cellar/users/aklie/data/ref/genomes/hg38/chrombpnet/splits/fold_0.json \
--project_name sc-islet-differentiation_10X-Multiome \
--run_name finetuned_v2 \
--backend tensorflow \
--gini_std_threshold 1.0 \
--learning_rate 1e-5 \
--epochs 30 \
--batch_size 64 \
--seq_len 2114 \
--top_k_percent 0.03 \
--seed 42 \
--lr_reduce_patience 3 \
--early_stopping_patience 6 \
--mixed_precision"
echo -e "Running command:\n$cmd\n"
eval $cmd

date
