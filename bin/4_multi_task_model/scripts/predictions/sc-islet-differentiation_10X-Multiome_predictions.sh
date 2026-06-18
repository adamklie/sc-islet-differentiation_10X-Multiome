#! /bin/bash

#####
# Script to generate prediction BigWigs with a CREsted peak regression model
# USAGE: sbatch \
#--job-name=crested_predictions \
#--partition=carter-gpu \
#--account=carter-gpu \
#--gpus=a30:1 \
#--output slurm_logs/%x.%A.%a.out \
#--mem=32G \
#-n 4 \
#-t 01-00:00:00 \
#--array=1-N%N \
#crested_predictions.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Set-up env
source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
# TF needs pip-installed nvidia libs on LD_LIBRARY_PATH
NVIDIA_LIBS=$VIRTUAL_ENV/lib/python3.11/site-packages/nvidia
export LD_LIBRARY_PATH=$NVIDIA_LIBS/cublas/lib:$NVIDIA_LIBS/cuda_runtime/lib:$NVIDIA_LIBS/cudnn/lib:$NVIDIA_LIBS/cufft/lib:$NVIDIA_LIBS/curand/lib:$NVIDIA_LIBS/cusolver/lib:$NVIDIA_LIBS/cusparse/lib:$NVIDIA_LIBS/nccl/lib:$NVIDIA_LIBS/nvjitlink/lib:$LD_LIBRARY_PATH

# Cell type arrays
celltypes=(
	DE
	ENP_phase1
	FB_FLT1
	PFG1
	PFG2
	PGT1
	PGT2
	PGT3
	PP1
	PP2
	SC_delta_GHRL
	early_ENP
	early_SC_EC
	early_SC_alpha
	early_SC_beta
	exocrine
	late_ENP
	late_SC_EC
	late_SC_alpha
	late_SC_beta
	liver
	proliferating_endocrine
)
output_dirs=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/DE
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/ENP_phase1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/FB_FLT1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/PFG1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/PFG2
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/PGT1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/PGT2
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/PGT3
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/PP1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/PP2
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/SC_delta_GHRL
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/early_ENP
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/early_SC_EC
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/early_SC_alpha
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/early_SC_beta
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/exocrine
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/late_ENP
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/late_SC_EC
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/late_SC_alpha
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/late_SC_beta
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/liver
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/predictions/proliferating_endocrine
)

# Grab each for this SLURM task
celltype=${celltypes[$SLURM_ARRAY_TASK_ID - 1]}
output_dir=${output_dirs[$SLURM_ARRAY_TASK_ID - 1]}

# Echo info
echo -e "Cell type: $celltype"
echo -e "Output directory: $output_dir\n"

# Make the output directory
mkdir -p $output_dir

# Run cmd
cmd="python /cellar/users/aklie/opt/crested/scripts/compute_predictions.py \
--checkpoint $(cat /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/evaluation/best_checkpoint.txt) \
--bigwigs_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/bigwigs \
--regions_file /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/consensus_peaks.bed \
--genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
--chrom_sizes /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
--output_dir $output_dir \
--splits_file /cellar/users/aklie/data/ref/genomes/hg38/chrombpnet/splits/fold_0.json \
--cell_types $celltype \
--backend tensorflow \
--batch_size 64"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Date
date
