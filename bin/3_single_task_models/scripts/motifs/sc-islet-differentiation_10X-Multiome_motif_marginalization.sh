#!/bin/bash

#####
# Script to run marginalization analysis with chromBPNet models
# USAGE:
# sbatch \
# --job-name=marginalize \
# --account carter-gpu \
# --partition carter-gpu \
# --gpus=a30:1 \
# --output slurm_logs/%x.%A.%a.out \
# --mem=32G \
# -n 4 \
# -t 02-00:00:00 \
# --array=1-20%10 \
# marginalize.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Set-up env
source activate eugene_tools
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/lib
python -c "import torch; print(torch.cuda.is_available())"

# File arrays
models=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/DE/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/ENP_phase1/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/FB_FLT1/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG1/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG2/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT1/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT2/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT3/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP1/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP2/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/SC_delta_GHRL/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_ENP/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_EC/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_alpha/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/exocrine/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_ENP/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_EC/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_alpha/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/liver/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/proliferating_endocrine/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
)
peaks=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/DE.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/ENP_phase1.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/FB_FLT1.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/PFG1.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/PFG2.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/PGT1.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/PGT2.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/PGT3.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/PP1.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/PP2.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/SC_delta_GHRL.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/early_ENP.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/early_SC_EC.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/early_SC_alpha.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/early_SC_beta.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/exocrine.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/late_ENP.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/late_SC_EC.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/late_SC_alpha.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/late_SC_beta.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/liver.narrowPeak
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/proliferating_endocrine.narrowPeak
)
output_dirs=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/DE/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/ENP_phase1/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/FB_FLT1/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG1/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG2/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT1/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT2/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT3/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP1/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP2/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/SC_delta_GHRL/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_ENP/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_EC/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_alpha/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/exocrine/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_ENP/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_EC/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_alpha/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/liver/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/proliferating_endocrine/fold_0/chrombpnet/0.5/marginalization/
)
n_seqs=100

# SLURM task-specific values
model=${models[$SLURM_ARRAY_TASK_ID - 1]}
peak=${peaks[$SLURM_ARRAY_TASK_ID - 1]}
output_dir=${output_dirs[$SLURM_ARRAY_TASK_ID - 1]}

# Constants
genome_fasta="/cellar/users/aklie/data/ref/genomes/hg38/hg38.fa"
motif_file="/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/motifs/gimme_cluster/meme/combined.meme"

echo -e "Model: $model"
echo -e "Peaks: $peak"
echo -e "Output: $output_dir\n"

mkdir -p $output_dir

# Run marginalization
cmd="python /cellar/users/aklie/projects/ML4GLand/chrombpnet/scripts/marginalization.py \
  -m $model \
  -f $motif_file \
  -g $genome_fasta \
  -p $peak \
  -o $output_dir \
  -n $n_seqs"

echo -e "Running command:\n$cmd\n"
eval $cmd

date
