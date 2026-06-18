#! /bin/bash

#####
# Script to run motif hit calling
# USAGE: sbatch \
#--job-name=motif_hits \
#--partition=carter-gpu \
#--account=carter-gpu \
#--gpus=a30:1 \
#--output=slurm_logs/%x.%A.%a.out \
#--mem=128G \
#-n 4 \
#-t 02-00:00:00 \
#--array=1-12%12 \
#motif_hits.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configure environment
source activate finemo_gpu
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/envs/finemo_gpu/lib/

# Define input files
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
input_h5s=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/DE/average/contributions/DE_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/ENP_phase1/average/contributions/ENP_phase1_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/FB_FLT1/average/contributions/FB_FLT1_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG1/average/contributions/PFG1_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG2/average/contributions/PFG2_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT1/average/contributions/PGT1_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT2/average/contributions/PGT2_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT3/average/contributions/PGT3_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP1/average/contributions/PP1_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP2/average/contributions/PP2_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/SC_delta_GHRL/average/contributions/SC_delta_GHRL_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_ENP/average/contributions/early_ENP_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_EC/average/contributions/early_SC_EC_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_alpha/average/contributions/early_SC_alpha_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_beta/average/contributions/early_SC_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/exocrine/average/contributions/exocrine_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_ENP/average/contributions/late_ENP_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_EC/average/contributions/late_SC_EC_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_alpha/average/contributions/late_SC_alpha_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_beta/average/contributions/late_SC_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/liver/average/contributions/liver_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/proliferating_endocrine/average/contributions/proliferating_endocrine_counts.h5
)
modiscos=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/DE/average/motifs/DE.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/ENP_phase1/average/motifs/ENP_phase1.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/FB_FLT1/average/motifs/FB_FLT1.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG1/average/motifs/PFG1.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG2/average/motifs/PFG2.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT1/average/motifs/PGT1.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT2/average/motifs/PGT2.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT3/average/motifs/PGT3.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP1/average/motifs/PP1.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP2/average/motifs/PP2.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/SC_delta_GHRL/average/motifs/SC_delta_GHRL.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_ENP/average/motifs/early_ENP.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_EC/average/motifs/early_SC_EC.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_alpha/average/motifs/early_SC_alpha.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_beta/average/motifs/early_SC_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/exocrine/average/motifs/exocrine.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_ENP/average/motifs/late_ENP.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_EC/average/motifs/late_SC_EC.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_alpha/average/motifs/late_SC_alpha.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_beta/average/motifs/late_SC_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/liver/average/motifs/liver.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/proliferating_endocrine/average/motifs/proliferating_endocrine.counts.modisco.h5
)
types=(
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
	counts
)
output_dirs=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/DE/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/ENP_phase1/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/FB_FLT1/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG1/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG2/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT1/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT2/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT3/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP1/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP2/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/SC_delta_GHRL/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_ENP/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_EC/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_alpha/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/exocrine/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_ENP/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_EC/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_alpha/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/liver/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/proliferating_endocrine/average/motifs/hits/counts
)
window=400

# Make the fold the slurm task id - 1
celltype=${celltypes[$SLURM_ARRAY_TASK_ID - 1]}
peak=${peaks[$SLURM_ARRAY_TASK_ID - 1]}
input_h5=${input_h5s[$SLURM_ARRAY_TASK_ID - 1]}
modisco=${modiscos[$SLURM_ARRAY_TASK_ID - 1]}
type=${types[$SLURM_ARRAY_TASK_ID - 1]}
output_dir=${output_dirs[$SLURM_ARRAY_TASK_ID - 1]}

# Create output directory
mkdir -p $output_dir

# Run cmd to preprocess
cmd="finemo extract-regions-chrombpnet-h5 \
-c $input_h5 \
-o $output_dir/${celltype}.${type}.finemo.npz \
-w $window"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Run cmd to call hits
cmd="finemo call-hits \
-r $output_dir/${celltype}.${type}.finemo.npz \
-m $modisco \
-o $output_dir \
-p $peak"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Run cmd for report
cmd="finemo report \
-r $output_dir/${celltype}.${type}.finemo.npz \
-p $peak \
-H $output_dir/hits.tsv \
-m $modisco \
-o $output_dir/report \
-W $window"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Completion message
date