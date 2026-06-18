#! /bin/bash

#####
# Script to run de novo motif discovery on contributions
# USAGE: sbatch \
#--job-name=motif_discovery \
#--partition=carter-compute \
#--output=slurm_logs/%x.%A.%a.out \
#--array=1-12%12 \
#--mem=128G \
#-n 4 \
#-t 02-00:00:00 \
#motif_discovery.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Set up environment
source activate eugene_tools
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/
export LD_LIBRARY_PATH=/cm/shared/apps/slurm/current/lib64
NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK

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
output_dirs=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/DE/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/ENP_phase1/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/FB_FLT1/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG1/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG2/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT1/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT2/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT3/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP1/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP2/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/SC_delta_GHRL/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_ENP/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_EC/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_alpha/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/exocrine/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_ENP/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_EC/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_alpha/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/liver/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/proliferating_endocrine/average/motifs
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
n_seqlets=100000
leiden_res=2
window=400

# Make the fold the slurm task id - 1
celltype=${celltypes[$SLURM_ARRAY_TASK_ID - 1]}
input_h5=${input_h5s[$SLURM_ARRAY_TASK_ID - 1]}
output_dir=${output_dirs[$SLURM_ARRAY_TASK_ID - 1]}
type=${types[$SLURM_ARRAY_TASK_ID - 1]}

# Create output directory
mkdir -p $output_dir

# Run motifs cmd
cmd="modisco motifs \
-i $input_h5 \
-n $n_seqlets \
-o $output_dir/${celltype}.${type}.modisco.h5 \
-l $leiden_res \
-w $window \
-v"
echo -e $cmd
eval $cmd

# Run report cmd
cmd="modisco report \
-i $output_dir/${celltype}.${type}.modisco.h5 \
-o $output_dir/${celltype}.${type}.modisco.report/ \
-s $output_dir/${celltype}.${type}.modisco.report/ \
-m /cellar/users/aklie/data/ref/motifs/jvierstra/motif-clustering-v2.0beta/motifs.meme"
echo -e $cmd
eval $cmd

# Run pfm generation script
mkdir -p $output_dir/pfms/$type
path_pfm_script=/cellar/users/aklie/projects/ML4GLand/chrombpnet/scripts/modisco_to_pfm.py
cmd="python $path_pfm_script \
-m $output_dir/${celltype}.${type}.modisco.h5 \
-o $output_dir/pfms/$type \
-op ${celltype}.${type} \
-f 2"
echo -e $cmd
eval $cmd

# Consolidate pfms
temp_file="$output_dir/pfms/$type/${celltype}.${type}.tmp.pfm"
final_file="$output_dir/pfms/$type/${celltype}.${type}.pfm"

# Ensure final file does not exist
rm -f "$final_file"

# Concatenate files with better formatting
for x in "$output_dir/pfms/$type"/${celltype}.${type}*.pfm; do
  echo ">$(basename "$x")" >> "$temp_file"
  cat "$x" >> "$temp_file"
done

# Move to final destination safely
mv "$temp_file" "$final_file"

echo "Final consolidated PFM file: $final_file"

# Completion message
date
