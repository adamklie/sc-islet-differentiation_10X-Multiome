#! /bin/bash

#####
# Script to run de novo motif discovery on contribution scores
# USAGE: sbatch \
#--job-name=crested_motif_discovery \
#--partition=carter-compute \
#--output=slurm_logs/%x.%A.%a.out \
#--array=1-N%N \
#--mem=128G \
#-n 4 \
#-t 02-00:00:00 \
#crested_motif_discovery.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Set up environment
source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
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
input_dirs=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/DE/DE
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/ENP_phase1/ENP_phase1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/FB_FLT1/FB_FLT1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/PFG1/PFG1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/PFG2/PFG2
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/PGT1/PGT1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/PGT2/PGT2
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/PGT3/PGT3
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/PP1/PP1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/PP2/PP2
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/SC_delta_GHRL/SC_delta_GHRL
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/early_ENP/early_ENP
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/early_SC_EC/early_SC_EC
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/early_SC_alpha/early_SC_alpha
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/early_SC_beta/early_SC_beta
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/exocrine/exocrine
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/late_ENP/late_ENP
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/late_SC_EC/late_SC_EC
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/late_SC_alpha/late_SC_alpha
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/late_SC_beta/late_SC_beta
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/liver/liver
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/contributions/proliferating_endocrine/proliferating_endocrine
)
output_dirs=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/DE
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/ENP_phase1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/FB_FLT1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PFG1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PFG2
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PGT1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PGT2
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PGT3
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PP1
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PP2
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/SC_delta_GHRL
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/early_ENP
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/early_SC_EC
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/early_SC_alpha
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/early_SC_beta
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/exocrine
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/late_ENP
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/late_SC_EC
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/late_SC_alpha
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/late_SC_beta
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/liver
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/proliferating_endocrine
)
n_seqlets=50000
leiden_res=2
window=400

# Grab each for this SLURM task
celltype=${celltypes[$SLURM_ARRAY_TASK_ID - 1]}
input_dir=${input_dirs[$SLURM_ARRAY_TASK_ID - 1]}
output_dir=${output_dirs[$SLURM_ARRAY_TASK_ID - 1]}

# Echo info
echo -e "Cell type: $celltype"
echo -e "Input directory: $input_dir"
echo -e "Output directory: $output_dir\n"

# Create output directory
mkdir -p $output_dir

# Run modisco motifs
cmd="modisco motifs \
-s $input_dir/one_hot_seqs.npz \
-a $input_dir/projected_contribs.npz \
-n $n_seqlets \
-o $output_dir/${celltype}.modisco.h5 \
-l $leiden_res \
-w $window \
-v"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Run report
cmd="modisco report \
-i $output_dir/${celltype}.modisco.h5 \
-o $output_dir/${celltype}.modisco.report/ \
-s $output_dir/${celltype}.modisco.report/ \
-m /cellar/users/aklie/data/ref/motifs/jvierstra/motif-clustering-v2.0beta/motifs.meme"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Run pfm generation
mkdir -p $output_dir/pfms
path_pfm_script=/cellar/users/aklie/projects/ML4GLand/chrombpnet/scripts/modisco_to_pfm.py
cmd="python $path_pfm_script \
-m $output_dir/${celltype}.modisco.h5 \
-o $output_dir/pfms \
-op $celltype \
-f 2"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Consolidate pfms
temp_file="$output_dir/pfms/${celltype}.tmp.pfm"
final_file="$output_dir/pfms/${celltype}.pfm"
rm -f "$final_file"
for x in "$output_dir/pfms"/${celltype}*.pfm; do
  echo ">$(basename "$x")" >> "$temp_file"
  cat "$x" >> "$temp_file"
done
mv "$temp_file" "$final_file"
echo "Final consolidated PFM file: $final_file"

# Completion message
date
