#! /bin/bash

#####
# Script to run motif clustering on de novo motifs
# USAGE: sbatch \
#--job-name=motif_clustering \
#--partition=carter-compute \
#--output=slurm_logs/%x.%A.out \
#--mem=128G \
#-n 4 \
#-t 02-00:00:00 \
#motif_clustering.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configure environment
source activate test_celloracle
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/envs/test_celloracle/lib/

# Define input files
paths_pfm=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/DE/average/motifs/pfms/counts/DE.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/ENP_phase1/average/motifs/pfms/counts/ENP_phase1.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/FB_FLT1/average/motifs/pfms/counts/FB_FLT1.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG1/average/motifs/pfms/counts/PFG1.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PFG2/average/motifs/pfms/counts/PFG2.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT1/average/motifs/pfms/counts/PGT1.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT2/average/motifs/pfms/counts/PGT2.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PGT3/average/motifs/pfms/counts/PGT3.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP1/average/motifs/pfms/counts/PP1.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/PP2/average/motifs/pfms/counts/PP2.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/SC_delta_GHRL/average/motifs/pfms/counts/SC_delta_GHRL.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_ENP/average/motifs/pfms/counts/early_ENP.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_EC/average/motifs/pfms/counts/early_SC_EC.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_alpha/average/motifs/pfms/counts/early_SC_alpha.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/early_SC_beta/average/motifs/pfms/counts/early_SC_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/exocrine/average/motifs/pfms/counts/exocrine.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_ENP/average/motifs/pfms/counts/late_ENP.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_EC/average/motifs/pfms/counts/late_SC_EC.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_alpha/average/motifs/pfms/counts/late_SC_alpha.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/late_SC_beta/average/motifs/pfms/counts/late_SC_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/liver/average/motifs/pfms/counts/liver.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/proliferating_endocrine/average/motifs/pfms/counts/proliferating_endocrine.counts.pfm
)
path_out=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/motifs/gimme_cluster
t=0.999

# Make the directory
mkdir -p $path_out/cluster
mkdir -p $path_out/meme
mkdir -p $path_out/tomtom

cat ${paths_pfm[@]} > $path_out/input_motifs.pfm

# Run cluster
cmd="gimme cluster \
$path_out/input_motifs.pfm \
$path_out/cluster \
-t $t"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Run conversion to meme
path_meme_script=/cellar/users/aklie/projects/ML4GLand/chrombpnet/scripts/pfm_to_meme.py
cmd="python $path_meme_script \
-i $path_out/cluster/clustered_motifs.pfm \
-o $path_out/meme"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Run TomTom for motif annotation
path_motif_database=/cellar/users/aklie/data/ref/motifs/jvierstra/motif-clustering-v2.0beta/motifs.meme
base_path=$path_out/meme
meme_files=($(find $base_path -name "*.meme" | sort))
for meme_file in "${meme_files[@]}"; do
  meme_id=$(basename $meme_file .meme)
  cmd="tomtom \
-no-ssc \
-oc $path_out/tomtom/$meme_id \
-verbosity 1 \
-text \
-min-overlap 5 \
-mi 1 \
-dist pearson \
-evalue \
-thresh 10.0 \
$meme_file \
$path_motif_database > $path_out/tomtom/${meme_id}.tomtom.txt"
  echo -e "Running command:\n$cmd\n"
  eval $cmd
done

# extract top TF matches
for x in $path_out/tomtom/*.tomtom.txt ; do cat $x | head -2 | tail -1 ; done | cut -f1,2,5 | sort -k3g > $path_out/tfs_initial.txt

# combine all meme files into one
combined_meme_path="$path_out/meme/combined.meme"
rm -f "$combined_meme_path"
(
    awk '/^MOTIF/ {exit} {print}' "$(find $path_out/meme -name '*.meme' ! -name 'combined.meme' | head -n1)"
    for f in $path_out/meme/*.meme; do
        [ "$(basename "$f")" = "combined.meme" ] && continue
        awk '/^MOTIF/ {p=1} p' "$f"
    done
) > "$combined_meme_path"

# Rename motifs in combined MEME using best TomTom hits
awk '
NR==FNR {split($0, a, "\t"); map[a[1]]=a[2]; next}
/^MOTIF/ {split($0, a, " "); if (a[2] in map) {print "MOTIF " map[a[2]]; next}}
{print}
' $path_out/tfs_initial.txt $combined_meme_path > $path_out/meme/combined_mapped.meme

# Completion message
date
