#! /bin/bash

#####
# Script to run motif clustering on de novo motifs
# USAGE: sbatch \
#--job-name=crested_motif_clustering \
#--partition=carter-compute \
#--output=slurm_logs/%x.%A.out \
#--mem=128G \
#-n 4 \
#-t 02-00:00:00 \
#crested_motif_clustering.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configure environment. gimme lives in test_celloracle; tomtom ships with the
# chrombpnet env — invoke tomtom via absolute path below.
source /cellar/users/aklie/opt/miniconda3/etc/profile.d/conda.sh
conda activate test_celloracle
TOMTOM=/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/bin/tomtom

# Define input files
paths_pfm=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/DE/pfms/DE.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/ENP_phase1/pfms/ENP_phase1.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/FB_FLT1/pfms/FB_FLT1.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PFG1/pfms/PFG1.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PFG2/pfms/PFG2.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PGT1/pfms/PGT1.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PGT2/pfms/PGT2.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PGT3/pfms/PGT3.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PP1/pfms/PP1.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/PP2/pfms/PP2.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/SC_delta_GHRL/pfms/SC_delta_GHRL.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/early_ENP/pfms/early_ENP.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/early_SC_EC/pfms/early_SC_EC.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/early_SC_alpha/pfms/early_SC_alpha.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/early_SC_beta/pfms/early_SC_beta.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/exocrine/pfms/exocrine.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/late_ENP/pfms/late_ENP.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/late_SC_EC/pfms/late_SC_EC.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/late_SC_alpha/pfms/late_SC_alpha.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/late_SC_beta/pfms/late_SC_beta.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/liver/pfms/liver.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/modisco/proliferating_endocrine/pfms/proliferating_endocrine.pfm
)
path_out=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model/crested/motifs
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
  cmd="$TOMTOM \
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
' $path_out/tfs_initial.txt $combined_meme_path > $path_out/combined_mapped.meme

# Completion message
date
