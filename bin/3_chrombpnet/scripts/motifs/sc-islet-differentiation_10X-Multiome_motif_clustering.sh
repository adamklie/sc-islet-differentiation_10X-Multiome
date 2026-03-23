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
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_delta/average/motifs/pfms/counts/D45_SC_delta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_LEFTY1+/average/motifs/pfms/counts/D4_DE_LEFTY1+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC_beta/average/motifs/pfms/counts/D45_early_SC_EC_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_DE_SOX17+/average/motifs/pfms/counts/D7_DE_SOX17+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_POLG2+/average/motifs/pfms/counts/D4_DE_POLG2+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_beta/average/motifs/pfms/counts/D22_SC_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_MitoHi/average/motifs/pfms/counts/D4_DE_MitoHi.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_alpha/average/motifs/pfms/counts/D22_late_SC_alpha.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_ONECUT1+/average/motifs/pfms/counts/D7_GT_ONECUT1+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_liver/average/motifs/pfms/counts/D12_liver.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC/average/motifs/pfms/counts/D15_pre_EC.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC/average/motifs/pfms/counts/D45_SC_EC.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_delta/average/motifs/pfms/counts/D15_pre_delta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ONECUT2+/average/motifs/pfms/counts/D4_DE_ONECUT2+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_FLT1+/average/motifs/pfms/counts/D7_FLT1+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_NECTIN3-AS1+/average/motifs/pfms/counts/D7_NECTIN3-AS1+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_CXCR4+/average/motifs/pfms/counts/D7_GT_CXCR4+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ENP/average/motifs/pfms/counts/D12_PP_ENP.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_others/average/motifs/pfms/counts/D45_early_others.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ERBB4+/average/motifs/pfms/counts/D15_PP_ERBB4+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_ENP/average/motifs/pfms/counts/D22_late_ENP.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC_beta/average/motifs/pfms/counts/D22_late_SC_EC_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_ENP_OCA2+/average/motifs/pfms/counts/D22_ENP_OCA2+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_pre_alpha/average/motifs/pfms/counts/D22_pre_alpha.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_ENP/average/motifs/pfms/counts/D9_ENP.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_beta/average/motifs/pfms/counts/D45_early_SC_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC_beta/average/motifs/pfms/counts/D15_pre_EC_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_liver/average/motifs/pfms/counts/D45_liver.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_FLT1+/average/motifs/pfms/counts/D15_FLT1+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_SCG2+/average/motifs/pfms/counts/D12_ENP_SCG2+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ERBB4+/average/motifs/pfms/counts/D4_DE_ERBB4+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ERBB4+/average/motifs/pfms/counts/D12_PP_ERBB4+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC/average/motifs/pfms/counts/D22_late_SC_EC.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_NECTIN3-AS1+/average/motifs/pfms/counts/D9_NECTIN3-AS1+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC_beta/average/motifs/pfms/counts/D45_SC_EC_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_LINC01924+/average/motifs/pfms/counts/D9_LINC01924+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_CDH8+/average/motifs/pfms/counts/D4_DE_CDH8+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_others/average/motifs/pfms/counts/D22_late_others.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_liver/average/motifs/pfms/counts/D7_liver.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha1/average/motifs/pfms/counts/D45_SC_alpha1.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_CREB5+/average/motifs/pfms/counts/D15_PP_CREB5+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_PLOG2+/average/motifs/pfms/counts/D9_PFG_PLOG2+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_CREB5+/average/motifs/pfms/counts/D12_PP_CREB5+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_proliferating/average/motifs/pfms/counts/D9_PFG_proliferating.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC/average/motifs/pfms/counts/D22_SC_EC.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_POLG2/average/motifs/pfms/counts/D7_GT_POLG2.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_ENP_EC/average/motifs/pfms/counts/D45_ENP_EC.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_ARX+/average/motifs/pfms/counts/D12_ENP_ARX+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_DE_OTX2+/average/motifs/pfms/counts/D9_DE_OTX2+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_TTR+/average/motifs/pfms/counts/D9_PFG_TTR+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_proliferating/average/motifs/pfms/counts/D7_GT_proliferating.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_alpha/average/motifs/pfms/counts/D45_early_SC_alpha.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_beta/average/motifs/pfms/counts/D22_late_SC_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_liver/average/motifs/pfms/counts/D15_liver.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_alpha/average/motifs/pfms/counts/D15_pre_alpha.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_beta/average/motifs/pfms/counts/D15_pre_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC_beta/average/motifs/pfms/counts/D22_SC_EC_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CREB5+/average/motifs/pfms/counts/D9_PFG_CREB5+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_liver/average/motifs/pfms/counts/D22_liver.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha2/average/motifs/pfms/counts/D45_SC_alpha2.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_ENP_OCA2+/average/motifs/pfms/counts/D15_ENP_OCA2+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_liver/average/motifs/pfms/counts/D9_liver.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_PP_ERBB4+/average/motifs/pfms/counts/D22_PP_ERBB4+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_alpha/average/motifs/pfms/counts/D22_SC_alpha.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_delta+/average/motifs/pfms/counts/D22_SC_delta+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_SOX4+/average/motifs/pfms/counts/D4_DE_SOX4+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC/average/motifs/pfms/counts/D45_early_SC_EC.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_HSPA6+/average/motifs/pfms/counts/D45_HSPA6+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_OCA2+/average/motifs/pfms/counts/D12_ENP_OCA2+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_NRG3+/average/motifs/pfms/counts/D9_PFG_NRG3+.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ENP/average/motifs/pfms/counts/D15_PP_ENP.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_pre_EC_beta/average/motifs/pfms/counts/D12_pre_EC_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_beta/average/motifs/pfms/counts/D45_SC_beta.counts.pfm
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CDH18+/average/motifs/pfms/counts/D9_PFG_CDH18+.counts.pfm
)
path_out=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/motifs
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

#
for x in $path_out/tomtom/*.tomtom.txt ; do cat $x | head -2 | tail -1 ; done | cut -f1,2,5 | sort -k3g > $path_out/tfs_initial.txt

# Completion message
date
