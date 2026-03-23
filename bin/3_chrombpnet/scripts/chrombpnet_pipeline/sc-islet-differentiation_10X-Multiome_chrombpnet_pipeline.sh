#! /bin/bash

#####
# Script to run chrombpnet model pipeline across multiple fragment files as a SLURM array job
# USAGE: sbatch \
#--job-name=chrombpnet_pipeline \
#--account carter-gpu \
#--partition carter-gpu \
#--gpus=a30:1 \
#--output slurm_logs/%x.%A.%a.out \
#--mem=128G \
#-n 4 \
#-t 02-00:00:00 \
#--array=1-12%12 \
#chrombpnet_pipeline.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Set-up env
source activate chrombpnet
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/lib
python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"

# file lists
beta=0.5
celltypes=(
	D45_SC_delta
	D4_DE_LEFTY1+
	D45_early_SC_EC_beta
	D7_DE_SOX17+
	D4_DE_POLG2+
	D22_SC_beta
	D4_DE_MitoHi
	D22_late_SC_alpha
	D7_GT_ONECUT1+
	D12_liver
	D15_pre_EC
	D45_SC_EC
	D15_pre_delta
	D4_DE_ONECUT2+
	D7_FLT1+
	D7_NECTIN3-AS1+
	D7_GT_CXCR4+
	D12_PP_ENP
	D45_early_others
	D15_PP_ERBB4+
	D22_late_ENP
	D22_late_SC_EC_beta
	D22_ENP_OCA2+
	D22_pre_alpha
	D9_ENP
	D45_early_SC_beta
	D15_pre_EC_beta
	D45_liver
	D15_FLT1+
	D12_ENP_SCG2+
	D4_DE_ERBB4+
	D12_PP_ERBB4+
	D22_late_SC_EC
	D9_NECTIN3-AS1+
	D45_SC_EC_beta
	D9_LINC01924+
	D4_DE_CDH8+
	D22_late_others
	D7_liver
	D45_SC_alpha1
	D15_PP_CREB5+
	D9_PFG_PLOG2+
	D12_PP_CREB5+
	D9_PFG_proliferating
	D22_SC_EC
	D7_GT_POLG2
	D45_ENP_EC
	D12_ENP_ARX+
	D9_DE_OTX2+
	D9_PFG_TTR+
	D7_GT_proliferating
	D45_early_SC_alpha
	D22_late_SC_beta
	D15_liver
	D15_pre_alpha
	D15_pre_beta
	D22_SC_EC_beta
	D9_PFG_CREB5+
	D22_liver
	D45_SC_alpha2
	D15_ENP_OCA2+
	D9_liver
	D22_PP_ERBB4+
	D22_SC_alpha
	D22_SC_delta+
	D4_DE_SOX4+
	D45_early_SC_EC
	D45_HSPA6+
	D12_ENP_OCA2+
	D9_PFG_NRG3+
	D15_PP_ENP
	D12_pre_EC_beta
	D45_SC_beta
	D9_PFG_CDH18+
)
fragments=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_SC_delta/merged_fragments_D45_SC_delta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D4_DE_LEFTY1+/merged_fragments_D4_DE_LEFTY1+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_early_SC_EC_beta/merged_fragments_D45_early_SC_EC_beta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D7_DE_SOX17+/merged_fragments_D7_DE_SOX17+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D4_DE_POLG2+/merged_fragments_D4_DE_POLG2+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_SC_beta/merged_fragments_D22_SC_beta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D4_DE_MitoHi/merged_fragments_D4_DE_MitoHi_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_late_SC_alpha/merged_fragments_D22_late_SC_alpha_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D7_GT_ONECUT1+/merged_fragments_D7_GT_ONECUT1+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D12_liver/merged_fragments_D12_liver_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D15_pre_EC/merged_fragments_D15_pre_EC_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_SC_EC/merged_fragments_D45_SC_EC_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D15_pre_delta/merged_fragments_D15_pre_delta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D4_DE_ONECUT2+/merged_fragments_D4_DE_ONECUT2+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D7_FLT1+/merged_fragments_D7_FLT1+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D7_NECTIN3-AS1+/merged_fragments_D7_NECTIN3-AS1+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D7_GT_CXCR4+/merged_fragments_D7_GT_CXCR4+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D12_PP_ENP/merged_fragments_D12_PP_ENP_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_early_others/merged_fragments_D45_early_others_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D15_PP_ERBB4+/merged_fragments_D15_PP_ERBB4+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_late_ENP/merged_fragments_D22_late_ENP_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_late_SC_EC_beta/merged_fragments_D22_late_SC_EC_beta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_ENP_OCA2+/merged_fragments_D22_ENP_OCA2+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_pre_alpha/merged_fragments_D22_pre_alpha_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D9_ENP/merged_fragments_D9_ENP_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_early_SC_beta/merged_fragments_D45_early_SC_beta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D15_pre_EC_beta/merged_fragments_D15_pre_EC_beta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_liver/merged_fragments_D45_liver_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D15_FLT1+/merged_fragments_D15_FLT1+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D12_ENP_SCG2+/merged_fragments_D12_ENP_SCG2+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D4_DE_ERBB4+/merged_fragments_D4_DE_ERBB4+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D12_PP_ERBB4+/merged_fragments_D12_PP_ERBB4+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_late_SC_EC/merged_fragments_D22_late_SC_EC_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D9_NECTIN3-AS1+/merged_fragments_D9_NECTIN3-AS1+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_SC_EC_beta/merged_fragments_D45_SC_EC_beta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D9_LINC01924+/merged_fragments_D9_LINC01924+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D4_DE_CDH8+/merged_fragments_D4_DE_CDH8+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_late_others/merged_fragments_D22_late_others_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D7_liver/merged_fragments_D7_liver_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_SC_alpha1/merged_fragments_D45_SC_alpha1_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D15_PP_CREB5+/merged_fragments_D15_PP_CREB5+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D9_PFG_PLOG2+/merged_fragments_D9_PFG_PLOG2+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D12_PP_CREB5+/merged_fragments_D12_PP_CREB5+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D9_PFG_proliferating/merged_fragments_D9_PFG_proliferating_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_SC_EC/merged_fragments_D22_SC_EC_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D7_GT_POLG2/merged_fragments_D7_GT_POLG2_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_ENP_EC/merged_fragments_D45_ENP_EC_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D12_ENP_ARX+/merged_fragments_D12_ENP_ARX+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D9_DE_OTX2+/merged_fragments_D9_DE_OTX2+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D9_PFG_TTR+/merged_fragments_D9_PFG_TTR+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D7_GT_proliferating/merged_fragments_D7_GT_proliferating_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_early_SC_alpha/merged_fragments_D45_early_SC_alpha_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_late_SC_beta/merged_fragments_D22_late_SC_beta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D15_liver/merged_fragments_D15_liver_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D15_pre_alpha/merged_fragments_D15_pre_alpha_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D15_pre_beta/merged_fragments_D15_pre_beta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_SC_EC_beta/merged_fragments_D22_SC_EC_beta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D9_PFG_CREB5+/merged_fragments_D9_PFG_CREB5+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_liver/merged_fragments_D22_liver_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_SC_alpha2/merged_fragments_D45_SC_alpha2_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D15_ENP_OCA2+/merged_fragments_D15_ENP_OCA2+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D9_liver/merged_fragments_D9_liver_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_PP_ERBB4+/merged_fragments_D22_PP_ERBB4+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_SC_alpha/merged_fragments_D22_SC_alpha_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D22_SC_delta+/merged_fragments_D22_SC_delta+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D4_DE_SOX4+/merged_fragments_D4_DE_SOX4+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_early_SC_EC/merged_fragments_D45_early_SC_EC_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_HSPA6+/merged_fragments_D45_HSPA6+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D12_ENP_OCA2+/merged_fragments_D12_ENP_OCA2+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D9_PFG_NRG3+/merged_fragments_D9_PFG_NRG3+_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D15_PP_ENP/merged_fragments_D15_PP_ENP_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D12_pre_EC_beta/merged_fragments_D12_pre_EC_beta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D45_SC_beta/merged_fragments_D45_SC_beta_sorted.tsv
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/fragments/D9_PFG_CDH18+/merged_fragments_D9_PFG_CDH18+_sorted.tsv
)
peaks=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_SC_delta/D45_SC_delta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D4_DE_LEFTY1+/D4_DE_LEFTY1+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_early_SC_EC_beta/D45_early_SC_EC_beta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D7_DE_SOX17+/D7_DE_SOX17+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D4_DE_POLG2+/D4_DE_POLG2+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_SC_beta/D22_SC_beta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D4_DE_MitoHi/D4_DE_MitoHi_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_late_SC_alpha/D22_late_SC_alpha_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D7_GT_ONECUT1+/D7_GT_ONECUT1+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D12_liver/D12_liver_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D15_pre_EC/D15_pre_EC_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_SC_EC/D45_SC_EC_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D15_pre_delta/D15_pre_delta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D4_DE_ONECUT2+/D4_DE_ONECUT2+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D7_FLT1+/D7_FLT1+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D7_NECTIN3-AS1+/D7_NECTIN3-AS1+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D7_GT_CXCR4+/D7_GT_CXCR4+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D12_PP_ENP/D12_PP_ENP_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_early_others/D45_early_others_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D15_PP_ERBB4+/D15_PP_ERBB4+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_late_ENP/D22_late_ENP_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_late_SC_EC_beta/D22_late_SC_EC_beta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_ENP_OCA2+/D22_ENP_OCA2+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_pre_alpha/D22_pre_alpha_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D9_ENP/D9_ENP_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_early_SC_beta/D45_early_SC_beta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D15_pre_EC_beta/D15_pre_EC_beta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_liver/D45_liver_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D15_FLT1+/D15_FLT1+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D12_ENP_SCG2+/D12_ENP_SCG2+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D4_DE_ERBB4+/D4_DE_ERBB4+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D12_PP_ERBB4+/D12_PP_ERBB4+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_late_SC_EC/D22_late_SC_EC_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D9_NECTIN3-AS1+/D9_NECTIN3-AS1+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_SC_EC_beta/D45_SC_EC_beta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D9_LINC01924+/D9_LINC01924+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D4_DE_CDH8+/D4_DE_CDH8+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_late_others/D22_late_others_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D7_liver/D7_liver_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_SC_alpha1/D45_SC_alpha1_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D15_PP_CREB5+/D15_PP_CREB5+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D9_PFG_PLOG2+/D9_PFG_PLOG2+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D12_PP_CREB5+/D12_PP_CREB5+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D9_PFG_proliferating/D9_PFG_proliferating_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_SC_EC/D22_SC_EC_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D7_GT_POLG2/D7_GT_POLG2_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_ENP_EC/D45_ENP_EC_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D12_ENP_ARX+/D12_ENP_ARX+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D9_DE_OTX2+/D9_DE_OTX2+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D9_PFG_TTR+/D9_PFG_TTR+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D7_GT_proliferating/D7_GT_proliferating_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_early_SC_alpha/D45_early_SC_alpha_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_late_SC_beta/D22_late_SC_beta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D15_liver/D15_liver_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D15_pre_alpha/D15_pre_alpha_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D15_pre_beta/D15_pre_beta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_SC_EC_beta/D22_SC_EC_beta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D9_PFG_CREB5+/D9_PFG_CREB5+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_liver/D22_liver_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_SC_alpha2/D45_SC_alpha2_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D15_ENP_OCA2+/D15_ENP_OCA2+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D9_liver/D9_liver_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_PP_ERBB4+/D22_PP_ERBB4+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_SC_alpha/D22_SC_alpha_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D22_SC_delta+/D22_SC_delta+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D4_DE_SOX4+/D4_DE_SOX4+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_early_SC_EC/D45_early_SC_EC_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_HSPA6+/D45_HSPA6+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D12_ENP_OCA2+/D12_ENP_OCA2+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D9_PFG_NRG3+/D9_PFG_NRG3+_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D15_PP_ENP/D15_PP_ENP_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D12_pre_EC_beta/D12_pre_EC_beta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D45_SC_beta/D45_SC_beta_peaks.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/peaks/D9_PFG_CDH18+/D9_PFG_CDH18+_peaks.bed
)
folds=(
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
	0
)
negatives=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_delta/fold_0/negatives/D45_SC_delta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_LEFTY1+/fold_0/negatives/D4_DE_LEFTY1+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC_beta/fold_0/negatives/D45_early_SC_EC_beta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_DE_SOX17+/fold_0/negatives/D7_DE_SOX17+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_POLG2+/fold_0/negatives/D4_DE_POLG2+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_beta/fold_0/negatives/D22_SC_beta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_MitoHi/fold_0/negatives/D4_DE_MitoHi_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_alpha/fold_0/negatives/D22_late_SC_alpha_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_ONECUT1+/fold_0/negatives/D7_GT_ONECUT1+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_liver/fold_0/negatives/D12_liver_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC/fold_0/negatives/D15_pre_EC_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC/fold_0/negatives/D45_SC_EC_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_delta/fold_0/negatives/D15_pre_delta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ONECUT2+/fold_0/negatives/D4_DE_ONECUT2+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_FLT1+/fold_0/negatives/D7_FLT1+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_NECTIN3-AS1+/fold_0/negatives/D7_NECTIN3-AS1+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_CXCR4+/fold_0/negatives/D7_GT_CXCR4+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ENP/fold_0/negatives/D12_PP_ENP_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_others/fold_0/negatives/D45_early_others_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ERBB4+/fold_0/negatives/D15_PP_ERBB4+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_ENP/fold_0/negatives/D22_late_ENP_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC_beta/fold_0/negatives/D22_late_SC_EC_beta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_ENP_OCA2+/fold_0/negatives/D22_ENP_OCA2+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_pre_alpha/fold_0/negatives/D22_pre_alpha_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_ENP/fold_0/negatives/D9_ENP_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_beta/fold_0/negatives/D45_early_SC_beta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC_beta/fold_0/negatives/D15_pre_EC_beta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_liver/fold_0/negatives/D45_liver_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_FLT1+/fold_0/negatives/D15_FLT1+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_SCG2+/fold_0/negatives/D12_ENP_SCG2+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ERBB4+/fold_0/negatives/D4_DE_ERBB4+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ERBB4+/fold_0/negatives/D12_PP_ERBB4+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC/fold_0/negatives/D22_late_SC_EC_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_NECTIN3-AS1+/fold_0/negatives/D9_NECTIN3-AS1+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC_beta/fold_0/negatives/D45_SC_EC_beta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_LINC01924+/fold_0/negatives/D9_LINC01924+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_CDH8+/fold_0/negatives/D4_DE_CDH8+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_others/fold_0/negatives/D22_late_others_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_liver/fold_0/negatives/D7_liver_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha1/fold_0/negatives/D45_SC_alpha1_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_CREB5+/fold_0/negatives/D15_PP_CREB5+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_PLOG2+/fold_0/negatives/D9_PFG_PLOG2+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_CREB5+/fold_0/negatives/D12_PP_CREB5+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_proliferating/fold_0/negatives/D9_PFG_proliferating_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC/fold_0/negatives/D22_SC_EC_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_POLG2/fold_0/negatives/D7_GT_POLG2_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_ENP_EC/fold_0/negatives/D45_ENP_EC_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_ARX+/fold_0/negatives/D12_ENP_ARX+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_DE_OTX2+/fold_0/negatives/D9_DE_OTX2+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_TTR+/fold_0/negatives/D9_PFG_TTR+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_proliferating/fold_0/negatives/D7_GT_proliferating_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_alpha/fold_0/negatives/D45_early_SC_alpha_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_beta/fold_0/negatives/D22_late_SC_beta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_liver/fold_0/negatives/D15_liver_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_alpha/fold_0/negatives/D15_pre_alpha_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_beta/fold_0/negatives/D15_pre_beta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC_beta/fold_0/negatives/D22_SC_EC_beta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CREB5+/fold_0/negatives/D9_PFG_CREB5+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_liver/fold_0/negatives/D22_liver_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha2/fold_0/negatives/D45_SC_alpha2_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_ENP_OCA2+/fold_0/negatives/D15_ENP_OCA2+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_liver/fold_0/negatives/D9_liver_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_PP_ERBB4+/fold_0/negatives/D22_PP_ERBB4+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_alpha/fold_0/negatives/D22_SC_alpha_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_delta+/fold_0/negatives/D22_SC_delta+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_SOX4+/fold_0/negatives/D4_DE_SOX4+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC/fold_0/negatives/D45_early_SC_EC_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_HSPA6+/fold_0/negatives/D45_HSPA6+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_OCA2+/fold_0/negatives/D12_ENP_OCA2+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_NRG3+/fold_0/negatives/D9_PFG_NRG3+_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ENP/fold_0/negatives/D15_PP_ENP_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_pre_EC_beta/fold_0/negatives/D12_pre_EC_beta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_beta/fold_0/negatives/D45_SC_beta_negatives.bed
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CDH18+/fold_0/negatives/D9_PFG_CDH18+_negatives.bed
)
bias_models=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_delta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_LEFTY1+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC_beta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_DE_SOX17+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_POLG2+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_beta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_MitoHi/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_alpha/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_ONECUT1+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_liver/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_delta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ONECUT2+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_FLT1+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_NECTIN3-AS1+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_CXCR4+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ENP/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_others/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ERBB4+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_ENP/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC_beta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_ENP_OCA2+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_pre_alpha/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_ENP/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_beta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC_beta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_liver/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_FLT1+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_SCG2+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ERBB4+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ERBB4+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_NECTIN3-AS1+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC_beta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_LINC01924+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_CDH8+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_others/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_liver/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha1/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_CREB5+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_PLOG2+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_CREB5+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_proliferating/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_POLG2/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_ENP_EC/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_ARX+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_DE_OTX2+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_TTR+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_proliferating/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_alpha/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_beta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_liver/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_alpha/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_beta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC_beta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CREB5+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_liver/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha2/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_ENP_OCA2+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_liver/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_PP_ERBB4+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_alpha/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_delta+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_SOX4+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_HSPA6+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_OCA2+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_NRG3+/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ENP/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_pre_EC_beta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_beta/fold_0/bias_model/0.5/models/bias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CDH18+/fold_0/bias_model/0.5/models/bias.h5
)
output_dirs=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_delta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_LEFTY1+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC_beta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_DE_SOX17+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_POLG2+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_beta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_MitoHi/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_alpha/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_ONECUT1+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_liver/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_delta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ONECUT2+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_FLT1+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_NECTIN3-AS1+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_CXCR4+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ENP/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_others/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ERBB4+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_ENP/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC_beta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_ENP_OCA2+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_pre_alpha/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_ENP/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_beta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC_beta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_liver/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_FLT1+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_SCG2+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ERBB4+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ERBB4+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_NECTIN3-AS1+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC_beta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_LINC01924+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_CDH8+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_others/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_liver/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha1/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_CREB5+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_PLOG2+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_CREB5+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_proliferating/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_POLG2/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_ENP_EC/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_ARX+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_DE_OTX2+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_TTR+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_proliferating/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_alpha/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_beta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_liver/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_alpha/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_beta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC_beta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CREB5+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_liver/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha2/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_ENP_OCA2+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_liver/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_PP_ERBB4+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_alpha/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_delta+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_SOX4+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_HSPA6+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_OCA2+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_NRG3+/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ENP/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_pre_EC_beta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_beta/fold_0/chrombpnet/0.5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CDH18+/fold_0/chrombpnet/0.5
)

# Grab each for this SLURM task
celltype=${celltypes[$SLURM_ARRAY_TASK_ID - 1]}
fragment=${fragments[$SLURM_ARRAY_TASK_ID - 1]}
peak=${peaks[$SLURM_ARRAY_TASK_ID - 1]}
fold=${folds[$SLURM_ARRAY_TASK_ID - 1]}
negative=${negatives[$SLURM_ARRAY_TASK_ID - 1]}
bias_model=${bias_models[$SLURM_ARRAY_TASK_ID - 1]}
output_dir=${output_dirs[$SLURM_ARRAY_TASK_ID - 1]}

# echo the celltype and peak
echo -e "Celltype: $celltype"
echo -e "Fragment: $fragment"
echo -e "Peakset: $peak"
echo -e "Fold: $fold"
echo -e "Negatives: $negative"
echo -e "Bias model: $bias_model"
echo -e "Output directory: $output_dir\n"

# make the output directory
mkdir -p $output_dir

# Run cmd
cmd="chrombpnet pipeline \
-ifrag $fragment \
-d "ATAC" \
-g /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
-c /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
-p $peak \
-n $negative \
-fl /cellar/users/aklie/data/ref/genomes/hg38/chrombpnet/splits/fold_${fold}.json \
-b $bias_model \
-o $output_dir"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Date
date
