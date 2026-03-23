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
input_h5s=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_delta/average/contributions/D45_SC_delta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_LEFTY1+/average/contributions/D4_DE_LEFTY1+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC_beta/average/contributions/D45_early_SC_EC_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_DE_SOX17+/average/contributions/D7_DE_SOX17+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_POLG2+/average/contributions/D4_DE_POLG2+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_beta/average/contributions/D22_SC_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_MitoHi/average/contributions/D4_DE_MitoHi_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_alpha/average/contributions/D22_late_SC_alpha_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_ONECUT1+/average/contributions/D7_GT_ONECUT1+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_liver/average/contributions/D12_liver_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC/average/contributions/D15_pre_EC_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC/average/contributions/D45_SC_EC_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_delta/average/contributions/D15_pre_delta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ONECUT2+/average/contributions/D4_DE_ONECUT2+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_FLT1+/average/contributions/D7_FLT1+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_NECTIN3-AS1+/average/contributions/D7_NECTIN3-AS1+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_CXCR4+/average/contributions/D7_GT_CXCR4+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ENP/average/contributions/D12_PP_ENP_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_others/average/contributions/D45_early_others_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ERBB4+/average/contributions/D15_PP_ERBB4+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_ENP/average/contributions/D22_late_ENP_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC_beta/average/contributions/D22_late_SC_EC_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_ENP_OCA2+/average/contributions/D22_ENP_OCA2+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_pre_alpha/average/contributions/D22_pre_alpha_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_ENP/average/contributions/D9_ENP_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_beta/average/contributions/D45_early_SC_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC_beta/average/contributions/D15_pre_EC_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_liver/average/contributions/D45_liver_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_FLT1+/average/contributions/D15_FLT1+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_SCG2+/average/contributions/D12_ENP_SCG2+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ERBB4+/average/contributions/D4_DE_ERBB4+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ERBB4+/average/contributions/D12_PP_ERBB4+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC/average/contributions/D22_late_SC_EC_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_NECTIN3-AS1+/average/contributions/D9_NECTIN3-AS1+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC_beta/average/contributions/D45_SC_EC_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_LINC01924+/average/contributions/D9_LINC01924+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_CDH8+/average/contributions/D4_DE_CDH8+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_others/average/contributions/D22_late_others_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_liver/average/contributions/D7_liver_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha1/average/contributions/D45_SC_alpha1_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_CREB5+/average/contributions/D15_PP_CREB5+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_PLOG2+/average/contributions/D9_PFG_PLOG2+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_CREB5+/average/contributions/D12_PP_CREB5+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_proliferating/average/contributions/D9_PFG_proliferating_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC/average/contributions/D22_SC_EC_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_POLG2/average/contributions/D7_GT_POLG2_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_ENP_EC/average/contributions/D45_ENP_EC_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_ARX+/average/contributions/D12_ENP_ARX+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_DE_OTX2+/average/contributions/D9_DE_OTX2+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_TTR+/average/contributions/D9_PFG_TTR+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_proliferating/average/contributions/D7_GT_proliferating_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_alpha/average/contributions/D45_early_SC_alpha_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_beta/average/contributions/D22_late_SC_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_liver/average/contributions/D15_liver_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_alpha/average/contributions/D15_pre_alpha_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_beta/average/contributions/D15_pre_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC_beta/average/contributions/D22_SC_EC_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CREB5+/average/contributions/D9_PFG_CREB5+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_liver/average/contributions/D22_liver_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha2/average/contributions/D45_SC_alpha2_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_ENP_OCA2+/average/contributions/D15_ENP_OCA2+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_liver/average/contributions/D9_liver_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_PP_ERBB4+/average/contributions/D22_PP_ERBB4+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_alpha/average/contributions/D22_SC_alpha_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_delta+/average/contributions/D22_SC_delta+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_SOX4+/average/contributions/D4_DE_SOX4+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC/average/contributions/D45_early_SC_EC_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_HSPA6+/average/contributions/D45_HSPA6+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_OCA2+/average/contributions/D12_ENP_OCA2+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_NRG3+/average/contributions/D9_PFG_NRG3+_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ENP/average/contributions/D15_PP_ENP_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_pre_EC_beta/average/contributions/D12_pre_EC_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_beta/average/contributions/D45_SC_beta_counts.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CDH18+/average/contributions/D9_PFG_CDH18+_counts.h5
)
modiscos=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_delta/average/motifs/D45_SC_delta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_LEFTY1+/average/motifs/D4_DE_LEFTY1+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC_beta/average/motifs/D45_early_SC_EC_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_DE_SOX17+/average/motifs/D7_DE_SOX17+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_POLG2+/average/motifs/D4_DE_POLG2+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_beta/average/motifs/D22_SC_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_MitoHi/average/motifs/D4_DE_MitoHi.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_alpha/average/motifs/D22_late_SC_alpha.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_ONECUT1+/average/motifs/D7_GT_ONECUT1+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_liver/average/motifs/D12_liver.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC/average/motifs/D15_pre_EC.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC/average/motifs/D45_SC_EC.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_delta/average/motifs/D15_pre_delta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ONECUT2+/average/motifs/D4_DE_ONECUT2+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_FLT1+/average/motifs/D7_FLT1+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_NECTIN3-AS1+/average/motifs/D7_NECTIN3-AS1+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_CXCR4+/average/motifs/D7_GT_CXCR4+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ENP/average/motifs/D12_PP_ENP.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_others/average/motifs/D45_early_others.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ERBB4+/average/motifs/D15_PP_ERBB4+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_ENP/average/motifs/D22_late_ENP.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC_beta/average/motifs/D22_late_SC_EC_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_ENP_OCA2+/average/motifs/D22_ENP_OCA2+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_pre_alpha/average/motifs/D22_pre_alpha.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_ENP/average/motifs/D9_ENP.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_beta/average/motifs/D45_early_SC_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC_beta/average/motifs/D15_pre_EC_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_liver/average/motifs/D45_liver.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_FLT1+/average/motifs/D15_FLT1+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_SCG2+/average/motifs/D12_ENP_SCG2+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ERBB4+/average/motifs/D4_DE_ERBB4+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ERBB4+/average/motifs/D12_PP_ERBB4+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC/average/motifs/D22_late_SC_EC.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_NECTIN3-AS1+/average/motifs/D9_NECTIN3-AS1+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC_beta/average/motifs/D45_SC_EC_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_LINC01924+/average/motifs/D9_LINC01924+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_CDH8+/average/motifs/D4_DE_CDH8+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_others/average/motifs/D22_late_others.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_liver/average/motifs/D7_liver.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha1/average/motifs/D45_SC_alpha1.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_CREB5+/average/motifs/D15_PP_CREB5+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_PLOG2+/average/motifs/D9_PFG_PLOG2+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_CREB5+/average/motifs/D12_PP_CREB5+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_proliferating/average/motifs/D9_PFG_proliferating.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC/average/motifs/D22_SC_EC.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_POLG2/average/motifs/D7_GT_POLG2.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_ENP_EC/average/motifs/D45_ENP_EC.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_ARX+/average/motifs/D12_ENP_ARX+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_DE_OTX2+/average/motifs/D9_DE_OTX2+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_TTR+/average/motifs/D9_PFG_TTR+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_proliferating/average/motifs/D7_GT_proliferating.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_alpha/average/motifs/D45_early_SC_alpha.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_beta/average/motifs/D22_late_SC_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_liver/average/motifs/D15_liver.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_alpha/average/motifs/D15_pre_alpha.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_beta/average/motifs/D15_pre_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC_beta/average/motifs/D22_SC_EC_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CREB5+/average/motifs/D9_PFG_CREB5+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_liver/average/motifs/D22_liver.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha2/average/motifs/D45_SC_alpha2.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_ENP_OCA2+/average/motifs/D15_ENP_OCA2+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_liver/average/motifs/D9_liver.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_PP_ERBB4+/average/motifs/D22_PP_ERBB4+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_alpha/average/motifs/D22_SC_alpha.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_delta+/average/motifs/D22_SC_delta+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_SOX4+/average/motifs/D4_DE_SOX4+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC/average/motifs/D45_early_SC_EC.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_HSPA6+/average/motifs/D45_HSPA6+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_OCA2+/average/motifs/D12_ENP_OCA2+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_NRG3+/average/motifs/D9_PFG_NRG3+.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ENP/average/motifs/D15_PP_ENP.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_pre_EC_beta/average/motifs/D12_pre_EC_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_beta/average/motifs/D45_SC_beta.counts.modisco.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CDH18+/average/motifs/D9_PFG_CDH18+.counts.modisco.h5
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
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_delta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_LEFTY1+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_DE_SOX17+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_POLG2+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_MitoHi/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_alpha/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_ONECUT1+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_liver/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_delta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ONECUT2+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_FLT1+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_NECTIN3-AS1+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_CXCR4+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ENP/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_others/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ERBB4+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_ENP/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_ENP_OCA2+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_pre_alpha/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_ENP/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_liver/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_FLT1+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_SCG2+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ERBB4+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ERBB4+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_NECTIN3-AS1+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_LINC01924+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_CDH8+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_others/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_liver/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha1/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_CREB5+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_PLOG2+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_CREB5+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_proliferating/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_POLG2/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_ENP_EC/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_ARX+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_DE_OTX2+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_TTR+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_proliferating/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_alpha/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_liver/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_alpha/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CREB5+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_liver/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha2/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_ENP_OCA2+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_liver/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_PP_ERBB4+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_alpha/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_delta+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_SOX4+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_HSPA6+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_OCA2+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_NRG3+/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ENP/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_pre_EC_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_beta/average/motifs/hits/counts
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CDH18+/average/motifs/hits/counts
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