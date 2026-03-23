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
output_dirs=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_delta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_LEFTY1+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_DE_SOX17+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_POLG2+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_MitoHi/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_alpha/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_ONECUT1+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_liver/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_delta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ONECUT2+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_FLT1+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_NECTIN3-AS1+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_CXCR4+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ENP/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_others/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ERBB4+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_ENP/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_ENP_OCA2+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_pre_alpha/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_ENP/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_liver/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_FLT1+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_SCG2+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ERBB4+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ERBB4+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_NECTIN3-AS1+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_LINC01924+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_CDH8+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_others/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_liver/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha1/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_CREB5+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_PLOG2+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_CREB5+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_proliferating/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_POLG2/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_ENP_EC/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_ARX+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_DE_OTX2+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_TTR+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_proliferating/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_alpha/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_liver/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_alpha/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CREB5+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_liver/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha2/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_ENP_OCA2+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_liver/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_PP_ERBB4+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_alpha/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_delta+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_SOX4+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_HSPA6+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_OCA2+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_NRG3+/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ENP/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_pre_EC_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_beta/average/motifs
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CDH18+/average/motifs
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
