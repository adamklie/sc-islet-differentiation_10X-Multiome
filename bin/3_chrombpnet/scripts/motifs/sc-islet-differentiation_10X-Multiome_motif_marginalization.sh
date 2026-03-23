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
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_delta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_LEFTY1+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_DE_SOX17+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_POLG2+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_MitoHi/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_alpha/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_ONECUT1+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_liver/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_delta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ONECUT2+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_FLT1+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_NECTIN3-AS1+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_CXCR4+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ENP/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_others/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ERBB4+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_ENP/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_ENP_OCA2+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_pre_alpha/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_ENP/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_liver/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_FLT1+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_SCG2+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ERBB4+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ERBB4+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_NECTIN3-AS1+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_LINC01924+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_CDH8+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_others/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_liver/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha1/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_CREB5+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_PLOG2+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_CREB5+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_proliferating/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_POLG2/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_ENP_EC/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_ARX+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_DE_OTX2+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_TTR+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_proliferating/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_alpha/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_liver/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_alpha/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CREB5+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_liver/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha2/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_ENP_OCA2+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_liver/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_PP_ERBB4+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_alpha/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_delta+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_SOX4+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_HSPA6+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_OCA2+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_NRG3+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ENP/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_pre_EC_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_beta/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CDH18+/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
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
output_dirs=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_delta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_LEFTY1+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_DE_SOX17+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_POLG2+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_MitoHi/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_alpha/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_ONECUT1+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_liver/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_delta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ONECUT2+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_FLT1+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_NECTIN3-AS1+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_CXCR4+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ENP/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_others/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ERBB4+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_ENP/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_ENP_OCA2+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_pre_alpha/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_ENP/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_liver/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_FLT1+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_SCG2+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ERBB4+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ERBB4+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_NECTIN3-AS1+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_LINC01924+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_CDH8+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_others/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_liver/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha1/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_CREB5+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_PLOG2+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_CREB5+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_proliferating/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_POLG2/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_ENP_EC/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_ARX+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_DE_OTX2+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_TTR+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_proliferating/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_alpha/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_liver/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_alpha/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CREB5+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_liver/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha2/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_ENP_OCA2+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_liver/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_PP_ERBB4+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_alpha/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_delta+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_SOX4+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_HSPA6+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_OCA2+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_NRG3+/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ENP/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_pre_EC_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_beta/fold_0/chrombpnet/0.5/marginalization/
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CDH18+/fold_0/chrombpnet/0.5/marginalization/
)
n_seqs=100

# SLURM task-specific values
model=${models[$SLURM_ARRAY_TASK_ID - 1]}
peak=${peaks[$SLURM_ARRAY_TASK_ID - 1]}
output_dir=${output_dirs[$SLURM_ARRAY_TASK_ID - 1]}

# Constants
genome_fasta="/cellar/users/aklie/data/ref/genomes/hg38/hg38.fa"
motif_file="/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/motifs/meme/combined.meme"

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
