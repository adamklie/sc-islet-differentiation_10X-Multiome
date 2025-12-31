#! /bin/bash

#####
# Script to get variant scores with a trained chrombpnet model
# USAGE: sbatch \
#--job-name=variant_scoring \
#--account carter-gpu \
#--partition carter-gpu \
#--gpus=a30:1 \
#--output slurm_logs/%x.%A.%a.out \
#--mem=32G \
#-n 4 \
#-t 02-00:00:00 \
#--array=1-12%12 \
#variant_scoring.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Set-up env
source activate chrombpnet
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/lib
python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"

# file lists
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
chrombpnet_nobias_models=(
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
output_dirs=(
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_delta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_LEFTY1+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC_beta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_DE_SOX17+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_POLG2+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_beta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_MitoHi/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_alpha/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_ONECUT1+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_liver/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_delta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ONECUT2+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_FLT1+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_NECTIN3-AS1+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_CXCR4+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ENP/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_others/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ERBB4+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_ENP/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC_beta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_ENP_OCA2+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_pre_alpha/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_ENP/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_beta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_EC_beta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_liver/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_FLT1+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_SCG2+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_ERBB4+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_ERBB4+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_EC/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_NECTIN3-AS1+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_EC_beta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_LINC01924+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_CDH8+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_others/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_liver/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha1/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_CREB5+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_PLOG2+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_PP_CREB5+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_proliferating/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_POLG2/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_ENP_EC/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_ARX+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_DE_OTX2+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_TTR+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D7_GT_proliferating/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_alpha/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_late_SC_beta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_liver/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_alpha/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_pre_beta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_EC_beta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CREB5+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_liver/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_alpha2/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_ENP_OCA2+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_liver/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_PP_ERBB4+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_alpha/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D22_SC_delta+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D4_DE_SOX4+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_early_SC_EC/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_HSPA6+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_ENP_OCA2+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_NRG3+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_PP_ENP/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D12_pre_EC_beta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D45_SC_beta/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
	/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D9_PFG_CDH18+/fold_0/chrombpnet/0.5/variant_scoring/T1D_Chiou_2021
)
list=/cellar/users/aklie/data/ref/variants/T1D/T1D_Chiou_2021_cred_set.wGWAS_info.minimal.snps.chrombpnet.bed
random_seed=1234

# Grab each for this SLURM task
celltype=${celltypes[$SLURM_ARRAY_TASK_ID - 1]}
chrombpnet_nobias_model=${chrombpnet_nobias_models[$SLURM_ARRAY_TASK_ID - 1]}
output_dir=${output_dirs[$SLURM_ARRAY_TASK_ID - 1]}

# echo the celltype and peak
echo -e "Celltype: $celltype"
echo -e "Model path: $chrombpnet_nobias_model"
echo -e "Output directory: $output_dir\n"
echo -e "Variant list: $list"
echo -e "Random seed: $random_seed\n"

# make the output directory
mkdir -p $output_dir

# Run cmd
variant_scoring_script=/cellar/users/aklie/opt/variant-scorer/src/variant_scoring.py
cmd="python $variant_scoring_script \
-l $list \
-g /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
-m $chrombpnet_nobias_model \
-o $output_dir/${celltype} \
-s /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
-r $random_seed \
-sc bed"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Date
date
