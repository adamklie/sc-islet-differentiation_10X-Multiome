#!/bin/bash
# Track C (cross-model): marginalize the MULTITASK catalog motifs (motifs_v3, 57)
# through each of the 22 single-task ChromBPNet models, to compare each motif's
# effect under the single-task models vs. the multitask model.
# Mirrors the constituent-analysis marginalize_global.sh, repointed to:
#   - the multitask motif MEME as the query
#   - the post-reorg models/<ct>/ ST model layout
#   - output under results/3_single_task_models/marginalization/mt_catalog/<ct>/
#
# Submit:
#   sbatch --partition=carter-gpu --account=carter-gpu --gpus=a30:1 --mem=32G -n 4 \
#     -t 12:00:00 --array=1-22%6 --output=slurm_logs/%x.%A.%a.out marginalize_mt_through_st.sh
set -e
date; echo -e "Job $SLURM_JOB_ID Task $SLURM_ARRAY_TASK_ID\n"
PY=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/lib
$PY -c "import torch; print('cuda', torch.cuda.is_available())"

DATA=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results
PEAK=$DATA/2_process_data/peak_calls/all
MROOT=$DATA/3_single_task_models/models
OUT=$DATA/3_single_task_models/marginalization/mt_catalog
MOTIF=$DATA/4_multi_task_model/crested/motifs_v3/combined_mapped.meme
GENOME=/cellar/users/aklie/data/ref/genomes/hg38/hg38.fa

celltypes=(DE ENP_phase1 FB_FLT1 PFG1 PFG2 PGT1 PGT2 PGT3 PP1 PP2 SC_delta_GHRL early_ENP early_SC_EC early_SC_alpha early_SC_beta exocrine late_ENP late_SC_EC late_SC_alpha late_SC_beta liver proliferating_endocrine)
ct=${celltypes[$SLURM_ARRAY_TASK_ID - 1]}
model=$MROOT/$ct/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5
peak=$PEAK/$ct.narrowPeak
outdir=$OUT/$ct
echo -e "CT $ct\nmodel $model\npeak $peak\nout $outdir\n"; mkdir -p $outdir
$PY /cellar/users/aklie/projects/ML4GLand/chrombpnet/scripts/marginalization.py \
  -m $model -f $MOTIF -g $GENOME -p $peak -o $outdir -n 100
echo "DONE $ct"; date
