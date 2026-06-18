#!/bin/bash
#####
# SLURM wrapper for the parallel sign-aware motif compendium build.
# Pairs with build_motif_compendium_all_peaks.py.
#
# Submit:
#   sbatch --job-name=pr_motif_compendium_all_peaks \
#     --partition=carter-gpu --account=carter-gpu --gres=gpu:a30:1 \
#     --cpus-per-task=8 --mem=64G -t 04:00:00 \
#     --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_multi_task_model/%x.%A.out \
#     build_motif_compendium_all_peaks.sh
#####
set -o pipefail
date
echo "Job ID: $SLURM_JOB_ID"

PY=/cellar/users/aklie/opt/miniconda3/envs/motifcompendium-gpu/bin/python
SCRIPT=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_multi_task_model/scripts/motifs/build_motif_compendium_all_peaks.py

${PY} ${SCRIPT}
date
