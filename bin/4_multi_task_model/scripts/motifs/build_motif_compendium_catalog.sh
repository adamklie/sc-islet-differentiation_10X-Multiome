#!/bin/bash
#####
# SLURM wrapper for build_motif_compendium_catalog.py — sign-aware MC catalog
# replacing gimme cluster. Output: results/4_multi_task_model/crested/motifs_compendium/
#
# Submit:
#   sbatch --job-name=pr_build_mc_catalog \
#     --partition=carter-gpu --account=carter-gpu --gres=gpu:a30:1 \
#     --cpus-per-task=8 --mem=64G -t 02:00:00 \
#     --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_multi_task_model/%x.%A.out \
#     build_motif_compendium_catalog.sh
#####
set -o pipefail
date
echo "Job ID: $SLURM_JOB_ID"

PY=/cellar/users/aklie/opt/miniconda3/envs/motifcompendium-gpu/bin/python
SCRIPT=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_multi_task_model/scripts/motifs/build_motif_compendium_catalog.py

${PY} ${SCRIPT}
date
