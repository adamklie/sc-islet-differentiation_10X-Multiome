#!/bin/bash
# SLURM wrapper for render_gimme_catalog_walls.py — renders pos + neg track
# motif walls (logo grids) for the dual-track gimme cluster catalog.
#
# Submit:
#   sbatch --job-name=pr_gimme_walls \
#     --partition=carter-gpu --account=carter-gpu --gres=gpu:a30:1 \
#     --cpus-per-task=4 --mem=16G -t 00:30:00 \
#     --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_multi_task_model/%x.%A.out \
#     render_gimme_catalog_walls.sh
set -o pipefail
date
echo "Job ID: $SLURM_JOB_ID"

PY=/cellar/users/aklie/opt/miniconda3/envs/motifcompendium-gpu/bin/python
SCRIPT=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_multi_task_model/scripts/motifs/render_gimme_catalog_walls.py

${PY} ${SCRIPT}
date
