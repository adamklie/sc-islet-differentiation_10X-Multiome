#!/bin/bash
# SLURM wrapper for finalize_mc_catalog.py — exports MEME, builds walls + name_map
# from the saved motif_compendium_clustered.mc.
set -o pipefail
date
echo "Job ID: $SLURM_JOB_ID"
PY=/cellar/users/aklie/opt/miniconda3/envs/motifcompendium-gpu/bin/python
SCRIPT=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_multi_task_model/scripts/motifs/finalize_mc_catalog.py
${PY} ${SCRIPT}
date
