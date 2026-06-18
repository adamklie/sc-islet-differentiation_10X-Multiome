#!/bin/bash
# SLURM wrapper — commit Scheme C as the canonical MC catalog.
set -o pipefail
date
echo "Job ID: $SLURM_JOB_ID"
PY=/cellar/users/aklie/opt/miniconda3/envs/motifcompendium-gpu/bin/python
SCRIPT=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_multi_task_model/scripts/motifs/build_scheme_C_catalog.py
${PY} ${SCRIPT}
date
