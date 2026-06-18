#!/bin/bash
# Run MoDISco on peak regression contributions (v3 run, one task per cell type).
# Differences vs v2:
#   - reads NPZs from contributions_specific_v3/
#   - writes h5/report/PFMs to modisco_v3/{CT}/
#   - modisco motifs uses -w 1000 (tutorial fidelity)
#
# USAGE:
#   sbatch --job-name=pr_modisco_v3 \
#     --partition=carter-compute --mem=128G --cpus-per-task=16 -t 14-00:00:00 \
#     --array=1-22%4 \
#     --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_multi_task_model/%x.%A_%a.out \
#     run_modisco_all_celltypes_v3.sh

date
echo -e "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID\n"

PYTHON=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python
MODISCO_BIN=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/modisco
PFM_SCRIPT=/cellar/users/aklie/projects/ML4GLand/chrombpnet/scripts/modisco_to_pfm.py
MEME_DB=/cellar/users/aklie/data/ref/motifs/jvierstra/motif-clustering-v2.0beta/motifs.meme
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK

celltypes=(
    DE ENP_phase1 FB_FLT1 PFG1 PFG2 PGT1 PGT2 PGT3 PP1 PP2
    SC_delta_GHRL early_ENP early_SC_EC early_SC_alpha early_SC_beta
    exocrine late_ENP late_SC_EC late_SC_alpha late_SC_beta liver
    proliferating_endocrine
)
CT=${celltypes[$SLURM_ARRAY_TASK_ID - 1]}
echo "Cell type: $CT"

BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model
CONTRIBS=${BASE}/crested/contributions_specific_v3
OUT=${BASE}/crested/modisco_v3/${CT}

if [ ! -f "${CONTRIBS}/${CT}_contrib.npz" ]; then
    echo "Contributions not found for ${CT}, skipping."; exit 0
fi
if [ -f "${OUT}/${CT}.modisco.h5" ]; then
    echo "MoDISco already done for ${CT}, skipping."; exit 0
fi
mkdir -p ${OUT}

# Step 1: build h5 from the CREsted per-class NPZs
echo "Step 1: building ${CT}_contribs.h5"
$PYTHON -c "
import numpy as np, h5py, sys
ct, contribs_dir, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
hypo = np.load(f'{contribs_dir}/{ct}_contrib.npz')['arr_0']   # (N, 4, L)
oh   = np.load(f'{contribs_dir}/{ct}_oh.npz')['arr_0']        # (N, 4, L)
proj = (hypo * oh).astype(np.float16)
h5_path = f'{out_dir}/{ct}_contribs.h5'
with h5py.File(h5_path, 'w') as f:
    f.create_dataset('projected_shap/seq', data=proj)
    f.create_dataset('shap/seq', data=hypo.astype(np.float16))
    f.create_dataset('raw/seq', data=oh.astype(np.int8))
print(f'  Saved: {h5_path}  shape={hypo.shape}')
" ${CT} ${CONTRIBS} ${OUT}

# Step 2: modisco motifs (window=1000, tutorial fidelity)
echo "Step 2: modisco motifs ${CT}"
$MODISCO_BIN motifs -i ${OUT}/${CT}_contribs.h5 -n 50000 -o ${OUT}/${CT}.modisco.h5 -l 2 -w 1000 -v

# Step 3: modisco report
echo "Step 3: modisco report ${CT}"
$MODISCO_BIN report -i ${OUT}/${CT}.modisco.h5 -o ${OUT}/${CT}.modisco.report/ -s ${OUT}/${CT}.modisco.report/ -m ${MEME_DB}

# Step 4: extract PFMs
echo "Step 4: extract PFMs ${CT}"
mkdir -p ${OUT}/pfms
$PYTHON ${PFM_SCRIPT} -m ${OUT}/${CT}.modisco.h5 -o ${OUT}/pfms -op ${CT} -f 2

# Consolidate PFMs
tmp=${OUT}/pfms/${CT}.tmp.pfm; final=${OUT}/pfms/${CT}.pfm
rm -f "$final" "$tmp"
for x in ${OUT}/pfms/${CT}*.pfm; do
    [ "$x" = "$final" ] && continue
    [ "$x" = "$tmp" ] && continue
    echo ">$(basename "$x")" >> "$tmp"; cat "$x" >> "$tmp"
done
[ -f "$tmp" ] && mv "$tmp" "$final"

echo "Done with ${CT}"
date
