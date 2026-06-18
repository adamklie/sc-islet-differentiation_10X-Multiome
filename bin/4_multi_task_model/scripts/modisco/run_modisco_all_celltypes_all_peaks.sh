#!/bin/bash
# Run MoDISco on peak-regression ALL-PEAK contributions (one task per cell type).
# Differences vs v3:
#   - reads NPZs from contributions_all_peaks/ (FULL per-CT peak set, 100-400K
#     regions @ 2114bp) instead of the top-K contributions_specific_v3 (~2000)
#   - -n 100000 (max_seqlets) to MATCH the ChromBPNet single-task pipeline
#     (v3 used 50000); the deeper region pool lets MoDISco draw its best 100K
#     seqlets, same rationale as ChromBPNet
#   - -w 400 to MATCH ChromBPNet (v3 used 1000); these are the same 2114bp
#     regions as the single-task models, so 400 gives cross-oracle consistency
#   - writes h5/report/PFMs to modisco_all_peaks/{CT}/ (leaves modisco_v3/ intact)
#
# Idempotent: skips CTs whose contributions aren't ready yet (exit 0) and CTs
# already done (modisco.h5 exists). Safe to re-submit as the last contributions
# land.
#
# Memory: the conversion-peak estimate (~40-50 GB for the biggest CT, DE @
# 399K regions) plus modisco motifs overhead fits in 128G per task. Cluster
# budget cap is 750 GB total, so we run at most 4 concurrent (4 × 128G = 512G,
# leaving headroom for the still-running contribution job).
#
# USAGE:
#   sbatch --job-name=pr_modisco_allpeak \
#     --partition=carter-compute --mem=128G --cpus-per-task=16 -t 14-00:00:00 \
#     --array=1-22%4 \
#     --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_multi_task_model/%x.%A_%a.out \
#     run_modisco_all_celltypes_all_peaks.sh

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
CONTRIBS=${BASE}/crested/contributions_all_peaks
OUT=${BASE}/crested/modisco_all_peaks/${CT}

if [ ! -f "${CONTRIBS}/${CT}_contrib.npz" ]; then
    echo "Contributions not found for ${CT}, skipping (will pick up on re-submit)."; exit 0
fi
if [ -f "${OUT}/${CT}.modisco.h5" ]; then
    echo "MoDISco already done for ${CT}, skipping."; exit 0
fi
mkdir -p ${OUT}

# Step 1: build h5 from the CREsted per-class NPZs (hypo contribs + one-hot)
echo "Step 1: building ${CT}_contribs.h5"
$PYTHON -c "
import numpy as np, h5py, sys
ct, contribs_dir, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
hypo = np.load(f'{contribs_dir}/{ct}_contrib.npz')['arr_0']   # (N, 4, L) float32
oh   = np.load(f'{contribs_dir}/{ct}_oh.npz')['arr_0']        # (N, 4, L)
proj = (hypo * oh).astype(np.float16)
h5_path = f'{out_dir}/{ct}_contribs.h5'
with h5py.File(h5_path, 'w') as f:
    f.create_dataset('projected_shap/seq', data=proj)
    f.create_dataset('shap/seq', data=hypo.astype(np.float16))
    f.create_dataset('raw/seq', data=oh.astype(np.int8))
print(f'  Saved: {h5_path}  shape={hypo.shape}')
" ${CT} ${CONTRIBS} ${OUT}

# Step 2: modisco motifs (n=100000 + w=400, ChromBPNet-matched)
echo "Step 2: modisco motifs ${CT}"
$MODISCO_BIN motifs -i ${OUT}/${CT}_contribs.h5 -n 100000 -o ${OUT}/${CT}.modisco.h5 -l 2 -w 400 -v

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
