#!/bin/bash
#####
# Motif clustering for the all-peak run — dual track (pos + neg in parallel).
#
# Steps:
#   1. Per-CT: extract pos+neg PFMs from modisco_all_peaks/{CT}/{CT}.modisco.h5
#      via extract_pfms_with_signs.py → per-CT pfms/{pos,neg}/{CT}.{pos,neg}_pattern_N.pfm
#   2. Per-CT concat → modisco_all_peaks/{CT}/pfms/{CT}.{pos,neg}.pfm
#   3. Per-sign aggregation across all 22 CTs → motifs_all_peaks/{pos,neg}/input_motifs.pfm
#   4. gimme cluster (-t 0.999, identical to v3) on each track separately
#   5. pfm_to_meme + TomTom (Vierstra) per sign
#   6. Rename motifs by best TomTom hit, then PREFIX with pos_/neg_ so downstream
#      can split tracks: pos_<TF_name> / neg_<TF_name>
#   7. Concatenate pos+neg mapped MEMEs into one combined_mapped.meme (still
#      one catalog file, but every motif is sign-tagged in the name)
#   8. build_unified_modisco_h5.py → synthetic h5 + name_map for finemo
#
# Output: motifs_all_peaks/
#   pos/, neg/                       per-sign intermediates (input.pfm, cluster/,
#                                     meme/, tomtom/, tfs_initial.txt, mapped.meme)
#   combined_mapped.meme             both tracks, sign-prefixed motif names
#   clustered_motifs.modisco.h5      finemo input
#   clustered_motifs_name_map.tsv    finemo name lookup
#
# Submit (depends on the all-peak MoDISco array finishing):
#   sbatch --dependency=afterok:11086563 \
#     --job-name=pr_motif_clustering_all_peaks \
#     --partition=carter-compute \
#     --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_multi_task_model/%x.%A.out \
#     --mem=64G --cpus-per-task=4 -t 04:00:00 \
#     sc-islet-differentiation_10X-Multiome_motif_clustering_all_peaks.sh
#####
# NOTE: do NOT use 'set -u' here — the test_celloracle conda activation hook
# references an unset $ADDR2LINE and would fatally exit before our pipeline runs.
set -o pipefail

date
echo "Job ID: $SLURM_JOB_ID"

source /cellar/users/aklie/opt/miniconda3/etc/profile.d/conda.sh
conda activate test_celloracle
TOMTOM=/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/bin/tomtom
EUGENE_PY=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python
PFM_TO_MEME=/cellar/users/aklie/projects/ML4GLand/chrombpnet/scripts/pfm_to_meme.py
VIERSTRA_MEME=/cellar/users/aklie/data/ref/motifs/jvierstra/motif-clustering-v2.0beta/motifs.meme
EXTRACT_SCRIPT=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_multi_task_model/scripts/modisco/extract_pfms_with_signs.py
BUILD_H5=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_multi_task_model/scripts/hits/build_unified_modisco_h5.py

BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model
MODISCO_ROOT=${BASE}/crested/modisco_all_peaks
OUT=${BASE}/crested/motifs_all_peaks
THRESH=0.999

celltypes=(
    DE ENP_phase1 FB_FLT1 PFG1 PFG2 PGT1 PGT2 PGT3 PP1 PP2
    SC_delta_GHRL early_ENP early_SC_EC early_SC_alpha early_SC_beta
    exocrine late_ENP late_SC_EC late_SC_alpha late_SC_beta liver
    proliferating_endocrine
)

mkdir -p ${OUT}/pos/cluster ${OUT}/pos/meme ${OUT}/pos/tomtom
mkdir -p ${OUT}/neg/cluster ${OUT}/neg/meme ${OUT}/neg/tomtom

############################################################################
# Step 1+2: extract pos+neg PFMs per CT, then concat per-CT-per-sign
############################################################################
echo
echo "=== Step 1+2: extract + concat PFMs per CT ==="
for ct in "${celltypes[@]}"; do
    h5=${MODISCO_ROOT}/${ct}/${ct}.modisco.h5
    if [ ! -f "$h5" ]; then
        echo "  WARN: missing $h5 — skipping $ct"
        continue
    fi
    pfm_dir=${MODISCO_ROOT}/${ct}/pfms
    mkdir -p ${pfm_dir}/pos ${pfm_dir}/neg
    # Extract (writes pos/, neg/ subdirs of pfm_dir, one .pfm per pattern)
    ${EUGENE_PY} ${EXTRACT_SCRIPT} -m ${h5} -o ${pfm_dir} -op ${ct} -f 2
    # Concat each side into a single per-CT-per-sign .pfm with named blocks
    for sign in pos neg; do
        merged=${pfm_dir}/${ct}.${sign}.pfm
        rm -f "$merged"
        for x in ${pfm_dir}/${sign}/${ct}.${sign}_pattern_*.pfm; do
            [ -f "$x" ] || continue
            echo ">$(basename "$x" .pfm)" >> "$merged"
            cat "$x" >> "$merged"
        done
        [ -s "$merged" ] && echo "  $(basename "$merged"): $(grep -c '^>' "$merged") patterns" || echo "  (no $sign patterns for $ct)"
    done
done

############################################################################
# Step 3-7: per-sign clustering pipeline (gimme cluster → pfm_to_meme → TomTom → mapped.meme)
############################################################################
for sign in pos neg; do
    echo
    echo "=== ${sign^^} TRACK ==="
    SDIR=${OUT}/${sign}
    # 3. aggregate across CTs
    : > ${SDIR}/input_motifs.pfm
    n_blocks=0
    for ct in "${celltypes[@]}"; do
        f=${MODISCO_ROOT}/${ct}/pfms/${ct}.${sign}.pfm
        if [ -s "$f" ]; then
            cat "$f" >> ${SDIR}/input_motifs.pfm
            n=$(grep -c '^>' "$f")
            n_blocks=$((n_blocks + n))
        fi
    done
    echo "Aggregated ${n_blocks} ${sign} PFMs across CTs → ${SDIR}/input_motifs.pfm"
    if [ $n_blocks -eq 0 ]; then
        echo "  (no ${sign} motifs anywhere — skipping ${sign} track)"
        continue
    fi

    # 4. gimme cluster
    echo "Step 4 (${sign}): gimme cluster -t ${THRESH}"
    gimme cluster ${SDIR}/input_motifs.pfm ${SDIR}/cluster -t ${THRESH}

    # 5. pfm_to_meme
    echo "Step 5 (${sign}): pfm_to_meme"
    python ${PFM_TO_MEME} -i ${SDIR}/cluster/clustered_motifs.pfm -o ${SDIR}/meme

    # 6. TomTom per cluster MEME
    echo "Step 6 (${sign}): TomTom per cluster"
    for meme_file in $(find ${SDIR}/meme -name "*.meme" | sort); do
        meme_id=$(basename ${meme_file} .meme)
        ${TOMTOM} -no-ssc -oc ${SDIR}/tomtom/${meme_id} -verbosity 1 -text \
            -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10.0 \
            ${meme_file} ${VIERSTRA_MEME} > ${SDIR}/tomtom/${meme_id}.tomtom.txt
    done

    # 7. extract top TomTom hits → tfs_initial.txt
    for x in ${SDIR}/tomtom/*.tomtom.txt; do
        cat "$x" | head -2 | tail -1
    done | cut -f1,2,5 | sort -k3g > ${SDIR}/tfs_initial.txt

    # combine all per-cluster MEMEs
    combined=${SDIR}/meme/combined.meme
    rm -f "${combined}"
    (
        awk '/^MOTIF/ {exit} {print}' "$(find ${SDIR}/meme -name '*.meme' ! -name 'combined.meme' | head -n1)"
        for f in ${SDIR}/meme/*.meme; do
            [ "$(basename "$f")" = "combined.meme" ] && continue
            awk '/^MOTIF/ {p=1} p' "$f"
        done
    ) > "${combined}"

    # 8. Rename motifs by best TomTom hit, then prefix with sign_ (pos_/neg_)
    awk -v sign="${sign}" '
        NR==FNR {split($0, a, "\t"); map[a[1]]=a[2]; next}
        /^MOTIF/ {
            split($0, a, " ");
            if (a[2] in map) {
                print "MOTIF " sign "_" map[a[2]]; next
            } else {
                print "MOTIF " sign "_" a[2]; next
            }
        }
        {print}
    ' ${SDIR}/tfs_initial.txt "${combined}" > ${SDIR}/mapped.meme

    n_motifs=$(grep -c "^MOTIF" ${SDIR}/mapped.meme || echo 0)
    echo "  ${sign} track: ${n_motifs} clustered motifs"
done

############################################################################
# Step 9: concatenate pos + neg mapped MEMEs → combined_mapped.meme
############################################################################
echo
echo "=== Step 9: concatenate pos + neg into one catalog ==="
COMBINED=${OUT}/combined_mapped.meme
rm -f ${COMBINED}
if [ -f ${OUT}/pos/mapped.meme ]; then
    # Header from pos
    awk '/^MOTIF/ {exit} {print}' ${OUT}/pos/mapped.meme > ${COMBINED}
    awk '/^MOTIF/ {p=1} p' ${OUT}/pos/mapped.meme >> ${COMBINED}
fi
if [ -f ${OUT}/neg/mapped.meme ]; then
    # If pos didn't write a header, take it from neg
    [ ! -s ${COMBINED} ] && awk '/^MOTIF/ {exit} {print}' ${OUT}/neg/mapped.meme > ${COMBINED}
    awk '/^MOTIF/ {p=1} p' ${OUT}/neg/mapped.meme >> ${COMBINED}
fi
n_total=$(grep -c "^MOTIF" ${COMBINED} || echo 0)
n_pos=$(grep -cE "^MOTIF pos_" ${COMBINED} || echo 0)
n_neg=$(grep -cE "^MOTIF neg_" ${COMBINED} || echo 0)
echo "  combined: ${n_total} motifs  (${n_pos} pos, ${n_neg} neg)"

############################################################################
# Step 10: build synthetic modisco h5 + name_map for finemo
############################################################################
echo
echo "=== Step 10: synthetic modisco h5 for finemo ==="
${EUGENE_PY} ${BUILD_H5} \
    --meme ${COMBINED} \
    --out_h5 ${OUT}/clustered_motifs.modisco.h5 \
    --out_name_map ${OUT}/clustered_motifs_name_map.tsv

# Ensure name_map has no duplicate names (we hit this with v3; finemo errors).
# If there are dups, write a _unique variant with positional suffixes.
UNIQ=${OUT}/clustered_motifs_name_map_unique.tsv
${EUGENE_PY} -c "
import pandas as pd
df = pd.read_csv('${OUT}/clustered_motifs_name_map.tsv', sep='\t', header=None, names=['key','name'])
dups = df['name'].duplicated(keep=False)
if dups.any():
    df.loc[dups, 'name'] = df.loc[dups, 'name'] + '_#' + df.loc[dups].groupby('name').cumcount().add(1).astype(str)
    print(f'  resolved {dups.sum()} duplicate names → _unique map')
else:
    print('  no duplicates — _unique = original')
df.to_csv('${UNIQ}', sep='\t', header=False, index=False)
"

echo
echo "DONE."
date
