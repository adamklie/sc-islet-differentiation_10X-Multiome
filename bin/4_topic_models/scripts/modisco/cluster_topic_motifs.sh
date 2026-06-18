#!/bin/bash

#####
# Cluster de novo motifs from CREsted MoDISco across all topics
# Mirrors the ChromBPNet motif clustering pipeline exactly:
#   1. Concatenate PFMs from all topics
#   2. gimme cluster (CellOracle)
#   3. Convert to MEME format
#   4. TomTom annotation against Vierstra database
#   5. Extract top TF matches + combined MEME
#
# USAGE: sbatch \
#   --job-name=cluster_topic_motifs \
#   --partition=carter-compute \
#   --output=slurm_logs/4_topic_models/modisco/%x.%A.out \
#   --mem=128G \
#   -n 4 \
#   -t 02-00:00:00 \
#   cluster_topic_motifs.sh <subset>
#
# Example: sbatch ... cluster_topic_motifs.sh endocrine_replicate
#
# Prerequisites: run_modisco_all_topics.sh must have completed for this subset
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Config — gimme cluster needs test_celloracle env
CELLORACLE_PYTHON=/cellar/users/aklie/opt/miniconda3/envs/test_celloracle/bin/python
EUGENE_PYTHON=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/envs/test_celloracle/lib/

PFM_TO_MEME=/cellar/users/aklie/projects/ML4GLand/chrombpnet/scripts/pfm_to_meme.py
MEME_DB=/cellar/users/aklie/data/ref/motifs/jvierstra/motif-clustering-v2.0beta/motifs.meme

# Subset from command line arg
SUBSET=${1:-endocrine_replicate}
BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/${SUBSET}
MODISCO_DIR=${BASE}/crested_model/modisco
MOTIFS_OUT=${BASE}/crested_model/motifs

mkdir -p ${MOTIFS_OUT}/cluster
mkdir -p ${MOTIFS_OUT}/meme
mkdir -p ${MOTIFS_OUT}/tomtom

# Step 1: Concatenate all topic PFMs
echo "Step 1: Concatenating PFMs from all topics..."
rm -f ${MOTIFS_OUT}/input_motifs.pfm
for topic_dir in ${MODISCO_DIR}/Topic*/; do
    topic=$(basename $topic_dir)
    pfm_file="${topic_dir}/pfms/${topic}.pfm"
    if [ -f "$pfm_file" ]; then
        echo "  Adding ${topic}"
        cat "$pfm_file" >> ${MOTIFS_OUT}/input_motifs.pfm
    fi
done
n_motifs=$(grep -c "^>" ${MOTIFS_OUT}/input_motifs.pfm 2>/dev/null || echo 0)
echo "  Total motifs: ${n_motifs}"

# Step 2: gimme cluster
echo "Step 2: Running gimme cluster (t=0.999)..."
# Need celloracle env for gimme
export PATH=/cellar/users/aklie/opt/miniconda3/envs/test_celloracle/bin:$PATH
gimme cluster \
    ${MOTIFS_OUT}/input_motifs.pfm \
    ${MOTIFS_OUT}/cluster \
    -t 0.999

# Step 3: Convert to MEME format
echo "Step 3: Converting to MEME format..."
${EUGENE_PYTHON} ${PFM_TO_MEME} \
    -i ${MOTIFS_OUT}/cluster/clustered_motifs.pfm \
    -o ${MOTIFS_OUT}/meme

# Step 4: TomTom annotation
echo "Step 4: Running TomTom..."
meme_files=($(find ${MOTIFS_OUT}/meme -name "*.meme" ! -name "combined.meme" | sort))
for meme_file in "${meme_files[@]}"; do
    meme_id=$(basename $meme_file .meme)
    tomtom \
        -no-ssc \
        -oc ${MOTIFS_OUT}/tomtom/$meme_id \
        -verbosity 1 \
        -text \
        -min-overlap 5 \
        -mi 1 \
        -dist pearson \
        -evalue \
        -thresh 10.0 \
        $meme_file \
        ${MEME_DB} > ${MOTIFS_OUT}/tomtom/${meme_id}.tomtom.txt
done

# Step 5: Extract top TF matches
echo "Step 5: Extracting top TF matches..."
for x in ${MOTIFS_OUT}/tomtom/*.tomtom.txt; do
    cat $x | head -2 | tail -1
done | cut -f1,2,5 | sort -k3g > ${MOTIFS_OUT}/tfs_initial.txt

# Step 6: Combined MEME file
echo "Step 6: Building combined MEME..."
combined_meme_path="${MOTIFS_OUT}/meme/combined.meme"
rm -f "$combined_meme_path"
(
    awk '/^MOTIF/ {exit} {print}' "$(find ${MOTIFS_OUT}/meme -name '*.meme' ! -name 'combined.meme' | head -n1)"
    for f in ${MOTIFS_OUT}/meme/*.meme; do
        [ "$(basename "$f")" = "combined.meme" ] && continue
        awk '/^MOTIF/ {p=1} p' "$f"
    done
) > "$combined_meme_path"

# Step 7: Rename motifs using best TomTom hits
echo "Step 7: Creating mapped MEME file..."
awk '
NR==FNR {split($0, a, "\t"); map[a[1]]=a[2]; next}
/^MOTIF/ {split($0, a, " "); if (a[2] in map) {print "MOTIF " map[a[2]]; next}}
{print}
' ${MOTIFS_OUT}/tfs_initial.txt ${combined_meme_path} > ${MOTIFS_OUT}/combined_mapped.meme

n_clustered=$(grep -c "^MOTIF" ${combined_meme_path} 2>/dev/null || echo 0)
echo ""
echo "=== Summary ==="
echo "  Input motifs: ${n_motifs}"
echo "  Clustered motifs: ${n_clustered}"
echo "  Output: ${MOTIFS_OUT}/"
echo ""

date
