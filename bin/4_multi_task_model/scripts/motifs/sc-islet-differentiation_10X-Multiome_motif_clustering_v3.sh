#! /bin/bash

#####
# Motif clustering for v3 run.
# Differences vs v2:
#   - reads PFMs from modisco_v3/*/pfms/*.pfm
#   - writes catalog to motifs_v3/{cluster,meme,tomtom,combined_mapped.meme,...}
# Keeps gimme cluster -t 0.999 (intentional, for ChromBPNet comparability).
# Also builds the synthetic modisco h5 + name_map needed by finemo at the end.
#
# USAGE: sbatch \
#--job-name=pr_motif_clustering_v3 \
#--partition=carter-compute \
#--output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_multi_task_model/%x.%A.out \
#--mem=128G \
#-n 4 \
#-t 02-00:00:00 \
#sc-islet-differentiation_10X-Multiome_motif_clustering_v3.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configure environment. gimme lives in test_celloracle; tomtom ships with the
# chrombpnet env — invoke tomtom via absolute path below.
source /cellar/users/aklie/opt/miniconda3/etc/profile.d/conda.sh
conda activate test_celloracle
TOMTOM=/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/bin/tomtom
EUGENE_PY=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python

BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model
MODISCO_ROOT=${BASE}/crested/modisco_v3
path_out=${BASE}/crested/motifs_v3
t=0.999

# Define input files (auto-discover PFMs from modisco_v3/*/pfms/<CT>.pfm)
celltypes=(
    DE ENP_phase1 FB_FLT1 PFG1 PFG2 PGT1 PGT2 PGT3 PP1 PP2
    SC_delta_GHRL early_ENP early_SC_EC early_SC_alpha early_SC_beta
    exocrine late_ENP late_SC_EC late_SC_alpha late_SC_beta liver
    proliferating_endocrine
)
paths_pfm=()
for ct in "${celltypes[@]}"; do
    pfm=${MODISCO_ROOT}/${ct}/pfms/${ct}.pfm
    if [ -f "$pfm" ]; then
        paths_pfm+=("$pfm")
    else
        echo "WARN: missing PFM $pfm (CT may have produced no patterns; skipping)"
    fi
done
echo "Combining ${#paths_pfm[@]} PFM files."

# Make the directory
mkdir -p $path_out/cluster
mkdir -p $path_out/meme
mkdir -p $path_out/tomtom

cat ${paths_pfm[@]} > $path_out/input_motifs.pfm

# Run cluster
cmd="gimme cluster \
$path_out/input_motifs.pfm \
$path_out/cluster \
-t $t"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Run conversion to meme
path_meme_script=/cellar/users/aklie/projects/ML4GLand/chrombpnet/scripts/pfm_to_meme.py
cmd="python $path_meme_script \
-i $path_out/cluster/clustered_motifs.pfm \
-o $path_out/meme"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Run TomTom for motif annotation
path_motif_database=/cellar/users/aklie/data/ref/motifs/jvierstra/motif-clustering-v2.0beta/motifs.meme
base_path=$path_out/meme
meme_files=($(find $base_path -name "*.meme" | sort))
for meme_file in "${meme_files[@]}"; do
  meme_id=$(basename $meme_file .meme)
  cmd="$TOMTOM \
-no-ssc \
-oc $path_out/tomtom/$meme_id \
-verbosity 1 \
-text \
-min-overlap 5 \
-mi 1 \
-dist pearson \
-evalue \
-thresh 10.0 \
$meme_file \
$path_motif_database > $path_out/tomtom/${meme_id}.tomtom.txt"
  echo -e "Running command:\n$cmd\n"
  eval $cmd
done

# extract top TF matches
for x in $path_out/tomtom/*.tomtom.txt ; do cat $x | head -2 | tail -1 ; done | cut -f1,2,5 | sort -k3g > $path_out/tfs_initial.txt

# combine all meme files into one
combined_meme_path="$path_out/meme/combined.meme"
rm -f "$combined_meme_path"
(
    awk '/^MOTIF/ {exit} {print}' "$(find $path_out/meme -name '*.meme' ! -name 'combined.meme' | head -n1)"
    for f in $path_out/meme/*.meme; do
        [ "$(basename "$f")" = "combined.meme" ] && continue
        awk '/^MOTIF/ {p=1} p' "$f"
    done
) > "$combined_meme_path"

# Rename motifs in combined MEME using best TomTom hits
awk '
NR==FNR {split($0, a, "\t"); map[a[1]]=a[2]; next}
/^MOTIF/ {split($0, a, " "); if (a[2] in map) {print "MOTIF " map[a[2]]; next}}
{print}
' $path_out/tfs_initial.txt $combined_meme_path > $path_out/combined_mapped.meme

# Build synthetic modisco h5 + name_map for finemo hit calling against this catalog
BUILD_SCRIPT=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_multi_task_model/scripts/hits/build_unified_modisco_h5.py
$EUGENE_PY $BUILD_SCRIPT \
    --meme $path_out/combined_mapped.meme \
    --out_h5 $path_out/clustered_motifs.modisco.h5 \
    --out_name_map $path_out/clustered_motifs_name_map.tsv

# Report final motif count
n_motifs=$(grep -c "^MOTIF" $path_out/combined_mapped.meme || echo 0)
echo "Final motif count in combined_mapped.meme: ${n_motifs}"

# Completion message
date
