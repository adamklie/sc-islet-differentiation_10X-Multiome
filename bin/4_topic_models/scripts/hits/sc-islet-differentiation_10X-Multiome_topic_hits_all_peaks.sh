#!/bin/bash
#####
# Motif hit calling for the n100 CREsted topic CNN (oracle) over ALL peaks
# (not the per-topic top-K). Topic-side analog of the MT all-peaks finemo run.
# Uses:
#   - contributions: contributions_all_peaks/<Topic>_{oh,contrib}.npz
#   - unified topic motif catalog:
#       motifs/clustered_motifs.modisco.h5 (built from combined_mapped.meme via
#       bin/4_multi_task_model/scripts/hits/build_unified_modisco_h5.py)
#     + motifs/clustered_motifs_name_map.tsv
# Output:
#   results/4_topic_models/all/n100/crested_model/motifs/hits_all_peaks/<Topic>/
#     {hits.tsv, finemo.npz, report}
#
# Submit via /cellar/users/aklie/opt/SLURM/gpu_array.sh:
#   gpu_array.sh -s <this>.sh -j topic_allpeak_hits \
#     -g a30:1 -c 4 -m 32G -t 14-00:00:00 -n 26 -x 4 \
#     -o bin/slurm_logs/4_topic_models/hits_all_peaks
#####

date
echo -e "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID\n"

# Configure environment
source activate finemo_gpu
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/envs/finemo_gpu/lib/

# 26 confident topics (one-to-one + lineage). Order matches
# topic_categories.tsv filtered by category in {one-to-one,lineage}.
topics=(
    Topic1 Topic9 Topic13 Topic15 Topic16 Topic17 Topic20 Topic29
    Topic31 Topic32 Topic33 Topic36 Topic38 Topic39 Topic44 Topic45
    Topic47 Topic52 Topic67 Topic76 Topic79 Topic80 Topic86 Topic94
    Topic96 Topic99
)
topic=${topics[$SLURM_ARRAY_TASK_ID - 1]}
echo "Topic: $topic"

BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome
N100=${BASE}/results/4_topic_models/all/n100

# Single consensus peak set for the topic all-peaks pipeline (narrowPeak with
# region center as summit), built from binarized/consensus_regions.bed.
PEAK=${N100}/binarized/consensus_regions.narrowPeak
OH=${N100}/crested_model/contributions_all_peaks/${topic}_oh.npz
CONTRIB=${N100}/crested_model/contributions_all_peaks/${topic}_contrib.npz

# Unified topic motif catalog (built once from combined_mapped.meme)
modisco_h5=${N100}/crested_model/motifs/clustered_motifs.modisco.h5
name_map=${N100}/crested_model/motifs/clustered_motifs_name_map.tsv
window=1000

output_dir=${N100}/crested_model/motifs/hits_all_peaks/${topic}
mkdir -p $output_dir
echo "Output: $output_dir"

# Step 1: Extract regions from oh/contrib npz pair
npz_path=$output_dir/${topic}.finemo.npz
if [ ! -f "$npz_path" ]; then
    cmd="finemo extract-regions-modisco-fmt \
    -s $OH \
    -a $CONTRIB \
    -o $npz_path \
    -w $window"
    echo -e "Running extract:\n$cmd\n"
    eval $cmd
else
    echo "Reusing existing npz: $npz_path"
fi

# Step 2: Call hits against unified topic catalog
cmd="finemo call-hits \
-r $npz_path \
-m $modisco_h5 \
-o $output_dir \
-p $PEAK \
-N $name_map \
-l 0.8"
echo -e "Running call-hits:\n$cmd\n"
eval $cmd

# Step 3: Report (best-effort - finemo report can fail on synthetic modisco h5;
# downstream aggregation only needs hits.tsv).
cmd="finemo report \
-r $npz_path \
-p $PEAK \
-H $output_dir/hits.tsv \
-m $modisco_h5 \
-o $output_dir/report \
-W $window \
-N $name_map"
echo -e "Running report (best-effort):\n$cmd\n"
eval $cmd || echo "WARN: finemo report failed (expected for synthetic modisco h5). Continuing."

date
