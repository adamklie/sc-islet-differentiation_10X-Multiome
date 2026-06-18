#!/bin/bash
#####
# Per-topic contribution scoring on the n100 CREsted topic CNN (oracle),
# computed over ALL peaks (no specificity filter) — restricted to the 26
# confident topics (categories one-to-one + lineage per
# results/4_topic_models/all/n100/summary/topic_categories.tsv).
#
# Writes per-topic NPZ pairs at
#   contributions_all_peaks/<Topic>_{oh,contrib}.npz
# which the topic finemo hits stage consumes (mirrors MT all-peaks layout).
#
# Submit via /cellar/users/aklie/opt/SLURM/gpu_array.sh:
#   gpu_array.sh -s <this>.sh -j topic_allpeak_contrib \
#     -g a30:1 -c 4 -m 64G -t 14-00:00:00 -n 26 -x 4 \
#     -o bin/slurm_logs/4_topic_models/contributions_all_peaks
#####

date
echo -e "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
export KERAS_BACKEND=torch

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

SCRIPT_DIR=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_topic_models/scripts/contributions
BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/n100

mkdir -p ${BASE}/crested_model/contributions_all_peaks

python ${SCRIPT_DIR}/compute_contributions_topic_all_peaks.py \
    --input_h5ad ${BASE}/crested_model/contributions_specific/adata_with_predictions.h5ad \
    --checkpoint ${BASE}/crested_model/checkpoints/60.keras \
    --genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --chrom_sizes /cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes \
    --output_dir ${BASE}/crested_model/contributions_all_peaks \
    --topics ${topic} \
    --method integrated_grad \
    --batch_size 128 \
    --chunk_size 20000 \
    --backend torch \
    --skip_if_exists

date
