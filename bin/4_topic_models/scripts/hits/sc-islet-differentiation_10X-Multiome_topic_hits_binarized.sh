#!/bin/bash
#####
# Per-topic finemo hit calling on the existing per-topic binarized contributions
# (the same 1000 regions/topic that MoDISco was run on). This closes the
# methodological gap between topic models and ChromBPNet/MT — topic hits had
# previously been "MoDISco-seqlets-as-hits", which is hits-by-construction.
# This pass is an independent finemo rescore against the unified n100 catalog.
#
# Inputs (per topic):
#   - one_hot_seqs.npz, hypothetical_contribs.npz under contributions/Topic{N}/
#   - clustered_motifs.modisco.h5 + clustered_motifs_name_map.tsv (n100 catalog)
#   - Topic{N}.narrowPeak (1000 rows, written by replay_topic_binarized_regions.py)
# Outputs:
#   - motifs/hits_finemo/Topic{N}/{hits.tsv, finemo.npz, report/}
#
# Submit (100 topics, concurrency 4, GPU finemo):
#   sbatch --partition=carter-gpu --account=carter-gpu --gres=gpu:a30:1 \
#     --cpus-per-task=4 --mem=16G --time=4:00:00 \
#     --array=1-100%4 --job-name=topic_hits_binarized \
#     --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_topic_models/hits_binarized/%x.%A.%a.out \
#     <this>.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID\n"

source activate finemo_gpu
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/envs/finemo_gpu/lib/

TOPIC="Topic${SLURM_ARRAY_TASK_ID}"
echo "Topic: $TOPIC"

BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/n100/crested_model

OH=${BASE}/contributions/${TOPIC}/one_hot_seqs.npz
CONTRIB=${BASE}/contributions/${TOPIC}/hypothetical_contribs.npz
PEAK=${BASE}/motifs/hits_finemo/peaks/${TOPIC}.narrowPeak
MODISCO_H5=${BASE}/motifs/clustered_motifs.modisco.h5
NAME_MAP=${BASE}/motifs/clustered_motifs_name_map_unique.tsv  # 5 dup TF names suffixed with #<idx> for finemo
WINDOW=400

OUT_DIR=${BASE}/motifs/hits_finemo/${TOPIC}
mkdir -p ${OUT_DIR}

# Step 1: extract per-region windows from oh/contrib pair
NPZ_PATH=${OUT_DIR}/${TOPIC}.finemo.npz
if [ ! -f "${NPZ_PATH}" ]; then
    finemo extract-regions-modisco-fmt -s ${OH} -a ${CONTRIB} -o ${NPZ_PATH} -w ${WINDOW}
else
    echo "Reusing existing ${NPZ_PATH}"
fi

# Step 2: call hits against unified n100 catalog
finemo call-hits \
    -r ${NPZ_PATH} \
    -m ${MODISCO_H5} \
    -o ${OUT_DIR} \
    -p ${PEAK} \
    -N ${NAME_MAP} \
    -l 0.8

# Step 3: report (best-effort — synthetic modisco h5 sometimes fails report)
finemo report \
    -r ${NPZ_PATH} \
    -p ${PEAK} \
    -H ${OUT_DIR}/hits.tsv \
    -m ${MODISCO_H5} \
    -o ${OUT_DIR}/report \
    -W ${WINDOW} \
    -N ${NAME_MAP} || echo "WARN: finemo report failed (expected on synthetic h5)."

date
