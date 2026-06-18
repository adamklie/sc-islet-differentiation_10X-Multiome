#!/bin/bash
# Run selected model analysis + topic QC for a given subset
# Usage: run_analysis_pipeline.sh <subset> <n_topics>
# Submit: cpu.sh -s <this>.sh -j analysis_<subset> -m 64G -t 04:00:00 -c 1

SUBSET=$1
N_TOPICS=${2:-50}
DATASET="sc-islet-differentiation_10X-Multiome"
BASE="/cellar/users/aklie/data/datasets/${DATASET}/results/4_topic_models/${SUBSET}"
SCRIPTS="/cellar/users/aklie/data/datasets/${DATASET}/bin/4_topic_models/scripts/analysis"
COLORS="/cellar/users/aklie/data/datasets/${DATASET}/config/cell_type_metadata.tsv"

date
echo "Subset: $SUBSET, n_topics: $N_TOPICS"

# Activate keras venv (has CPU-only torch — avoids CUDA lib errors on CPU nodes)
conda deactivate 2>/dev/null
source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate

# Step 1: Selected model analysis
echo "=== Step 1: Selected model analysis ==="
python $SCRIPTS/run_selected_model_analysis.py \
    --cistopic_obj ${BASE}/objects/${DATASET}_${SUBSET}.pkl \
    --models ${BASE}/models_mallet/merged \
    --output_dir ${BASE}/analysis \
    --n_topics $N_TOPICS \
    --groupby cell_type \
    --colors_tsv $COLORS

# Step 2: Topic QC metrics
echo "=== Step 2: Topic QC metrics ==="
mkdir -p ${BASE}/annotated/topic_qc
python $SCRIPTS/compute_topic_qc.py \
    --cistopic_obj ${BASE}/analysis/cistopic_obj_with_model.pkl \
    --annot_var cell_type \
    --output_dir ${BASE}/annotated/topic_qc

date
