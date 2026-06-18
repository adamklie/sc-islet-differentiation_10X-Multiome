#!/bin/bash
# Analysis + QC for all/n100 (100 topics on "all" subset)
# Submit: cpu.sh -s <this>.sh -j analysis_all100 -m 64G -t 04:00:00 -c 1

date
echo "Running analysis pipeline for all/n100 (100 topics)"

conda deactivate 2>/dev/null
source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate

DATASET="sc-islet-differentiation_10X-Multiome"
BASE="/cellar/users/aklie/data/datasets/${DATASET}/results/4_topic_models"
SCRIPTS="/cellar/users/aklie/data/datasets/${DATASET}/bin/4_topic_models/scripts/analysis"
COLORS="/cellar/users/aklie/data/datasets/${DATASET}/config/cell_type_metadata.tsv"

# Step 1: Selected model analysis
echo "=== Step 1: Selected model analysis ==="
python $SCRIPTS/run_selected_model_analysis.py \
    --cistopic_obj ${BASE}/all/objects/${DATASET}_all.pkl \
    --models ${BASE}/all/models_mallet/merged \
    --output_dir ${BASE}/all/n100/analysis \
    --n_topics 100 \
    --groupby cell_type \
    --colors_tsv $COLORS

# Step 2: Topic QC metrics
echo "=== Step 2: Topic QC metrics ==="
python $SCRIPTS/compute_topic_qc.py \
    --cistopic_obj ${BASE}/all/n100/analysis/cistopic_obj_with_model.pkl \
    --annot_var cell_type \
    --output_dir ${BASE}/all/n100/annotated/topic_qc

date
