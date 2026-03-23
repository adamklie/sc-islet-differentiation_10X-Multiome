#!/bin/bash
# Merge individual MALLET model .pkl files into a single directory for evaluate_models.py
# Submit: /cellar/users/aklie/opt/SLURM/cpu.sh -s <this>.sh -j cistopic_merge -m 4G -t 00:30:00 -c 1
# Tip: use --dependency=afterok:<array_job_id> when submitting after the MALLET array job

date
echo -e "Job ID: $SLURM_JOB_ID\n"

MALLET_OUTPUT_DIR=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_deeptopic/non-endocrine_replicate/models_mallet
MERGED_DIR=$MALLET_OUTPUT_DIR/merged

echo -e "MALLET output dir: $MALLET_OUTPUT_DIR"
echo -e "Merged dir: $MERGED_DIR\n"

mkdir -p "$MERGED_DIR"

count=0
for model_dir in "$MALLET_OUTPUT_DIR"/topic_*/; do
    topic_name=$(basename "$model_dir")
    src="$model_dir/models.pkl"
    if [ -f "$src" ]; then
        cp "$src" "$MERGED_DIR/${topic_name}.pkl"
        echo "Copied: $src -> $MERGED_DIR/${topic_name}.pkl"
        count=$((count + 1))
    else
        echo "WARNING: $src not found, skipping"
    fi
done

echo -e "\nMerged $count model files into $MERGED_DIR"
echo "Use with: evaluate_models.py --models $MERGED_DIR"

date
