#!/bin/bash
# HOMER motif enrichment — array job (one topic per task, 0-indexed)
# Submit: cpu_array.sh -s <this>.sh -j motif_enrich -m 8G -t 04:00:00 -c 4 -n 50 -x 10

date
echo -e "Job ID: $SLURM_JOB_ID, Array Task ID: $SLURM_ARRAY_TASK_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv/bin/activate

BEDS_DIR=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all_replicate/binarized/beds

# Get the Nth BED file (0-indexed from array)
BED_FILES=($(ls -1 ${BEDS_DIR}/*.bed | sort))
BED=${BED_FILES[$SLURM_ARRAY_TASK_ID]}
NAME=$(basename ${BED} .bed)

echo "Processing: ${BED}"

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/annotate/motif_enrichment.py \
    --foreground ${BED} \
    --background /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all_replicate/binarized/consensus_regions.bed \
    --subtract_foreground \
    --max_background auto \
    --genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all_replicate/annotated/motifs/${NAME} \
    --homer_path /cellar/users/aklie/opt/homer/bin/findMotifsGenome.pl \
    --size 200 \
    --n_cpu 4 \
    --mask"
echo -e "Running:\n$cmd\n"
eval $cmd

date
