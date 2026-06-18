#!/bin/bash
# Topic-metadata association analysis
# Submit: cpu.sh -s <this>.sh -j meta_enrich -m 64G -t 01:00:00 -c 1

date
echo -e "Job ID: $SLURM_JOB_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv/bin/activate

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/annotate/metadata_enrichment.py \
    --cistopic_obj /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/objects/sc-islet-differentiation_10X-Multiome_all.pkl \
    --models /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/models_mallet/merged \
    --n_topics 50 \
    --output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/annotated/metadata \
    --columns cell_type"
echo -e "Running:\n$cmd\n"
eval $cmd

date
