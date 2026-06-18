#!/bin/bash
# Metadata enrichment for all_100
# Submit: cpu.sh -s <this>.sh -j meta_enrich_100 -m 64G -t 01:00:00 -c 1
# Note: depends on analysis step completing first (needs cistopic_obj_with_model.pkl)

date
source /cellar/users/aklie/opt/deeptopic/.venv/bin/activate

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/annotate/metadata_enrichment.py \
    --cistopic_obj /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/objects/sc-islet-differentiation_10X-Multiome_all.pkl \
    --models /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/models_mallet/merged \
    --n_topics 100 \
    --output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/n100/annotated/metadata \
    --columns cell_type"
echo -e "Running:\n$cmd\n"
eval $cmd

date
