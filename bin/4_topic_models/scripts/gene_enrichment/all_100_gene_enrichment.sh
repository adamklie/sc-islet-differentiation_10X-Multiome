#!/bin/bash
# Gene enrichment for all_100
# Submit: cpu.sh -s <this>.sh -j gene_enrich_100 -m 8G -t 04:00:00 -c 1

date
source /cellar/users/aklie/opt/deeptopic/.venv/bin/activate

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/annotate/gene_enrichment.py \
    --beds_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/n100/binarized/beds \
    --genome hg38 \
    --organism hsapiens \
    --output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/n100/annotated/genes \
    --homer_path /cellar/users/aklie/opt/homer/bin/annotatePeaks.pl"
echo -e "Running:\n$cmd\n"
eval $cmd

date
