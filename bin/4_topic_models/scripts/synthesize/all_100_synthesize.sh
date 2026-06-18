#!/bin/bash
# Synthesize all topic annotations for all_100
# Submit: cpu.sh -s <this>.sh -j synthesize_100 -m 4G -t 00:30:00 -c 1

date
source /cellar/users/aklie/opt/deeptopic/.venv/bin/activate

ANNOT_DIR=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/n100/annotated

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/annotate/synthesize_annotations.py \
    --motif_summary ${ANNOT_DIR}/motifs/motif_enrichment_summary.tsv \
    --gene_summary ${ANNOT_DIR}/genes/gene_enrichment_summary.tsv \
    --metadata_top ${ANNOT_DIR}/metadata/cell_type_top.tsv \
    --metadata_column cell_type \
    --output_dir ${ANNOT_DIR}/synthesis \
    --use_llm"
echo -e "Running:\n$cmd\n"
eval $cmd

date
