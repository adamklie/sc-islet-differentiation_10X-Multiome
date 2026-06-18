#!/bin/bash
#####
# One task per cluster: run TomTom (per motif) for that cluster's constituents+consensus.
# USAGE:
# sbatch \
# --job-name=tomtom_constituent \
# --partition carter-compute \
# --mem=8G \
# -n 1 \
# -t 02-00:00:00 \
# --array=1-76%8 \
# --output /carter/users/aklie/projects/islet_organoid_differentiation/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/3_single_task_models/constituent_sheets/%x.%A.%a.out \
# tomtom_array.sh
#####
set -e
date; echo "Task $SLURM_ARRAY_TASK_ID"
TT=/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/bin/tomtom
DB=/cellar/users/aklie/data/ref/motifs/jvierstra/motif-clustering-v2.0beta/motifs.meme
CONST=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/motifs/gimme_cluster/marginalization/constituent
PM=$CONST/per_motif
OUT=$CONST/tomtom; mkdir -p $OUT/.oc
mapfile -t clusters < $CONST/clusters.txt
cl=${clusters[$SLURM_ARRAY_TASK_ID - 1]}
echo "cluster $cl"
n=0
while read -r name; do
  [ -z "$name" ] && continue
  $TT -no-ssc -oc $OUT/.oc/$name -verbosity 1 -text -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10.0 \
     $PM/$name.meme $DB > $OUT/$name.tomtom.txt 2>>$OUT/.oc/tomtom.$SLURM_ARRAY_TASK_ID.log
  n=$((n+1))
done < $CONST/cluster_motifs/$cl.txt
echo "ran $n tomtoms for $cl"; date
