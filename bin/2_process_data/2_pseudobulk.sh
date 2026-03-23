#!/bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/2_process_data/%x.%A_%a.out
#SBATCH --cpus-per-task=12
#SBATCH --mem=512G
#SBATCH --time=14-00:00:00

#####
# USAGE:
# sbatch --job-name=2_pseudobulk.sh 2_pseudobulk.sh
#####

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configure env
source activate /cellar/users/aklie/opt/miniconda3/envs/scverse-lite-py39
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/lib/

# --- File definitions ---
SCRIPT=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/2_process_data/scripts/generate_bws.py
ANNDATASET=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/subset/_dataset.h5ads
OUTDIR=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/subset
BEDTOOLS=/cellar/users/aklie/opt/miniconda3/envs/chrombpnet/bin/bedtools
BIGWIG=/cellar/users/aklie/opt/ucsc/bedGraphToBigWig

# --- Run pipeline ---
cmd="python $SCRIPT \
--anndataset $ANNDATASET \
--outdir $OUTDIR \
--bedtools $BEDTOOLS \
--bedgraph_to_bw $BIGWIG \
--skip-export \
--threads 8"
echo $cmd
eval $cmd

# Date
date
