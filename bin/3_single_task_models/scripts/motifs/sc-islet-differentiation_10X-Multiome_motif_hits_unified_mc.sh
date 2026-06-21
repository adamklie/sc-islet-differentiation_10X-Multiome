#!/bin/bash

#####
# Motif hit calling with the UNIFIED MotifCompendium ("mc") clustered catalog
# (motif_compendium_cluster, 76 motifs named cluster_final#N).
#
# Mirrors ..._motif_hits_unified.sh but:
#   - points finemo at the mc modisco catalog instead of the gimme one
#   - reads per-CT contributions from the post-reorg models/<ct>/ layout
#   - writes to counts_unified_mc/ so it does NOT clobber the gimme hits
#
# The mc modisco h5 pattern keys ARE the final names (cluster_final#N), so no
# name_map is passed — finemo emits motif_name = cluster_final#N directly, which
# matches the catalog display_name used by Phase 6/7.
#
# Submit (one CT to test, then the full array):
#   sbatch --partition=carter-gpu --account=carter-gpu --gpus=a30:1 \
#     --mem=128G -n 4 -t 14-00:00:00 --array=1%1   ...  (DE only, test)
#   sbatch ... --array=1-22%8   ...                     (all 22)
#####

date
echo -e "Job ID: $SLURM_JOB_ID  Task: $SLURM_ARRAY_TASK_ID\n"

# Configure environment
source activate finemo_gpu
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/envs/finemo_gpu/lib/

BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome
ST=$BASE/results/3_single_task_models
PEAKDIR=$BASE/results/2_process_data/peak_calls/all

celltypes=(
	DE ENP_phase1 FB_FLT1 PFG1 PFG2 PGT1 PGT2 PGT3 PP1 PP2 SC_delta_GHRL
	early_ENP early_SC_EC early_SC_alpha early_SC_beta exocrine
	late_ENP late_SC_EC late_SC_alpha late_SC_beta liver proliferating_endocrine
)

# Unified mc motif catalog (76 motifs, keys = cluster_final#N)
modisco_h5=$ST/motifs/motif_compendium_cluster/clustered_motifs.modisco.h5
window=400

# SLURM task-specific values
celltype=${celltypes[$SLURM_ARRAY_TASK_ID - 1]}
peak=$PEAKDIR/${celltype}.narrowPeak
input_h5=$ST/models/${celltype}/average/contributions/${celltype}_counts.h5
output_dir=$ST/models/${celltype}/average/motifs/hits/counts_unified_mc

mkdir -p $output_dir
echo "Cell type: $celltype"
echo "Contrib:   $input_h5"
echo "Peaks:     $peak"
echo "Catalog:   $modisco_h5"
echo "Output:    $output_dir"

# Step 1: Extract regions (reuse the gimme run's npz if present — same contributions/window)
npz_path=$output_dir/${celltype}.counts.finemo.npz
gimme_npz=$ST/models/${celltype}/average/motifs/hits/counts_unified/${celltype}.counts.finemo.npz
if [ -f "$npz_path" ]; then
    echo "Reusing existing npz: $npz_path"
elif [ -f "$gimme_npz" ]; then
    echo "Linking gimme npz (same contributions+window): $gimme_npz"
    ln -s "$gimme_npz" "$npz_path"
else
    cmd="finemo extract-regions-chrombpnet-h5 -c $input_h5 -o $npz_path -w $window"
    echo -e "Running extract:\n$cmd\n"; eval $cmd
fi

# Step 2: Call hits against the mc catalog
cmd="finemo call-hits -r $npz_path -m $modisco_h5 -o $output_dir -p $peak -l 0.8"
echo -e "Running call-hits:\n$cmd\n"; eval $cmd

# Step 3: Report
cmd="finemo report -r $npz_path -p $peak -H $output_dir/hits.tsv -m $modisco_h5 -o $output_dir/report -W $window"
echo -e "Running report:\n$cmd\n"; eval $cmd

date
