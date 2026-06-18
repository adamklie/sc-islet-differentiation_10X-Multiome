#!/bin/bash

#####
# v3 motif hit calling with the v3 clustered motif catalog for the
# multi-task (peak regression) CREsted model (finetuned_v2).
#
# Inputs:
#   - per-CT NPZ pairs at contributions_specific_v3/<CT>_{oh,contrib}.npz
#   - unified v3 catalog at motifs_v3/clustered_motifs.modisco.h5 + name_map.tsv
# Outputs:
#   - motifs_v3/hits/<CT>/{hits.tsv, finemo.npz, report}
#
# Submit (per-CT array, 22 cell types, 4 concurrent):
#   sbatch --partition=carter-gpu --account=carter-gpu --gres=gpu:a30:1 \
#     --cpus-per-task=4 --mem=32G --time=14-00:00:00 \
#     --array=1-22%4 --job-name=pr_hits_v3 \
#     --output=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_multi_task_model/%x.%A.%a.out \
#     <this>.sh
#####

date
echo -e "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID\n"

# Configure environment
source activate finemo_gpu
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/opt/miniconda3/envs/finemo_gpu/lib/

celltypes=(
    DE
    ENP_phase1
    FB_FLT1
    PFG1
    PFG2
    PGT1
    PGT2
    PGT3
    PP1
    PP2
    SC_delta_GHRL
    early_ENP
    early_SC_EC
    early_SC_alpha
    early_SC_beta
    exocrine
    late_ENP
    late_SC_EC
    late_SC_alpha
    late_SC_beta
    liver
    proliferating_endocrine
)

# Resolve per-task values
celltype=${celltypes[$SLURM_ARRAY_TASK_ID - 1]}
echo "Cell type: $celltype"

BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome
PEAK=${BASE}/results/2_process_data/peak_calls/all/${celltype}.narrowPeak
OH=${BASE}/results/4_multi_task_model/crested/contributions_specific_v3/${celltype}_oh.npz
CONTRIB=${BASE}/results/4_multi_task_model/crested/contributions_specific_v3/${celltype}_contrib.npz

# Unified v3 motif catalog
modisco_h5=${BASE}/results/4_multi_task_model/crested/motifs_v3/clustered_motifs.modisco.h5
name_map=${BASE}/results/4_multi_task_model/crested/motifs_v3/clustered_motifs_name_map_unique.tsv  # 5 dup TF names suffixed with #<idx> for finemo
window=1000

output_dir=${BASE}/results/4_multi_task_model/crested/motifs_v3/hits/${celltype}
mkdir -p $output_dir
echo "Output: $output_dir"

# Step 1: Extract regions from oh/contrib npz pair
npz_path=$output_dir/${celltype}.finemo.npz
if [ ! -f "$npz_path" ]; then
    cmd="finemo extract-regions-modisco-fmt \
    -s $OH \
    -a $CONTRIB \
    -o $npz_path \
    -w $window"
    echo -e "Running extract:\n$cmd\n"
    eval $cmd
else
    echo "Reusing existing npz: $npz_path"
fi

# Step 2: Call hits against unified catalog
cmd="finemo call-hits \
-r $npz_path \
-m $modisco_h5 \
-o $output_dir \
-p $PEAK \
-N $name_map \
-l 0.8"
echo -e "Running call-hits:\n$cmd\n"
eval $cmd

# Step 3: Report (best-effort — finemo report can fail on synthetic modisco h5,
# downstream aggregation reads hits.tsv directly)
cmd="finemo report \
-r $npz_path \
-p $PEAK \
-H $output_dir/hits.tsv \
-m $modisco_h5 \
-o $output_dir/report \
-W $window \
-N $name_map"
echo -e "Running report (best-effort):\n$cmd\n"
eval $cmd || echo "WARN: finemo report failed (expected for synthetic modisco h5). Continuing."

date
