#!/bin/bash

#####
# Motif hit calling with UNIFIED clustered motif catalog for the
# multi-task (peak regression) CREsted model (finetuned_v2).
#
# Uses the 15-motif catalog at
#   results/4_multi_task_model/crested/motifs/clustered_motifs.modisco.h5
# and the per-CT contributions stored as .npz pairs (oh + contrib) at
#   results/4_multi_task_model/crested/contributions_specific_v2/<CT>_{oh,contrib}.npz
#
# Mirrors bin/3_single_task_models/scripts/motifs/sc-islet-differentiation_10X-Multiome_motif_hits_unified.sh
# but switches the extract step to `finemo extract-regions-modisco-fmt`
# since the multi-task model writes .npz pairs (not a single ChromBPNet h5).
#
# Submit with gpu_array.sh: -n 22 -x 4 -t 14-00:00:00
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
OH=${BASE}/results/4_multi_task_model/crested/contributions_specific_v2/${celltype}_oh.npz
CONTRIB=${BASE}/results/4_multi_task_model/crested/contributions_specific_v2/${celltype}_contrib.npz

# Unified motif catalog (15 motifs)
modisco_h5=${BASE}/results/4_multi_task_model/crested/motifs/clustered_motifs.modisco.h5
name_map=${BASE}/results/4_multi_task_model/crested/motifs/clustered_motifs_name_map.tsv
window=1000

output_dir=${BASE}/results/4_multi_task_model/crested/motifs/hits/${celltype}
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

# Step 3: Report (optional; will fail on synthetic modisco h5 lacking seqlets/* — that's fine,
# downstream aggregation reads hits.tsv directly).
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
