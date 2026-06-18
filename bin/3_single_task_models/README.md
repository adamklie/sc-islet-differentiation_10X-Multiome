# ChromBPNet Single-Task Models Pipeline

Pipeline for training and analyzing single-task ChromBPNet models for each cell type.

## Dependencies

See the [chrombpnet runner README](../../tools/chrombpnet/README.md) for full dependency setup (ChromBPNet, wiggletools, wigToBigWig, gimme, etc.).

Key conda environments:
- `chrombpnet` — model training, predictions, contributions (GPU)
- `eugene_tools` — averaging, motif discovery, clustering (CPU)

## Dataset specifics

- **Cell types:** 22 (DE, PGT1-3, PFG1-2, PP1-2, ENP_phase1, early/late_ENP, early/late_SC_beta, early/late_SC_alpha, early/late_SC_EC, SC_delta_GHRL, proliferating_endocrine, exocrine, liver, FB_FLT1)
- **Folds:** 1 (single fold, fold_0)
- **Bias model:** `early_SC_EC` with `beta=0.5`
- **Peaks:** `results/2_process_data/peak_calls/all/{ct}.narrowPeak`
- **Fragments:** `results/2_process_data/subset/`

## How scripts are generated

All bash scripts in `scripts/` are generated from templates in `tools/chrombpnet/templates/` using the `1_make_scripts_from_templates.ipynb` notebook. **Do not hand-edit scripts** — modify the notebook and re-run to regenerate.

See `tools/chrombpnet/README.md` for detailed step-by-step instructions and SLURM submission patterns.

## Pipeline status

| Step | Description | Status | Notes |
|------|-------------|--------|-------|
| 1. Negatives | Matched negative region sets | Done | Run by Han on external cluster |
| 2. Bias pipeline | Train bias model (early_SC_EC, beta=0.5) | Done | Run by Han on external cluster |
| 3. ChromBPNet pipeline | Train models for 22 cell types | Done | Run by Han on external cluster |
| 4. Predictions | Generate prediction bigwigs | Done | Run by Han on external cluster |
| 5. Average predictions | Average predictions across folds | Done | Single fold → cp |
| 6. Contributions | Generate contribution scores | Done | Run by Han on external cluster |
| 7. Average contributions | Average contribution scores across folds | Done | Single fold → cp for bw, h5 |
| 8. Motif discovery | TF-MoDISco on averaged contributions | Done | 22/22 cell types (PFG2 HTML report failed, h5 OK) |
| 9. Motif clustering | Cluster motifs across cell types (gimme) | Done | 70 clustered motifs, TomTom annotated |
| 10. Motif hit calling | finemo hit calling at peaks | Done | Unified catalog (70 motifs), lambda=0.8 |
| 11. Motif marginalization | Marginalize motifs through models | Done | 22/22 cell types, 70 motifs each |

### Downstream analysis (notebooks 4–7)

| Step | Description | Status | Notes |
|------|-------------|--------|-------|
| 12. Marginalization analysis | Clustermap, ANOVA, profile grids | Done | Notebook 4 |
| 13. Hit analysis | Count/importance heatmaps, QC, composition | Done | Notebook 5, two versions: `hit_analysis/` (per-celltype motifs) and `hit_analysis_unified/` (shared catalog) |
| 14. Motif curation | Quality tiers, TF resolution via expression | Done | Notebook 6, 70 → 59 curated (11 flagged poor/redundant) |
| 15. Co-occurrence | Jaccard pairwise matrices, peak complexity | Done | Notebook 7, 22 per-celltype CSVs |

## Output structure

```
results/3_single_task_models/
├── _flat/                          # Original flat structure from Han's cluster
├── summary/                        # Model performance metrics (CSV, PDF)
├── {ct}/
│   ├── fold_0/chrombpnet/0.5/
│   │   ├── models/                 # chrombpnet.h5, chrombpnet_nobias.h5, bias_model_scaled.h5
│   │   ├── auxiliary/              # peaks, nonpeaks, footprints, interpret_subsample
│   │   ├── evaluation/             # metrics, reports, plots
│   │   ├── logs/                   # training logs
│   │   ├── predictions/            # {ct}_bias.bw, {ct}_chrombpnet.bw, {ct}_chrombpnet_nobias.bw, .bed
│   │   ├── contributions/          # {ct}.counts_scores.{h5,bw}, {ct}.profile_scores.{h5,bw}
│   │   └── marginalization/        # Per-motif marginalization PNGs + data (npz, html)
│   └── average/
│       ├── predictions/            # {ct}_bias.bw, {ct}_chrombpnet.bw, {ct}_chrombpnet_nobias.bw
│       ├── contributions/          # {ct}_counts.{h5,bw}, {ct}_profile.{h5,bw}
│       └── motifs/                 # modisco h5, pfms, hits (per-celltype + unified)
└── motifs/                         # Clustered motifs across all cell types
    ├── cluster/                    # gimme cluster output
    ├── meme/                       # MEME format motifs
    ├── tomtom/                     # TomTom annotation results
    ├── analysis/                   # Marginalization summary (clustermap, ANOVA, TF resolution)
    ├── curation/                   # Curated annotations (70→59), quality tiers, catalog PDF
    ├── cooccurrence/               # Jaccard matrices, peak complexity, per-celltype CSVs
    ├── hit_analysis/               # Hit analysis with per-celltype motifs (original)
    └── hit_analysis_unified/       # Hit analysis with shared 70-motif catalog (canonical)
```

Note: `fold_0/` contents are symlinks to `_flat/` (original Han structure). `average/` is generated by steps 5 and 7. `hit_analysis_unified/` is the canonical hit analysis using the shared motif catalog.

## Running steps

### Step 5: Average predictions

```bash
SLURM_HELPER=/cellar/users/aklie/projects/ML4GLand/chrombpnet/SLURM/cpu.sh
SCRIPT_DIR=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/3_single_task_models/scripts
LOG_DIR=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/3_single_task_models/average_predictions
mkdir -p $LOG_DIR

celltypes=(DE PGT1 PGT2 PGT3 PFG1 PFG2 PP1 PP2 ENP_phase1 early_ENP late_ENP early_SC_beta late_SC_beta early_SC_alpha late_SC_alpha early_SC_EC late_SC_EC SC_delta_GHRL proliferating_endocrine exocrine liver FB_FLT1)

for celltype in ${celltypes[@]}; do
    cmd="$SLURM_HELPER \
-s $SCRIPT_DIR/average_predictions/sc-islet-differentiation_10X-Multiome_${celltype}_average_predictions.sh \
-j sc-islet-diff_${celltype}_avg_pred \
-p carter-compute \
-a carter-compute \
-c 1 \
-m 4G \
-o $LOG_DIR"
    echo -e "Running command:\n$cmd"
    eval $cmd
done
```

### Step 7: Average contributions

```bash
LOG_DIR=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/3_single_task_models/average_contributions
mkdir -p $LOG_DIR

for celltype in ${celltypes[@]}; do
    cmd="$SLURM_HELPER \
-s $SCRIPT_DIR/average_contributions/sc-islet-differentiation_10X-Multiome_${celltype}_average_contributions.sh \
-j sc-islet-diff_${celltype}_avg_contrib \
-p carter-compute \
-a carter-compute \
-c 1 \
-m 4G \
-o $LOG_DIR"
    echo -e "Running command:\n$cmd"
    eval $cmd
done
```

### Step 8: Motif discovery

```bash
cmd="/cellar/users/aklie/projects/ML4GLand/chrombpnet/SLURM/cpu_array.sh \
-s $SCRIPT_DIR/motifs/sc-islet-differentiation_10X-Multiome_motif_discovery.sh \
-j sc-islet-diff_motif_discovery \
-p carter-compute \
-a carter-compute \
-c 32 \
-m 32G \
-o /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/3_single_task_models/motif_discovery \
-n 22 \
-x 5"
echo -e "Running command:\n$cmd"
eval $cmd
```

### Step 9: Motif clustering

```bash
cmd="/cellar/users/aklie/projects/ML4GLand/chrombpnet/SLURM/cpu.sh \
-s $SCRIPT_DIR/motifs/sc-islet-differentiation_10X-Multiome_motif_clustering.sh \
-j sc-islet-diff_motif_clustering \
-p carter-compute \
-a carter-compute \
-c 4 \
-m 16G \
-o /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/3_single_task_models/motif_clustering"
echo -e "Running command:\n$cmd"
eval $cmd
```

### Step 10: Motif hit calling

```bash
cmd="/cellar/users/aklie/projects/ML4GLand/chrombpnet/SLURM/cpu_array.sh \
-s $SCRIPT_DIR/motifs/sc-islet-differentiation_10X-Multiome_motif_hits.sh \
-j sc-islet-diff_motif_hits \
-p carter-compute \
-a carter-compute \
-c 4 \
-m 16G \
-o /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/3_single_task_models/motif_hits \
-n 22 \
-x 5"
echo -e "Running command:\n$cmd"
eval $cmd
```

### Step 11: Motif marginalization

```bash
cmd="/cellar/users/aklie/projects/ML4GLand/chrombpnet/SLURM/gpu_array.sh \
-s $SCRIPT_DIR/motifs/sc-islet-differentiation_10X-Multiome_motif_marginalization.sh \
-j sc-islet-diff_motif_marginalization \
-p carter-gpu \
-a carter-gpu \
-g a30:1 \
-c 4 \
-m 16G \
-t 14-00:00:00 \
-o /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/3_single_task_models/motif_marginalization \
-n 22 \
-x 5"
echo -e "Running command:\n$cmd"
eval $cmd
```
