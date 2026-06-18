# deeptopic

# Dependencies

## pycisTopic
You will need a virtualenv with `pycisTopic` installed. Follow the [pycisTopic installation guide](https://pycistopic.readthedocs.io/en/latest/installation.html) or set up the environment as needed.

The current environment used is:
```bash
/cellar/users/aklie/projects/igvf/single_cell_utilities/pycistopic/.venv
```

## Auxiliary scripts
Several Python scripts in the `single_cell_utilities/pycistopic` directory are necessary for running the pipeline steps:
```bash
/cellar/users/aklie/projects/igvf/single_cell_utilities/pycistopic
```

## SLURM submission scripts
SLURM submission wrappers are located at:
```bash
/cellar/users/aklie/opt/SLURM/
```

## MALLET (optional)
If you want to run MALLET-based topic models in addition to CGS, you'll need to install the [MALLET](https://mimno.github.io/Mallet/) toolkit and provide the path to the binary.

# Overview
The workflow is broken down into the following steps:

0. [Make scripts](#0-make-scripts)
1. [Build cisTopic object](#1-build-cistopic-object)
2. [Run topic models (CGS)](#2-run-topic-models-cgs)
3. [Run topic models (MALLET)](#3-run-topic-models-mallet-optional) (optional)
4. [Evaluate models](#4-evaluate-models)
5. [Select model and binarize](#5-select-model-and-binarize)
6. [Train CREsted model](#6-train-crested-model)

The `1_make_scripts_from_templates.ipynb` notebook generates scripts for **all subsets** at once. It reads templates from the `deeptopic/templates/` directory and fills in per-subset paths and shared parameters.

Scripts are organized by step, with one script per subset:
```
scripts/
├── build/
│   ├── endocrine_replicate_build.sh
│   ├── all_build.sh
│   └── ...
├── cgs/
├── mallet/
├── merge_mallet/
├── eval/
└── binarize/
```

These scripts are submitted to the cluster using two SLURM wrappers:
1. `cpu.sh` — Single job. Used for build, CGS, merge, and evaluation steps.
2. `cpu_array.sh` — Array job. Used for MALLET, where each topic count runs as a separate parallel job.

# Steps

## 0. Make scripts
Edit and run `1_make_scripts_from_templates.ipynb` in the dataset's `bin/4_deeptopic/` directory. Configure the subsets list, topic range, and paths, then run all cells.

## 1. Build cisTopic object
Build the cisTopic object from the peak matrix h5ad, optionally filtering blacklisted regions.

Actual usage: ~53G RAM, ~9 min (for 48k cells x 234k regions).

```bash
SUBSET=endocrine_replicate
SCRIPTS=scripts
LOGS=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_deeptopic

cmd="/cellar/users/aklie/opt/SLURM/cpu.sh \
-s $SCRIPTS/build/${SUBSET}_build.sh \
-j build_${SUBSET} \
-c 1 \
-m 96G \
-t 01:00:00 \
-o $LOGS/build"
echo -e "Running command:\n$cmd"
eval $cmd
```

## 2. Run topic models (CGS)
Train LDA topic models using Collapsed Gibbs Sampling across all topic counts in a single job. CGS handles parallelization internally via `n_cpu`.

```bash
cmd="/cellar/users/aklie/opt/SLURM/cpu.sh \
-s $SCRIPTS/cgs/${SUBSET}_cgs.sh \
-j cgs_${SUBSET} \
-c 8 \
-m 128G \
-t 14-00:00:00 \
-o $LOGS/cgs"
echo -e "Running command:\n$cmd"
eval $cmd
```

## 3. Run topic models (MALLET) (optional)
MALLET parallelization happens *within* each model, so we use `cpu_array.sh` to submit one SLURM job per topic count. After all array tasks complete, run the merge script to collect individual model files.

Actual usage per job: 20-31G RAM, 50-100 min wall time (8 CPUs).

```bash
# Step 3a: MALLET array job
# -n = number of topic values, -x = max concurrent tasks
cmd="/cellar/users/aklie/opt/SLURM/cpu_array.sh \
-s $SCRIPTS/mallet/${SUBSET}_mallet.sh \
-j mallet_${SUBSET} \
-c 8 \
-m 64G \
-t 04:00:00 \
-n 17 \
-x 4 \
-o $LOGS/mallet"
echo -e "Running command:\n$cmd"
eval $cmd

# Step 3b: Merge (run after array completes)
cmd="/cellar/users/aklie/opt/SLURM/cpu.sh \
-s $SCRIPTS/merge_mallet/${SUBSET}_merge_mallet.sh \
-j merge_${SUBSET} \
-c 1 \
-m 4G \
-t 00:30:00 \
-o $LOGS/merge_mallet"
echo -e "Running command:\n$cmd"
eval $cmd
```

## 4. Evaluate models
Compute metrics (Arun_2010, Cao_Juan_2009, Minmo_2011, loglikelihood) across topic counts. Optionally generates UMAPs colored by a metadata column.

```bash
cmd="/cellar/users/aklie/opt/SLURM/cpu.sh \
-s $SCRIPTS/eval/${SUBSET}_eval.sh \
-j eval_${SUBSET} \
-c 1 \
-m 32G \
-t 04:00:00 \
-o $LOGS/eval"
echo -e "Running command:\n$cmd"
eval $cmd
```

## 5. Select model and binarize
Select a model (by `--n_topics` or auto-select) and binarize topic-region distributions into per-topic BED files using Otsu thresholding. These BED files are the input for CREsted model training.

```bash
cmd="/cellar/users/aklie/opt/SLURM/cpu.sh \
-s $SCRIPTS/binarize/${SUBSET}_binarize.sh \
-j binarize_${SUBSET} \
-c 1 \
-m 32G \
-t 01:00:00 \
-o $LOGS/binarize"
echo -e "Running command:\n$cmd"
eval $cmd
```

## 6. Train CREsted model
Train a DeepTopic CNN to predict topic assignments from DNA sequences. Uses the binarized BED files from step 5.

```bash
cmd="/cellar/users/aklie/opt/SLURM/gpu.sh \
-s $SCRIPTS/train/${SUBSET}_train.sh \
-j train_${SUBSET} \
-c 4 \
-m 32G \
-t 12:00:00 \
-o $LOGS/train"
echo -e "Running command:\n$cmd"
eval $cmd
```

# Preflight check
Before running the pipeline, you can validate your input h5ad file using the preflight script:

```bash
python preflight.py --h5ad_path <path_to_h5ad>
```

This checks for proper region formatting, duplicate barcodes/regions, NaN/Inf values, and other potential issues.

# Output structure
```
results/4_deeptopic/<subset>/
├── objects/
│   └── <dataset>_<subset>.pkl          # cisTopic object
├── models_cgs/
│   ├── models.pkl                      # Trained CGS models (all topics)
│   └── models_params.yaml
├── models_mallet/                      # (if MALLET was run)
│   ├── topic_2/                        # one subdir per n_topics value
│   │   ├── models.pkl
│   │   └── models_params.yaml
│   ├── topic_5/
│   │   └── ...
│   └── merged/                         # flat dir for evaluate_models.py
│       ├── topic_2.pkl
│       ├── topic_5.pkl
│       └── ...
├── evaluation/
│   ├── model_evaluation_metrics.pdf
│   ├── model_evaluation_metrics.csv
│   ├── umap_per_topic_count.pdf        # (if groupby specified)
│   └── evaluation_params.yaml
└── binarized/
    ├── beds/
    │   ├── Topic_1.bed
    │   ├── Topic_2.bed
    │   └── ...
    ├── consensus_regions.bed
    ├── binarization_thresholds.pdf
    └── binarization_params.yaml
```
