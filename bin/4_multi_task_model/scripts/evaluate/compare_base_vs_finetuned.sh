#!/bin/bash
# Compare base vs finetuned on full and cell-type-specific test sets
# Submit: gpu.sh -s <this>.sh -j compare_models -m 64G -t 04:00:00 -c 4 -g a30:1

date
echo -e "Job ID: $SLURM_JOB_ID\n"

source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
NVIDIA_LIBS=$VIRTUAL_ENV/lib/python3.11/site-packages/nvidia
export LD_LIBRARY_PATH=$NVIDIA_LIBS/cublas/lib:$NVIDIA_LIBS/cuda_runtime/lib:$NVIDIA_LIBS/cudnn/lib:$NVIDIA_LIBS/cufft/lib:$NVIDIA_LIBS/curand/lib:$NVIDIA_LIBS/cusolver/lib:$NVIDIA_LIBS/cusparse/lib:$NVIDIA_LIBS/nccl/lib:$NVIDIA_LIBS/nvjitlink/lib:$LD_LIBRARY_PATH
export KERAS_BACKEND=tensorflow

python << 'PYEOF'
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import crested
import keras
import pysam

base_dir = '/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_multi_task_model'
genome_fa = '/cellar/users/aklie/data/ref/genomes/hg38/hg38.fa'
chrom_sizes = '/cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes'

genome = crested.Genome(genome_fa, chrom_sizes)
crested.register_genome(genome)

print("Importing BigWigs...")
adata = crested.import_bigwigs(
    bigwigs_folder=f'{base_dir}/bigwigs',
    regions_file='/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/peak_calls/all/consensus_peaks.bed',
    chromsizes_file=chrom_sizes,
    target='count',
)
print(f"  {adata.shape}")

crested.pp.train_val_test_split(adata, val_chroms=['chr8', 'chr20'], test_chroms=['chr1', 'chr3', 'chr6'])
crested.pp.normalize_peaks(adata, top_k_percent=0.03)

test_mask = adata.var['split'] == 'test'
X_test = adata[:, test_mask].X
test_regions = adata.var_names[test_mask].tolist()
cell_types = adata.obs_names.tolist()
n_test = test_mask.sum()
print(f"Test: {n_test} regions, {len(cell_types)} cell types")

# Gini for cell-type-specific filtering
def gini(x):
    x = np.sort(np.abs(x))
    n = len(x)
    if x.sum() == 0: return 0
    idx = np.arange(1, n + 1)
    return (2 * np.sum(idx * x) - (n + 1) * np.sum(x)) / (n * np.sum(x))

gini_scores = np.array([gini(X_test[:, i]) for i in range(n_test)])
threshold = np.mean(gini_scores) + 1.0 * np.std(gini_scores)
specific_mask = gini_scores > threshold
print(f"Specific: {specific_mask.sum()}/{n_test} ({100*specific_mask.mean():.1f}%)")

# Load models
print("Loading base model...")
base_model = keras.saving.load_model(f'{base_dir}/crested/basemodel/checkpoints/15.keras', compile=False)
print("Loading finetuned model...")
fine_model = keras.saving.load_model(f'{base_dir}/crested/finetuned/checkpoints/23.keras', compile=False)

# Predict
fasta = pysam.FastaFile(genome_fa)
seq_len = 2114

def one_hot(seq):
    m = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    arr = np.zeros((len(seq), 4), dtype=np.float32)
    for i, b in enumerate(seq.upper()):
        if b in m: arr[i, m[b]] = 1
    return arr

def parse_region(r):
    p = r.replace('-', ':').split(':')
    c = (int(p[1]) + int(p[2])) // 2
    return p[0], c - seq_len // 2, c + seq_len // 2

batch_size = 512
base_preds, fine_preds = [], []
n_batches = (n_test + batch_size - 1) // batch_size
for b in range(0, n_test, batch_size):
    e = min(b + batch_size, n_test)
    seqs = []
    for i in range(b, e):
        ch, s, en = parse_region(test_regions[i])
        seq = fasta.fetch(ch, max(0, s), en)
        if len(seq) < seq_len: seq += 'N' * (seq_len - len(seq))
        seqs.append(one_hot(seq[:seq_len]))
    X = np.stack(seqs)
    base_preds.append(base_model.predict(X, verbose=0))
    fine_preds.append(fine_model.predict(X, verbose=0))
    batch_num = b // batch_size
    if batch_num % 25 == 0:
        print(f"  Batch {batch_num}/{n_batches}")

fasta.close()
base_preds = np.concatenate(base_preds)
fine_preds = np.concatenate(fine_preds)
print(f"Predictions: {base_preds.shape}")

# Compare
results = []
for i, ct in enumerate(cell_types):
    results.append({
        'cell_type': ct,
        'base_r_full': pearsonr(X_test[i], base_preds[:, i])[0],
        'fine_r_full': pearsonr(X_test[i], fine_preds[:, i])[0],
        'base_rho_full': spearmanr(X_test[i], base_preds[:, i])[0],
        'fine_rho_full': spearmanr(X_test[i], fine_preds[:, i])[0],
        'base_r_specific': pearsonr(X_test[i, specific_mask], base_preds[specific_mask, i])[0],
        'fine_r_specific': pearsonr(X_test[i, specific_mask], fine_preds[specific_mask, i])[0],
        'base_rho_specific': spearmanr(X_test[i, specific_mask], base_preds[specific_mask, i])[0],
        'fine_rho_specific': spearmanr(X_test[i, specific_mask], fine_preds[specific_mask, i])[0],
    })

df = pd.DataFrame(results)
df['diff_full'] = df['fine_r_full'] - df['base_r_full']
df['diff_specific'] = df['fine_r_specific'] - df['base_r_specific']
df = df.sort_values('diff_specific', ascending=False)
df.to_csv(f'{base_dir}/crested/base_vs_finetuned_comparison.tsv', sep='\t', index=False)

print(f"\n{'cell_type':30s}  {'base_full':>9s} {'fine_full':>9s} {'diff':>7s}  ||  {'base_spec':>9s} {'fine_spec':>9s} {'diff':>7s}")
print("-" * 105)
for _, r in df.iterrows():
    m = " ***" if r['cell_type'].startswith('late_SC') else ""
    print(f"{r['cell_type']:30s}  {r['base_r_full']:9.3f} {r['fine_r_full']:9.3f} {r['diff_full']:+7.3f}  ||  {r['base_r_specific']:9.3f} {r['fine_r_specific']:9.3f} {r['diff_specific']:+7.3f}{m}")

print(f"\nMean full:     base={df['base_r_full'].mean():.3f}  fine={df['fine_r_full'].mean():.3f}  diff={df['diff_full'].mean():+.3f}")
print(f"Mean specific: base={df['base_r_specific'].mean():.3f}  fine={df['fine_r_specific'].mean():.3f}  diff={df['diff_specific'].mean():+.3f}")
print(f"\nImproved on full:     {(df['diff_full'] > 0).sum()}/{len(df)}")
print(f"Improved on specific: {(df['diff_specific'] > 0).sum()}/{len(df)}")
PYEOF

date
