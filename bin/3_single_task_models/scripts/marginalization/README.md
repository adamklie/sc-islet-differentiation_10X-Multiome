# Marginalization pipeline (single-task models)

In-silico marginalization of motif sets against the 22 per-cell-type single-task
ChromBPNet models, producing a `motif × 22 cell-type` Δcount summary matrix per motif set.

## Motif sets

| set            | # motifs | query source                                    | build script              |
|----------------|----------|-------------------------------------------------|---------------------------|
| `gimme_average`| 70       | original gimme-cluster pipeline (pre-existing)  | (legacy, not in this dir) |
| `consensus`    | 76       | one PPM per final cluster (`08_final_averages.mc`) | `build_consensus_meme.py` |
| `constituent`  | 751      | consensus + constituents of every multi-constituent cluster | `build_inputs_from_mc.py` |

## Data layout

- Query MEMEs (one per set):
  `results/3_single_task_models/motifs/gimme_cluster/marginalization/<set>/query.meme`
- Per-cell-type marginalization NPZs:
  `results/3_single_task_models/<CT>/fold_0/chrombpnet/0.5/marginalization/<set>/marginalization_data.{counts,motifs,profiles,X}.npz`
- Summary matrices:
  `results/3_single_task_models/motifs/gimme_cluster/marginalization/<set>/marginalization_matrix.csv`
  (+ `marginalization_matrix_zscore.csv`)

`results` root = `/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results`.

## Pipeline

1. **Build the query MEME** (env `motifcompendium-gpu` — needs `MotifCompendium`):
   ```bash
   MCPY=/cellar/users/aklie/opt/miniconda3/envs/motifcompendium-gpu/bin/python
   $MCPY build_consensus_meme.py     # -> consensus/query.meme   (76 motifs)
   $MCPY build_inputs_from_mc.py     # -> constituent/query.meme (751) + manifest.json
   ```

2. **Marginalize across the 22 models** (GPU SLURM array; env `eugene_tools`).
   Reads `<set>/query.meme`, writes per-CT NPZs to each model's `marginalization/<set>/`:
   ```bash
   sbatch --array=1-22%10 marginalize_consensus.sh   # consensus
   sbatch --array=1-22%10 marginalize_global.sh      # constituent
   ```
   Full sbatch flags (account/partition/gpus/output) are in each script's header.
   Logs: `bin/slurm_logs/3_single_task_models/marginalization/<set>/`.

3. **Aggregate** to the summary matrix (cheap CPU; env `eugene_tools`):
   ```bash
   PY=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python
   $PY aggregate_global.py --motif-set consensus     # -> consensus/marginalization_matrix.csv  (76 x 22)
   $PY aggregate_global.py --motif-set constituent   # -> constituent/marginalization_matrix.csv (751 x 22)
   ```
   `score = mean_over_seqs(c_after - c_before)`. A non-fatal orthogonal check compares
   shared rows against the `gimme_average` matrix (skipped when motif names don't overlap,
   which is expected: new `cl<cf>_*` names vs old `Average_<NNN>`).

## Scripts

- `build_consensus_meme.py` — write the 76-motif `consensus` query MEME.
- `build_inputs_from_mc.py` — write the 751-motif `constituent` query MEME + manifest.json
  (provides the shared `parse_ppm_meme` / `block` / `HEADER` / `MC` / `MARG` helpers).
- `marginalize_consensus.sh` / `marginalize_global.sh` — GPU SLURM array marginalization.
- `aggregate_global.py` — `--motif-set {consensus,constituent}` summary matrix + zscore.
