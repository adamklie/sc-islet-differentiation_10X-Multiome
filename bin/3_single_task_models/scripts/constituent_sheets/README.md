# Constituent sheets

Per-cluster constituent-motif report sheets ("is this cluster chimeric?") for the
`constituent` motif set. For each multi-constituent cluster, renders the consensus +
constituent motifs with marginalization profiles, TomTom hits, Vierstra archetype
candidates, and marg-vs-expression scatters; plus a searchable chimerism index.

## Workspace

Self-contained under the `constituent` marginalization summary dir:
`results/3_single_task_models/motifs/gimme_cluster/marginalization/constituent/`

Inputs already produced by the marginalization pipeline (see `../marginalization/README.md`):
- `query.meme`, `manifest.json`
- `marginalization_matrix.csv`, `marginalization_matrix_zscore.csv`

Produced here: `per_motif/`, `cluster_motifs/`, `clusters.txt`, `tomtom/`, `sheets/`,
`index.html`, `chimerism_summary.tsv`.

## Steps

```bash
PY=/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/python

# 1. Split query.meme into per-motif MEMEs + per-cluster motif lists
$PY split_meme.py

# 2. TomTom each motif vs Vierstra v2.0beta (CPU SLURM array, one task per cluster)
sbatch --array=1-76%8 tomtom_array.sh      # see header for full flags

# 3. Render one sheet per cluster (matplotlib; --index or --cluster)
$PY render_all_sheets.py --index 1         # run per-cluster via array in practice

# 4. Build the searchable chimerism index across all sheets
$PY build_index.py
```

`render_all_sheets.py` imports helper libs from
`scratch/2026_05_25_motif_review/lib/` (logos, TomTom parsing, expression heatmaps).

## Scripts

- `split_meme.py` — query.meme -> per-motif MEMEs + cluster motif lists.
- `tomtom_array.sh` — per-cluster TomTom vs Vierstra archetype DB.
- `render_all_sheets.py` — render one cluster's constituent sheet.
- `build_index.py` — chimerism index (`index.html` + `chimerism_summary.tsv`).
