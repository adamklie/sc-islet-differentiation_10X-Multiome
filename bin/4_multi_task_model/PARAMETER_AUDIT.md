# Peak regression (CREsted multi-task) — parameter audit vs CREsted tutorial

**Date:** 2026-05-26
**Tutorial reference:** [enhancer_code_analysis](https://crested.readthedocs.io/en/stable/tutorials/enhancer_code_analysis.html)
**Audit triggered by:** only 15 motifs in the gimme-clustered catalog vs the 5–50/CT range expected from per-CT MoDISco. Root cause confirmed: per-CT MoDISco produces 5–19 patterns each, totaling 241 across 22 CTs and 27K seqlets — but downstream filters (`modisco_to_pfm min_seqlets=100`) cull 68% before clustering.

## Side-by-side audit

| step | parameter | tutorial | first-run (v2) | new run (v3) |
|---|---|---|---|---|
| **sort_and_filter** | `top_k` | 2000 | 1000 | **2000** ✓ |
| **sort_and_filter** | `method` | gini | gini | gini |
| **sort_and_filter** | `model_name` | "combined" | "combined" | "combined" |
| **contribution_scores_specific** | `method` | integrated_grad | integrated_grad | integrated_grad |
| **contribution_scores_specific** | `target_idx` | None (all classes) | per-CT (filtered to CT's channel) | per-CT |
| **tfmodisco** | `window` | 1000 | 400 | **1000** ✓ |
| **tfmodisco** | `max_seqlets` | 20000 | 50000 | 50000 (≥ tutorial) |
| **clustering** | tool | `crested.tl.modisco.process_patterns` | gimme cluster (`-t 0.999`) | **gimme cluster** (intentional — keeps comparability with ChromBPNet single-task catalog) |
| **PFM extraction** | filter | n/a (process_patterns uses IC filters) | `modisco_to_pfm min_seqlets=100` | same — kept for ChromBPNet comparability |
| **model** | checkpoint | the trained CREsted model | `finetuned_v2/checkpoints/28.keras` | same |

## Decisions for the new (v3) run

1. **Adopt tutorial fidelity on upstream params:** `top_k=2000` and MoDISco `window=1000`. These two changes restore tutorial parity at the contribution + seqlet-extraction layer.

2. **Keep gimme cluster + `modisco_to_pfm` downstream:** even though the tutorial uses `crested.tl.modisco.process_patterns` and `create_pattern_matrix`, we keep our existing clustering pipeline so the peak-regression motif catalog can be compared apples-to-apples with the ChromBPNet single-task catalog (which also uses gimme cluster). This is a deliberate deviation from tutorial-pure for cross-oracle comparability.

3. **Marginalization + finemo hit calling:** not covered in the CREsted tutorial — we keep our existing approach (plant motif into k=2 shuffled per-CT background peaks → predict via multi-task model → delta on log-counts; finemo on the per-CT `contributions_specific` arrays with the unified catalog at λ=0.8). Same recipe as ChromBPNet so the per-CT effect-size + hit-count matrices are directly comparable.

## Expected impact

- Doubling `top_k` (1000→2000) → ~2× contribution input → ~2× seqlets per CT
- Widening `window` (400→1000) → 2.5× wider sequence context per region → MoDISco can find more patterns and richer/longer motifs
- Combined ≈ 4–5× more seqlets → many more patterns pass the `min_seqlets=100` cutoff → larger final catalog (target: 30–60 motifs vs current 15)

## Run plan (5-stage SLURM chain)

```
v3-contributions  (GPU, 22-CT array, top_k=2000)
       │
       ▼ afterok
v3-modisco        (CPU, 22-CT array, window=1000, max_seqlets=50000)
       │
       ▼ afterok
v3-cluster        (CPU, single, gimme cluster -t 0.999)
       │
       ▼ afterok
v3-hits           (GPU, 22-CT array, finemo λ=0.8)
       │
       ▼ afterok
v3-marg           (GPU, 22-CT array)
       │
       ▼ afterok
v3-aggregate      (CPU, small)
```

Total realistic wall-time: ~3-4h with `--dependency=afterok` chaining. Outputs land at `results/4_multi_task_model/crested/{contributions_specific_v3, modisco_v3, motifs_v3}/`.

## What stays the same (v2 results not deleted)

- `contributions_specific_v2/` keeps v2 outputs for diff
- `modisco/` keeps v2 per-CT h5s
- `motifs/{combined_mapped.meme, cluster, tomtom, analysis, hits}` keeps v2 catalog + matrices

After v3 completes we'll point the motif-review build script at the v3 paths.
