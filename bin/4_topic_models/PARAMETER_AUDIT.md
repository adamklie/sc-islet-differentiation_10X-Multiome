# Topic classification (CREsted DeepTopic n100) — parameter audit vs CREsted tutorials

**Date:** 2026-05-26
**Tutorial references:**
- [topic_classification](https://crested.readthedocs.io/en/stable/tutorials/topic_classification.html) — training + architecture
- [enhancer_code_analysis](https://crested.readthedocs.io/en/stable/tutorials/enhancer_code_analysis.html) — downstream contributions / MoDISco / pattern processing (CREsted has no topic-specific downstream tutorial; both tutorials share the same `contribution_scores_specific` + `tfmodisco` recipe so the EC tutorial is the de facto reference)
**Trained model:** `results/4_topic_models/all/n100/crested_model/checkpoints/60.keras` (best per `evaluation/best_checkpoint.txt`)
**Audit triggered by:** parity check against the MT pipeline audit (`bin/4_multi_task_model/PARAMETER_AUDIT.md`) before adding topic-level marginalization. Headline result: n100 pipeline is **largely tutorial-faithful on training**, but drifts on the downstream contributions/MoDISco params in the same direction the MT v2 pipeline did — `top_k=1000` (tut: 2000) and MoDISco `window=400` (tut: 1000). Final motif catalog is 51 (vs 15 for MT v2), so the cumulative drift is less severe than MT v2 (likely because n100 has many more output channels and broader seqlet pools per topic).

## CNN head activation — confirmed

`crested.tl.zoo.deeptopic_cnn(..., output_activation="sigmoid")` is the default and is what we use (we pass no override). Loss is `BinaryCrossentropy(from_logits=False)` from `default_configs("topic_classification")`. **Outputs are independent per-topic probabilities in [0,1]** (multi-label). Marginalization `c_after - c_before` is a delta in sigmoid probability points per topic and is directly comparable across topics for the same motif.

(Not softmax — no zero-sum constraint across topics, no need for special interpretation.)

## Side-by-side audit

### Training
| param | tutorial | ours (n100) | drift |
|---|---|---|---|
| `seq_len` | 500 | **501** | OK (501 is odd, matches our consensus regions; CREsted accepts both) |
| `num_classes` | 80 | 100 | dataset choice |
| `output_activation` | sigmoid | sigmoid | OK |
| `loss` | BinaryCrossentropy | BinaryCrossentropy | OK |
| `optimizer` | Adam(lr=1e-3) | Adam(lr=1e-3) | OK |
| `batch_size` | 128 | 128 | OK |
| `epochs` | 100 | 100 | OK |
| `max_stochastic_shift` | 3 | 3 | OK |
| `always_reverse_complement` | True | True | OK |
| val/test chromosomes | chr8,chr10 / chr9,chr18 | same | OK |
| best epoch | n/a | 60 | — |
| test auROC / auPR | n/a | 0.826 / 0.306 | (catacc 0.07 is meaningless for multi-label) |

### Contributions (per-topic)
| param | tutorial (enhancer_code_analysis) | ours `all_n100_contributions_specific.sh` | drift |
|---|---|---|---|
| `sort_and_filter_regions_on_specificity.top_k` | **2000** | **1000** | DRIFT (same as MT v2) |
| `sort_and_filter_regions_on_specificity.method` | gini | gini | OK |
| `sort_and_filter_regions_on_specificity.model_name` | "combined" | "combined" | OK |
| `contribution_scores_specific.method` | integrated_grad | integrated_grad | OK |
| `contribution_scores_specific.target_idx` | None (all classes) | per-topic (filtered to that topic's channel) | drift in spirit but identical numerically — we compute per-topic specifically rather than all-at-once for memory + per-topic NPZ output |
| `batch_size` | unspecified (notebook default) | 128 | OK |
| backend | torch | torch | OK |

There is **also a second contribution path** in this pipeline:
`all_n100_contributions_array.sh` → `compute_contribution_scores.py` (uses `--method expected_integrated_grad`, `--max_regions 1000`). This is the *array* job that feeds MoDISco via the projected/hypothetical NPZ format. The two scripts coexist:

- `contributions_array.sh` (expected_integrated_grad, max_regions=1000) → MoDISco
- `contributions_specific.sh` (integrated_grad, top_k=1000) → only writes adata predictions + per-topic NPZs not used downstream

The MoDISco-feeding path is the EIG/1000-region one, **which deviates from the CREsted tutorial's `contribution_scores_specific` recipe** but matches the original DeepTopic recipe and the ChromBPNet single-task pipeline.

### MoDISco
| param | tutorial | ours `run_modisco_all_topics.sh` | drift |
|---|---|---|---|
| `-n` / max_seqlets | 20000 | **50000** | OK (we sample more aggressively; tutorial floor) |
| `-w` / window | **1000** | **400** | DRIFT (same as MT v2) |
| `-l` / n_leiden | 2 (CLI default) | 2 | OK |
| input h5 schema | CREsted-internal | ChromBPNet-style h5 (`projected_shap/seq`, `shap/seq`, `raw/seq`, float16/int8) | intentional — we convert to the ChromBPNet schema before invoking `modisco motifs` so the downstream PFM extraction tooling is shared with the ChromBPNet pipeline |

### Pattern → PFM → catalog
| step | tutorial | ours | drift |
|---|---|---|---|
| pattern processing | `crested.tl.modisco.process_patterns` (sim_threshold=6.5, trim_ic=0.05, discard_ic=0.2) | `modisco_to_pfm.py` (`-t 0.45`, `-s 100` min_seqlets, `-f 2` flank) | intentional swap — keeps comparability with the ChromBPNet + MT catalogs which use the same script |
| pattern matrix | `crested.tl.modisco.create_pattern_matrix` | gimme `cluster -t 0.999` → `pfm_to_meme.py` → TomTom against Vierstra | intentional — gimme cluster + Vierstra TomTom is the shared catalog-building stack |
| final catalog | (in-memory pattern matrix) | `motifs/combined_mapped.meme` (**51 motifs**) | — |
| min_seqlets cutoff | 100–300 (visualization-dependent) | 100 (`modisco_to_pfm -s 100`) | OK (matches tutorial floor) |

51 motifs for 100 topics ≈ 0.5 motifs/topic catalog density. For reference, MT v2 was 15 motifs / 22 CTs ≈ 0.7/CT, and MT v3 (with `top_k=2000` + `window=1000`) is expected to land in 30–60. The n100 catalog is already richer than MT v2 in absolute terms, likely because (a) 100 topics → more independent contribution signals, (b) `max_seqlets=50000` was applied here from the start.

### Cooccurrence + downstream
| step | tutorial | ours | drift |
|---|---|---|---|
| finemo hit calling | not covered by tutorial | finemo with combined_mapped.meme at λ=0.8 | inherited from ChromBPNet pipeline; matches MT |
| motif activity / cooccurrence | not covered | bespoke notebooks `5_*_n100`, `6_*_n100`, `7_*_n100` | — |
| marginalization | not covered by tutorial | **does not yet exist for topic CNN** — Task B in this audit cycle | — |

## Drift summary

The drift profile is **identical to MT v2**:
1. `top_k=1000` instead of tutorial 2000 (halves the contribution input to MoDISco)
2. MoDISco `window=400` instead of tutorial 1000 (2.5× narrower context)

Because n100 still produced 51 catalog motifs (~3× MT v2's 15), the practical cost of drift looks smaller — but it is the same drift, and an n100 v3-style rerun (top_k=2000, window=1000) would likely push the catalog into the 80–150 motif range and reveal motifs with weaker / longer-range support that the current run misses.

## Decision: **no v3 rerun for n100 at this time**

Rationale:
- The 51-motif catalog already meets downstream needs (hits, cooccurrence, activity matrices were all produced and looked clean in `5/6/7_*_n100.ipynb`).
- Adam has not flagged a missing-motif problem for n100 — the impetus for the MT v3 rerun was the 15-motif sparsity; n100 does not have that symptom.
- A v3 rerun would be CPU+GPU expensive (100-topic arrays, MoDISco compute scales linearly with `max_seqlets * window`). Not worth the cost without a concrete motif-coverage gap to close.
- The audit is recorded here so we can revisit if a missing-motif problem appears in topic-level downstream analyses (e.g. marginalization or interpretation notebooks).

If/when we do an n100 v3:
- `top_k=2000` in `contributions_array.sh` and `contributions_specific.sh`
- `-w 1000` in `run_modisco_all_topics.sh`
- Keep `-n 50000` (already ≥ tutorial)
- Keep gimme cluster + `modisco_to_pfm` (cross-oracle comparability)
- Land outputs at `crested_model/{contributions_v3, contributions_specific_v3, modisco_v3, motifs_v3}/`

## Marginalization (Task B addition)

Topic-level marginalization scripts have been added at `bin/4_topic_models/scripts/marginalization/`:
- `marginalize_topic.py` — plant motif into k=2-shuffled per-topic background sequences (sourced from `binarized/beds/Topic{N}.bed`), predict sigmoid probability for the target topic before/after.
- `sc-islet-differentiation_10X-Multiome_topic_marginalization.sh` — 100-topic SLURM array
- `aggregate_marginalization.py` — collapses per-topic NPZs into `motifs/analysis/marginalization_matrix{,_zscore}.csv` (motif × 100 topics)
- `run_aggregate.sh` — SLURM wrapper for the aggregator (chained `afterok` on the array)

Output keys match the MT marginalization schema (`c_before`, `c_after`, `names`) so the renderer/notebook code is reusable, with the caveat that the *units* are sigmoid probability points per topic (not log-counts as in MT/ChromBPNet). See "CNN head activation" section above.

## v2 pilot (2026-05-28) — region count + method, on 5 topics (1, 9, 16, 44, 80)

**Goal:** vet whether a full n100 v3 rebuild at `max_regions=2000` (vs current 1000) is worth the compute, before committing 100-topic arrays. Pilot scripts:
- `scripts/contributions/sc-islet-differentiation_10X-Multiome_topic_pilot_v2.sh` — 2000 regions, **integrated_grad** → `crested_model/{contributions_pilot_v2, modisco_pilot_v2}/`
- `scripts/contributions/sc-islet-differentiation_10X-Multiome_topic_baseline_ig_1000.sh` — matched 1000-region **integrated_grad** baseline → `crested_model/{contributions_baseline_ig_1000, modisco_baseline_ig_1000}/`

**Confound discovered mid-pilot:** the v1 *bundle* catalog (`modisco/` → `motifs/combined_mapped.meme`) was fed by `contributions_array.sh`, which uses **`expected_integrated_grad` @ 1000** (see "second contribution path" note above). The pilot script had been seeded from that array script but was first run with `expected_integrated_grad` and **hung at `Model: 0%` for the full 6 h cap** — expected-gradients samples many baselines per region (~100× a forward pass), so 2000 regions is intractable in 6 h (≈16 h+/topic). Switching the pilot to plain `integrated_grad` made it complete in ~28–32 min/topic (CPU-bound; the contribution step runs on CPU here regardless of GPU allocation). So the pilot tests `integrated_grad`, NOT the v1 method.

**Clean region-count comparison (both `integrated_grad`, 1000 vs 2000):**

| Topic | 1000-IG patterns | 2000-IG patterns | Δ |
|---|---|---|---|
| Topic1 | 11 | 12 | +1 |
| Topic9 | 20 | 27 | +7 |
| Topic16 | 10 | 11 | +1 |
| Topic44 | 6 | 6 | +0 |
| Topic80 | 6 | 8 | +2 |
| **total** | **53** | **64** | **+21%** |

Seqlet support roughly doubles (mechanical from 2× regions). Net: **+21% patterns + ~2× seqlet support** — a modest, real sensitivity gain concentrated in already-rich topics (Topic9); Topic44 gained nothing.

**Method effect (the bigger finding) — `expected_IG` vs `IG`, both @ 1000 regions:**

| Topic | expected-IG (v1) pos/neg | IG pos/neg |
|---|---|---|
| Topic1 | 9 / **5** | 11 / **0** |
| Topic9 | 5 / **4** | 20 / **0** |
| Topic16 | 8 / **2** | 10 / **0** |
| Topic80 | 8 / **6** | 6 / **0** |

- `integrated_grad` finds **more *positive* (activating)** patterns (Topic9: 5→20 pos).
- `expected_integrated_grad` is the **only** method that surfaces **negative (repressive)** patterns — IG against a zero baseline yields zero negatives across all topics. Expected-gradients averages over data-sampled baselines, which is what exposes repressive signal.

The negative-pattern loss is **100% a method effect, not region count** (confirmed by the matched 1000-region IG run also having zero negatives).

## OPEN DECISION (deferred 2026-05-28 — revisit) — topic contribution method + v3

Region count is a minor knob (+21%). The real choice is the **contribution method**, and it is scientific, not just compute:
- **`integrated_grad`** — fast (~28 min/topic, CPU), more activating motifs, full v3 rebuild tractable (~hours). **No repressive motifs.**
- **`expected_integrated_grad`** — surfaces the repressive (negative) motifs that are in the *current* v1 bundle, but ~16 h+/topic → a full v3 is multi-day.

Not yet decided. When revisiting, also fold in the standing v3 region-count/window plan in the "If/when we do an n100 v3" section above (`top_k/max_regions=2000`, `-w 1000`). Pilot artifacts are preserved under `crested_model/{contributions_pilot_v2, modisco_pilot_v2, contributions_baseline_ig_1000, modisco_baseline_ig_1000}/` for whoever picks this up.
