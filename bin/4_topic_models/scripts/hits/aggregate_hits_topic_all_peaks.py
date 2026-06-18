#!/usr/bin/env python
"""Aggregate per-topic finemo hits.tsv from the all-peaks topic run into
topic x motif matrices.

Topic-side analog of bin/4_multi_task_model/scripts/hits/aggregate_hits.py with
two differences:
  - Rows are topics (e.g. Topic1, Topic9, ...) instead of cell types.
  - Output rows are ordered by the canonical 26-confident-topic list so the
    matrix has a stable orientation: rows=26 topics, cols=motifs.

Inputs (one per topic, written by finemo call-hits):
    results/4_topic_models/all/n100/crested_model/motifs/hits_all_peaks/
        <Topic>/hits.tsv

Outputs (under args.output_dir):
    hit_count_matrix.csv         (topic rows x motif cols, int)
    hit_importance_matrix.csv    (topic rows x motif cols, sum hit_coefficient_global)
    hit_importance_fraction.csv  (rows sum to 1)
    hit_count_zscore.csv         (z-scored across topics per motif)
    hit_importance_zscore.csv
    hit_quality_summary.csv
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd


# Canonical confident-topic ordering: matches
# results/4_topic_models/all/n100/summary/topic_categories.tsv filtered to
# category in {one-to-one, lineage}.
CONFIDENT_TOPICS_DEFAULT = [
    "Topic1", "Topic9", "Topic13", "Topic15", "Topic16", "Topic17", "Topic20",
    "Topic29", "Topic31", "Topic32", "Topic33", "Topic36", "Topic38",
    "Topic39", "Topic44", "Topic45", "Topic47", "Topic52", "Topic67",
    "Topic76", "Topic79", "Topic80", "Topic86", "Topic94", "Topic96",
    "Topic99",
]


def main(args):
    hit_files = sorted(glob.glob(os.path.join(args.hits_root, "*", "hits.tsv")))
    print(f"Found {len(hit_files)} topic hits.tsv files under {args.hits_root}")
    if not hit_files:
        raise SystemExit(f"No hits.tsv under {args.hits_root}/*/")

    count_rows = []
    importance_rows = []
    qc_rows = []
    for f in hit_files:
        topic = Path(f).parent.name
        df = pd.read_csv(f, sep="\t")
        qc_rows.append({
            "topic": topic,
            "n_hits": len(df),
            "n_peaks": df["peak_id"].nunique() if "peak_id" in df.columns else np.nan,
            "median_similarity": df["hit_similarity"].median(),
            "mean_coefficient": df["hit_coefficient"].mean(),
            "mean_coefficient_global": df["hit_coefficient_global"].mean(),
            "median_importance": df["hit_importance"].median(),
        })
        grouped = df.groupby("motif_name")
        count_rows.append(grouped.size().rename(topic))
        importance_rows.append(grouped["hit_coefficient_global"].sum().rename(topic))
        print(f"  {topic}: {len(df):,} hits, "
              f"{df['peak_id'].nunique() if 'peak_id' in df.columns else 'NA'} peaks")
        del df

    count_matrix = pd.DataFrame(count_rows).fillna(0).astype(int)
    importance_matrix = pd.DataFrame(importance_rows).fillna(0.0)
    qc = pd.DataFrame(qc_rows).set_index("topic")

    # Reorder rows to canonical confident-topic order (drop topics absent from
    # the run; warn if any are missing).
    confident_topics = (args.topic_order.split(",")
                        if args.topic_order else CONFIDENT_TOPICS_DEFAULT)
    present = [t for t in confident_topics if t in count_matrix.index]
    missing = [t for t in confident_topics if t not in count_matrix.index]
    extra = [t for t in count_matrix.index if t not in confident_topics]
    if missing:
        print(f"WARNING: confident topics missing from hits_all_peaks: {missing}")
    if extra:
        print(f"NOTE: extra topics present in hits_all_peaks (not in confident list): {extra}")
    # Final row order: canonical confident topics first (those present), then extras
    row_order = present + sorted(extra)
    count_matrix = count_matrix.reindex(row_order).fillna(0).astype(int)
    importance_matrix = importance_matrix.reindex(row_order).fillna(0.0)
    qc = qc.reindex(row_order)
    qc["hits_per_peak"] = qc["n_hits"] / qc["n_peaks"]

    # Derived
    importance_frac = importance_matrix.div(importance_matrix.sum(axis=1), axis=0)
    count_z = (count_matrix - count_matrix.mean(axis=0)) / count_matrix.std(axis=0)
    count_z = count_z.fillna(0)
    imp_z = (importance_matrix - importance_matrix.mean(axis=0)) / importance_matrix.std(axis=0)
    imp_z = imp_z.fillna(0)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    count_matrix.to_csv(out_dir / "hit_count_matrix.csv")
    importance_matrix.to_csv(out_dir / "hit_importance_matrix.csv")
    importance_frac.to_csv(out_dir / "hit_importance_fraction.csv")
    count_z.to_csv(out_dir / "hit_count_zscore.csv")
    imp_z.to_csv(out_dir / "hit_importance_zscore.csv")
    qc.to_csv(out_dir / "hit_quality_summary.csv")
    print(f"\nWrote 5 matrices + QC under {out_dir}")
    print(f"  count_matrix       shape={count_matrix.shape}  (rows=topics, cols=motifs)")
    print(f"  importance_matrix  shape={importance_matrix.shape}")
    print(f"  Total motif instances: {count_matrix.sum().sum():,}")
    print(f"  Top motifs by total importance:")
    for m, v in importance_matrix.sum(axis=0).sort_values(ascending=False).head(15).items():
        print(f"    {m:35s}  {v:12.4f}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--hits_root", required=True,
                   help="Dir containing per-topic subdirs with hits.tsv")
    p.add_argument("--output_dir", required=True)
    p.add_argument("--topic_order", default=None,
                   help="Optional comma-separated topic order. Defaults to the "
                        "26 confident topics (one-to-one + lineage).")
    main(p.parse_args())
