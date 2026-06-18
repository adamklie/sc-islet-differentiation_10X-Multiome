#!/usr/bin/env python
"""Aggregate per-CT finemo hits into motif × CT matrices.

Mirrors the schema produced by ChromBPNet's
``bin/3_single_task_models/5_motif_hit_analysis.ipynb`` (saved at
``results/3_single_task_models/motifs/hit_analysis_unified/``):

  - hit_count_matrix.csv         (CT rows × motif cols, int)
  - hit_importance_matrix.csv    (CT rows × motif cols, sum of hit_coefficient_global)
  - hit_importance_fraction.csv  (rows sum to 1)
  - hit_count_zscore.csv         (z-scored across CTs per motif)
  - hit_importance_zscore.csv

Inputs (one per CT, written by finemo call-hits):
    results/4_multi_task_model/crested/motifs/hits/<CT>/hits.tsv

Cell-type ordering uses cell_type_metadata.tsv if available.
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd


def main(args):
    hit_files = sorted(glob.glob(os.path.join(args.hits_root, "*", "hits.tsv")))
    print(f"Found {len(hit_files)} CT hits.tsv files")
    if not hit_files:
        raise SystemExit(f"No hits.tsv under {args.hits_root}/*/")

    count_rows = []
    importance_rows = []
    qc_rows = []
    for f in hit_files:
        celltype = Path(f).parent.name
        df = pd.read_csv(f, sep="\t")
        qc_rows.append({
            "celltype": celltype,
            "n_hits": len(df),
            "n_peaks": df["peak_id"].nunique() if "peak_id" in df.columns else np.nan,
            "median_similarity": df["hit_similarity"].median(),
            "mean_coefficient": df["hit_coefficient"].mean(),
            "mean_coefficient_global": df["hit_coefficient_global"].mean(),
            "median_importance": df["hit_importance"].median(),
        })
        grouped = df.groupby("motif_name")
        count_rows.append(grouped.size().rename(celltype))
        importance_rows.append(grouped["hit_coefficient_global"].sum().rename(celltype))
        print(f"  {celltype}: {len(df):,} hits, {df['peak_id'].nunique():,} peaks")
        del df

    count_matrix = pd.DataFrame(count_rows).fillna(0).astype(int)
    importance_matrix = pd.DataFrame(importance_rows).fillna(0.0)
    qc = pd.DataFrame(qc_rows).set_index("celltype")

    # Optional CT reorder
    if args.metadata is not None and Path(args.metadata).exists():
        meta = pd.read_csv(args.metadata, sep="\t")
        if "cell_type" in meta.columns:
            ct_order = [c for c in meta["cell_type"].tolist() if c in count_matrix.index]
            count_matrix = count_matrix.reindex(ct_order)
            importance_matrix = importance_matrix.reindex(ct_order)
            qc = qc.reindex(ct_order)
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
    print(f"  count_matrix       shape={count_matrix.shape}")
    print(f"  importance_matrix  shape={importance_matrix.shape}")
    print(f"  Total motif instances: {count_matrix.sum().sum():,}")
    print(f"  Top motifs by total importance:")
    for m, v in importance_matrix.sum(axis=0).sort_values(ascending=False).head(15).items():
        print(f"    {m:35s}  {v:12.4f}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--hits_root", required=True,
                   help="Dir containing per-CT subdirs with hits.tsv")
    p.add_argument("--output_dir", required=True)
    p.add_argument("--metadata", default=None,
                   help="Optional cell_type_metadata.tsv for row ordering")
    main(p.parse_args())
