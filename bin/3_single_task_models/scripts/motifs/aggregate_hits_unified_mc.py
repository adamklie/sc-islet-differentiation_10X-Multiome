#!/usr/bin/env python
"""Aggregate per-CT finemo hits (MC unified catalog) into motif x CT matrices.

Same schema as bin/4_multi_task_model/scripts/hits/aggregate_hits.py, but:
  - reads the single-task MC hits at
      models/<CT>/average/motifs/hits/counts_unified_mc/hits.tsv
  - strips the finemo group prefix from motif_name
      (pos_patterns.cluster_final#N / neg_patterns.cluster_final#N -> cluster_final#N)
    so motif columns match the catalog display_name (the .mc `name`).

Writes to results/3_single_task_models/motifs/motif_compendium_cluster/hits/unified/:
  hit_count_matrix.csv, hit_importance_matrix.csv, hit_importance_fraction.csv,
  hit_count_zscore.csv, hit_importance_zscore.csv, hit_quality_summary.csv
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd

ST = "/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models"
PREFIX_RE = re.compile(r"^(pos_patterns|neg_patterns)\.")


def main(args):
    pattern = f"{args.models_root}/*/average/motifs/hits/counts_unified_mc/hits.tsv"
    hit_files = sorted(glob.glob(pattern))
    print(f"Found {len(hit_files)} CT hits.tsv under {pattern}")
    if not hit_files:
        raise SystemExit("No MC hits.tsv found — has the finemo array finished?")

    count_rows, importance_rows, qc_rows = [], [], []
    for f in hit_files:
        # .../models/<CT>/average/motifs/hits/counts_unified_mc/hits.tsv
        # parents: [0]=counts_unified_mc [1]=hits [2]=motifs [3]=average [4]=<CT>
        celltype = Path(f).parents[4].name
        df = pd.read_csv(f, sep="\t")
        df["motif_name"] = df["motif_name"].str.replace(PREFIX_RE, "", regex=True)
        qc_rows.append({
            "celltype": celltype,
            "n_hits": len(df),
            "n_peaks": df["peak_id"].nunique() if "peak_id" in df.columns else np.nan,
            "median_similarity": df["hit_similarity"].median(),
            "mean_coefficient_global": df["hit_coefficient_global"].mean(),
            "median_importance": df["hit_importance"].median(),
        })
        grouped = df.groupby("motif_name")
        count_rows.append(grouped.size().rename(celltype))
        importance_rows.append(grouped["hit_coefficient_global"].sum().rename(celltype))
        print(f"  {celltype}: {len(df):,} hits, {df['peak_id'].nunique():,} peaks, "
              f"{df['motif_name'].nunique()} motifs")
        del df

    count_matrix = pd.DataFrame(count_rows).fillna(0).astype(int)
    importance_matrix = pd.DataFrame(importance_rows).fillna(0.0)
    qc = pd.DataFrame(qc_rows).set_index("celltype")

    if args.metadata and Path(args.metadata).exists():
        meta = pd.read_csv(args.metadata, sep="\t")
        if "cell_type" in meta.columns:
            ct_order = [c for c in meta["cell_type"].tolist() if c in count_matrix.index]
            count_matrix = count_matrix.reindex(ct_order)
            importance_matrix = importance_matrix.reindex(ct_order)
            qc = qc.reindex(ct_order)
    qc["hits_per_peak"] = qc["n_hits"] / qc["n_peaks"]

    importance_frac = importance_matrix.div(importance_matrix.sum(axis=1), axis=0)
    count_z = ((count_matrix - count_matrix.mean(axis=0)) / count_matrix.std(axis=0)).fillna(0)
    imp_z = ((importance_matrix - importance_matrix.mean(axis=0)) / importance_matrix.std(axis=0)).fillna(0)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    count_matrix.to_csv(out_dir / "hit_count_matrix.csv")
    importance_matrix.to_csv(out_dir / "hit_importance_matrix.csv")
    importance_frac.to_csv(out_dir / "hit_importance_fraction.csv")
    count_z.to_csv(out_dir / "hit_count_zscore.csv")
    imp_z.to_csv(out_dir / "hit_importance_zscore.csv")
    qc.to_csv(out_dir / "hit_quality_summary.csv")
    print(f"\nWrote 5 matrices + QC under {out_dir}")
    print(f"  count_matrix shape={count_matrix.shape}  total hits={count_matrix.sum().sum():,}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--models_root", default=f"{ST}/models")
    p.add_argument("--output_dir", default=f"{ST}/motifs/motif_compendium_cluster/hits/unified")
    p.add_argument("--metadata",
                   default="/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/config/cell_type_metadata.tsv")
    main(p.parse_args())
