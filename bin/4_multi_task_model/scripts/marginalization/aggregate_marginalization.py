#!/usr/bin/env python
"""Aggregate per-CT multi-task marginalization NPZs into motif × CT matrices.

Mirrors the schema produced for ChromBPNet in
``results/3_single_task_models/motifs/analysis/marginalization_matrix{,_zscore}.csv``.

Inputs (one per CT, written by marginalize_multitask.py):
    results/4_multi_task_model/crested/marginalization/<CT>/marginalization_data.counts.npz

Outputs:
    results/4_multi_task_model/crested/motifs/analysis/marginalization_matrix.csv
    results/4_multi_task_model/crested/motifs/analysis/marginalization_matrix_zscore.csv
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd


def main(args):
    in_paths = sorted(glob.glob(os.path.join(args.marginalization_root, "*",
                                             "marginalization_data.counts.npz")))
    print(f"Found {len(in_paths)} CT marginalization NPZs")
    if not in_paths:
        raise SystemExit("No NPZs found.")

    # Synthetic MoDISco h5 catalogs may carry duplicate MEME names (e.g. when
    # gimme cluster lands two distinct PWMs on the same TomTom-mapped TF).
    # Suffix the second+ occurrences with a positional index so the pivot below
    # has unique row labels — same fix as aggregate_marginalization.py for topics.
    def _dedupe(names):
        seen = {}
        out = []
        for i, n in enumerate(names):
            if n in seen:
                seen[n] += 1
                out.append(f"{n}#{i}")
            else:
                seen[n] = 1
                out.append(n)
        return np.array(out)

    dfs = []
    motif_names = None
    for p in in_paths:
        ct = Path(p).parent.name
        data = np.load(p, allow_pickle=True)
        c_before = data["c_before"]  # (n_motifs, n_seqs)
        c_after = data["c_after"]
        names = _dedupe(data["names"].astype(str))
        c_diff_mean = (c_after - c_before).mean(axis=1)
        if motif_names is None:
            motif_names = names
        dfs.append(pd.DataFrame({
            "name": names,
            "celltype": ct,
            "c_diff_mean": c_diff_mean,
        }))
    df = pd.concat(dfs, ignore_index=True)
    print(f"  {df['celltype'].nunique()} CTs x {df['name'].nunique()} motifs")

    # Optional cell-type ordering using metadata
    ct_order = None
    if args.metadata is not None and Path(args.metadata).exists():
        meta = pd.read_csv(args.metadata, sep="\t")
        if "cell_type" in meta.columns:
            ct_order = [c for c in meta["cell_type"].tolist() if c in df["celltype"].unique()]

    matrix = df.pivot(index="name", columns="celltype", values="c_diff_mean")
    if ct_order is not None:
        matrix = matrix[ct_order]

    # Z-score across cell types (per-motif)
    mu = matrix.mean(axis=1).values[:, np.newaxis]
    sd = matrix.std(axis=1).values[:, np.newaxis]
    sd[sd == 0] = 1.0
    matrix_z = pd.DataFrame((matrix.values - mu) / sd, index=matrix.index, columns=matrix.columns)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_marg = out_dir / "marginalization_matrix.csv"
    out_z = out_dir / "marginalization_matrix_zscore.csv"
    matrix.to_csv(out_marg, index_label="name")
    matrix_z.to_csv(out_z, index_label="name")
    print(f"Wrote {out_marg}  shape={matrix.shape}")
    print(f"Wrote {out_z}    shape={matrix_z.shape}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--marginalization_root", required=True,
                   help="Directory containing per-CT subdirs with marginalization_data.counts.npz")
    p.add_argument("--output_dir", required=True)
    p.add_argument("--metadata", default=None,
                   help="Optional cell_type_metadata.tsv for column ordering")
    main(p.parse_args())
