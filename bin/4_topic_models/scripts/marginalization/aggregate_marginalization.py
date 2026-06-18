#!/usr/bin/env python
"""Aggregate per-topic DeepTopic marginalization NPZs into motif × topic matrices.

Mirrors the schema produced for ChromBPNet and the multi-task model
(``results/4_multi_task_model/crested/motifs/analysis/marginalization_matrix{,_zscore}.csv``).

Inputs (one per topic, written by marginalize_topic.py):
    results/4_topic_models/all/n100/crested_model/marginalization/Topic{N}/marginalization_data.counts.npz

Outputs:
    <output_dir>/marginalization_matrix.csv          motif × 100 topics (mean delta, sigmoid prob)
    <output_dir>/marginalization_matrix_zscore.csv   z-scored per motif across topics

Columns are sorted **numerically by topic index** (Topic1 ... Topic100), not
lexicographically — this fixes the natural sort issue from import_beds.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


_TOPIC_NUM_RE = re.compile(r"Topic(\d+)$")


def _topic_sort_key(name: str) -> int:
    m = _TOPIC_NUM_RE.match(name)
    return int(m.group(1)) if m else 10**9


def main(args):
    in_paths = sorted(glob.glob(os.path.join(args.marginalization_root, "Topic*",
                                             "marginalization_data.counts.npz")))
    print(f"Found {len(in_paths)} topic marginalization NPZs")
    if not in_paths:
        raise SystemExit("No NPZs found.")

    # Each NPZ has the same motif catalog (in the same order). Resolve duplicate
    # MEME names (e.g. CTCF_MA0139.1 appears twice across distinct PWMs) by
    # suffixing the second+ occurrences with a positional index so the pivot
    # below has unique row labels.
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
        topic = Path(p).parent.name
        data = np.load(p, allow_pickle=True)
        c_before = data["c_before"]  # (n_motifs, n_seqs)
        c_after = data["c_after"]
        names = _dedupe(data["names"].astype(str))
        c_diff_mean = (c_after - c_before).mean(axis=1)
        if motif_names is None:
            motif_names = names
        dfs.append(pd.DataFrame({
            "name": names,
            "topic": topic,
            "c_diff_mean": c_diff_mean,
        }))
    df = pd.concat(dfs, ignore_index=True)
    print(f"  {df['topic'].nunique()} topics x {df['name'].nunique()} motifs")

    matrix = df.pivot(index="name", columns="topic", values="c_diff_mean")

    # Sort columns numerically: Topic1, Topic2, ..., Topic100
    topic_order = sorted(matrix.columns.tolist(), key=_topic_sort_key)
    matrix = matrix[topic_order]

    # Z-score across topics (per-motif)
    mu = matrix.mean(axis=1).values[:, np.newaxis]
    sd = matrix.std(axis=1).values[:, np.newaxis]
    sd[sd == 0] = 1.0
    matrix_z = pd.DataFrame((matrix.values - mu) / sd,
                            index=matrix.index, columns=matrix.columns)

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
                   help="Directory containing per-topic subdirs with marginalization_data.counts.npz")
    p.add_argument("--output_dir", required=True)
    main(p.parse_args())
