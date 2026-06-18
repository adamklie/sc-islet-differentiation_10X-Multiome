#!/usr/bin/env python
"""Replay the per-topic region selection used by
`compute_contribution_scores.py` (deeptopic) so we can run finemo on the
existing per-topic contributions NPZs.

The deeptopic script selects, per topic:
    positives = (adata.X[topic_idx] > 0)           # binarized labels
    if |positives| > max_regions:
        rank positives by model predictions for this topic,
        take top-K, then sort by index (genomic order).

We replay that logic by reading `adata_with_predictions.h5ad` (which already
has both the binarized matrix in X and stored model predictions in
layers/predictions), and emit a per-topic narrowPeak file in the SAME row
order as the contributions NPZs (1000 rows each).

The narrowPeak we write has 10 columns so finemo's `call-hits -p` accepts it
(only chr/start/end matter downstream, the rest are placeholders).

Sanity check: aborts if the resulting region count doesn't match the NPZ.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import h5py


def main(args):
    adata_path = args.adata
    contribs_root = Path(args.contribs_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(adata_path, "r") as f:
        obs_names = [b.decode() for b in f["obs/_index"][:]]
        chr_codes = f["var/chr/codes"][:]
        chr_cats = [b.decode() for b in f["var/chr/categories"][:]]
        chrs = np.array([chr_cats[c] for c in chr_codes])
        starts = f["var/start"][:]
        ends = f["var/end"][:]
        n_topics, n_vars = f["X"].shape
        print(f"adata: {n_topics} topics x {n_vars} regions")

        topics = args.topics or obs_names
        for topic in topics:
            if topic not in obs_names:
                print(f"  SKIP {topic}: not in obs_names")
                continue
            ti = obs_names.index(topic)

            x_row = f["X"][ti, :]
            pos_idx = np.where(x_row > 0)[0]
            n_pos = len(pos_idx)

            if n_pos > args.max_regions:
                preds = f["layers/predictions"][ti, :][pos_idx]
                top_local = np.argsort(preds)[::-1][: args.max_regions]
                top_global = pos_idx[top_local]
                top_global.sort()  # genomic order
                kept = top_global
            else:
                kept = pos_idx

            # Sanity: count must match contributions NPZ
            npz_path = contribs_root / topic / "one_hot_seqs.npz"
            if npz_path.exists():
                with np.load(npz_path) as z:
                    npz_n = z[z.files[0]].shape[0]
                if npz_n != len(kept):
                    raise SystemExit(
                        f"{topic}: NPZ has {npz_n} regions but replay gave {len(kept)} "
                        f"(positives={n_pos}, max={args.max_regions}). Selection drifted; "
                        f"refusing to write misaligned BED."
                    )
            else:
                print(f"  WARN {topic}: no NPZ found at {npz_path}, will still write BED")

            out_bed = out_dir / f"{topic}.narrowPeak"
            with open(out_bed, "w") as o:
                for i in kept:
                    c = chrs[i]
                    s = int(starts[i])
                    e = int(ends[i])
                    summit = (e - s) // 2
                    # chr start end name score strand signal pvalue qvalue summit
                    o.write(f"{c}\t{s}\t{e}\t{c}:{s}-{e}\t1000\t.\t0\t0\t0\t{summit}\n")
            print(f"  {topic}: wrote {len(kept)} regions -> {out_bed}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--adata", required=True,
                   help="adata_with_predictions.h5ad")
    p.add_argument("--contribs_dir", required=True,
                   help="Directory containing per-topic contribution NPZ subdirs (for sanity check)")
    p.add_argument("--output_dir", required=True)
    p.add_argument("--max_regions", type=int, default=1000)
    p.add_argument("--topics", nargs="*", default=None,
                   help="Specific topics to process (default: all)")
    main(p.parse_args())
