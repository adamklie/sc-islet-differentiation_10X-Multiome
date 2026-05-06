#!/usr/bin/env python
"""Convert compute_contributions_specific.py output (approach B, one NPZ per
topic, modisco-lite-ready) into the directory layout expected by
run_modisco_all_topics.sh and the validated endo_rep characterization
notebooks (approach A).

Approach B layout (input):
  {input_dir}/Topic{N}_contrib.npz   ['arr_0'] -> (k, 4, L) float64
  {input_dir}/Topic{N}_oh.npz        ['arr_0'] -> (k, 4, L) float32

Approach A layout (output):
  {output_dir}/Topic{N}/one_hot_seqs.npz         ['arr_0']  ← Topic{N}_oh.npz
  {output_dir}/Topic{N}/hypothetical_contribs.npz['arr_0']  ← Topic{N}_contrib.npz
  {output_dir}/Topic{N}/projected_contribs.npz   ['arr_0']  ← contrib * one_hot

(region_names.npy is not needed by run_modisco_all_topics.sh; downstream
notebooks that need it can recover region indices from the AnnData.)
"""

import argparse
import glob
import os
import re
import sys

import numpy as np


def main(args):
    pat = re.compile(r"^(Topic\d+)_contrib\.npz$")
    contrib_files = sorted(
        f for f in glob.glob(os.path.join(args.input_dir, "Topic*_contrib.npz"))
        if pat.match(os.path.basename(f))
    )
    if not contrib_files:
        sys.exit(f"No Topic*_contrib.npz under {args.input_dir}")
    print(f"Converting {len(contrib_files)} topics: {args.input_dir} -> {args.output_dir}")
    for cf in contrib_files:
        topic = pat.match(os.path.basename(cf)).group(1)
        oh = os.path.join(args.input_dir, f"{topic}_oh.npz")
        if not os.path.exists(oh):
            print(f"  SKIP {topic}: missing {topic}_oh.npz"); continue
        out = os.path.join(args.output_dir, topic)
        os.makedirs(out, exist_ok=True)
        if (os.path.exists(os.path.join(out, "projected_contribs.npz"))
                and not args.force):
            print(f"  skip {topic} (already converted, use --force to overwrite)"); continue
        contrib = np.load(cf)["arr_0"]
        one_hot = np.load(oh)["arr_0"]
        if contrib.shape != one_hot.shape:
            print(f"  WARN {topic}: shape mismatch {contrib.shape} vs {one_hot.shape}")
        np.savez_compressed(os.path.join(out, "hypothetical_contribs.npz"), arr_0=contrib)
        np.savez_compressed(os.path.join(out, "one_hot_seqs.npz"), arr_0=one_hot)
        np.savez_compressed(os.path.join(out, "projected_contribs.npz"), arr_0=contrib * one_hot)
        print(f"  {topic}: ({contrib.shape[0]}, {contrib.shape[1]}, {contrib.shape[2]}) -> {out}/")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input_dir", required=True,
                   help="Dir with Topic{N}_{contrib,oh}.npz from compute_contributions_specific.py")
    p.add_argument("--output_dir", required=True,
                   help="Dir to write Topic{N}/{one_hot_seqs,hypothetical_contribs,projected_contribs}.npz")
    p.add_argument("--force", action="store_true",
                   help="Overwrite existing converted topics")
    main(p.parse_args())
