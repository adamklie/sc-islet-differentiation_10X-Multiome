"""Extract PFMs from a tfmodisco-lite h5, keeping pos and neg patterns separate.

Mirrors chrombpnet's ``modisco_to_pfm.py`` (same trim threshold, min_seqlets,
flank semantics) with two changes:

1. Processes BOTH ``pos_patterns`` AND ``neg_patterns`` (the original script
   silently dropped neg). Output filenames carry a sign tag so the downstream
   clustering pipeline can keep the two tracks separate.
2. Fixes the original's off-by-one: ``num_patterns = len(group) - 1`` dropped
   the last pattern in every CT. We iterate all patterns here. The bug is
   benign (one motif dropped per CT) but means this catalog is slightly larger
   than v3 even before counting neg patterns.

Output layout:
  <output_dir>/pos/<prefix>.pos_pattern_<i>.pfm
  <output_dir>/neg/<prefix>.neg_pattern_<i>.pfm

The per-CT concat (one big .pfm per CT-per-sign) is done by the shell wrapper.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import h5py
import numpy as np


def trim_ppm(ppm, t=0.45, min_length=3, flank=0):
    """Trim PPM to first/last position with max>=t; preserve original behavior."""
    maxes = np.max(ppm, -1)
    maxes = np.where(maxes >= t)
    if (len(maxes[0]) == 0) or (maxes[0][-1] + 1 - maxes[0][0] < min_length):
        return None
    return ppm[max(maxes[0][0] - flank, 0):maxes[0][-1] + 1 + flank]


def extract_side(f, group_key: str, sign: str, out_dir: Path, prefix: str,
                 threshold: float, min_length: int, flank: int,
                 min_seqlets: int, normalize: bool) -> int:
    """Extract PFMs for one side (pos or neg). Returns count written."""
    if group_key not in f:
        return 0
    group = f[group_key]
    # iterate ALL patterns (no off-by-one); sort by index to match historical order
    pattern_keys = sorted(group.keys(), key=lambda x: int(x.split("_")[-1]))
    side_dir = out_dir / sign
    side_dir.mkdir(parents=True, exist_ok=True)
    n_written = 0
    for pname in pattern_keys:
        idx = int(pname.split("_")[-1])
        ppm = np.asarray(group[pname]["sequence"])
        trimmed = trim_ppm(ppm, t=threshold, min_length=min_length, flank=flank)
        if trimmed is None:
            continue
        n_seqlets = len(group[pname]["seqlets"]["sequence"])
        if n_seqlets < min_seqlets:
            continue
        pfm = (trimmed * n_seqlets).astype(int)
        if normalize:
            pfm = (pfm + 1) / np.sum(pfm, axis=1, keepdims=True)
            fmt = "%.2f"
        else:
            fmt = "%d"
        out_path = side_dir / f"{prefix}.{sign}_pattern_{idx}.pfm"
        np.savetxt(out_path, pfm, fmt=fmt, delimiter="\t")
        n_written += 1
    return n_written


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-m", "--modisco_h5py", required=True, type=str)
    ap.add_argument("-o", "--output_dir", required=True, type=str)
    ap.add_argument("-t", "--threshold", default=0.45, type=float)
    ap.add_argument("-ml", "--min_length", default=3, type=int)
    ap.add_argument("-f", "--flank_add", default=0, type=int)
    ap.add_argument("-s", "--min_seqlets", default=100, type=int)
    ap.add_argument("-n", "--normalize", default=False, action="store_true")
    ap.add_argument("-op", "--output_prefix", default="", type=str,
                    help="Prefix for output .pfm files (typically the cell type)")
    args = ap.parse_args()

    out_dir = Path(args.output_dir)
    with h5py.File(args.modisco_h5py, "r") as f:
        n_pos = extract_side(f, "pos_patterns", "pos", out_dir, args.output_prefix,
                              args.threshold, args.min_length, args.flank_add,
                              args.min_seqlets, args.normalize)
        n_neg = extract_side(f, "neg_patterns", "neg", out_dir, args.output_prefix,
                              args.threshold, args.min_length, args.flank_add,
                              args.min_seqlets, args.normalize)
    print(f"  {args.output_prefix}: wrote {n_pos} pos, {n_neg} neg PFMs to {out_dir}/{{pos,neg}}/")


if __name__ == "__main__":
    main()
