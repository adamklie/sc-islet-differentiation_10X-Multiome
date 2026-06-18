#!/usr/bin/env python
"""Build a synthetic tfmodisco-lite h5 + name_map from a MEME catalog.

Mirrors what the ChromBPNet pipeline does for finemo's unified-catalog hit
calling (see ``results/3_single_task_models/motifs/clustered_motifs.modisco.h5``
and its sibling ``clustered_motifs_name_map.tsv``).

Each MEME motif becomes one ``pos_patterns/pattern_<i>`` group with
``contrib_scores`` / ``hypothetical_contribs`` set to ``PPM * IC_per_position``
(IC = ``sum(PPM * max(0, log2(PPM/0.25)))`` per row) and ``sequence`` = PPM.
All patterns are padded to a common width so finemo can stack them.

The name_map TSV maps ``pos_patterns.pattern_<i>`` -> the original MEME name
(``Average_XXX`` etc.) so finemo's hit table reports human-readable IDs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def parse_meme(meme_path: str | Path):
    motifs = []
    with open(meme_path) as fh:
        cur_name = None
        cur_rows = []
        in_matrix = False
        for raw in fh:
            line = raw.rstrip()
            if line.startswith("MOTIF"):
                if cur_name is not None and cur_rows:
                    motifs.append((cur_name, np.asarray(cur_rows, dtype=np.float64)))
                cur_name = line.split()[1]
                cur_rows = []
                in_matrix = False
            elif line.startswith("letter-probability matrix"):
                in_matrix = True
            elif in_matrix:
                if not line.strip():
                    in_matrix = False
                    continue
                try:
                    parts = [float(x) for x in line.split()]
                except ValueError:
                    in_matrix = False
                    continue
                if len(parts) >= 4:
                    cur_rows.append(parts[:4])
                else:
                    in_matrix = False
        if cur_name is not None and cur_rows:
            motifs.append((cur_name, np.asarray(cur_rows, dtype=np.float64)))
    return motifs


def ic_per_pos(ppm: np.ndarray) -> np.ndarray:
    eps = 1e-9
    log_ratio = np.maximum(np.log2((ppm + eps) / 0.25), 0.0)
    return (ppm * log_ratio).sum(axis=1)


def pad_centered(matrix: np.ndarray, target_len: int) -> np.ndarray:
    L = matrix.shape[0]
    if L >= target_len:
        return matrix[:target_len]
    pad_left = (target_len - L) // 2
    pad_right = target_len - L - pad_left
    out = np.zeros((target_len, matrix.shape[1]), dtype=matrix.dtype)
    out[pad_left:pad_left + L] = matrix
    return out


def main(args):
    motifs = parse_meme(args.meme)
    if not motifs:
        raise SystemExit(f"No motifs parsed from {args.meme}")
    target_len = max(m[1].shape[0] for m in motifs)
    print(f"Parsed {len(motifs)} motifs from {args.meme}; padding to length {target_len}")

    name_map = {}
    out_h5 = Path(args.out_h5)
    out_h5.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_h5, "w") as f:
        grp = f.create_group("pos_patterns")
        for i, (motif_name, ppm) in enumerate(motifs):
            ppm = ppm + 1e-9
            ppm = ppm / ppm.sum(axis=1, keepdims=True)
            ic = ic_per_pos(ppm)
            cs = ppm * ic[:, None]
            seq = ppm.copy()
            cs_padded = pad_centered(cs, target_len)
            seq_padded = pad_centered(seq, target_len)
            pname = f"pattern_{i}"
            sub = grp.create_group(pname)
            sub.create_dataset("contrib_scores", data=cs_padded)
            sub.create_dataset("hypothetical_contribs", data=cs_padded)
            sub.create_dataset("sequence", data=seq_padded)
            name_map[f"pos_patterns.{pname}"] = motif_name
    print(f"Wrote {out_h5}")

    name_map_tsv = Path(args.out_name_map)
    name_map_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(name_map_tsv, "w") as fh:
        for k, v in name_map.items():
            fh.write(f"{k}\t{v}\n")
    print(f"Wrote {name_map_tsv}")

    # Optional: also write JSON
    json_path = name_map_tsv.with_suffix(".json")
    with open(json_path, "w") as fh:
        json.dump(name_map, fh, indent=2)
    print(f"Wrote {json_path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--meme", required=True, help="MEME file with the unified motif catalog")
    p.add_argument("--out_h5", required=True)
    p.add_argument("--out_name_map", required=True, help="TSV mapping pos_patterns.pattern_<i> -> motif name")
    main(p.parse_args())
