#!/usr/bin/env python
"""Aggregate per-cell-type marginalization NPZs -> motif x 22 matrix for a given
motif set (consensus | constituent).
score = mean_over_seqs(c_after - c_before).

Reads per-model NPZs from:
  results/3_single_task_models/<CT>/fold_0/chrombpnet/0.5/marginalization/<motif_set>/
Writes the summary matrix + zscore to:
  results/3_single_task_models/motifs/gimme_cluster/marginalization/<motif_set>/

Orthogonal check (non-fatal): rows shared by name with the existing gimme_average
pipeline matrix. The new pipeline names motifs cl<cf>_consensus / cl<cf>_<ct>_<n>_<posneg>
while the old gimme pipeline uses Average_<NNN>, so direct name overlap is empty unless a
cl<cf>->Average map is applied; the guard keeps this non-fatal.
"""
from __future__ import annotations
import argparse
import os
import numpy as np
import pandas as pd

RESULTS = "/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models"
MARG_SUMMARY = os.path.join(RESULTS, "motifs", "gimme_cluster", "marginalization")
EXIST = os.path.join(MARG_SUMMARY, "gimme_average", "marginalization_matrix.csv")
CT_ORDER = ["DE","PGT1","PGT2","PGT3","PFG1","PFG2","PP1","PP2","ENP_phase1","early_ENP","late_ENP",
            "early_SC_beta","late_SC_beta","early_SC_alpha","late_SC_alpha","early_SC_EC","late_SC_EC",
            "SC_delta_GHRL","proliferating_endocrine","exocrine","liver","FB_FLT1"]


def per_ct_dir(ct, motif_set):
    return os.path.join(RESULTS, ct, "fold_0", "chrombpnet", "0.5", "marginalization", motif_set)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--motif-set", required=True, choices=["consensus", "constituent"])
    args = ap.parse_args()
    motif_set = args.motif_set

    out_dir = os.path.join(MARG_SUMMARY, motif_set)
    os.makedirs(out_dir, exist_ok=True)

    cols = {}; names = None
    for ct in CT_ORDER:
        d_dir = per_ct_dir(ct, motif_set)
        p = os.path.join(d_dir, "marginalization_data.counts.npz")
        mp = os.path.join(d_dir, "marginalization_data.motifs.npz")
        if not (os.path.exists(p) and os.path.exists(mp)):
            print(f"  WARN missing {ct}"); continue
        d = np.load(p, allow_pickle=True); m = np.load(mp, allow_pickle=True)
        cur = m["motifs"].astype(str)
        if names is None:
            names = cur
        else:
            assert list(cur) == list(names), f"{ct}: order mismatch"
        cols[ct] = (d["c_after"] - d["c_before"]).mean(axis=1).reshape(len(cur))
    mat = pd.DataFrame(cols, index=names)[[c for c in CT_ORDER if c in cols]]
    mat.to_csv(os.path.join(out_dir, "marginalization_matrix.csv"), index_label="motif")
    mu = mat.mean(1).values[:, None]; sd = mat.std(1).values[:, None]; sd[sd == 0] = 1
    pd.DataFrame((mat.values-mu)/sd, index=mat.index, columns=mat.columns).to_csv(
        os.path.join(out_dir, "marginalization_matrix_zscore.csv"), index_label="motif")
    print(f"[{motif_set}] matrix {mat.shape}; {len([c for c in cols])} cell types")
    print(f"  wrote {os.path.join(out_dir, 'marginalization_matrix.csv')}")

    # Orthogonal check vs the existing gimme_average pipeline matrix (non-fatal).
    if not os.path.exists(EXIST):
        print(f"[ORTHOGONAL CHECK] skipped — {EXIST} not found.")
        return
    old = pd.read_csv(EXIST, index_col=0)
    common = [c for c in mat.columns if c in old.columns]
    shared = [n for n in mat.index if n in old.index]
    if not shared:
        print(f"[ORTHOGONAL CHECK] skipped — 0 rows share a name with {os.path.basename(EXIST)} "
              f"(new naming cl<cf>_* vs old Average_<NNN>; needs a cl->Average map to compare).")
        return
    rs = np.array([np.corrcoef(mat.loc[c, common].astype(float).values,
                               old.loc[c, common].astype(float).values)[0, 1] for c in shared])
    print(f"[ORTHOGONAL CHECK] {len(shared)} shared rows vs existing matrix:")
    print(f"  Pearson r: min={rs.min():.4f} mean={rs.mean():.4f} (n with r<0.999: {(rs<0.999).sum()})")


if __name__ == "__main__":
    main()
