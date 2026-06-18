"""Aggregate per-CT motif_genomic_context.tsv into long + wide tables for Phase 7.3.

Inputs:
  {ANN_DIR}/{CT}/motif_genomic_context.tsv  for each of 22 CTs.

Outputs (under {ANN_DIR}/):
  per_motif_per_ct.tsv     — long: motif_name, celltype, n_hits, region pcts,
                              median_tss_distance_abs/signed, median_summit_distance_abs
                              (used for box plots & compartment-aware aggregation in 7.3)

  per_motif_summary.tsv    — wide: one row per motif, aggregated:
                              total_hits, n_ct_with_hits, hit-weighted region pcts,
                              hit-weighted median TSS / summit distances
                              (the "wide" form that drives the Fig 3e roll-up)
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


CTS = [
    "DE", "ENP_phase1", "FB_FLT1", "PFG1", "PFG2", "PGT1", "PGT2", "PGT3", "PP1", "PP2",
    "SC_delta_GHRL", "early_ENP", "early_SC_EC", "early_SC_alpha", "early_SC_beta",
    "exocrine", "late_ENP", "late_SC_EC", "late_SC_alpha", "late_SC_beta", "liver",
    "proliferating_endocrine",
]
REGIONS = ["promoter", "5UTR", "3UTR", "exon", "intron", "distal"]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ann_dir", required=True)
    args = p.parse_args()

    ann = Path(args.ann_dir)
    long_rows = []
    for ct in CTS:
        f = ann / ct / "motif_genomic_context.tsv"
        if not f.exists():
            print(f"  WARN missing {f}")
            continue
        df = pd.read_csv(f, sep="\t")
        df["celltype"] = ct
        long_rows.append(df)
    long_df = pd.concat(long_rows, ignore_index=True)
    long_out = ann / "per_motif_per_ct.tsv"
    long_df.to_csv(long_out, sep="\t", index=False)
    print(f"wrote {long_out} ({len(long_df)} motif*CT rows)")

    # --- wide per-motif summary: hit-weighted aggregates across CTs ---
    # Region pct: weighted average of pct_<region> by n_hits per CT (i.e. equivalent
    # to fraction across the pooled hit set).
    grp = long_df.groupby("motif_name", as_index=False)

    rows = []
    for name, sub in grp:
        nh = sub["n_hits"].astype(float)
        total_hits = int(nh.sum())
        n_ct = int((nh > 0).sum())
        # weighted region pcts
        pcts = {}
        for r in REGIONS:
            col = f"pct_{r}"
            if col not in sub.columns or total_hits == 0:
                pcts[col] = float("nan")
                continue
            pcts[col] = float((sub[col] * nh).sum() / nh.sum()) if nh.sum() else float("nan")
        # weighted-median proxy: use hit-weighted mean of per-CT medians.
        # Reproducing the paper's box plot exactly requires the per-CT median
        # per motif (already in long_df) — keep this as a summary; long_df is
        # the source for the box plot.
        wm_tss_abs = float((sub["median_tss_distance_abs"] * nh).sum() / nh.sum()) if nh.sum() else float("nan")
        wm_tss_sig = float((sub["median_tss_distance_signed"] * nh).sum() / nh.sum()) if nh.sum() else float("nan")
        wm_sum_abs = float((sub["median_summit_distance_abs"] * nh).sum() / nh.sum()) if nh.sum() else float("nan")

        row = {"motif_name": name, "total_hits": total_hits, "n_ct_with_hits": n_ct,
               "wm_median_tss_distance_abs": wm_tss_abs,
               "wm_median_tss_distance_signed": wm_tss_sig,
               "wm_median_summit_distance_abs": wm_sum_abs}
        row.update(pcts)
        rows.append(row)

    wide_df = pd.DataFrame(rows).sort_values("total_hits", ascending=False)
    wide_out = ann / "per_motif_summary.tsv"
    wide_df.to_csv(wide_out, sep="\t", index=False)
    print(f"wrote {wide_out} ({len(wide_df)} motifs)")
    print("\nhead of summary:")
    print(wide_df.head(8).to_string(index=False))


if __name__ == "__main__":
    main()
