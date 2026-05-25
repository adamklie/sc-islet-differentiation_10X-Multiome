#!/usr/bin/env python
"""Generate v2 topic annotations + topic_categories.tsv using the skip-structural framework.

This is the framework that was applied to endo_rep/n50 (see issue #2 comment 2026-04-10).
Tiers: exclusive / strong / moderate / broad / promoter_structural / batch_artifact
Confidence is based on zscore_gap (cell-type specificity) and eta2_batch (batch confounding).

The output topic_categories.tsv uses the v2 tiers but writes the column layout that the
characterization notebooks (3-7) expect: topic, category, manual_override, cell_types_list,
n_cell_types, batch_flag.  Notebooks key off `category` for column coloring; we map
v2 tiers -> notebook category palette (one-to-one, lineage, shared, structural,
batch-driven, unannotated) so the existing plots stay readable while preserving the
v2 informational columns in topic_annotations_v2.tsv.

Inputs (all read from <base>/summary/topic_master_summary.tsv built by
generate_topic_report.py, plus a few raw tables):
  - summary/topic_master_summary.tsv  (one row per topic, all annotations joined)
  - annotated/topic_qc/topic_batch_correlation.tsv  (raw eta2)
  - annotated/synthesis/topic_annotations.tsv       (top_motifs string per topic)
  - annotated/metadata/cell_type_top.tsv            (per-rank z-scores)

Outputs (in <base>/summary/):
  - topic_annotations_v2.tsv  -- full v2 schema with tier/confidence/structural_fraction/etc.
  - topic_categories.tsv      -- overwritten with v2-derived categories (notebook-compatible)
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd


STRUCTURAL_MOTIFS = {
    "CTCF", "BORIS", "CTCFL", "NFY", "NFYA", "NFYB", "NFYC",
    "Sp1", "SP1", "Sp2", "SP2", "SP3", "SP5", "Sp5",
    "NRF1", "NRF", "YY1", "TYY1", "ZN143", "ZNF143",
    "KLF1", "KLF3", "KLF4", "KLF5", "KLF10", "KLF15",
    "MAZ", "Maz",
    # SeqBias-style promoter "motifs"
    "SeqBias:", "Promoters", "promoter",
}


def _is_structural(motif_str):
    """Return True if the motif's TF family is in the structural set."""
    if pd.isna(motif_str):
        return False
    s = str(motif_str).strip()
    base = s.split("(")[0].split("/")[0].strip()
    # Try full token + leading word
    if base in STRUCTURAL_MOTIFS:
        return True
    head = base.split()[0]
    if head in STRUCTURAL_MOTIFS:
        return True
    upper = base.upper()
    return any(sm.upper() in upper for sm in STRUCTURAL_MOTIFS)


def first_nonstructural(motif_str):
    """Given a 'Cux2; HNF6; CTCF; ...' string, return the first non-structural TF."""
    if pd.isna(motif_str) or not str(motif_str).strip():
        return None
    parts = [m.strip() for m in str(motif_str).split(";") if m.strip()]
    for m in parts:
        if not _is_structural(m):
            return m
    return parts[0] if parts else None


def structural_fraction(motif_str):
    """Fraction of top motifs that are structural."""
    if pd.isna(motif_str) or not str(motif_str).strip():
        return 0.0
    parts = [m.strip() for m in str(motif_str).split(";") if m.strip()]
    if not parts:
        return 0.0
    return sum(1 for m in parts if _is_structural(m)) / len(parts)


def categorise_v2(row):
    """Return (tier, primary_cell_type, informative_motif, structural_fraction, description).

    Rule order (first match wins):
      1. batch_artifact     : eta2_batch > eta2_celltype  OR  batch_flag=='BATCH-DRIVEN'
                             AND eta2_batch > 0.10
      2. promoter_structural: structural_fraction >= 0.6  (CTCF/Sp1/NFY dominate top motifs)
      3. exclusive          : n_cell_types == 1
      4. strong             : zscore_gap >= 1.0 AND n_cell_types <= 3
      5. moderate           : 0.5 <= zscore_gap < 1.0 AND n_cell_types <= 4
      6. broad              : otherwise
    """
    motifs = row.get("top_motifs_syn", "")
    sf = structural_fraction(motifs)
    informative = first_nonstructural(motifs)

    eta2_b = row.get("eta2_batch", 0) or 0
    eta2_c = row.get("eta2_celltype", 0) or 0
    batch_flag = row.get("batch_flag", "OK")
    gap = row.get("zscore_gap", 0) or 0
    n_ct = int(row.get("n_cell_types", 0) or 0)
    primary = row.get("top_cell_type", "")

    # Rule 1: batch_artifact
    if (eta2_b > eta2_c or batch_flag == "BATCH-DRIVEN") and eta2_b > 0.10:
        desc = (f"Likely batch artifact (eta2_batch={eta2_b:.2f} vs eta2_celltype={eta2_c:.2f}). "
                f"Assigned to {primary} but confounded. "
                + (f"Informative motif: {informative}." if informative else ""))
        return "batch_artifact", primary, informative, sf, desc

    # Rule 2: promoter_structural
    if sf >= 0.6:
        desc = (f"Structural/promoter program. CTCF/insulator-dominated "
                f"({int(sf*100)}% structural motifs in top 5). "
                f"Weakly assigned to {primary}. n_regions captured via master table.")
        return "promoter_structural", primary, informative, sf, desc

    # Rule 3: exclusive
    if n_ct == 1:
        desc = (f"Exclusive to {primary}. "
                + (f"Driven by {informative} (gap={gap:.2f})." if informative
                   else f"gap={gap:.2f}."))
        return "exclusive", primary, informative, sf, desc

    # Rule 4: strong
    if gap >= 1.0 and n_ct <= 3:
        cts_str = row.get("cell_types_list", "")
        other = ", ".join(c.strip() for c in cts_str.split(",")[1:3])
        desc = (f"Strongly enriched in {primary} (gap={gap:.2f})"
                + (f", also in {other}" if other else "") + ". "
                + (f"Key TF: {informative}." if informative else ""))
        return "strong", primary, informative, sf, desc

    # Rule 5: moderate
    if gap >= 0.5 and n_ct <= 4:
        cts_str = row.get("cell_types_list", "")
        desc = (f"Moderately enriched in {primary} (gap={gap:.2f}). "
                f"Shared across {cts_str}. "
                + (f"Key TF: {informative}." if informative else ""))
        return "moderate", primary, informative, sf, desc

    # Rule 6: broad
    cts_str = row.get("cell_types_list", "")
    desc = (f"Broadly active across {n_ct} cell types ({cts_str[:80]}...). "
            f"Weak cell type specificity (gap={gap:.2f}). "
            + (f"Top TF: {informative}." if informative else ""))
    return "broad", primary, informative, sf, desc


def confidence_score(tier, gap, eta2_batch):
    """Three-level confidence: high / medium / low."""
    if tier == "batch_artifact":
        return "low"
    if tier == "promoter_structural":
        return "low"
    if tier in ("exclusive", "strong"):
        if gap >= 1.5 and eta2_batch < 0.10:
            return "high"
        if gap >= 1.0 and eta2_batch < 0.20:
            return "high"
        return "medium"
    if tier == "moderate":
        if gap >= 0.7 and eta2_batch < 0.10:
            return "medium"
        return "low"
    # broad
    return "low"


# Map v2 tiers -> v1 categories that the notebook col_colors palette knows about.
# This keeps the existing motif characterization notebook visuals readable while we
# still ship the full v2 schema in topic_annotations_v2.tsv.
TIER_TO_NB_CATEGORY = {
    "exclusive": "one-to-one",
    "strong": "lineage",
    "moderate": "lineage",
    "broad": "shared",
    "promoter_structural": "structural",
    "batch_artifact": "batch-driven",
}


def main():
    ap = argparse.ArgumentParser(description="Generate v2 topic annotations + categories")
    ap.add_argument("--base", required=True, help="Topic model subset results dir")
    args = ap.parse_args()

    base = args.base
    summary_dir = os.path.join(base, "summary")
    master_path = os.path.join(summary_dir, "topic_master_summary.tsv")
    if not os.path.exists(master_path):
        sys.exit(f"ERROR: topic_master_summary.tsv not found at {master_path}. "
                 "Run generate_topic_report.py first.")

    print(f"Loading master summary: {master_path}")
    master = pd.read_csv(master_path, sep="\t")
    print(f"  {len(master)} topics")

    # Apply v2 categorization
    print("Applying v2 categorisation framework (skip-structural)...")
    out_rows = []
    for _, row in master.iterrows():
        tier, primary, informative, sf, desc = categorise_v2(row)
        motifs = row.get("top_motifs_syn", "")
        top3 = "; ".join([m.strip() for m in str(motifs).split(";")[:3]])
        conf = confidence_score(tier, row.get("zscore_gap", 0) or 0,
                                row.get("eta2_batch", 0) or 0)
        out_rows.append({
            "topic": row["topic"],
            "tier": tier,
            "confidence": conf,
            "primary_cell_type": primary,
            "cell_types": row.get("cell_types_list", ""),
            "n_cell_types": row.get("n_cell_types", 0),
            "zscore_gap": round(float(row.get("zscore_gap", 0) or 0), 3),
            "informative_motif": informative if informative else "",
            "informative_motifs_top3": top3,
            "structural_fraction": round(sf, 2),
            "Gini": round(float(row.get("Gini", 0) or 0), 3),
            "n_regions": int(row.get("n_regions", 0) or 0),
            "eta2_batch": round(float(row.get("eta2_batch", 0) or 0), 4),
            "eta2_celltype": round(float(row.get("eta2_celltype", 0) or 0), 4),
            "batch_flag": row.get("batch_flag", "OK"),
            "description": desc,
        })

    v2_df = pd.DataFrame(out_rows)
    # Order by topic number
    v2_df["_tnum"] = v2_df["topic"].str.replace("Topic", "").astype(int)
    v2_df = v2_df.sort_values("_tnum").drop(columns="_tnum").reset_index(drop=True)

    # Save v2 annotations
    v2_path = os.path.join(summary_dir, "topic_annotations_v2.tsv")
    v2_df.to_csv(v2_path, sep="\t", index=False)
    print(f"  Saved: {v2_path}")

    # Tier distribution
    print("\n  Tier distribution:")
    for tier in ["exclusive", "strong", "moderate", "broad",
                 "promoter_structural", "batch_artifact"]:
        n = (v2_df["tier"] == tier).sum()
        if n > 0:
            print(f"    {tier:25s} {n:>3d}")

    confident = (v2_df["confidence"].isin(["high", "medium"])
                 & v2_df["tier"].isin(["exclusive", "strong", "moderate"])).sum()
    print(f"\n  Confident cell-type-specific topics: {confident}/{len(v2_df)}")

    # Re-derive notebook-compatible topic_categories.tsv with v2 -> v1 mapping
    cat_df = pd.DataFrame({
        "topic": v2_df["topic"],
        "category": v2_df["tier"].map(TIER_TO_NB_CATEGORY).fillna("unannotated"),
        "manual_override": "",
        "cell_types_list": v2_df["cell_types"],
        "n_cell_types": v2_df["n_cell_types"],
        "batch_flag": v2_df["batch_flag"],
    })
    cat_path = os.path.join(summary_dir, "topic_categories.tsv")
    cat_df.to_csv(cat_path, sep="\t", index=False)
    print(f"  Saved: {cat_path}  (v2 tiers mapped to notebook palette)")

    print("\n  Notebook-palette distribution:")
    for cat in ["one-to-one", "lineage", "shared", "structural",
                "batch-driven", "unannotated"]:
        n = (cat_df["category"] == cat).sum()
        if n > 0:
            print(f"    {cat:25s} {n:>3d}")


if __name__ == "__main__":
    main()
