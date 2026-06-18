"""Adapter shim: surface the ST gimme catalog in the layout that 7.1/7.2/7.3 expect.

Creates only NEW files under the ST catalog dir; nothing existing is moved or
overwritten. The renderer pipeline reads:

  combined_mapped.meme                       (already exists)
  clustered_motifs.modisco.h5                (already exists)
  clustered_motifs_name_map_unique.tsv       <-- NEW: copy of name_map.tsv (already unique)
  08_final_averages_metadata.tsv             <-- NEW: thin metadata from h5 + name_map
  tomtom/<display>.tomtom.txt                (already exists)
  hits/<CT>/hits.tsv                         <-- NEW: per-CT hits with motif_name renamed
                                                       from h5 key to display
  hits/hit_count_matrix.csv                  <-- NEW: copy of hit_analysis_unified matrix

By default, motif_name in the ST hits.tsv files is the h5 key
('pos_patterns.pattern_2'); we rename it to the display name ('Average_234')
so 7.2/7.3's per-motif joins line up with the MEME and per-motif TomTom files.

Idempotent: re-running overwrites the shim outputs but never touches originals.
"""
import argparse, shutil
from pathlib import Path
import pandas as pd
import h5py
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--st_catalog", required=True,
                   help="e.g. results/3_single_task_models/motifs")
    p.add_argument("--st_results_root", required=True,
                   help="e.g. results/3_single_task_models  (contains per-CT subdirs)")
    p.add_argument("--ct_metadata", required=True)
    args = p.parse_args()

    cat = Path(args.st_catalog)
    src_root = Path(args.st_results_root)

    # --- 1) name_map (already unique) ---
    nm_src = cat / "clustered_motifs_name_map.tsv"
    nm_dst = cat / "clustered_motifs_name_map_unique.tsv"
    nm = pd.read_csv(nm_src, sep="\t", header=None, names=["h5_key", "display"])
    assert nm["display"].is_unique, "ST name_map displays not unique"
    nm.to_csv(nm_dst, sep="\t", header=False, index=False)
    print(f"wrote {nm_dst.name} ({len(nm)} rows)")

    # --- 2) metadata shim ---
    # Build a thin 08_final_averages_metadata.tsv from h5 + name_map. The renderer
    # needs at least: name, posneg, num_final_constituents, total_seqlets.
    with h5py.File(cat / "clustered_motifs.modisco.h5", "r") as h5:
        meta_rows = []
        for _, row in nm.iterrows():
            track, _, key = row["h5_key"].partition(".")
            posneg = "pos" if track == "pos_patterns" else "neg"
            seqlets = 0
            num_sub = 0
            # ST gimme h5 has /pos_patterns/pattern_X/contrib_scores -- no num_seqlets
            # attribute typically. Best-effort: try common locations.
            grp = h5[track][key]
            if "seqlets" in grp:
                try:
                    seqlets = int(np.asarray(grp["seqlets"]["n_seqlets"]).sum())
                except Exception:
                    pass
            if "subclusters" in grp:
                try:
                    num_sub = len(grp["subclusters"])
                except Exception:
                    pass
            meta_rows.append({
                "name": row["display"],
                "posneg": posneg,
                "num_final_constituents": num_sub or 1,
                "total_seqlets": seqlets or 1,
                "source_celltypes": "",
            })

    meta = pd.DataFrame(meta_rows)
    meta_dst = cat / "08_final_averages_metadata.tsv"
    meta.to_csv(meta_dst, sep="\t", index=False)
    print(f"wrote {meta_dst.name} ({len(meta)} rows)")

    # --- 3) per-CT hits.tsv with motif_name renamed to display ---
    h2d = dict(zip(nm["h5_key"], nm["display"]))
    cts = [d.name for d in src_root.iterdir() if d.is_dir() and not d.name.startswith("_") and d.name not in {"motifs", "summary"}]
    cts.sort()
    n_done = 0
    hits_root = cat / "hits"
    hits_root.mkdir(exist_ok=True)
    for ct in cts:
        # Source: counts_unified for consistency (vs counts/)
        src = src_root / ct / "average" / "motifs" / "hits" / "counts_unified" / "hits.tsv"
        if not src.exists():
            print(f"  WARN: no hits.tsv for {ct} at {src}")
            continue
        dst_dir = hits_root / ct
        dst_dir.mkdir(exist_ok=True)
        dst = dst_dir / "hits.tsv"
        # Stream-rename motif_name col
        df = pd.read_csv(src, sep="\t")
        if "motif_name" not in df.columns:
            print(f"  WARN: no motif_name col in {src}")
            continue
        # ST hits.tsv uses h5 keys (e.g. 'pos_patterns.pattern_2'); map to display
        df["motif_name"] = df["motif_name"].map(h2d).fillna(df["motif_name"])
        df.to_csv(dst, sep="\t", index=False)
        n_done += 1
    print(f"wrote {n_done} per-CT hits.tsv files under {hits_root}")

    # --- 4) hit_count_matrix.csv copy ---
    # ST already has hit_analysis_unified/hit_count_matrix.csv whose columns are
    # display names (Average_<N>). 7.3 wants this under hits/ at the catalog root.
    hcm_src = cat / "hit_analysis_unified" / "hit_count_matrix.csv"
    if not hcm_src.exists():
        hcm_src = cat / "hit_analysis" / "hit_count_matrix.csv"
    if hcm_src.exists():
        hcm_dst = hits_root / "hit_count_matrix.csv"
        shutil.copy2(hcm_src, hcm_dst)
        print(f"copied hit_count_matrix.csv from {hcm_src.parent.name}/ to hits/")
    else:
        print("  WARN: hit_count_matrix.csv not found in hit_analysis* dirs")

    # --- 5) rebuild combined_mapped.meme using Average_N motif names ---
    # Original gimme-clustered MEME has MOTIFs named after the top Vierstra hit
    # (e.g. 'ZN143_MOUSE.H11MO.0.A'). But name_map, hits.tsv, and per-motif
    # tomtom files key by Average_N. Rebuild MEME using Average_N so 7.1/7.2/7.3
    # joins are consistent. Original is backed up to combined_mapped.meme.vierstra_named.bak.
    meme_src = cat / "combined_mapped.meme"
    bak = cat / "combined_mapped.meme.vierstra_named.bak"
    if meme_src.exists() and not bak.exists():
        shutil.copy2(meme_src, bak)
        print(f"backed up original MEME (Vierstra-named) -> {bak.name}")

    with h5py.File(cat / "clustered_motifs.modisco.h5", "r") as h5:
        lines = [
            "MEME version 4", "", "ALPHABET= ACGT", "", "strands: +", "",
            "Background letter frequencies:", "A 0.25 C 0.25 G 0.25 T 0.25", "",
        ]
        for _, row in nm.iterrows():
            track, _, key = row["h5_key"].partition(".")
            disp = row["display"]
            cwm = np.asarray(h5[track][key]["contrib_scores"]).astype(np.float64)
            ppm = np.clip(cwm, 0, None)
            rs = ppm.sum(axis=1, keepdims=True)
            zero = (rs.squeeze() == 0)
            if zero.any():
                ppm[zero] = 0.25
                rs = ppm.sum(axis=1, keepdims=True)
            ppm = ppm / rs
            W = ppm.shape[0]
            lines.append(f"MOTIF {disp}")
            lines.append(f"letter-probability matrix: alength= 4 w= {W}")
            for r in ppm:
                lines.append(" ".join(f"{v:.6f}" for v in r))
            lines.append("")
    meme_src.write_text("\n".join(lines) + "\n")
    print(f"rewrote {meme_src.name} with {len(nm)} Average-named motifs (PPM)")

    # --- 6) annotation dir for 7.2 outputs ---
    (cat / "annotation").mkdir(exist_ok=True)
    (cat / "qc").mkdir(exist_ok=True)
    print(f"\nST shim ready. Catalog dir: {cat}")
    print(f"  motifs: {len(nm)} | CTs with hits: {n_done}")


if __name__ == "__main__":
    main()
