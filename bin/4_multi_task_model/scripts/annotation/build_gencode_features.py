"""One-shot: build GENCODE v43 TSS bed + non-overlapping feature bed for hg38.

Outputs (under {ANN_DIR}/):
  gencode_v43_tss.bed          — chr, start, end (=start+1), gene_id, gene_name, strand
  gencode_v43_features.bed     — chr, start, end, feature, gene_id (sorted, NON-overlapping
                                  by collapse priority: promoter > 5UTR > 3UTR > exon > intron > distal)

Feature priority follows ArchR's createGeneAnnotation convention (matches the
Liu 2026 paper's annotation source). Promoter window is TSS ± 2kb. 'distal'
is not written as a feature; absence-from-features-bed = distal.
"""
import argparse, gzip, subprocess, tempfile
from pathlib import Path
import pandas as pd

GTF = "/cellar/users/aklie/data/ref/hg38-ref-data-2026/gencode.v43.primary_assembly.annotation.gtf.gz"
PROMOTER_HALF = 2000  # TSS +/- 2kb


STD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}


def parse_gtf_records(gtf_path):
    """Yield (chr, start, end, strand, feature, gene_id, gene_name) per GTF row.
    Filters to STD_CHROMS only (drops GL/KI scaffold contigs).
    1-based GTF -> 0-based BED is handled at write time.
    """
    op = gzip.open if gtf_path.endswith(".gz") else open
    with op(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, src, feat, start, end, _, strand, _, attrs = parts
            if chrom not in STD_CHROMS:
                continue
            if feat not in {"gene", "exon", "UTR"}:
                continue
            ad = {}
            for f_attr in attrs.split(";"):
                f_attr = f_attr.strip()
                if not f_attr:
                    continue
                k, _, v = f_attr.partition(" ")
                ad[k] = v.strip('"')
            yield chrom, int(start) - 1, int(end), strand, feat, ad.get("gene_id", ""), ad.get("gene_name", ""), ad.get("gene_type", ""), ad.get("transcript_id", "")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", default=GTF)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--chrom_sizes", default="/cellar/users/aklie/data/ref/genomes/hg38/hg38.chrom.sizes")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- 1) Parse GTF into per-feature lists ---
    print(f"[1/4] Parsing {args.gtf}", flush=True)
    gene_rows, exon_rows, utr_rows = [], [], []
    for chrom, s, e, strand, feat, gid, gname, gtype, tid in parse_gtf_records(args.gtf):
        if gtype != "protein_coding":
            continue
        if feat == "gene":
            gene_rows.append((chrom, s, e, strand, gid, gname))
        elif feat == "exon":
            exon_rows.append((chrom, s, e, strand, gid))
        elif feat == "UTR":
            utr_rows.append((chrom, s, e, strand, gid, tid))
    print(f"  genes={len(gene_rows)} exons={len(exon_rows)} UTRs={len(utr_rows)}", flush=True)

    # --- 2) TSS bed (one per gene) ---
    print("[2/4] TSS bed", flush=True)
    tss_rows = []
    for chrom, s, e, strand, gid, gname in gene_rows:
        tss = s if strand == "+" else e - 1
        tss_rows.append((chrom, tss, tss + 1, gid, gname, strand))
    tss = pd.DataFrame(tss_rows, columns=["chr", "start", "end", "gene_id", "gene_name", "strand"])
    tss = tss.sort_values(["chr", "start"]).reset_index(drop=True)
    tss_path = out_dir / "gencode_v43_tss.bed"
    tss.to_csv(tss_path, sep="\t", index=False, header=False)
    print(f"  wrote {tss_path} ({len(tss)} TSSs)", flush=True)

    # --- 3) Per-category beds for promoter/5UTR/3UTR/exon/intron ---
    # Promoter = TSS +/- PROMOTER_HALF
    print("[3/4] Per-category bed assembly", flush=True)
    promoters = tss.assign(
        start=(tss["start"] - PROMOTER_HALF).clip(lower=0),
        end=tss["start"] + PROMOTER_HALF,
    )[["chr", "start", "end", "gene_id"]].assign(feature="promoter")

    # 5UTR / 3UTR: GTF's UTR rows. To split into 5' vs 3', we need transcript orientation.
    # For each UTR, decide by position relative to its parent gene's TSS:
    # +strand gene: UTR with end <= gene CDS start (approx using gene start) = 5UTR.
    # Simpler/canonical: 5UTR = the UTR whose end coincides with start of CDS (or whose
    # start coincides with gene's start codon for + strand). We don't have CDS rows
    # loaded; approximate via gene start/strand:
    gene_df = pd.DataFrame(gene_rows, columns=["chr", "g_start", "g_end", "strand", "gene_id", "gene_name"])
    utr_df = pd.DataFrame(utr_rows, columns=["chr", "start", "end", "strand", "gene_id", "transcript_id"])
    utr_df = utr_df.merge(gene_df[["gene_id", "g_start", "g_end", "strand"]].rename(columns={"strand": "g_strand"}), on="gene_id", how="left")
    # 5UTR: on + strand, UTR start near g_start (i.e. mid < gene midpoint).
    utr_df["mid"] = (utr_df["start"] + utr_df["end"]) // 2
    utr_df["g_mid"] = (utr_df["g_start"] + utr_df["g_end"]) // 2
    is_5utr = ((utr_df["g_strand"] == "+") & (utr_df["mid"] <= utr_df["g_mid"])) | \
              ((utr_df["g_strand"] == "-") & (utr_df["mid"] >= utr_df["g_mid"]))
    utr_df["feature"] = "3UTR"
    utr_df.loc[is_5utr, "feature"] = "5UTR"

    utr_5 = utr_df.loc[utr_df["feature"] == "5UTR", ["chr", "start", "end", "gene_id"]].assign(feature="5UTR")
    utr_3 = utr_df.loc[utr_df["feature"] == "3UTR", ["chr", "start", "end", "gene_id"]].assign(feature="3UTR")

    exon_df = pd.DataFrame(exon_rows, columns=["chr", "start", "end", "strand", "gene_id"])[["chr", "start", "end", "gene_id"]].assign(feature="exon")

    # Intron = gene span minus exons. Build per-gene; precompute via bedtools subtract.
    gene_bed = gene_df[["chr", "g_start", "g_end", "gene_id"]].rename(columns={"g_start": "start", "g_end": "end"})

    bt = "/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/bedtools"

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        # exon sort/merge for subtract
        exon_sorted = td / "exons.sorted.bed"
        exon_df.sort_values(["chr", "start"])[["chr", "start", "end"]].to_csv(exon_sorted, sep="\t", index=False, header=False)
        exon_merged = td / "exons.merged.bed"
        subprocess.run(f"{bt} merge -i {exon_sorted} > {exon_merged}", shell=True, check=True)
        # gene sorted
        gene_sorted = td / "genes.sorted.bed"
        gene_bed.sort_values(["chr", "start"]).to_csv(gene_sorted, sep="\t", index=False, header=False)
        # intron = gene - merged_exons
        intron_path = td / "introns.bed"
        subprocess.run(f"{bt} subtract -a {gene_sorted} -b {exon_merged} > {intron_path}", shell=True, check=True)
        intron_df = pd.read_csv(intron_path, sep="\t", header=None, names=["chr", "start", "end", "gene_id"])
        intron_df["feature"] = "intron"

    # --- 4) Assemble priority-resolved features bed ---
    # Priority: promoter > 5UTR > 3UTR > exon > intron
    # Use 3-col BEDs throughout the subtract chain (avoids bedtools type-check
    # complaints on heterogeneous column counts); reattach feature label per
    # category after subtraction.
    print("[4/4] Priority-resolved features bed", flush=True)
    bt = "/cellar/users/aklie/opt/miniconda3/envs/eugene_tools/bin/bedtools"

    priorities = [
        ("promoter", promoters[["chr", "start", "end"]]),
        ("5UTR",     utr_5[["chr", "start", "end"]]),
        ("3UTR",     utr_3[["chr", "start", "end"]]),
        ("exon",     exon_df[["chr", "start", "end"]]),
        ("intron",   intron_df[["chr", "start", "end"]]),
    ]
    kept_dfs = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        cumulative_mask = None  # path to merged union of all kept spans so far

        for label, df in priorities:
            df = df[df["chr"].isin(STD_CHROMS)].copy()
            df["start"] = df["start"].astype(int)
            df["end"]   = df["end"].astype(int)
            df = df[df["end"] > df["start"]]
            if df.empty:
                continue

            # write + sort + merge category to a clean 3-col bed
            raw = td / f"{label}.raw.bed"
            srt = td / f"{label}.merged.bed"
            df.sort_values(["chr", "start"]).to_csv(raw, sep="\t", index=False, header=False)
            subprocess.run(f"sort -k1,1 -k2,2n {raw} | {bt} merge -i - > {srt}",
                           shell=True, check=True)

            # trim by cumulative mask (higher-priority spans)
            if cumulative_mask is None:
                trimmed_path = srt
            else:
                trimmed_path = td / f"{label}.trim.bed"
                subprocess.run(f"{bt} subtract -a {srt} -b {cumulative_mask} > {trimmed_path}",
                               shell=True, check=True)

            trimmed_df = pd.read_csv(trimmed_path, sep="\t", header=None,
                                     names=["chr", "start", "end"])
            trimmed_df["feature"] = label
            kept_dfs.append(trimmed_df)

            # update cumulative mask
            new_mask = td / f"mask_after_{label}.bed"
            if cumulative_mask is None:
                subprocess.run(f"cp {srt} {new_mask}", shell=True, check=True)
            else:
                subprocess.run(
                    f"cat {cumulative_mask} {trimmed_path} | sort -k1,1 -k2,2n | {bt} merge -i - > {new_mask}",
                    shell=True, check=True,
                )
            cumulative_mask = new_mask

        feat_df = pd.concat(kept_dfs, ignore_index=True).sort_values(["chr", "start"])
        feat_path = out_dir / "gencode_v43_features.bed"
        feat_df[["chr", "start", "end", "feature"]].to_csv(feat_path, sep="\t", index=False, header=False)
        print(f"  wrote {feat_path} ({len(feat_df)} non-overlapping intervals)", flush=True)
        print(feat_df["feature"].value_counts().to_string(), flush=True)


if __name__ == "__main__":
    main()
