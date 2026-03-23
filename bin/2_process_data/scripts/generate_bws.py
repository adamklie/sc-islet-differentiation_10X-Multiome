#!/usr/bin/env python

import argparse
import glob
import gzip
import logging
import os
import subprocess

import snapatac2 as snap
from concurrent.futures import ProcessPoolExecutor, as_completed


def setup_logger(outdir, tag):
    os.makedirs(outdir, exist_ok=True)
    log_file = os.path.join(outdir, f"{tag}.log")

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        filename=log_file,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console.setFormatter(formatter)
    logging.getLogger("").addHandler(console)

    return log_file


def process_celltype(celltype, frag_bed, chr_sizes_path, bedtools, bedgraph_to_bw, keep_bdg=False):
    try:
        scale_out = frag_bed.replace(".bed.gz", ".scale_factor.txt")
        bw_out = frag_bed.replace(".bed.gz", ".bw")
        bdg_out = frag_bed.replace(".bed.gz", ".bdg")

        if os.path.exists(bw_out) and os.path.exists(scale_out):
            return f"{celltype}: SKIPPED (already has {bw_out} + {scale_out})"

        # count fragments
        with gzip.open(frag_bed, 'rt') as f:
            n_fragments = sum(1 for _ in f)
        scale_factor = 1_000_000 / (2 * n_fragments)

        # write scale factor
        with open(scale_out, "w") as f:
            f.write(str(scale_factor) + "\n")

        # make .bdg
        if not os.path.exists(bdg_out):
            with open(bdg_out, 'w') as f_bdg:
                sort_proc = subprocess.Popen([bedtools, "sort", "-i", frag_bed], stdout=subprocess.PIPE)
                cov_proc = subprocess.Popen(
                    [bedtools, "genomecov", "-bg", "-i", "/dev/stdin",
                     "-g", chr_sizes_path, "-scale", str(scale_factor)],
                    stdin=sort_proc.stdout,
                    stdout=f_bdg,
                )
                sort_proc.stdout.close()
                sort_proc.wait()
                cov_proc.wait()

        # convert to bigWig
        if not os.path.exists(bw_out):
            subprocess.run([bedgraph_to_bw, bdg_out, chr_sizes_path, bw_out], check=True)

        if not keep_bdg and os.path.exists(bdg_out):
            os.remove(bdg_out)

        return f"{celltype}: {n_fragments} fragments, scale={scale_factor:.6g}, wrote {bw_out}"

    except Exception as e:
        logging.error(f"Error in {celltype}: {e}")
        return f"ERROR processing {celltype}: {str(e)}"


def run_pipeline(args):

    # --- Setup ---
    os.makedirs(args.outdir, exist_ok=True)
    log_file = setup_logger(args.outdir, "generate_bws")

    # --- Export fragments or discover existing ones ---
    if args.skip_export:
        # Verify chr_sizes.tsv exists
        chr_sizes_path = os.path.join(args.outdir, "chr_sizes.tsv")
        if not os.path.exists(chr_sizes_path):
            raise FileNotFoundError(
                f"{chr_sizes_path} not found. Cannot use --skip-export without existing chr_sizes.tsv"
            )

        # Discover existing .bed.gz files
        exported = {}
        for bed in glob.glob(os.path.join(args.outdir, "*.bed.gz")):
            ct = os.path.basename(bed).replace(".bed.gz", "")
            exported[ct] = bed
        logging.info(f"Skipping export, found {len(exported)} existing .bed.gz files")

    else:
        logging.info(f"=== Processing AnnDataset {args.anndataset} ===")
        dataset = snap.read_dataset(args.anndataset)
        logging.info(f"Loaded AnnDataset with {dataset.n_obs} cells")

        # Export chr sizes
        chr_sizes_path = os.path.join(args.outdir, "chr_sizes.tsv")
        dataset.uns['reference_sequences'].to_pandas().to_csv(
            chr_sizes_path, sep="\t", header=False, index=False
        )

        # Export fragments by cell type
        exported = snap.ex.export_fragments(
            dataset,
            groupby="cell_type",
            out_dir=args.outdir,
            suffix=".bed.gz"
        )

    # --- Parallel process each cell type ---
    results = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(process_celltype, ct, bed, chr_sizes_path,
                            args.bedtools, args.bedgraph_to_bw, args.keep_bdg): ct
            for ct, bed in exported.items()
        }
        for fut in as_completed(futures):
            results.append(fut.result())

    for r in results:
        logging.info(r)

    logging.info("=== Finished ===")
    logging.info(f"Log written to {log_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export normalized bigWigs directly from processed AnnDataset")
    parser.add_argument("--anndataset", required=True, help="Path to AnnDataset .h5ads")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--bedtools", default="bedtools", help="Path to bedtools executable")
    parser.add_argument("--bedgraph_to_bw", default="bedGraphToBigWig", help="Path to UCSC bedGraphToBigWig")
    parser.add_argument("--keep-bdg", action="store_true", help="Keep intermediate .bed files instead of deleting them")
    parser.add_argument("--skip-export", action="store_true", help="Skip fragment export and AnnDataset loading; reuse existing .bed.gz and chr_sizes.tsv")
    parser.add_argument("--threads", type=int, default=4, help="Parallel workers per dataset")
    args = parser.parse_args()
    run_pipeline(args)
    