#!/usr/bin/env python
"""Predict accessibility with a trained CREsted model and save the AnnData with
predictions + a "combined" layer (mean of ground truth and predictions), plus
the recommended top-k QC plot.

Run ONCE per model. Per-cell-type contribution jobs then load this h5ad,
run crested.pp.sort_and_filter_regions_on_specificity for the desired top_k,
and call crested.tl.contribution_scores_specific.
"""

import argparse
import json
import logging
import os
import sys

import numpy as np

logging.basicConfig(format="%(asctime)s %(levelname)-8s %(message)s",
                    level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)


def main(args):
    os.environ["KERAS_BACKEND"] = args.backend
    import keras
    import crested

    crested.register_genome(crested.Genome(args.genome, args.chrom_sizes))

    logger.info("Importing BigWigs from %s", args.bigwigs_dir)
    adata = crested.import_bigwigs(
        bigwigs_folder=args.bigwigs_dir,
        regions_file=args.regions_file,
        target_region_width=args.target_region_width,
        target=args.target,
    )
    logger.info("  %d cell types x %d regions", adata.n_obs, adata.n_vars)

    if args.splits_file:
        splits = json.load(open(args.splits_file))
        val_chroms = splits.get("valid", splits.get("val", args.val_chroms))
        test_chroms = splits.get("test", args.test_chroms)
    else:
        val_chroms, test_chroms = args.val_chroms, args.test_chroms

    crested.pp.train_val_test_split(adata, strategy="chr", val_chroms=val_chroms, test_chroms=test_chroms)
    crested.pp.change_regions_width(adata, args.seq_len)
    crested.pp.normalize_peaks(adata, top_k_percent=args.top_k_percent)

    checkpoint_path = args.checkpoint
    if checkpoint_path is None:
        ckpt_file = os.path.join(args.training_dir, "evaluation", "best_checkpoint.txt")
        if not os.path.exists(ckpt_file):
            logger.error("No best_checkpoint.txt; pass --checkpoint."); sys.exit(1)
        checkpoint_path = open(ckpt_file).read().strip()
    logger.info("Loading model from %s", checkpoint_path)
    model = keras.models.load_model(checkpoint_path, compile=False)

    logger.info("Predicting on %d regions...", adata.n_vars)
    preds = crested.tl.predict(input=adata, model=model, batch_size=args.predict_batch_size)
    adata.layers["predictions"] = preds.T

    X_dense = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    adata.layers["combined"] = (X_dense + adata.layers["predictions"]) / 2.0

    if args.qc_plot:
        logger.info("Plotting sort_and_filter_cutoff QC → %s", args.qc_plot)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        crested.pl.qc.sort_and_filter_cutoff(
            adata, model_name="combined", cutoffs=args.qc_cutoffs, max_k=args.qc_max_k,
        )
        os.makedirs(os.path.dirname(args.qc_plot) or ".", exist_ok=True)
        plt.savefig(args.qc_plot, bbox_inches="tight"); plt.close("all")

    os.makedirs(os.path.dirname(args.output_h5ad) or ".", exist_ok=True)
    adata.write_h5ad(args.output_h5ad)
    logger.info("Saved → %s", args.output_h5ad)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--training_dir", default=None)
    p.add_argument("--checkpoint", default=None)
    p.add_argument("--bigwigs_dir", required=True)
    p.add_argument("--regions_file", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--chrom_sizes", required=True)
    p.add_argument("--output_h5ad", required=True)
    p.add_argument("--qc_plot", default=None, help="Optional sort_and_filter_cutoff QC PDF path.")
    p.add_argument("--qc_cutoffs", type=int, nargs="+", default=[500, 1000, 2000])
    p.add_argument("--qc_max_k", type=int, default=5000)
    p.add_argument("--splits_file", default=None)
    p.add_argument("--target", default="count")
    p.add_argument("--target_region_width", type=int, default=None)
    p.add_argument("--top_k_percent", type=float, default=0.03)
    p.add_argument("--seq_len", type=int, default=2114)
    p.add_argument("--val_chroms", nargs="+", default=["chr8", "chr20"])
    p.add_argument("--test_chroms", nargs="+", default=["chr1", "chr3", "chr6"])
    p.add_argument("--predict_batch_size", type=int, default=128)
    p.add_argument("--backend", default="tensorflow")
    main(p.parse_args())
