#!/usr/bin/env python
"""Predict topic probabilities with a trained CREsted topic model and save an
AnnData with predictions + "combined" layer, mirroring the CREsted
enhancer_code_analysis recipe. Topic-model variant of predict_and_save_adata.py.

Differences vs peak-regression variant:
  - Input is per-topic BED directory (crested.import_beds) instead of BigWigs.
  - No train_val_test_split / change_regions_width / normalize_peaks needed
    (topic models are trained directly on binarized membership).
"""

import argparse
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

    logger.info("Importing BED files from %s", args.beds_dir)
    kw = {"beds_folder": args.beds_dir}
    if args.regions_file:
        kw["regions_file"] = args.regions_file
    adata = crested.import_beds(**kw)
    logger.info("  %d topics x %d regions", adata.n_obs, adata.n_vars)

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
    adata.layers["predictions"] = preds.T  # (classes, regions)

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
    # Drop obs columns that anndata can't serialize (e.g. crested.import_beds
    # populates `file_path` with object-typed values that fail vlen-string
    # serialization). We don't need them for downstream contribution scoring.
    for col in ("file_path",):
        if col in adata.obs.columns:
            del adata.obs[col]
    adata.write_h5ad(args.output_h5ad)
    logger.info("Saved → %s", args.output_h5ad)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--training_dir", default=None)
    p.add_argument("--checkpoint", default=None)
    p.add_argument("--beds_dir", required=True)
    p.add_argument("--regions_file", default=None, help="Optional consensus BED to intersect with.")
    p.add_argument("--genome", required=True)
    p.add_argument("--chrom_sizes", required=True)
    p.add_argument("--output_h5ad", required=True)
    p.add_argument("--qc_plot", default=None)
    p.add_argument("--qc_cutoffs", type=int, nargs="+", default=[500, 1000, 2000])
    p.add_argument("--qc_max_k", type=int, default=5000)
    p.add_argument("--predict_batch_size", type=int, default=128)
    p.add_argument("--backend", default="torch")
    main(p.parse_args())
