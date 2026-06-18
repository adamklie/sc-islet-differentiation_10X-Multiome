#!/usr/bin/env python
"""Per-cell-type contribution scoring for a CREsted peak regression model,
following the CREsted enhancer-code analysis recipe.

Inputs:
  --input_h5ad : AnnData produced by predict_and_save_adata.py
                 (must have layers["combined"]).
  --checkpoint : trained model checkpoint (same one used to build the h5ad).
  --cell_types : one or more class names. For SLURM array jobs, pass one per task.

Per task:
  1. Load the AnnData (predictions + combined layer already computed).
  2. Run crested.pp.sort_and_filter_regions_on_specificity(top_k=..., model_name="combined").
  3. Subset to the requested cell type(s) via the "Class name" annotation.
  4. Call crested.tl.contribution_scores_specific → writes one NPZ per class
     into --output_dir, ready for tfmodisco-lite.
"""

import argparse
import logging
import os
import sys

logging.basicConfig(format="%(asctime)s %(levelname)-8s %(message)s",
                    level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)


def main(args):
    os.environ["KERAS_BACKEND"] = args.backend
    import keras
    import anndata as ad
    import crested

    crested.register_genome(crested.Genome(args.genome, args.chrom_sizes))

    logger.info("Loading prepared AnnData from %s", args.input_h5ad)
    adata = ad.read_h5ad(args.input_h5ad)
    logger.info("  %d cell types x %d regions; layers=%s",
                adata.n_obs, adata.n_vars, list(adata.layers.keys()))
    if "combined" not in adata.layers:
        logger.error("Input AnnData missing layers['combined']; rerun predict_and_save_adata.py")
        sys.exit(1)

    cell_type_names = list(adata.obs_names)
    targets = args.cell_types or cell_type_names
    missing = [ct for ct in targets if ct not in cell_type_names]
    if missing:
        logger.warning("Cell types not found, skipping: %s", missing)
    targets = [ct for ct in targets if ct in cell_type_names]
    if not targets:
        logger.error("No valid cell types"); sys.exit(1)

    logger.info("Loading model from %s", args.checkpoint)
    model = keras.models.load_model(args.checkpoint, compile=False)

    logger.info("sort_and_filter_regions_on_specificity(top_k=%d, model_name='combined', method='gini')",
                args.top_k_regions)
    adata_filtered = crested.pp.sort_and_filter_regions_on_specificity(
        adata, top_k=args.top_k_regions, model_name="combined", method="gini", inplace=False,
    )
    logger.info("  Filtered: %d regions across %d classes (~%d/class)",
                adata_filtered.n_vars, adata_filtered.n_obs, args.top_k_regions)

    target_idx = [cell_type_names.index(ct) for ct in targets]
    logger.info("Computing contributions for %s (target_idx=%s, method=%s) → %s",
                targets, target_idx, args.method, args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)
    crested.tl.contribution_scores_specific(
        input=adata_filtered,
        target_idx=target_idx,
        model=model,
        method=args.method,
        transpose=True,
        batch_size=args.batch_size,
        output_dir=args.output_dir,
        verbose=True,
    )
    logger.info("Done.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input_h5ad", required=True, help="From predict_and_save_adata.py")
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--chrom_sizes", required=True)
    p.add_argument("--output_dir", required=True)
    p.add_argument("--cell_types", nargs="+", default=None,
                   help="Subset (default: all). Use one per array task.")
    p.add_argument("--top_k_regions", type=int, default=1000)
    p.add_argument("--method", default="integrated_grad",
                   choices=["integrated_grad", "expected_integrated_grad", "saliency_map", "mutagenesis"])
    p.add_argument("--batch_size", type=int, default=128)
    p.add_argument("--backend", default="tensorflow")
    main(p.parse_args())
