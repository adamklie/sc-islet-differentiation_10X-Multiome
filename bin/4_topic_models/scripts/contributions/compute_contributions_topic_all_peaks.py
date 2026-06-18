#!/usr/bin/env python
"""Per-topic contribution scoring on ALL peaks for a CREsted topic CNN.

Topic-side analog of
``bin/4_multi_task_model/scripts/contributions/compute_contributions_all_peaks.py``:

  - Loads the predictions+combined AnnData built by
    ``all_n100_predict_and_save_adata.sh``  (obs = 100 topics; var = consensus
    peaks; layers = ``predictions``, ``combined``).
  - Does NOT call ``sort_and_filter_regions_on_specificity`` -> uses ALL input
    peaks.
  - Iterates regions in chunks to keep memory bounded.
  - Saves per-topic NPZ pair compatible with finemo:
        {Topic}_contrib.npz : (N_regions, 4, L)
        {Topic}_oh.npz      : (N_regions, 4, L)

Output layout matches the MT all-peaks run so the same finemo hits pipeline
shape applies (one subdir per topic under ``hits_all_peaks/``).
"""

import argparse
import logging
import os
import sys
import time

logging.basicConfig(format="%(asctime)s %(levelname)-8s %(message)s",
                    level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)


def main(args):
    os.environ["KERAS_BACKEND"] = args.backend
    import keras
    import numpy as np
    import anndata as ad
    import crested

    crested.register_genome(crested.Genome(args.genome, args.chrom_sizes))

    logger.info("Loading prepared AnnData from %s", args.input_h5ad)
    adata = ad.read_h5ad(args.input_h5ad)
    logger.info("  %d topics x %d regions; layers=%s",
                adata.n_obs, adata.n_vars, list(adata.layers.keys()))

    topic_names = list(adata.obs_names)
    targets = args.topics or topic_names
    missing = [t for t in targets if t not in topic_names]
    if missing:
        logger.warning("Topics not found, skipping: %s", missing)
    targets = [t for t in targets if t in topic_names]
    if not targets:
        logger.error("No valid topics"); sys.exit(1)

    logger.info("Loading model from %s", args.checkpoint)
    model = keras.models.load_model(args.checkpoint, compile=False)

    os.makedirs(args.output_dir, exist_ok=True)

    all_regions = list(adata.var_names)
    N = len(all_regions)
    chunk = args.chunk_size
    logger.info("Will compute contributions over ALL %d regions; chunk_size=%d", N, chunk)

    for topic in targets:
        target_id = topic_names.index(topic)
        contrib_out = os.path.join(args.output_dir, f"{topic}_contrib.npz")
        oh_out = os.path.join(args.output_dir, f"{topic}_oh.npz")
        if os.path.exists(contrib_out) and os.path.exists(oh_out):
            if args.skip_if_exists:
                logger.info("Skip %s (outputs exist).", topic)
                continue
            else:
                logger.info("Outputs exist for %s, but --skip_if_exists not set; will recompute.", topic)

        t0 = time.time()
        logger.info("=== %s (target_idx=%d) - %d regions in chunks of %d ===",
                    topic, target_id, N, chunk)

        # Per-chunk shape: (chunk, 1, 4, L) contrib + (chunk, 4, L) oh.
        # Total in-memory: N * 4 * L * 4B (contrib float32) + N * 4 * L * 4B (oh float32)
        # ~= 8 * N * 4 * L bytes. For N=605K, L=2114: contrib ~20 GB, oh ~20 GB -> need ~50GB.
        # 64G mem allocation should fit with model + per-chunk transients.
        contrib_all = None
        oh_all = None
        n_done = 0
        for s in range(0, N, chunk):
            e = min(N, s + chunk)
            regs = all_regions[s:e]
            t_chunk = time.time()
            scores, ohs = crested.tl.contribution_scores(
                input=regs,
                target_idx=[target_id],
                model=model,
                method=args.method,
                transpose=True,            # gives (N, C, 4, L) and (N, 4, L)
                batch_size=args.batch_size,
                output_dir=None,           # we save the concatenated result ourselves
                verbose=False,
            )
            # scores: (n, 1, 4, L) -> squeeze channel dim to (n, 4, L)
            scores_2d = scores[:, 0, :, :].astype(np.float32, copy=False)
            ohs_2d = ohs.astype(np.float32, copy=False)

            if contrib_all is None:
                _, ch, L = scores_2d.shape
                contrib_all = np.empty((N, ch, L), dtype=np.float32)
                oh_all = np.empty((N, ch, L), dtype=np.float32)

            contrib_all[s:e] = scores_2d
            oh_all[s:e] = ohs_2d
            n_done = e
            logger.info("  %s chunk %d-%d: %.1fs  (%.1f regions/s)",
                        topic, s, e, time.time() - t_chunk,
                        (e - s) / max(1e-6, time.time() - t_chunk))

        logger.info("  Writing %s ...", contrib_out)
        np.savez_compressed(contrib_out, contrib_all)
        logger.info("  Writing %s ...", oh_out)
        np.savez_compressed(oh_out, oh_all)
        logger.info("=== %s done in %.1f min (%d regions) ===",
                    topic, (time.time() - t0) / 60, n_done)
        del contrib_all, oh_all

    logger.info("All requested topics complete.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input_h5ad", required=True, help="From predict_and_save_adata_topics.py")
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--chrom_sizes", required=True)
    p.add_argument("--output_dir", required=True)
    p.add_argument("--topics", nargs="+", default=None,
                   help="Subset (default: all topics in AnnData). Use one per array task.")
    p.add_argument("--method", default="integrated_grad",
                   choices=["integrated_grad", "expected_integrated_grad", "saliency_map", "mutagenesis"])
    p.add_argument("--batch_size", type=int, default=128)
    p.add_argument("--chunk_size", type=int, default=20000,
                   help="Number of regions per CREsted call (bounds peak GPU/CPU memory).")
    p.add_argument("--skip_if_exists", action="store_true",
                   help="If both NPZs already exist for a topic, skip it.")
    p.add_argument("--backend", default="torch",
                   help="Topic model was trained with torch backend; matches existing topic pipeline.")
    main(p.parse_args())
