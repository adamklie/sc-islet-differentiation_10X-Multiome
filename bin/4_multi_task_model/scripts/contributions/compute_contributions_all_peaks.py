#!/usr/bin/env python
"""Per-cell-type contribution scoring on ALL peaks for a CREsted peak regression
model. Mirrors compute_contributions_specific_v2.py BUT:

  - Does NOT call sort_and_filter_regions_on_specificity (uses all input peaks)
  - Iterates regions in chunks to keep memory bounded
  - Saves per-CT NPZ pair compatible with finemo:
        {CT}_contrib.npz : (N_regions, 4, L)
        {CT}_oh.npz      : (N_regions, 4, L)

Output layout matches the v2/v3 specific runs so the same hits pipeline works.
"""

import argparse
import logging
import os
import sys
import time

logging.basicConfig(format="%(asctime)s %(levelname)-8s %(message)s",
                    level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)


def _load_narrowpeak_regions(bed_path: str, window: int):
    """Generate "chr:start-end" window strings centered on each peak's summit.

    narrowPeak format (10 columns): chrom, chromStart, chromEnd, name, score,
    strand, signalValue, pValue, qValue, peak (summit offset relative to start).
    Returns regions in the SAME order as rows in the BED so the resulting
    contributions NPZ row order matches the peak file (finemo's required
    alignment).
    """
    import pandas as pd
    half = window // 2
    df = pd.read_csv(bed_path, header=None, sep="\t",
                     usecols=[0, 1, 2, 9], names=["chrom", "start", "end", "summit"])
    df["center"] = df["start"] + df["summit"]
    df["wstart"] = df["center"] - half
    df["wend"] = df["center"] + (window - half)
    # Width check
    widths = df["wend"] - df["wstart"]
    if not (widths == window).all():
        raise SystemExit(f"derived window widths != {window}: min={widths.min()} max={widths.max()}")
    return [f"{c}:{s}-{e}" for c, s, e in zip(df["chrom"], df["wstart"], df["wend"])]


def main(args):
    os.environ["KERAS_BACKEND"] = args.backend
    import keras
    import numpy as np
    import anndata as ad
    import crested

    crested.register_genome(crested.Genome(args.genome, args.chrom_sizes))

    logger.info("Loading prepared AnnData from %s", args.input_h5ad)
    adata = ad.read_h5ad(args.input_h5ad)
    logger.info("  %d cell types x %d regions; layers=%s",
                adata.n_obs, adata.n_vars, list(adata.layers.keys()))

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
    in_len = model.input_shape[1]
    logger.info("Model input length: %d", in_len)

    os.makedirs(args.output_dir, exist_ok=True)

    # Region source: either per-CT narrowPeak (new) or full adata var_names (legacy)
    if args.peaks_bed:
        if len(targets) != 1:
            raise SystemExit("--peaks_bed requires exactly one --cell_types target "
                             "(per-CT narrowPeak is single-CT by construction).")
        all_regions = _load_narrowpeak_regions(args.peaks_bed, window=in_len)
        logger.info("Loaded %d %dbp summit-centered regions from %s",
                    len(all_regions), in_len, args.peaks_bed)
    else:
        all_regions = list(adata.var_names)
    N = len(all_regions)
    chunk = args.chunk_size
    logger.info("Will compute contributions over %d regions; chunk_size=%d", N, chunk)

    for ct in targets:
        target_id = cell_type_names.index(ct)
        contrib_out = os.path.join(args.output_dir, f"{ct}_contrib.npz")
        oh_out = os.path.join(args.output_dir, f"{ct}_oh.npz")
        if os.path.exists(contrib_out) and os.path.exists(oh_out):
            if args.skip_if_exists:
                logger.info("Skip %s (outputs exist).", ct)
                continue
            else:
                logger.info("Outputs exist for %s, but --skip_if_exists not set; will recompute.", ct)

        t0 = time.time()
        logger.info("=== %s (target_idx=%d) — %d regions in chunks of %d ===",
                    ct, target_id, N, chunk)

        # Pre-allocate full output arrays on disk-like in-memory float32 to keep memory bounded.
        # Per-chunk shape: (chunk, 1, 4, L) contrib + (chunk, 4, L) oh.
        # Total memory: N * 4 * L * 4B (contrib float32) + N * 4 * L * 1B (oh) ~= 5*N*4*L bytes.
        # For N=605K, L=2114: contrib ~20 GB, oh ~5 GB. Fits in 64G mem allocation.
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
                output_dir=None,           # we'll save the concatenated result ourselves
                verbose=False,
            )
            # scores: (n, 1, 4, L); we'll squeeze channel dim to match v2 NPZ shape (n, 4, L).
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
                        ct, s, e, time.time() - t_chunk,
                        (e - s) / max(1e-6, time.time() - t_chunk))

        logger.info("  Writing %s ...", contrib_out)
        np.savez_compressed(contrib_out, contrib_all)
        logger.info("  Writing %s ...", oh_out)
        np.savez_compressed(oh_out, oh_all)
        logger.info("=== %s done in %.1f min (%d regions) ===",
                    ct, (time.time() - t0) / 60, n_done)
        # Free
        del contrib_all, oh_all

    logger.info("All requested CTs complete.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input_h5ad", required=True, help="From predict_and_save_adata.py")
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--chrom_sizes", required=True)
    p.add_argument("--output_dir", required=True)
    p.add_argument("--cell_types", nargs="+", default=None,
                   help="Subset (default: all). Use one per array task.")
    p.add_argument("--method", default="integrated_grad",
                   choices=["integrated_grad", "expected_integrated_grad", "saliency_map", "mutagenesis"])
    p.add_argument("--batch_size", type=int, default=128)
    p.add_argument("--chunk_size", type=int, default=20000,
                   help="Number of regions per CREsted call (bounds peak GPU/CPU memory).")
    p.add_argument("--skip_if_exists", action="store_true",
                   help="If both NPZs already exist for a CT, skip it.")
    p.add_argument("--backend", default="tensorflow")
    p.add_argument("--peaks_bed", default=None,
                   help="Per-CT narrowPeak BED. If given, contributions are computed "
                        "on summit-centered windows of model.input_length (typically "
                        "2114bp) IN ROW ORDER. Required to align with finemo's "
                        "per-CT hit caller. Without this flag, falls back to "
                        "adata.var_names (consensus regions, e.g. 605k).")
    main(p.parse_args())
