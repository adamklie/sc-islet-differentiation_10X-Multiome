#!/usr/bin/env python
"""Compute per-cell-type contribution scores using CREsted for peak regression models.

For each cell type (or a subset), runs CREsted's contribution_scores on the
cell-type-specific regions (high accessibility). Outputs NPZ files compatible
with modisco-lite.

Outputs per cell type (in --output_dir/{cell_type}/):
    hypothetical_contribs.npz  — (N, 4, L) hypothetical attribution scores
    projected_contribs.npz     — (N, 4, L) projected scores (scores * one_hot)
    one_hot_seqs.npz           — (N, 4, L) one-hot encoded sequences
    region_names.npy           — region identifiers

These can be fed directly to modisco-lite:
    modisco motifs -s one_hot_seqs.npz -a hypothetical_contribs.npz -n 50000 -o modisco.h5
"""

import argparse
import json
import logging
import os
import sys

import numpy as np

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def main(args):
    os.environ["KERAS_BACKEND"] = args.backend

    import keras

    import crested

    # --- Register genome ---
    genome = crested.Genome(args.genome, args.chrom_sizes)
    crested.register_genome(genome)

    # --- Import BigWig data ---
    logger.info("Importing BigWig files from %s", args.bigwigs_dir)
    adata = crested.import_bigwigs(
        bigwigs_folder=args.bigwigs_dir,
        regions_file=args.regions_file,
        target_region_width=args.target_region_width,
        target=args.target,
    )
    cell_type_names = list(adata.obs_names)
    logger.info("  %d cell types x %d regions", adata.n_obs, adata.n_vars)

    # --- Preprocessing (must match training) ---
    if args.splits_file:
        splits = json.load(open(args.splits_file))
        val_chroms = splits.get("valid", splits.get("val", args.val_chroms))
        test_chroms = splits.get("test", args.test_chroms)
    else:
        val_chroms = args.val_chroms
        test_chroms = args.test_chroms

    crested.pp.train_val_test_split(
        adata, strategy="chr", val_chroms=val_chroms, test_chroms=test_chroms,
    )
    crested.pp.change_regions_width(adata, args.seq_len)
    crested.pp.normalize_peaks(adata, top_k_percent=args.top_k_percent)

    # --- Load model ---
    checkpoint_path = args.checkpoint
    if checkpoint_path is None:
        ckpt_file = os.path.join(args.training_dir, "evaluation", "best_checkpoint.txt")
        if os.path.exists(ckpt_file):
            with open(ckpt_file) as f:
                checkpoint_path = f.read().strip()
        else:
            logger.error("No best_checkpoint.txt found. Specify --checkpoint explicitly.")
            sys.exit(1)

    logger.info("Loading model from %s", checkpoint_path)
    model = keras.models.load_model(checkpoint_path, compile=False)

    # --- Determine cell types to process ---
    if args.cell_types:
        cell_types = args.cell_types
    else:
        cell_types = cell_type_names

    logger.info("Will compute contribution scores for %d cell types", len(cell_types))

    # --- Compute contribution scores per cell type ---
    for ct_name in cell_types:
        if ct_name not in cell_type_names:
            logger.warning("Cell type %s not found, skipping", ct_name)
            continue

        ct_idx = cell_type_names.index(ct_name)
        logger.info("Processing %s (idx=%d)...", ct_name, ct_idx)

        # Get accessibility matrix (all cell types x regions)
        X_all = adata.X
        if hasattr(X_all, "toarray"):
            X_all = X_all.toarray()
        else:
            X_all = np.array(X_all)

        ct_values = X_all[ct_idx, :]

        # Cell-type-specific region selection:
        # 1. Compute Gini per region (how variable across cell types)
        # 2. Filter to Gini > mean + 1*std (cell-type-variable regions)
        # 3. Z-score each region across cell types
        # 4. Take top K by z-score for this cell type
        def _gini(x):
            x = np.sort(np.abs(x))
            n = len(x)
            if x.sum() == 0:
                return 0
            idx = np.arange(1, n + 1)
            return (2 * np.sum(idx * x) - (n + 1) * np.sum(x)) / (n * np.sum(x))

        n_regions = X_all.shape[1]
        gini_scores = np.array([_gini(X_all[:, i]) for i in range(n_regions)])
        gini_threshold = np.mean(gini_scores) + 1.0 * np.std(gini_scores)
        gini_mask = gini_scores > gini_threshold
        logger.info(
            "  Gini filter: %d/%d cell-type-specific regions (%.1f%%)",
            gini_mask.sum(), n_regions, 100 * gini_mask.mean(),
        )

        # Z-score this cell type's values across cell types (per region)
        means = X_all[:, gini_mask].mean(axis=0)
        stds = X_all[:, gini_mask].std(axis=0)
        stds[stds == 0] = 1
        z_scores = (ct_values[gini_mask] - means) / stds

        # Take top K by z-score
        n_select = min(args.top_k_regions or 5000, len(z_scores))
        top_z_indices = np.argsort(z_scores)[::-1][:n_select]
        top_z_indices = np.sort(top_z_indices)

        # Map back to full adata indices
        gini_indices = np.where(gini_mask)[0]
        full_indices = gini_indices[top_z_indices]
        selected_regions = adata.var_names[full_indices]

        logger.info(
            "  Selected top %d regions by z-score (z range: %.2f - %.2f, count range: %.2f - %.2f)",
            n_select,
            z_scores[top_z_indices].min(),
            z_scores[top_z_indices].max(),
            ct_values[full_indices].min(),
            ct_values[full_indices].max(),
        )

        if False:  # Skip model-prediction subsampling — we already selected by z-score
            pass
            top_idx = np.argsort(ct_preds)[::-1][:args.max_regions]
            top_idx = np.sort(top_idx)
            selected_regions = selected_regions[top_idx]

        logger.info("  %d regions selected", len(selected_regions))

        if len(selected_regions) == 0:
            logger.warning("  No regions selected for %s, skipping", ct_name)
            continue

        # Compute contribution scores
        scores, one_hots = crested.tl.contribution_scores(
            input=list(selected_regions),
            target_idx=ct_idx,
            model=model,
            method=args.method,
            transpose=True,
            batch_size=args.batch_size,
            seed=args.seed,
            verbose=True,
        )
        logger.info("  Scores shape: %s, One-hots shape: %s", scores.shape, one_hots.shape)

        # Squeeze class dimension if present
        if scores.ndim == 4:
            scores = scores[:, 0, :, :]  # (N, 4, L)

        # Compute projected scores
        projected_scores = scores * one_hots

        # Save output
        ct_dir = os.path.join(args.output_dir, ct_name)
        os.makedirs(ct_dir, exist_ok=True)

        np.savez_compressed(
            os.path.join(ct_dir, "hypothetical_contribs.npz"),
            arr_0=scores,
        )
        np.savez_compressed(
            os.path.join(ct_dir, "projected_contribs.npz"),
            arr_0=projected_scores,
        )
        np.savez_compressed(
            os.path.join(ct_dir, "one_hot_seqs.npz"),
            arr_0=one_hots,
        )
        np.save(
            os.path.join(ct_dir, "region_names.npy"),
            np.array(list(selected_regions)),
        )

        logger.info("  Saved to %s", ct_dir)

    logger.info("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute per-cell-type contribution scores from a CREsted peak regression model.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # All cell types, top 10000 accessible regions each
  python compute_contribution_scores.py \\
      --training_dir ./crested/finetuned \\
      --bigwigs_dir ./bigwigs \\
      --regions_file ./consensus_peaks.bed \\
      --genome /path/to/mm10.fa \\
      --chrom_sizes /path/to/mm10.chrom.sizes \\
      --output_dir ./crested/contributions \\
      --top_k_regions 10000

  # Specific cell types
  python compute_contribution_scores.py \\
      --checkpoint ./crested/finetuned/checkpoints/best.keras \\
      --bigwigs_dir ./bigwigs \\
      --regions_file ./consensus_peaks.bed \\
      --genome /path/to/mm10.fa \\
      --chrom_sizes /path/to/mm10.chrom.sizes \\
      --output_dir ./crested/contributions \\
      --cell_types ATL ENDO PTS1 \\
      --top_k_regions 10000
""",
    )
    parser.add_argument("--training_dir", default=None, help="Training output dir")
    parser.add_argument("--checkpoint", default=None, help="Specific checkpoint path")
    parser.add_argument("--bigwigs_dir", required=True, help="Directory with per-cell-type BigWig files")
    parser.add_argument("--regions_file", required=True, help="Consensus peaks BED file")
    parser.add_argument("--genome", required=True, help="Path to genome FASTA")
    parser.add_argument("--chrom_sizes", required=True, help="Path to chromosome sizes file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--splits_file", default=None, help="JSON with chromosome splits")
    parser.add_argument("--cell_types", nargs="+", default=None, help="Specific cell types (default: all)")
    parser.add_argument("--method", default="expected_integrated_grad",
                        choices=["integrated_grad", "expected_integrated_grad", "saliency_map"])
    parser.add_argument("--top_k_regions", type=int, default=None,
                        help="Use top K most accessible regions per cell type")
    parser.add_argument("--accessibility_threshold", type=float, default=None,
                        help="Use regions above this accessibility threshold")
    parser.add_argument("--max_regions", type=int, default=None,
                        help="Max regions per cell type (subsample by prediction score)")
    parser.add_argument("--target", default="count", choices=["count", "mean"])
    parser.add_argument("--target_region_width", type=int, default=1000)
    parser.add_argument("--seq_len", type=int, default=2114)
    parser.add_argument("--top_k_percent", type=float, default=0.03)
    parser.add_argument("--val_chroms", nargs="+", default=["chr8", "chr10"])
    parser.add_argument("--test_chroms", nargs="+", default=["chr9", "chr18"])
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--backend", default="torch", choices=["torch", "tensorflow"])

    args = parser.parse_args()
    if args.training_dir is None and args.checkpoint is None:
        parser.error("Must specify either --training_dir or --checkpoint")
    main(args)
