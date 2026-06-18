"""Build a principled sign-aware motif catalog with MotifCompendium.

Replaces the gimme-cluster path with a CWM-aware pipeline:
  1. Load motif_compendium_initial.mc (1905 raw motifs across 22 CTs)
  2. Re-label `posneg` using actual avg_contrib sign — moves the ~200
     posneg_inverted motifs to the correct track instead of dropping them
  3. Apply standard QC filters (drop ~500 noisy/under-supported; KEEP both
     pos and neg tracks)
  4. Two-stage clustering (per MotifCompendium tutorial 3 best practice):
       stage 1: strict per-CT (sim=0.95, cluster_within=celltype)
       stage 2: across CTs WITHIN posneg (sim=0.85, cluster_on=stage1,
                cluster_within=posneg) — guarantees pos & neg never merge
  5. Export:
       - clustered_motifs.modisco.h5  (pos_patterns + neg_patterns groups,
         finemo-ready; this is the key advantage over gimme)
       - catalog.meme
       - cluster_to_name.tsv
       - catalog_walls.html  (CWM-based logos, so neg actually looks neg)
"""

from __future__ import annotations

import os
import logging
from pathlib import Path

import numpy as np
import pandas as pd

CUDA_ENV = "/cellar/users/aklie/opt/miniconda3/envs/motifcompendium-gpu"
os.environ["CUDA_PATH"] = CUDA_ENV
os.environ["CPLUS_INCLUDE_PATH"] = f"{CUDA_ENV}/include"
os.environ["CUPY_NVRTC_OPTIONS"] = f"-I{CUDA_ENV}/include"

import cupy  # noqa
_ = cupy.ones((10,), dtype=cupy.float16); del _

import MotifCompendium
import MotifCompendium.utils.analysis as utils_analysis
import MotifCompendium.utils.plotting as utils_plotting
from MotifCompendium.utils.plotting import LogoPlottingInput
from MotifCompendium.utils.visualization import motif_collection_html


BASE = Path("/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome")
PATH_MC_INIT  = BASE / "results/4_multi_task_model/crested/motif_compendium_all_peaks/motif_compendium_initial.mc"
PATH_OUT      = BASE / "results/4_multi_task_model/crested/motifs_compendium"
PATH_OUT.mkdir(parents=True, exist_ok=True)

# Clustering threshold — single-stage cross-CT clustering within posneg.
# (Two-stage strict-then-merge needs cluster_within_on=tuple in this MC version;
# single-stage with cluster_within="posneg" is simpler and still guarantees
# pos and neg never merge.)
CLUSTER_THRESHOLD = 0.95


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(str(PATH_OUT / "build_catalog.log")),
        ],
        force=True,
    )
    return logging.getLogger()


def main():
    log = setup_logging()

    # --- 1. Load saved initial MC ---
    log.info(f"Loading initial MC: {PATH_MC_INIT}")
    mc = MotifCompendium.load(str(PATH_MC_INIT))
    # Fix celltype parsing — the init step split on "." which left the "-pos"/"-neg"
    # suffix attached (e.g. "DE-pos" instead of "DE"). Re-derive cleanly here.
    mc["celltype"] = mc["name"].str.split("-").str[0]
    log.info(f"  loaded {len(mc)} raw motifs across {mc['celltype'].nunique()} CTs (re-parsed)")

    # --- 2. Re-label posneg by avg_contrib sign (correct the inverted ones) ---
    orig_posneg = mc.metadata["posneg"].copy()
    new_posneg = np.where(mc["avg_contrib"] >= 0, "pos", "neg")
    mc.metadata["posneg_orig"] = orig_posneg
    mc.metadata["posneg"] = new_posneg
    changed = (orig_posneg != new_posneg).sum()
    log.info(f"Re-labeled posneg by avg_contrib sign: {changed} motifs moved")
    log.info(f"  new posneg counts: {pd.Series(new_posneg).value_counts().to_dict()}")

    # --- 3. Apply QC filters (use existing flag_remove from init step) ---
    if "flag_remove" not in mc.metadata.columns:
        log.warning("flag_remove not present — initial QC may have been skipped")
    else:
        before = len(mc)
        mc = mc[~mc.metadata["flag_remove"]]
        log.info(f"Dropped {before - len(mc)} motifs flagged by initial QC filters")
        log.info(f"  retained: {len(mc)}  (posneg counts: {mc.metadata['posneg'].value_counts().to_dict()})")

    # --- 4. Single-stage clustering within posneg (pos and neg never merge) ---
    log.info(f"Clustering: sim>={CLUSTER_THRESHOLD}, within=posneg")
    mc.cluster(
        similarity_threshold=CLUSTER_THRESHOLD,
        save_name="final_clusters",
        cluster_within="posneg",
    )
    n_final = mc["final_clusters"].nunique()
    pos_clusters = mc.metadata[mc["posneg"] == "pos"]["final_clusters"].nunique()
    neg_clusters = mc.metadata[mc["posneg"] == "neg"]["final_clusters"].nunique()
    log.info(f"  -> {n_final} final clusters ({pos_clusters} pos + {neg_clusters} neg)")

    # Save the post-filter post-cluster MC
    mc.save(str(PATH_OUT / "motif_compendium_clustered.mc"))
    mc.metadata.to_csv(PATH_OUT / "motif_compendium_clustered_metadata.tsv", sep="\t", index=False)

    # --- 5. Export to MEME and synthetic modisco h5 (for finemo) ---
    log.info("Exporting catalog formats...")
    # MEME — uses cluster representative for each cluster
    utils_analysis.export_compendium_clustered_modisco(
        mc,
        "final_clusters",
        str(PATH_OUT / "clustered_motifs.modisco.h5"),
        export_subpatterns=True,
    )
    log.info(f"  wrote clustered_motifs.modisco.h5")

    utils_analysis.export_compendium_meme(
        mc, "final_clusters", str(PATH_OUT / "catalog.meme")
    )
    log.info(f"  wrote catalog.meme")

    # cluster → list-of-member-motifs mapping
    cluster_members = (
        mc.metadata.groupby("final_clusters")["name"]
        .apply(lambda s: ";".join(s.astype(str)))
        .reset_index()
    )
    cluster_members.to_csv(PATH_OUT / "cluster_to_members.tsv", sep="\t", index=False)

    # --- 6. Render walls via MC's built-in renderer ---
    # Split MC by posneg → two HTMLs, each grouped by cluster with average_motif on.
    # MC's renderer chooses the appropriate representation (CWM-based, sign-aware
    # when the underlying motif data is signed).
    log.info("Rendering walls via mc.motif_collection_html (one per track, grouped by cluster)")
    mc_pos = mc[mc["posneg"] == "pos"]
    mc_neg = mc[mc["posneg"] == "neg"]
    mc_pos.motif_collection_html(
        str(PATH_OUT / "catalog_walls_pos.html"),
        group_by="final_clusters",
    )
    mc_neg.motif_collection_html(
        str(PATH_OUT / "catalog_walls_neg.html"),
        group_by="final_clusters",
    )
    log.info(f"  wrote catalog_walls_pos.html and catalog_walls_neg.html")

    # --- 7. Headline summary ---
    log.info("=" * 60)
    log.info("PRINCIPLED MC CATALOG — HEADLINE SUMMARY")
    log.info("=" * 60)
    log.info(f"Final motif clusters:           {n_final}")
    log.info(f"  pos:                          {pos_clusters}")
    log.info(f"  neg:                          {neg_clusters}")
    log.info(f"For comparison, gimme catalog:  149 (62 pos + 87 neg)")
    log.info(f"Output: {PATH_OUT}/")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
