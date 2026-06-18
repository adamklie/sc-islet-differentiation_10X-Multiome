"""Finalize the MC catalog — runs the steps that failed in the main build:
  - MEME export (correct arg order: mc, save_loc, name_col)
  - Per-track walls via mc.motif_collection_html (group_by=cluster, average_motif)
  - Build name_map.tsv (synthetic-h5 idx → human-readable cluster annotation)
    so finemo can produce hits with readable names.

Reads the saved motif_compendium_clustered.mc from the prior build run; does
not re-do clustering.
"""

from __future__ import annotations

import os
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import h5py

CUDA_ENV = "/cellar/users/aklie/opt/miniconda3/envs/motifcompendium-gpu"
os.environ["CUDA_PATH"] = CUDA_ENV
os.environ["CPLUS_INCLUDE_PATH"] = f"{CUDA_ENV}/include"
os.environ["CUPY_NVRTC_OPTIONS"] = f"-I{CUDA_ENV}/include"

import cupy  # noqa
_ = cupy.ones((10,), dtype=cupy.float16); del _

import MotifCompendium
import MotifCompendium.utils.analysis as utils_analysis


BASE = Path("/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome")
PATH_OUT = BASE / "results/4_multi_task_model/crested/motifs_compendium"
PATH_MC = PATH_OUT / "motif_compendium_clustered.mc"
H5 = PATH_OUT / "clustered_motifs.modisco.h5"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(),
              logging.FileHandler(str(PATH_OUT / "finalize_catalog.log"))],
    force=True,
)
log = logging.getLogger()


def main():
    log.info(f"Loading clustered MC from {PATH_MC}")
    mc = MotifCompendium.load(str(PATH_MC))
    log.info(f"  {len(mc)} motifs in {mc['final_clusters'].nunique()} clusters")
    log.info(f"  pos clusters: {mc.metadata[mc['posneg']=='pos']['final_clusters'].nunique()}")
    log.info(f"  neg clusters: {mc.metadata[mc['posneg']=='neg']['final_clusters'].nunique()}")

    # 1. MEME export (correct arg order: mc, save_loc, name_col)
    meme_path = PATH_OUT / "catalog.meme"
    log.info(f"Exporting MEME → {meme_path}")
    utils_analysis.export_compendium_meme(mc, str(meme_path), "final_clusters")

    # 2. Walls per track
    log.info("Rendering walls — POS track")
    mc_pos = mc[mc["posneg"] == "pos"]
    mc_pos.motif_collection_html(
        str(PATH_OUT / "catalog_walls_pos.html"),
        "final_clusters",
    )
    log.info("Rendering walls — NEG track")
    mc_neg = mc[mc["posneg"] == "neg"]
    mc_neg.motif_collection_html(
        str(PATH_OUT / "catalog_walls_neg.html"),
        "final_clusters",
    )

    # 3. Build name_map.tsv for finemo
    # The synthetic h5 keys are integers (cluster ids). For finemo we map
    # {posneg}_patterns.pattern_<i> → human-readable.
    log.info("Building name_map.tsv for finemo")
    # Read keys from the synthetic h5 to get the actual indexing
    with h5py.File(H5, "r") as f:
        pos_keys = sorted(f["pos_patterns"].keys(), key=lambda x: int(x))
        neg_keys = sorted(f["neg_patterns"].keys(), key=lambda x: int(x))

    # Pick a representative name per cluster (highest num_seqlets motif's annotation)
    def cluster_label(grp):
        rep = grp.loc[grp["num_seqlets"].idxmax()]
        ann = rep.get("annotation_database_name0", "") or ""
        ct = rep.get("celltype", "")
        cid = rep["final_clusters"]
        # name format: cluster_<id>_<annotation>_<celltype_of_rep>
        parts = [f"cluster_{cid}"]
        if ann: parts.append(str(ann)[:30].replace(" ", "_"))
        if ct: parts.append(f"src={ct}")
        return "_".join(parts)

    cluster_to_name = mc.metadata.groupby("final_clusters").apply(cluster_label).to_dict()

    rows = []
    for k in pos_keys:
        cid = int(k)
        nm = cluster_to_name.get(cid, f"cluster_{cid}")
        rows.append((f"pos_patterns.pattern_{k}", f"pos_{nm}"))
    for k in neg_keys:
        cid = int(k)
        nm = cluster_to_name.get(cid, f"cluster_{cid}")
        rows.append((f"neg_patterns.pattern_{k}", f"neg_{nm}"))

    name_map_df = pd.DataFrame(rows, columns=["synthetic_key", "name"])
    name_map_df.to_csv(PATH_OUT / "clustered_motifs_name_map.tsv", sep="\t", index=False, header=False)

    # _unique variant: resolve duplicates with positional suffix
    dups = name_map_df["name"].duplicated(keep=False)
    if dups.any():
        name_map_df.loc[dups, "name"] = (
            name_map_df.loc[dups, "name"]
            + "_#"
            + name_map_df.loc[dups].groupby("name").cumcount().add(1).astype(str)
        )
        log.info(f"  resolved {int(dups.sum())} duplicate cluster names → _unique map")
    name_map_df.to_csv(PATH_OUT / "clustered_motifs_name_map_unique.tsv", sep="\t", index=False, header=False)

    log.info("DONE — full MC catalog ready for hit calling")
    log.info(f"  modisco h5:        {H5}")
    log.info(f"  MEME:              {PATH_OUT / 'catalog.meme'}")
    log.info(f"  name_map_unique:   {PATH_OUT / 'clustered_motifs_name_map_unique.tsv'}")
    log.info(f"  walls (pos/neg):   {PATH_OUT / 'catalog_walls_pos.html'}, _neg.html")


if __name__ == "__main__":
    main()
