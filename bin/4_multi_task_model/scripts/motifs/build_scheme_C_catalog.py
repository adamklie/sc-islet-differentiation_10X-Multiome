"""Build the Scheme C catalog from the QC'd MC.

The Phase A.5 QC saved qc/post_qc_clustered.mc with three cluster columns:
  final_clusters  — Scheme A (sim=0.95 within=posneg, 325 clusters, 52% singletons)
  clusters_B      — Scheme B (sim=0.97 within=posneg, 407 clusters, worse)
  clusters_C      — Scheme C (2-stage tutorial-3: strict per-CT then merge
                     across CTs within posneg; 150 clusters, 89 pos + 61 neg)

This script commits to Scheme C and overwrites the canonical paths
(clustered_motifs.modisco.h5, name_map_unique.tsv, catalog.meme, walls) used
by motif_hits_mc_catalog.sh. Scheme A artifacts get renamed with _schemeA
suffix so we don't lose them.
"""

from __future__ import annotations

import os
import shutil
import logging
from pathlib import Path

import h5py
import pandas as pd

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
PATH_MC = PATH_OUT / "qc" / "post_qc_clustered.mc"

CLUSTER_COL = "clusters_C"   # commit to Scheme C

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(),
              logging.FileHandler(str(PATH_OUT / "build_scheme_C.log"))],
    force=True,
)
log = logging.getLogger()


def backup_schemeA():
    """Rename existing canonical Scheme A artifacts so we don't clobber them."""
    pairs = [
        ("clustered_motifs.modisco.h5",           "clustered_motifs_schemeA.modisco.h5"),
        ("clustered_motifs_name_map.tsv",         "clustered_motifs_name_map_schemeA.tsv"),
        ("clustered_motifs_name_map_unique.tsv",  "clustered_motifs_name_map_schemeA_unique.tsv"),
        ("catalog.meme",                          "catalog_schemeA.meme"),
        ("catalog_walls_pos.html",                "catalog_walls_schemeA_pos.html"),
        ("catalog_walls_neg.html",                "catalog_walls_schemeA_neg.html"),
    ]
    for src, dst in pairs:
        src_p = PATH_OUT / src
        dst_p = PATH_OUT / dst
        if src_p.is_file():
            if dst_p.is_file():
                log.info(f"  {dst} already exists — leaving as is, removing live {src}")
                src_p.unlink()
            else:
                shutil.move(str(src_p), str(dst_p))
                log.info(f"  renamed {src} -> {dst}")
        else:
            log.info(f"  no {src} to backup")


def main():
    log.info(f"Loading post-QC MC from {PATH_MC}")
    mc = MotifCompendium.load(str(PATH_MC))
    log.info(f"  loaded {len(mc)} motifs")
    if CLUSTER_COL not in mc.metadata.columns:
        raise RuntimeError(f"{CLUSTER_COL} not present in saved MC")
    n = mc[CLUSTER_COL].nunique()
    n_pos = mc.metadata[mc["posneg"] == "pos"][CLUSTER_COL].nunique()
    n_neg = mc.metadata[mc["posneg"] == "neg"][CLUSTER_COL].nunique()
    log.info(f"Scheme C: {n} clusters ({n_pos} pos + {n_neg} neg)")

    log.info("Backing up Scheme A canonical artifacts")
    backup_schemeA()

    # --- export modisco h5 (overwrites canonical path) ---
    h5_out = PATH_OUT / "clustered_motifs.modisco.h5"
    log.info(f"Exporting Scheme C → {h5_out}")
    utils_analysis.export_compendium_clustered_modisco(
        mc, CLUSTER_COL, str(h5_out), export_subpatterns=True
    )

    # --- export MEME (canonical path) ---
    meme_out = PATH_OUT / "catalog.meme"
    log.info(f"Exporting MEME → {meme_out}")
    utils_analysis.export_compendium_meme(mc, str(meme_out), CLUSTER_COL)

    # --- walls per track via mc.motif_collection_html, grouped by clusters_C ---
    log.info("Rendering walls per track")
    mc_pos = mc[mc["posneg"] == "pos"]
    mc_neg = mc[mc["posneg"] == "neg"]
    mc_pos.motif_collection_html(str(PATH_OUT / "catalog_walls_pos.html"), CLUSTER_COL)
    mc_neg.motif_collection_html(str(PATH_OUT / "catalog_walls_neg.html"), CLUSTER_COL)

    # --- name_map.tsv (synthetic h5 key -> readable name) ---
    log.info("Building name_map_unique.tsv from new h5 + metadata")
    def cluster_label(grp):
        rep = grp.loc[grp["num_seqlets"].idxmax()]
        ann = str(rep.get("annotation_database_name0", "") or "")
        ct = str(rep.get("celltype", ""))
        cid = rep[CLUSTER_COL]
        parts = [f"cluster_{cid}"]
        if ann: parts.append(ann[:30].replace(" ", "_").replace("(", "").replace(")", ""))
        if ct: parts.append(f"src{ct}")
        return "_".join(parts)
    cluster_to_name = mc.metadata.groupby(CLUSTER_COL).apply(cluster_label).to_dict()

    with h5py.File(h5_out, "r") as f:
        pos_keys = sorted(f["pos_patterns"].keys(), key=lambda x: int(x)) if "pos_patterns" in f else []
        neg_keys = sorted(f["neg_patterns"].keys(), key=lambda x: int(x)) if "neg_patterns" in f else []
    log.info(f"  h5 has {len(pos_keys)} pos_patterns + {len(neg_keys)} neg_patterns groups")

    rows = []
    for k in pos_keys:
        cid = int(k)
        nm = cluster_to_name.get(cid, f"cluster_{cid}")
        rows.append((f"pos_patterns.pattern_{k}", f"pos_{nm}"))
    for k in neg_keys:
        cid = int(k)
        nm = cluster_to_name.get(cid, f"cluster_{cid}")
        rows.append((f"neg_patterns.pattern_{k}", f"neg_{nm}"))

    nm = pd.DataFrame(rows, columns=["key", "name"])
    nm.to_csv(PATH_OUT / "clustered_motifs_name_map.tsv", sep="\t", index=False, header=False)
    dup = nm["name"].duplicated(keep=False)
    if dup.any():
        nm.loc[dup, "name"] = (
            nm.loc[dup, "name"] + "_#"
            + nm.loc[dup].groupby("name").cumcount().add(1).astype(str)
        )
        log.info(f"  resolved {int(dup.sum())} duplicate names → _unique")
    nm.to_csv(PATH_OUT / "clustered_motifs_name_map_unique.tsv", sep="\t", index=False, header=False)

    log.info("Scheme C catalog committed to canonical paths:")
    log.info(f"  modisco h5:   {h5_out}")
    log.info(f"  name map:     {PATH_OUT / 'clustered_motifs_name_map_unique.tsv'}")
    log.info(f"  MEME:         {meme_out}")
    log.info(f"  walls:        catalog_walls_pos.html / _neg.html")
    log.info("Scheme A files preserved with _schemeA suffix")


if __name__ == "__main__":
    main()
