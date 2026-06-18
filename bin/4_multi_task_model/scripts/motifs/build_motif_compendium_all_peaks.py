"""Build an initial MotifCompendium from the all-peak MoDISco runs.

Purpose: parallel sign-aware clustering of the all-peak motifs, alongside the
gimme cluster output that lives at ``motifs_all_peaks/``. The gimme path
operates on PPMs (sign-agnostic), so neg-vs-pos can't be distinguished from
the values — only from the name prefix we put there. MotifCompendium is
purpose-built for this: ``build_from_modisco`` populates ``posneg`` and
``avg_contrib`` directly from the modisco h5, and ``posneg_inverted`` is a
built-in QC metric for the "modisco said neg but the contribution sign
actually says pos" mis-categorization case.

Mirrors ``stimulated_sc-islets/manuscripts/ChromBetaNet/bin/7_motif2cRE_linking/
1_create_initial_motif_compedium.ipynb`` but with celltype as the only metadata
axis (no condition / no model head split here).

Outputs land at ``results/4_multi_task_model/crested/motif_compendium_all_peaks/``:
    motif_compendium_initial.mc
    motif_compendium_initial_metadata.tsv
    motif_compendium_flagged_by_filters.tsv
    qc/                                 (entropy + similarity plots)
    examples/                           (per-archetype example logos)
"""

from __future__ import annotations

import os
import glob
import logging
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# CUDA / cupy setup — same pattern as the starter notebook.
CUDA_ENV = "/cellar/users/aklie/opt/miniconda3/envs/motifcompendium-gpu"
os.environ["CUDA_PATH"] = CUDA_ENV
os.environ["CPLUS_INCLUDE_PATH"] = f"{CUDA_ENV}/include"
os.environ["CUPY_NVRTC_OPTIONS"] = f"-I{CUDA_ENV}/include"

import cupy  # noqa: E402  forces the cupy header initialization
_warmup = cupy.ones((10,), dtype=cupy.float16)
del _warmup

import MotifCompendium  # noqa: E402
import MotifCompendium.utils.analysis as utils_analysis  # noqa: E402
import MotifCompendium.utils.plotting as utils_plotting  # noqa: E402


PATH_REF_MC_MEME = "/cellar/users/aklie/opt/MotifCompendium/pipeline/data/MotifCompendium-Database-Human.meme.txt"
PATH_REF_DATABASE_MEME = "/cellar/users/aklie/opt/MotifCompendium/pipeline/data/HUMAN-JASPAR2024-HOCOMOCOv13.meme.txt"

BASE = Path("/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome")
MODISCO_GLOB = str(BASE / "results/4_multi_task_model/crested/modisco_all_peaks/*/[A-z]*.modisco.h5")
PATH_OUT = BASE / "results/4_multi_task_model/crested/motif_compendium_all_peaks"
PATH_OUT.mkdir(parents=True, exist_ok=True)
(PATH_OUT / "qc").mkdir(exist_ok=True)
(PATH_OUT / "examples").mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(str(PATH_OUT / "motif_compendium_step1.log")),
    ],
    force=True,
)
log = logging.getLogger()


def main():
    # 1. Load all 22 MoDISco outputs into a dict
    path_modiscos = sorted(glob.glob(MODISCO_GLOB))
    log.info(f"Found {len(path_modiscos)} MoDISco files")
    if len(path_modiscos) != 22:
        log.warning(f"Expected 22 modisco h5 files, got {len(path_modiscos)}")
    modisco_dict = {os.path.basename(p).split(".modisco.h5")[0]: p for p in path_modiscos}

    # 2. Build the compendium. modisco_region_width defaults to 400 (matches our -w 400).
    log.info("Building MotifCompendium from modisco h5 files (this loads contrib_scores + PPMs)...")
    mc = MotifCompendium.build_from_modisco(modisco_dict)
    log.info(f"MotifCompendium built with {len(mc)} motifs total")

    # Metadata: parse celltype from name (single axis here — no condition / head)
    mc["celltype"] = mc["name"].str.split(".").str[0]
    log.info(f"  cell types: {mc['celltype'].nunique()} unique ({sorted(mc['celltype'].unique())})")

    # Sanity check: build_from_modisco should auto-populate posneg + avg_contrib
    cols = list(mc.metadata.columns)
    log.info(f"  metadata columns ({len(cols)}): {cols}")
    if "posneg" in cols:
        log.info(f"  posneg breakdown: {mc.metadata['posneg'].value_counts().to_dict()}")
    if "avg_contrib" in cols:
        log.info(f"  avg_contrib: min={mc['avg_contrib'].min():+.3f}  median={mc['avg_contrib'].median():+.3f}  max={mc['avg_contrib'].max():+.3f}")

    # 3. Annotate against MotifCompendium reference
    log.info(f"Annotating against {PATH_REF_MC_MEME}...")
    utils_analysis.assign_label_from_pfms(
        mc=mc,
        pfm_file=PATH_REF_MC_MEME,
        save_col_prefix="annotation_database",
        min_score=0.7,
        save_images=False,
        logo_trimming=0,
    )

    # 4. Calculate QC metrics — including the diagnostic the user asked about
    log.info("Calculating QC metrics (including posneg_inverted — the diagnostic flag)...")
    metric_list = [
        "motif_entropy",
        "weighted_base_entropy",
        "weighted_position_entropy",
        "posbase_entropy_score",
        "copair_entropy_score",
        "copair_composition",
        "dinuc_entropy_score",
        "dinuc_composition",
        "dinuc_score",
        "posneg_inverted",
        "truncated",
    ]
    utils_analysis.calculate_filters(mc=mc, metric_list=metric_list)

    # 5. Apply standard filter rules from the starter
    @dataclass
    class _FilterArgs:
        name: str
        metric: str
        operation: str
        threshold: Union[float, bool]

        def to_dict(self):
            return asdict(self)

    motif_filters = [
        _FilterArgs("1_singlepeak",       "motif_entropy",            "<",  0.35),
        _FilterArgs("2_noisemix",         "motif_entropy",            ">",  0.75),
        _FilterArgs("3_noisypeaks",       "weighted_base_entropy",    ">",  0.6),
        _FilterArgs("4_broadsingle_1",    "weighted_position_entropy", ">", 0.85),
        _FilterArgs("4_broadsingle_2",    "posbase_entropy_score",    ">",  0.5),
        _FilterArgs("5_gcbias_1",         "copair_entropy_score",     ">",  0.35),
        _FilterArgs("5_gcbias_2",         "copair_composition",       ">",  0.45),
        _FilterArgs("6_dinucrepeat_1",    "dinuc_entropy_score",      ">",  0.5),
        _FilterArgs("6_dinucrepeat_2",    "dinuc_composition",        ">",  0.875),
        _FilterArgs("6_dinucrepeat_3",    "dinuc_score",              ">",  0.45),
        _FilterArgs("7_minseqlets",       "num_seqlets",              "<",  100),
        _FilterArgs("8_posneg_inverted",  "posneg_inverted",          "==", True),
        _FilterArgs("9_truncated",        "truncated",                "==", True),
        _FilterArgs("10_no_annotation",   "annotation_database_name0","==", ""),
        # Note: starter has a "11_neg_pattern" filter that flags posneg==neg. We are
        # explicitly KEEPING neg patterns here so we can do the sign-aware analysis.
    ]

    def apply_filter_threshold(mc, flag_col, metric, operation, threshold):
        col = mc[metric] if metric in mc.metadata.columns else mc.metadata.get(metric)
        if col is None:
            log.warning(f"  metric '{metric}' not found — skipping")
            return
        if operation == "<":
            mc.metadata[flag_col] = mc.metadata[flag_col] | (col < threshold)
        elif operation == "<=":
            mc.metadata[flag_col] = mc.metadata[flag_col] | (col <= threshold)
        elif operation == ">":
            mc.metadata[flag_col] = mc.metadata[flag_col] | (col > threshold)
        elif operation == ">=":
            mc.metadata[flag_col] = mc.metadata[flag_col] | (col >= threshold)
        elif operation == "==":
            mc.metadata[flag_col] = mc.metadata[flag_col] | (col == threshold)

    flag_col = "flag_remove"
    mc.metadata[flag_col] = False
    mc.metadata["flagged_by"] = ""

    for f in motif_filters:
        before = int(mc.metadata[flag_col].sum())
        apply_filter_threshold(mc, flag_col, f.metric, f.operation, f.threshold)
        after = int(mc.metadata[flag_col].sum())
        newly_flagged = after - before
        newly_flagged_idx = mc.metadata.index[mc.metadata[flag_col] & (mc.metadata["flagged_by"] == "")]
        mc.metadata.loc[newly_flagged_idx, "flagged_by"] = f.name
        log.info(f"  {f.name}: flagged {newly_flagged} additional ({after} cumulative)")

    flagged_motifs_path = PATH_OUT / "motif_compendium_flagged_by_filters.tsv"
    mc.metadata.loc[mc.metadata[flag_col], ["name", "flagged_by"]].to_csv(
        flagged_motifs_path, sep="\t", index=False
    )
    log.info(f"Wrote flagged motif list with reasons → {flagged_motifs_path}")

    # 6. QC plots
    log.info("Writing QC plots...")
    plt.figure(figsize=(6, 4))
    sns.histplot(mc.metadata["num_seqlets"], bins=50)
    plt.title("Distribution of num_seqlets")
    plt.savefig(PATH_OUT / "qc" / "qc_num_seqlets.png", bbox_inches="tight"); plt.close()

    if "avg_contrib" in mc.metadata.columns:
        plt.figure(figsize=(6, 4))
        sns.histplot(mc.metadata["avg_contrib"], bins=50)
        plt.axvline(0, color="grey", linestyle="--", lw=0.7)
        plt.title("Distribution of avg_contrib (signed; should be bimodal pos vs neg)")
        plt.savefig(PATH_OUT / "qc" / "qc_avg_contrib.png", bbox_inches="tight"); plt.close()

    if "posneg" in mc.metadata.columns:
        plt.figure(figsize=(6, 4))
        sns.countplot(data=mc.metadata, x="posneg", hue="posneg_inverted")
        plt.title("posneg label vs posneg_inverted flag")
        plt.savefig(PATH_OUT / "qc" / "qc_posneg_inverted.png", bbox_inches="tight"); plt.close()

    # 7. Save objects
    mc.save(str(PATH_OUT / "motif_compendium_initial.mc"))
    mc.metadata.to_csv(PATH_OUT / "motif_compendium_initial_metadata.tsv", sep="\t", index=False)
    log.info(f"Saved compendium + metadata to {PATH_OUT}/")

    # 8. Headline diagnostics — what the user asked about
    log.info("=" * 60)
    log.info("HEADLINE DIAGNOSTICS")
    log.info("=" * 60)
    log.info(f"Total motifs in MC:                  {len(mc):>5}")
    if "posneg" in mc.metadata.columns:
        for v, n in mc.metadata["posneg"].value_counts().items():
            log.info(f"  posneg = {v!s:<5}                  {n:>5}")
    if "posneg_inverted" in mc.metadata.columns:
        n_inv = int(mc.metadata["posneg_inverted"].sum())
        log.info(f"posneg_inverted (modisco label mismatched with actual sign): {n_inv} / {len(mc)}")
    n_flag = int(mc.metadata[flag_col].sum())
    log.info(f"Total flagged for removal by all rules: {n_flag} / {len(mc)}")
    log.info(f"Would keep {len(mc) - n_flag} motifs after standard QC (still includes neg)")
    log.info(f"Compare: gimme cluster (separate-tracks) catalog has 149 motifs (62 pos + 87 neg)")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
