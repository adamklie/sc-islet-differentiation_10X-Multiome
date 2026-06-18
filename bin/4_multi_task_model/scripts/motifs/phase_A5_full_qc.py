"""Comprehensive Phase A.5 QC for the MC catalog.

Vets the catalog along the lines of stim_sc-islets/manuscripts/ChromBetaNet/
bin/7_motif2cRE_linking notebooks 1+2 + MotifCompendium tutorials 3, 4, 6, 7
before any hit calling proceeds.

Sections:
  1. Catalog snapshot
  2. Cluster-quality diagnostics for current scheme (A: sim=0.95, within=posneg)
       - judge_clustering, heatmap, similarity distribution, per-cluster sizes
  3. Alternative clustering schemes
       - Scheme B: sim=0.97 single-stage within posneg (stricter)
       - Scheme C: two-stage strict-then-merge (tutorial 3 recommended)
                     stage 1: sim=0.95 within=celltype  -> strict_clusters
                     stage 2: sim=0.85 cluster_within_on=("posneg", "strict_clusters")
       - Compare scheme cluster counts (pos vs neg) side by side
  4. Per-CT splits: save subset MC + judge_clustering + similarity heatmap per CT
  5. Annotation × celltype crosstab + heatmap
  6. Motif entropy by annotation status (boxplot)
  7. Archetype example logos (8 archetypes)
  8. Save mc_filtered + mc_removed objects separately

All outputs go to motifs_compendium/qc/.
"""

from __future__ import annotations

import os
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

CUDA_ENV = "/cellar/users/aklie/opt/miniconda3/envs/motifcompendium-gpu"
os.environ["CUDA_PATH"] = CUDA_ENV
os.environ["CPLUS_INCLUDE_PATH"] = f"{CUDA_ENV}/include"
os.environ["CUPY_NVRTC_OPTIONS"] = f"-I{CUDA_ENV}/include"

import cupy  # noqa
_ = cupy.ones((10,), dtype=cupy.float16); del _

import MotifCompendium
import MotifCompendium.utils.analysis as utils_analysis
import MotifCompendium.utils.plotting as utils_plotting


BASE = Path("/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome")
PATH_MC      = BASE / "results/4_multi_task_model/crested/motifs_compendium/motif_compendium_clustered.mc"
PATH_MC_INIT = BASE / "results/4_multi_task_model/crested/motif_compendium_all_peaks/motif_compendium_initial.mc"
PATH_OUT     = BASE / "results/4_multi_task_model/crested/motifs_compendium"
QC = PATH_OUT / "qc"
QC.mkdir(parents=True, exist_ok=True)
(QC / "schemes").mkdir(exist_ok=True)
(QC / "per_ct").mkdir(exist_ok=True)
(QC / "archetypes").mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(),
              logging.FileHandler(str(QC / "phase_A5_qc.log"))],
    force=True,
)
log = logging.getLogger()


def safe(section_name):
    """Decorator: log section, swallow exceptions so one bad step doesn't kill the rest."""
    def deco(fn):
        def wrapper(*args, **kwargs):
            log.info(f"=== {section_name} ===")
            try:
                return fn(*args, **kwargs)
            except Exception as e:
                log.exception(f"  FAILED: {section_name}: {e}")
        return wrapper
    return deco


@safe("Section 1: catalog snapshot")
def section_1_snapshot(mc):
    log.info(f"motifs: {len(mc)}")
    log.info(f"clusters (current 'final_clusters' scheme): {mc['final_clusters'].nunique()}")
    log.info(f"posneg counts: {mc.metadata['posneg'].value_counts().to_dict()}")
    log.info(f"per-posneg cluster count: pos={mc.metadata[mc['posneg']=='pos']['final_clusters'].nunique()}  neg={mc.metadata[mc['posneg']=='neg']['final_clusters'].nunique()}")
    log.info(f"celltypes: {sorted(mc['celltype'].unique())}")


@safe("Section 2: cluster quality diagnostics (Scheme A current)")
def section_2_diagnostics(mc):
    # judge_clustering — the canonical clustering quality plot
    log.info("judge_clustering on final_clusters")
    utils_analysis.judge_clustering(
        mc, "final_clusters",
        save_loc=str(QC / "schemes" / "schemeA_judge_clustering.png"),
    )

    # Similarity heatmap sorted by cluster
    log.info("similarity heatmap sorted by final_clusters")
    mc_sorted = mc.sort(by="final_clusters")
    mc_sorted.heatmap(save_loc=str(QC / "schemes" / "schemeA_heatmap_sorted.png"))
    mc_sorted.heatmap(similarity_threshold=0.95, save_loc=str(QC / "schemes" / "schemeA_heatmap_thresh095.png"))

    # Similarity distribution (global)
    log.info("global similarity distribution")
    utils_analysis.plot_similarity_distribution(
        mc, str(QC / "schemes" / "schemeA_similarity_distribution.html"), tolerance=0.001
    )

    # Per-cluster sizes
    log.info("per-cluster member-count distribution")
    sizes = mc.metadata.groupby("final_clusters").size().sort_values(ascending=False)
    sizes_per_posneg = mc.metadata.groupby(["posneg", "final_clusters"]).size()

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(sizes.values, bins=40)
    axes[0].set_title("Cluster sizes (Scheme A)")
    axes[0].set_xlabel("members per cluster"); axes[0].set_ylabel("n clusters")
    pos_sizes = sizes_per_posneg.xs("pos").values if "pos" in sizes_per_posneg.index.get_level_values(0) else []
    neg_sizes = sizes_per_posneg.xs("neg").values if "neg" in sizes_per_posneg.index.get_level_values(0) else []
    axes[1].hist([pos_sizes, neg_sizes], bins=30, label=["pos", "neg"], color=["green", "red"], alpha=0.6, stacked=False)
    axes[1].set_title("By posneg"); axes[1].set_xlabel("members per cluster"); axes[1].legend()
    fig.tight_layout(); fig.savefig(QC / "schemes" / "schemeA_cluster_sizes.png", dpi=110, bbox_inches="tight"); plt.close(fig)

    log.info(f"  cluster size summary: min={sizes.min()}, median={int(sizes.median())}, max={sizes.max()}")
    log.info(f"  singletons (size=1): {(sizes == 1).sum()} of {len(sizes)}")


@safe("Section 3: alternative clustering schemes")
def section_3_alt_schemes(mc):
    results = []
    # current = Scheme A
    A_n = mc["final_clusters"].nunique()
    A_pos = mc.metadata[mc["posneg"] == "pos"]["final_clusters"].nunique()
    A_neg = mc.metadata[mc["posneg"] == "neg"]["final_clusters"].nunique()
    results.append(("A: sim=0.95 within=posneg (current)", A_n, A_pos, A_neg))

    # Scheme B: stricter
    log.info("Scheme B: sim=0.97 within=posneg")
    mc.cluster(
        similarity_threshold=0.97,
        save_name="clusters_B",
        cluster_within="posneg",
    )
    B_n = mc["clusters_B"].nunique()
    B_pos = mc.metadata[mc["posneg"] == "pos"]["clusters_B"].nunique()
    B_neg = mc.metadata[mc["posneg"] == "neg"]["clusters_B"].nunique()
    results.append(("B: sim=0.97 within=posneg", B_n, B_pos, B_neg))
    utils_analysis.judge_clustering(mc, "clusters_B", save_loc=str(QC / "schemes" / "schemeB_judge_clustering.png"))

    # Scheme C: two-stage per tutorial 3 — strict per-CT, then merge within posneg
    log.info("Scheme C stage 1: sim=0.95 within=celltype  -> strict_clusters")
    mc.cluster(
        similarity_threshold=0.95,
        save_name="strict_clusters",
        cluster_within="celltype",
    )
    log.info(f"  strict_clusters: {mc['strict_clusters'].nunique()}")
    log.info("Scheme C stage 2: sim=0.85 cluster_within_on=('posneg', 'strict_clusters')")
    mc.cluster(
        similarity_threshold=0.85,
        save_name="clusters_C",
        cluster_within_on=("posneg", "strict_clusters"),
    )
    C_n = mc["clusters_C"].nunique()
    C_pos = mc.metadata[mc["posneg"] == "pos"]["clusters_C"].nunique()
    C_neg = mc.metadata[mc["posneg"] == "neg"]["clusters_C"].nunique()
    results.append(("C: 2-stage strict(0.95 within=ct) -> merge(0.85 within=posneg,on=strict)", C_n, C_pos, C_neg))
    utils_analysis.judge_clustering(mc, "clusters_C", save_loc=str(QC / "schemes" / "schemeC_judge_clustering.png"))

    # Comparison table
    df = pd.DataFrame(results, columns=["scheme", "n_total", "n_pos", "n_neg"])
    log.info(f"\n{df.to_string(index=False)}")
    df.to_csv(QC / "schemes" / "scheme_comparison.tsv", sep="\t", index=False)


@safe("Section 4: per-CT splits + judge_clustering + heatmaps")
def section_4_per_ct(mc):
    summary = []
    cts = sorted(mc["celltype"].unique())
    for ct in cts:
        sub = mc[mc["celltype"] == ct]
        n = len(sub)
        n_cl = sub["final_clusters"].nunique()
        n_pos = sub.metadata[sub["posneg"] == "pos"]["final_clusters"].nunique() if "pos" in sub["posneg"].values else 0
        n_neg = sub.metadata[sub["posneg"] == "neg"]["final_clusters"].nunique() if "neg" in sub["posneg"].values else 0
        summary.append({"celltype": ct, "n_motifs": n, "n_clusters": n_cl, "n_pos_clusters": n_pos, "n_neg_clusters": n_neg})
        log.info(f"  {ct}: {n} motifs, {n_cl} clusters ({n_pos} pos + {n_neg} neg)")

        # Save subset MC + metadata + judge_clustering + heatmap (if n>=2 patterns)
        ct_dir = QC / "per_ct" / ct
        ct_dir.mkdir(exist_ok=True)
        sub.metadata.to_csv(ct_dir / f"{ct}_metadata.tsv", sep="\t", index=False)
        if n_cl >= 2:
            try:
                utils_analysis.judge_clustering(sub, "final_clusters", save_loc=str(ct_dir / f"{ct}_judge_clustering.png"))
                sub_sorted = sub.sort(by="final_clusters")
                sub_sorted.heatmap(save_loc=str(ct_dir / f"{ct}_heatmap.png"))
            except Exception as e:
                log.warning(f"  per-CT diagnostics failed for {ct}: {e}")

    pd.DataFrame(summary).to_csv(QC / "per_ct" / "summary_table.tsv", sep="\t", index=False)


@safe("Section 5: annotation x celltype crosstab + heatmap")
def section_5_annotation_x_ct(mc):
    ann_col = "annotation_database_name0"
    if ann_col not in mc.metadata.columns:
        log.warning(f"  column {ann_col} not in metadata — skipping")
        return
    ct = pd.crosstab(mc["celltype"], mc[ann_col])
    # restrict to top-30 most common annotations for legibility
    top = ct.sum(axis=0).sort_values(ascending=False).head(30).index
    ct_top = ct[top]
    ct_top.to_csv(QC / "annotation_x_celltype_top30.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(12, 7))
    sns.heatmap(ct_top, cmap="viridis", annot=True, fmt="d", cbar_kws={"label": "n motifs"}, ax=ax)
    ax.set_title("Top-30 annotations × celltype")
    fig.tight_layout(); fig.savefig(QC / "annotation_x_celltype_heatmap.png", dpi=110, bbox_inches="tight"); plt.close(fig)


@safe("Section 6: motif_entropy by annotation status")
def section_6_entropy_by_annotation(mc):
    if "annotation_database_name0" not in mc.metadata.columns:
        log.warning("  no annotation column")
        return
    md = mc.metadata.copy()
    md["annotated"] = md["annotation_database_name0"].astype(str).ne("")
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.boxplot(data=md, x="annotated", y="motif_entropy", hue="posneg", ax=ax)
    ax.set_title("motif_entropy by annotation × posneg")
    fig.tight_layout(); fig.savefig(QC / "entropy_by_annotation_posneg.png", dpi=110, bbox_inches="tight"); plt.close(fig)


@safe("Section 7: archetype examples (top motif per archetype)")
def section_7_archetypes(mc):
    archetypes = {
        "1_sharp_single_peak":           ("motif_entropy",            True),
        "2_random_noise_mix":            ("motif_entropy",            False),
        "3_noisy_uncertain_peaks":       ("weighted_base_entropy",    False),
        "4_broad_single_nucleotide":     ("posbase_entropy_score",    False),
        "5_gc_at_bias":                  ("copair_entropy_score",     False),
        "6_dinuc_repeats":               ("dinuc_entropy_score",      False),
        "7_few_seqlets":                 ("num_seqlets",              True),
        "8_posneg_inverted":             ("posneg_inverted",          False),
    }
    for arch_name, (metric, ascending) in archetypes.items():
        try:
            if metric not in mc.metadata.columns:
                log.warning(f"  {arch_name}: missing metric {metric}, skip")
                continue
            mc_sorted = mc.sort(by=metric, ascending=ascending)
            stack = mc_sorted.get_standard_motif_stack()
            if not len(stack):
                log.warning(f"  {arch_name}: empty stack")
                continue
            utils_plotting.plot_motif(stack[0], show=False, save_loc=str(QC / "archetypes" / f"archetype{arch_name}.png"))
        except Exception as e:
            log.warning(f"  archetype {arch_name} failed: {e}")


@safe("Section 8: mc_filtered vs mc_removed saved separately")
def section_8_save_splits(mc_initial):
    if "flag_remove" not in mc_initial.metadata.columns:
        log.warning("  no flag_remove on initial MC — skipping")
        return
    mc_filtered = mc_initial[~mc_initial.metadata["flag_remove"]]
    mc_removed = mc_initial[mc_initial.metadata["flag_remove"]]
    mc_filtered.save(str(PATH_OUT / "motif_compendium_filtered.mc"))
    mc_removed.save(str(PATH_OUT / "motif_compendium_removed.mc"))
    mc_filtered.metadata.to_csv(PATH_OUT / "motif_compendium_filtered_metadata.tsv", sep="\t", index=False)
    mc_removed.metadata.to_csv(PATH_OUT / "motif_compendium_removed_metadata.tsv", sep="\t", index=False)
    log.info(f"  filtered: {len(mc_filtered)}  removed: {len(mc_removed)}")


def main():
    log.info("Loading clustered MC for QC")
    mc = MotifCompendium.load(str(PATH_MC))
    log.info(f"  loaded {len(mc)} motifs")

    section_1_snapshot(mc)
    section_2_diagnostics(mc)
    section_3_alt_schemes(mc)
    section_4_per_ct(mc)
    section_5_annotation_x_ct(mc)
    section_6_entropy_by_annotation(mc)
    section_7_archetypes(mc)

    # Section 8 needs the INIT MC (with all 1905 motifs + flag_remove); load fresh
    log.info("Loading initial MC for mc_filtered/removed splits")
    mc_initial = MotifCompendium.load(str(PATH_MC_INIT))
    section_8_save_splits(mc_initial)

    # Save the post-QC clustered MC (now has scheme A + B + C columns)
    mc.metadata.to_csv(QC / "post_qc_metadata.tsv", sep="\t", index=False)
    mc.save(str(QC / "post_qc_clustered.mc"))
    log.info(f"Saved post-QC artifacts to {QC}/")
    log.info("DONE — Phase A.5 full QC complete. Review the qc/ outputs before proceeding to hit calling.")


if __name__ == "__main__":
    main()
