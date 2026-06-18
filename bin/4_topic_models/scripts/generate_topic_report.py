#!/usr/bin/env python
"""Generate a Topic MultiQC HTML report for the DeepTopic/CREsted pipeline.

Integrates all annotation layers (16+ data sources) into a standalone HTML report
with interactive tables, per-topic drill-down, and auto-assigned topic categories.

Usage:
    python generate_topic_report.py --base <results_dir> [--config <config_dir>] [--output-dir <dir>]
"""

import argparse
import base64
import io
import os
import sys
import traceback
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore", category=FutureWarning)

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 16,
})

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DATASET = "sc-islet-differentiation_10X-Multiome"
DEFAULT_BASE = f"/cellar/users/aklie/data/datasets/{DATASET}/results/4_topic_models/endocrine_replicate"
DEFAULT_CONFIG = f"/cellar/users/aklie/data/datasets/{DATASET}/config"

CATEGORY_COLORS = {
    "one-to-one": "#2196F3",
    "lineage": "#4CAF50",
    "shared": "#FF9800",
    "structural": "#9C27B0",
    "batch-driven": "#F44336",
    "unannotated": "#9E9E9E",
}

CATEGORY_ORDER = ["one-to-one", "lineage", "shared", "structural", "batch-driven", "unannotated"]

# cCRE colors
CCRE_COLORS = {
    "PLS": "#d62728",
    "pELS": "#ff7f0e",
    "dELS": "#2ca02c",
    "CTCF-only": "#9467bd",
    "DNase-H3K4me3": "#8c564b",
    "no_cCRE": "#cccccc",
}

# Structural motif families (top-level match)
STRUCTURAL_MOTIFS = {"CTCF", "BORIS", "NFY", "NFYB", "NRF", "NRF1", "Sp1", "Sp2", "YY1", "TYY1"}

# Differentiation adjacency graph
# Each cell type maps to its set of neighbors (undirected)
ADJACENCY = {
    "ENP_phase1": {"early_ENP"},
    "early_ENP": {"ENP_phase1", "late_ENP", "early_SC_beta", "early_SC_alpha", "early_SC_EC", "SC_delta_GHRL", "proliferating_endocrine"},
    "late_ENP": {"early_ENP", "early_SC_beta", "early_SC_alpha", "early_SC_EC", "SC_delta_GHRL", "proliferating_endocrine"},
    "early_SC_beta": {"early_ENP", "late_ENP", "late_SC_beta", "early_SC_alpha", "early_SC_EC", "SC_delta_GHRL", "proliferating_endocrine"},
    "late_SC_beta": {"early_SC_beta", "late_SC_alpha", "late_SC_EC"},
    "early_SC_alpha": {"early_ENP", "late_ENP", "late_SC_alpha", "early_SC_beta", "early_SC_EC", "SC_delta_GHRL", "proliferating_endocrine"},
    "late_SC_alpha": {"early_SC_alpha", "late_SC_beta", "late_SC_EC"},
    "early_SC_EC": {"early_ENP", "late_ENP", "late_SC_EC", "early_SC_beta", "early_SC_alpha", "SC_delta_GHRL", "proliferating_endocrine"},
    "late_SC_EC": {"early_SC_EC", "late_SC_beta", "late_SC_alpha"},
    "SC_delta_GHRL": {"early_ENP", "late_ENP", "early_SC_beta", "early_SC_alpha", "early_SC_EC", "proliferating_endocrine"},
    "proliferating_endocrine": {"early_ENP", "late_ENP", "early_SC_beta", "early_SC_alpha", "early_SC_EC", "SC_delta_GHRL"},
}


def _are_lineage_adjacent(cell_types):
    """Check if all cell types are pairwise adjacent in the differentiation hierarchy."""
    cts = [ct.strip() for ct in cell_types if ct.strip() in ADJACENCY]
    if len(cts) < 2:
        return True  # single or unknown cell types
    for i in range(len(cts)):
        for j in range(i + 1, len(cts)):
            if cts[j] not in ADJACENCY.get(cts[i], set()):
                return False
    return True


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------
def fig_to_base64(fig, dpi=150):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return b64


def embed_img(fig, dpi=150, width="100%"):
    b64 = fig_to_base64(fig, dpi=dpi)
    return f'<img src="data:image/png;base64,{b64}" style="width:{width}; max-width:1400px;">'


def html_section(title, content, level=2):
    anchor = title.lower().replace(" ", "-").replace("/", "-").replace("(", "").replace(")", "").replace("&", "and")
    tag = f"h{level}"
    return f'<{tag} id="{anchor}">{title}</{tag}>\n{content}\n'


def html_subsection(title, content):
    return html_section(title, content, level=3)


def table_html(df, max_rows=200, table_id=None):
    display_df = df.head(max_rows).copy()
    for col in display_df.select_dtypes(include=["float64", "float32"]).columns:
        display_df[col] = display_df[col].apply(lambda x: f"{x:.3f}" if pd.notna(x) else "---")
    display_df = display_df.fillna("---")
    return display_df.to_html(
        classes="data-table display", index=False,
        border=0, escape=False, table_id=table_id,
    )


def info_box(text, color="#e8f4fd", border="#2196F3"):
    return (
        f'<div style="background:{color}; border-left:4px solid {border}; '
        f'padding:12px 16px; margin:16px 0; border-radius:4px; font-size:15px; line-height:1.6;">'
        f'{text}</div>\n'
    )


def note_box(text):
    return info_box(text, color="#fff8e1", border="#ffc107")


def category_badge(cat):
    color = CATEGORY_COLORS.get(cat, "#999")
    return f'<span style="background:{color}; color:white; padding:2px 8px; border-radius:4px; font-size:12px; font-weight:500;">{cat}</span>'


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_all_data(base, config_dir):
    """Load all pre-computed annotation data. Returns dict of DataFrames (None if missing)."""
    data = {}

    def _load(key, path, **kwargs):
        if os.path.exists(path):
            try:
                data[key] = pd.read_csv(path, **kwargs)
                print(f"  Loaded {key}: {data[key].shape}")
            except Exception as e:
                print(f"  WARNING: Failed to load {key}: {e}")
                data[key] = None
        else:
            print(f"  MISSING: {key} ({path})")
            data[key] = None

    ann = os.path.join(base, "annotated")
    crested = os.path.join(base, "crested_model")

    _load("qc_metrics", f"{ann}/topic_qc/topic_qc_metrics.tsv", sep="\t", index_col=0)
    _load("auto_annotations", f"{ann}/topic_qc/topic_annotations_auto.tsv", sep="\t", index_col=0)
    _load("batch_correlation", f"{ann}/topic_qc/topic_batch_correlation.tsv", sep="\t")
    _load("batch_within_ct", f"{ann}/topic_qc/topic_batch_within_celltype.tsv", sep="\t")
    _load("ct_zscores", f"{ann}/metadata/cell_type_zscores.tsv", sep="\t", index_col=0)
    _load("ct_top", f"{ann}/metadata/cell_type_top.tsv", sep="\t")
    _load("ct_kruskal", f"{ann}/metadata/cell_type_kruskal_wallis.tsv", sep="\t")
    _load("motif_enrichment", f"{ann}/motifs/motif_enrichment_summary.tsv", sep="\t")
    _load("gene_enrichment", f"{ann}/genes/gene_enrichment_summary.tsv", sep="\t")
    _load("ccre_overlap", f"{ann}/ccre_overlap/topic_ccre_overlap.tsv", sep="\t")
    _load("region_annotations", f"{ann}/region_annotation/topic_region_annotations.tsv", sep="\t")
    _load("chromhmm_E087", f"{ann}/region_annotation/topic_chromhmm_E087.tsv", sep="\t")
    _load("chromhmm_E003", f"{ann}/region_annotation/topic_chromhmm_E003.tsv", sep="\t")
    _load("chromhmm_E066", f"{ann}/region_annotation/topic_chromhmm_E066.tsv", sep="\t")
    _load("chromhmm_E115", f"{ann}/region_annotation/topic_chromhmm_E115.tsv", sep="\t")
    _load("synthesis", f"{ann}/synthesis/topic_annotations.tsv", sep="\t")
    _load("crested_metrics", f"{crested}/evaluation/per_topic_metrics.tsv", sep="\t")
    _load("training_history", f"{crested}/training_history.csv", sep=",")
    _load("ct_metadata", f"{config_dir}/cell_type_metadata.tsv", sep="\t")

    return data


# ---------------------------------------------------------------------------
# Master table
# ---------------------------------------------------------------------------
def build_master_table(data):
    """Build one-row-per-topic summary table from all data sources."""

    # Start from QC metrics index
    qc = data["qc_metrics"]
    topics = sorted(qc.index.tolist(), key=lambda x: int(x.replace("Topic", "")))
    master = pd.DataFrame({"topic": topics})

    # QC metrics
    master["Gini"] = master["topic"].map(qc["Gini_index"]).values
    master["coherence"] = master["topic"].map(qc["Coherence"]).values
    master["assignments"] = master["topic"].map(qc["Assignments"]).values

    # CREsted metrics
    cm = data["crested_metrics"]
    if cm is not None:
        cm_idx = cm.set_index("topic")
        master["auROC"] = master["topic"].map(cm_idx["auroc"]).values
        master["auPR"] = master["topic"].map(cm_idx["aupr"]).values
        n_pos = master["topic"].map(cm_idx["n_positive"]).values
        n_tot = master["topic"].map(cm_idx["n_total"]).values
        aupr_random = n_pos / n_tot
        master["aupr_fc"] = master["auPR"] / aupr_random

    # Batch correlation
    bc = data["batch_correlation"]
    if bc is not None:
        bc_idx = bc.set_index("topic")
        master["batch_flag"] = master["topic"].map(bc_idx["batch_flag"]).values
        master["eta2_batch"] = master["topic"].map(bc_idx["eta2_batch"]).values
        master["eta2_celltype"] = master["topic"].map(bc_idx["eta2_celltype"]).values

    # Cell type info from auto_annotations
    aa = data["auto_annotations"]
    if aa is not None:
        master["cell_types_list"] = master["topic"].map(aa["cell_type"]).values
        master["n_cell_types"] = master["cell_types_list"].apply(
            lambda x: len([c.strip() for c in str(x).split(",") if c.strip()]) if pd.notna(x) else 0
        )

    # Top cell type from ct_top (rank 1)
    ct_top = data["ct_top"]
    if ct_top is not None:
        rank1 = ct_top[ct_top["rank"] == 1].copy()
        rank1_idx = rank1.set_index("topic")
        master["top_cell_type"] = master["topic"].map(rank1_idx["category"]).values
        master["top_ct_zscore"] = master["topic"].map(rank1_idx["z_score"]).values
        master["top_ct_weight"] = master["topic"].map(rank1_idx["mean_weight"]).values

    # Z-score gap (rank1 - rank2)
    if ct_top is not None:
        rank2 = ct_top[ct_top["rank"] == 2].copy()
        rank2_idx = rank2.set_index("topic")
        rank1_z = master["topic"].map(rank1_idx["z_score"])
        rank2_z = master["topic"].map(rank2_idx["z_score"])
        master["zscore_gap"] = rank1_z.values - rank2_z.values

    # cCRE overlap
    ccre = data["ccre_overlap"]
    if ccre is not None:
        ccre_idx = ccre.set_index("topic")
        master["n_regions"] = master["topic"].map(ccre_idx["n_regions"]).values
        # Dominant cCRE class
        ccre_cols = ["frac_PLS", "frac_pELS", "frac_dELS", "frac_CTCF-only", "frac_no_cCRE"]
        ccre_labels = ["PLS", "pELS", "dELS", "CTCF-only", "no_cCRE"]
        def _dominant_ccre(topic):
            if topic not in ccre_idx.index:
                return "---"
            row = ccre_idx.loc[topic]
            vals = [row.get(c, 0) for c in ccre_cols]
            return ccre_labels[np.argmax(vals)]
        master["dominant_cCRE"] = master["topic"].apply(_dominant_ccre)

    # Top HOMER motif (rank 1)
    me = data["motif_enrichment"]
    if me is not None:
        rank1_motif = me[me["rank"] == 1].copy()
        rank1_motif_idx = rank1_motif.set_index("topic")
        master["top_homer_motif"] = master["topic"].map(rank1_motif_idx["motif_name"]).values
        master["top_motif_logp"] = master["topic"].map(rank1_motif_idx["log_p_value"]).values
        master["top_motif_fold"] = master["topic"].map(rank1_motif_idx["fold_enrichment"]).values

    # Top pathway (GO:BP, first per topic)
    ge = data["gene_enrichment"]
    if ge is not None:
        gobp = ge[ge["source"] == "GO:BP"].copy()
        first_per_topic = gobp.groupby("topic").first().reset_index()
        fpt_idx = first_per_topic.set_index("topic")
        master["top_pathway"] = master["topic"].map(fpt_idx["term_name"]).values

    # Synthesis label
    syn = data["synthesis"]
    if syn is not None:
        syn_idx = syn.set_index("topic")
        master["label"] = master["topic"].map(syn_idx["label"]).values
        master["top_motifs_syn"] = master["topic"].map(syn_idx["top_motifs"]).values

    return master


# ---------------------------------------------------------------------------
# Categorization
# ---------------------------------------------------------------------------
def assign_categories(master, data):
    """Assign topic categories using rule-based logic. First match wins."""
    categories = []

    for _, row in master.iterrows():
        # Rule 1: batch-driven
        if row.get("batch_flag") == "BATCH-DRIVEN":
            categories.append("batch-driven")
            continue

        # Rule 2: structural
        top_motif = str(row.get("top_homer_motif", ""))
        motif_family = top_motif.split("(")[0].split("/")[0].strip()
        is_structural_motif = motif_family in STRUCTURAL_MOTIFS
        gini = row.get("Gini", 1.0)
        n_ct = row.get("n_cell_types", 0)
        if is_structural_motif and gini < 0.5 and n_ct >= 4:
            categories.append("structural")
            continue

        # Rule 3: one-to-one
        if n_ct == 1:
            categories.append("one-to-one")
            continue

        # Rules 4-5: lineage vs shared (2-3 cell types)
        if 2 <= n_ct <= 3:
            cts = [c.strip() for c in str(row.get("cell_types_list", "")).split(",")]
            if _are_lineage_adjacent(cts):
                categories.append("lineage")
            else:
                categories.append("shared")
            continue

        # Rule 6: unannotated (>= 4 cell types, not structural)
        categories.append("unannotated")

    master["category"] = categories
    master["manual_override"] = ""
    return master


def apply_manual_overrides(master, overrides_path):
    """Apply manual category overrides from a TSV file."""
    if not os.path.exists(overrides_path):
        return master
    overrides = pd.read_csv(overrides_path, sep="\t")
    if "manual_override" not in overrides.columns or "topic" not in overrides.columns:
        return master
    override_map = dict(zip(overrides["topic"], overrides["manual_override"]))
    for topic, override in override_map.items():
        if pd.notna(override) and str(override).strip():
            mask = master["topic"] == topic
            master.loc[mask, "category"] = str(override).strip()
            master.loc[mask, "manual_override"] = str(override).strip()
    return master


# ---------------------------------------------------------------------------
# Section 1: Summary Table
# ---------------------------------------------------------------------------
def gen_section_summary(master, data, ct_colors):
    html = ""

    html += info_box(
        "Master summary table integrating all annotation layers. One row per topic. "
        "Click column headers to sort. Use the search box to filter."
    )

    # Prepare display table
    display_cols = [
        "topic", "category", "label", "top_cell_type", "top_ct_zscore", "zscore_gap",
        "n_cell_types", "batch_flag", "auROC", "auPR", "aupr_fc",
        "n_regions", "Gini", "coherence", "top_homer_motif", "dominant_cCRE", "top_pathway",
    ]
    display = master[[c for c in display_cols if c in master.columns]].copy()

    # Add category badges
    if "category" in display.columns:
        display["category"] = display["category"].apply(category_badge)

    html += table_html(display, max_rows=50, table_id="master-summary")

    return html_section("1. Summary Table", html)


# ---------------------------------------------------------------------------
# Section 2: Cell Type Identity
# ---------------------------------------------------------------------------
def gen_section_celltype(master, data, ct_colors):
    html = ""

    # 2a: Z-score heatmap
    ct_zscores = data["ct_zscores"]
    if ct_zscores is not None:
        html += html_subsection("2a. Cell Type Z-Score Heatmap", "")

        # Sort topics by category then by top_ct_zscore descending
        sort_order = master.sort_values(
            ["category", "top_ct_zscore"],
            key=lambda x: x.map({c: i for i, c in enumerate(CATEGORY_ORDER)}) if x.name == "category" else -x,
            ascending=[True, True],
        )["topic"].tolist()

        # Filter to topics present in z-score matrix
        sort_order = [t for t in sort_order if t in ct_zscores.index]
        plot_df = ct_zscores.loc[sort_order]

        fig = plt.figure(figsize=(14, max(10, len(sort_order) * 0.25)))
        gs = gridspec.GridSpec(1, 3, width_ratios=[0.5, 10, 2], wspace=0.05)

        # Left: category color bar
        ax_cat = fig.add_subplot(gs[0])
        cat_map = dict(zip(master["topic"], master["category"]))
        for i, topic in enumerate(sort_order):
            cat = cat_map.get(topic, "unannotated")
            ax_cat.add_patch(plt.Rectangle((0, i), 1, 1, color=CATEGORY_COLORS.get(cat, "#999")))
        ax_cat.set_xlim(0, 1)
        ax_cat.set_ylim(0, len(sort_order))
        ax_cat.set_xticks([])
        ax_cat.set_yticks([])
        ax_cat.set_xlabel("Cat.", fontsize=9)
        ax_cat.invert_yaxis()

        # Center: heatmap
        ax_heat = fig.add_subplot(gs[1])
        sns.heatmap(
            plot_df.values, ax=ax_heat, cmap="RdBu_r", center=0, vmin=-3, vmax=3,
            xticklabels=plot_df.columns, yticklabels=sort_order,
            cbar_kws={"shrink": 0.5, "label": "Z-score"},
        )
        ax_heat.set_xlabel("")
        ax_heat.set_ylabel("")
        ax_heat.tick_params(axis="y", labelsize=8)
        ax_heat.tick_params(axis="x", labelsize=9, rotation=45)

        # Right: auROC + n_regions bars
        ax_side = fig.add_subplot(gs[2])
        auroc_map = dict(zip(master["topic"], master.get("auROC", pd.Series(dtype=float))))
        aurocs = [auroc_map.get(t, 0) for t in sort_order]
        ax_side.barh(range(len(sort_order)), aurocs, color="#2196F3", alpha=0.7)
        ax_side.set_xlim(0.5, 1.0)
        ax_side.set_ylim(-0.5, len(sort_order) - 0.5)
        ax_side.set_yticks([])
        ax_side.set_xlabel("auROC", fontsize=9)
        ax_side.invert_yaxis()

        fig.suptitle("Topic x Cell Type Z-Scores (sorted by category)", fontsize=14, y=1.01)
        html += embed_img(fig, dpi=120)

        # Legend
        legend_html = '<div style="margin:8px 0;">'
        for cat in CATEGORY_ORDER:
            c = CATEGORY_COLORS[cat]
            legend_html += f'<span style="display:inline-block; width:14px; height:14px; background:{c}; margin-right:4px; vertical-align:middle; border-radius:2px;"></span>{cat} &nbsp;&nbsp;'
        legend_html += '</div>'
        html += legend_html

    # 2b: Z-score gap scatter
    if "zscore_gap" in master.columns and "auROC" in master.columns:
        html += html_subsection("2b. Cell Type Specificity vs Sequence Predictability", "")

        fig, ax = plt.subplots(figsize=(8, 6))
        for cat in CATEGORY_ORDER:
            sub = master[master["category"] == cat]
            if len(sub) == 0:
                continue
            ax.scatter(sub["zscore_gap"], sub["auROC"], c=CATEGORY_COLORS[cat],
                       label=cat, s=60, alpha=0.8, edgecolors="white", linewidth=0.5)

        # Label outliers
        for _, row in master.iterrows():
            if row.get("zscore_gap", 0) > 2.0 or row.get("auROC", 0) > 0.85:
                ax.annotate(row["topic"].replace("Topic", "T"), (row["zscore_gap"], row["auROC"]),
                            fontsize=7, ha="left", va="bottom")

        ax.set_xlabel("Z-score gap (rank 1 - rank 2)")
        ax.set_ylabel("CREsted auROC")
        ax.legend(fontsize=8, loc="lower right")
        ax.set_title("Specificity vs Predictability")
        ax.grid(True, alpha=0.3)
        html += embed_img(fig, dpi=120, width="70%")

    return html_section("2. Cell Type Identity", html)


# ---------------------------------------------------------------------------
# Section 3: Quality & Confidence
# ---------------------------------------------------------------------------
def gen_section_quality(master, data, ct_colors):
    html = ""

    # 3a: Report card
    html += html_subsection("3a. Topic Report Card", "")

    metrics = ["auROC", "aupr_fc", "top_ct_zscore", "Gini", "coherence"]
    metric_labels = ["auROC", "auPR FC", "Top CT Z-score", "Gini Index", "Coherence"]
    available = [m for m in metrics if m in master.columns]

    if available:
        sorted_master = master.sort_values("auROC", ascending=True) if "auROC" in master.columns else master
        n = len(sorted_master)

        fig, axes = plt.subplots(1, len(available), figsize=(3 * len(available), max(8, n * 0.22)),
                                 sharey=True)
        if len(available) == 1:
            axes = [axes]

        for ax, metric, label in zip(axes, available, [metric_labels[metrics.index(m)] for m in available]):
            vals = sorted_master[metric].values
            colors = [CATEGORY_COLORS.get(c, "#999") for c in sorted_master["category"]]
            ax.barh(range(n), vals, color=colors, alpha=0.8)
            ax.set_xlabel(label, fontsize=10)
            ax.tick_params(axis="y", labelsize=7)

        axes[0].set_yticks(range(n))
        axes[0].set_yticklabels(sorted_master["topic"].values, fontsize=7)
        fig.suptitle("Topic Report Card (sorted by auROC)", fontsize=13, y=1.01)
        fig.tight_layout()
        html += embed_img(fig, dpi=100)

    # 3b: Training curves
    th = data.get("training_history")
    if th is not None:
        html += html_subsection("3b. CREsted Training Curves", "")
        curve_cols = [c for c in th.columns if c not in ("epoch", "Unnamed: 0")]
        # Group by train/val
        train_cols = [c for c in curve_cols if not c.startswith("val_")]
        val_cols = [c for c in curve_cols if c.startswith("val_")]

        fig, axes = plt.subplots(1, len(train_cols), figsize=(5 * len(train_cols), 4))
        if len(train_cols) == 1:
            axes = [axes]
        for ax, tc in zip(axes, train_cols):
            vc = f"val_{tc}" if f"val_{tc}" in th.columns else None
            ax.plot(th[tc], label=f"train", alpha=0.8)
            if vc:
                ax.plot(th[vc], label=f"val", alpha=0.8)
            ax.set_xlabel("Epoch")
            ax.set_ylabel(tc)
            ax.set_title(tc)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
        fig.tight_layout()
        html += embed_img(fig, dpi=100, width="90%")

    # 3c: Batch effects scatter
    bc = data.get("batch_correlation")
    if bc is not None:
        html += html_subsection("3c. Batch Effects", "")

        fig, ax = plt.subplots(figsize=(7, 6))
        flag_colors = {"OK": "#4CAF50", "MIXED": "#FF9800", "BATCH-DRIVEN": "#F44336"}
        for flag, color in flag_colors.items():
            sub = bc[bc["batch_flag"] == flag]
            ax.scatter(sub["eta2_batch"], sub["eta2_celltype"], c=color, label=flag, s=60, alpha=0.8)
        # Label batch-driven
        for _, row in bc[bc["batch_flag"] == "BATCH-DRIVEN"].iterrows():
            ax.annotate(row["topic"].replace("Topic", "T"),
                        (row["eta2_batch"], row["eta2_celltype"]),
                        fontsize=8, ha="left")
        ax.plot([0, 0.6], [0, 0.6], "k--", alpha=0.3, label="equal")
        ax.set_xlabel("eta2 (batch)")
        ax.set_ylabel("eta2 (cell type)")
        ax.set_title("Batch vs Cell Type Variance Explained")
        ax.legend()
        ax.grid(True, alpha=0.3)
        html += embed_img(fig, dpi=120, width="60%")

        # Summary
        counts = bc["batch_flag"].value_counts()
        html += info_box(
            f"<b>Batch flag distribution:</b> "
            + ", ".join(f"{k}: {v}" for k, v in counts.items())
        )

    # 3d: Performance vs regions
    if "n_regions" in master.columns and "auROC" in master.columns:
        html += html_subsection("3d. Performance vs Region Count", "")
        fig, ax = plt.subplots(figsize=(7, 5))
        for cat in CATEGORY_ORDER:
            sub = master[master["category"] == cat]
            if len(sub) == 0:
                continue
            ax.scatter(sub["n_regions"], sub["auROC"], c=CATEGORY_COLORS[cat],
                       label=cat, s=60, alpha=0.8)
        ax.set_xlabel("Number of Regions")
        ax.set_ylabel("auROC")
        ax.set_title("CREsted Performance vs Region Count")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        html += embed_img(fig, dpi=120, width="60%")

    return html_section("3. Quality and Confidence", html)


# ---------------------------------------------------------------------------
# Section 4: Regulatory Composition
# ---------------------------------------------------------------------------
def gen_section_regulatory(master, data, ct_colors):
    html = ""

    # 4a: cCRE stacked bar
    ccre = data.get("ccre_overlap")
    if ccre is not None:
        html += html_subsection("4a. cCRE Class Composition", "")

        # Sort by category
        topic_order = master.sort_values(
            ["category"],
            key=lambda x: x.map({c: i for i, c in enumerate(CATEGORY_ORDER)}),
        )["topic"].tolist()

        ccre_idx = ccre.set_index("topic")
        ccre_cols = ["frac_PLS", "frac_pELS", "frac_dELS", "frac_CTCF-only", "frac_DNase-H3K4me3", "frac_no_cCRE"]
        ccre_labels = ["PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3", "no cCRE"]
        ccre_colors = [CCRE_COLORS.get(l, "#ccc") for l in ["PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3", "no_cCRE"]]

        available_topics = [t for t in topic_order if t in ccre_idx.index]
        plot_data = ccre_idx.loc[available_topics, [c for c in ccre_cols if c in ccre_idx.columns]]

        fig, ax = plt.subplots(figsize=(14, 6))
        bottom = np.zeros(len(available_topics))
        for i, (col, label, color) in enumerate(zip(
            [c for c in ccre_cols if c in plot_data.columns],
            ccre_labels[:len([c for c in ccre_cols if c in plot_data.columns])],
            ccre_colors,
        )):
            vals = plot_data[col].values
            ax.bar(range(len(available_topics)), vals, bottom=bottom, color=color, label=label, width=0.8)
            bottom += vals

        ax.set_xticks(range(len(available_topics)))
        ax.set_xticklabels([t.replace("Topic", "T") for t in available_topics], rotation=90, fontsize=7)
        ax.set_ylabel("Fraction")
        ax.set_title("cCRE Class Composition per Topic")
        ax.legend(fontsize=8, ncol=3, loc="upper right")
        fig.tight_layout()
        html += embed_img(fig, dpi=120)

    # 4b: ChromHMM heatmap (E087 pancreas)
    chromhmm = data.get("chromhmm_E087")
    if chromhmm is not None:
        html += html_subsection("4b. ChromHMM States (E087 - Pancreatic Islets)", "")

        ch_idx = chromhmm.set_index("topic") if "topic" in chromhmm.columns else chromhmm
        # Get fraction columns
        frac_cols = [c for c in ch_idx.columns if c.startswith("E087_frac_")]
        if frac_cols:
            # Clean column names
            clean_names = [c.replace("E087_frac_", "").replace("_", " ") for c in frac_cols]

            # Filter to informative states (max > 0.03)
            maxvals = ch_idx[frac_cols].max()
            keep = maxvals[maxvals > 0.03].index.tolist()
            keep_clean = [c.replace("E087_frac_", "").replace("_", " ") for c in keep]

            topic_order = master.sort_values(
                ["category"],
                key=lambda x: x.map({c: i for i, c in enumerate(CATEGORY_ORDER)}),
            )["topic"].tolist()
            avail = [t for t in topic_order if t in ch_idx.index]

            fig, ax = plt.subplots(figsize=(10, max(8, len(avail) * 0.22)))
            sns.heatmap(
                ch_idx.loc[avail, keep].values, ax=ax,
                xticklabels=keep_clean, yticklabels=avail,
                cmap="YlOrRd", vmin=0, vmax=0.5,
                cbar_kws={"shrink": 0.5, "label": "Fraction"},
            )
            ax.set_title("ChromHMM State Fractions (E087 Pancreatic Islets)")
            ax.tick_params(axis="y", labelsize=7)
            ax.tick_params(axis="x", labelsize=8, rotation=45)
            fig.tight_layout()
            html += embed_img(fig, dpi=100)

    return html_section("4. Regulatory Composition", html)


# ---------------------------------------------------------------------------
# Section 5: Motif & Gene Enrichment
# ---------------------------------------------------------------------------
def gen_section_enrichment(master, data, ct_colors):
    html = ""

    # 5a: Motif dotplot
    me = data.get("motif_enrichment")
    if me is not None:
        html += html_subsection("5a. Top HOMER Motifs per Topic", "")

        top3 = me[me["rank"] <= 3].copy()

        # Get unique motifs (by name) and topics
        topic_order = master.sort_values(
            ["category"],
            key=lambda x: x.map({c: i for i, c in enumerate(CATEGORY_ORDER)}),
        )["topic"].tolist()

        # Pivot: for each (topic, motif_name), store log_p_value and fold_enrichment
        motif_names = top3["motif_name"].unique()

        # Build grid
        fig, ax = plt.subplots(figsize=(16, max(6, len(motif_names) * 0.25)))

        for _, row in top3.iterrows():
            if row["topic"] not in topic_order:
                continue
            x = topic_order.index(row["topic"])
            y_labels = list(motif_names)
            if row["motif_name"] not in y_labels:
                continue
            y = y_labels.index(row["motif_name"])
            size = min(abs(row["log_p_value"]), 500) / 500 * 200 + 10
            color = min(row["fold_enrichment"], 5) / 5
            ax.scatter(x, y, s=size, c=[plt.cm.Reds(color)], alpha=0.8,
                       edgecolors="grey", linewidth=0.3)

        ax.set_xticks(range(len(topic_order)))
        ax.set_xticklabels([t.replace("Topic", "T") for t in topic_order], rotation=90, fontsize=7)
        ax.set_yticks(range(len(motif_names)))
        ax.set_yticklabels(motif_names, fontsize=7)
        ax.set_title("HOMER Motif Enrichment (top 3 per topic)")
        ax.set_xlabel("Topics (sorted by category)")
        ax.grid(True, alpha=0.15)
        fig.tight_layout()
        html += embed_img(fig, dpi=100)

        html += note_box(
            "Dot size proportional to -log10(p-value) (capped at 500). "
            "Color intensity proportional to fold enrichment (capped at 5x)."
        )

    # 5b: Top pathway table
    ge = data.get("gene_enrichment")
    if ge is not None:
        html += html_subsection("5b. Top GO:BP Pathway per Topic", "")
        gobp = ge[ge["source"] == "GO:BP"].copy()
        first_per_topic = gobp.groupby("topic").first().reset_index()
        display = first_per_topic[["topic", "term_name", "p_value", "intersection_size", "query_size"]].copy()
        display["coverage"] = (display["intersection_size"] / display["query_size"] * 100).round(1).astype(str) + "%"
        display = display[["topic", "term_name", "p_value", "coverage"]]
        display.columns = ["Topic", "Top GO:BP Term", "P-value", "Coverage"]
        # Sort by topic number
        display["_sort"] = display["Topic"].str.replace("Topic", "").astype(int)
        display = display.sort_values("_sort").drop("_sort", axis=1)
        html += table_html(display, table_id="pathway-table")

    # 5c: De novo motif wall placeholder
    html += html_subsection("5c. De Novo Motifs (CREsted MoDISco)", "")
    motifs_dir = os.path.join(os.path.dirname(data.get("_base", "")), "crested_model", "motifs")
    cluster_pfm = os.path.join(motifs_dir, "cluster", "clustered_motifs.pfm") if os.path.isdir(motifs_dir) else ""
    tfs_file = os.path.join(motifs_dir, "tfs_initial.txt") if os.path.isdir(motifs_dir) else ""

    if os.path.exists(cluster_pfm) and os.path.exists(tfs_file) and os.path.getsize(tfs_file) > 0:
        html += info_box("De novo motif clustering complete. Motif wall generation available in notebook 4.")
    else:
        html += note_box(
            "De novo motif clustering is still running or has not completed yet. "
            "Re-run this report after the clustering job finishes to populate this section."
        )

    return html_section("5. Motif and Gene Enrichment", html)


# ---------------------------------------------------------------------------
# Section 6: Topic Categories
# ---------------------------------------------------------------------------
def gen_section_categories(master, data, ct_colors):
    html = ""

    # 6a: Category definitions
    html += html_subsection("6a. Category Definitions", "")
    defs = """
    <table class="data-table" style="width:auto;">
    <thead><tr><th>Category</th><th>Definition</th><th>Rule</th></tr></thead>
    <tbody>
    <tr><td>{badge_oto}</td><td>Maps to exactly 1 cell type</td><td>1 cell type in auto annotation</td></tr>
    <tr><td>{badge_lin}</td><td>Maps to 2-3 adjacent cell types in differentiation</td><td>2-3 CTs, all pairwise adjacent in hierarchy</td></tr>
    <tr><td>{badge_sha}</td><td>Maps to 2-3 non-adjacent cell types</td><td>2-3 CTs, not all adjacent</td></tr>
    <tr><td>{badge_str}</td><td>Structural/housekeeping program</td><td>Top motif is CTCF/NFY/Sp1 + Gini &lt; 0.5 + n_CT &ge; 4</td></tr>
    <tr><td>{badge_bat}</td><td>Dominated by batch effects</td><td>batch_flag = BATCH-DRIVEN</td></tr>
    <tr><td>{badge_una}</td><td>No clear cell type assignment</td><td>&ge; 4 cell types, not structural</td></tr>
    </tbody></table>
    """.format(
        badge_oto=category_badge("one-to-one"),
        badge_lin=category_badge("lineage"),
        badge_sha=category_badge("shared"),
        badge_str=category_badge("structural"),
        badge_bat=category_badge("batch-driven"),
        badge_una=category_badge("unannotated"),
    )
    html += defs

    html += note_box(
        "Categories are assigned by rule (first match wins, in order listed above). "
        "The <code>topic_categories.tsv</code> file includes a <code>manual_override</code> column "
        "you can edit. Re-run with <code>--manual-overrides</code> to apply your edits."
    )

    # 6b: Distribution bar chart
    html += html_subsection("6b. Category Distribution", "")
    counts = master["category"].value_counts()
    fig, ax = plt.subplots(figsize=(8, 4))
    bars = [counts.get(c, 0) for c in CATEGORY_ORDER]
    colors = [CATEGORY_COLORS[c] for c in CATEGORY_ORDER]
    ax.bar(CATEGORY_ORDER, bars, color=colors)
    ax.set_ylabel("Number of Topics")
    ax.set_title(f"Topic Categories (n={len(master)})")
    for i, v in enumerate(bars):
        if v > 0:
            ax.text(i, v + 0.3, str(v), ha="center", fontsize=11, fontweight="bold")
    ax.set_ylim(0, max(bars) * 1.2 if bars else 10)
    plt.xticks(rotation=30, ha="right")
    fig.tight_layout()
    html += embed_img(fig, dpi=120, width="60%")

    # 6c: Category summary stats
    html += html_subsection("6c. Summary Statistics by Category", "")
    agg_cols = {"topic": "count"}
    for c in ["auROC", "auPR", "aupr_fc", "n_regions", "Gini", "coherence", "top_ct_zscore"]:
        if c in master.columns:
            agg_cols[c] = "mean"

    stats = master.groupby("category").agg(agg_cols).rename(columns={"topic": "n_topics"})
    stats = stats.reindex([c for c in CATEGORY_ORDER if c in stats.index])
    stats = stats.reset_index()
    html += table_html(stats, table_id="category-stats")

    return html_section("6. Topic Categories", html)


# ---------------------------------------------------------------------------
# Section 7: Per-Topic Drill-Down
# ---------------------------------------------------------------------------
def gen_section_drilldown(master, data, ct_colors):
    html = ""
    html += info_box("Click a topic tab to view its full annotation profile.")

    topics = sorted(master["topic"].tolist(), key=lambda x: int(x.replace("Topic", "")))

    # Tab buttons
    html += '<div class="tab-container">\n<div class="tab-buttons">\n'
    for i, topic in enumerate(topics):
        active = " active" if i == 0 else ""
        short = topic.replace("Topic", "T")
        html += f'<button class="tab-btn{active}" onclick="showDrillTab(\'{topic}\')" style="font-size:11px; padding:4px 8px;">{short}</button>\n'
    html += '</div>\n'

    # Tab content
    ct_zscores = data.get("ct_zscores")
    ccre = data.get("ccre_overlap")
    chromhmm = data.get("chromhmm_E087")
    me = data.get("motif_enrichment")
    ge = data.get("gene_enrichment")
    bc = data.get("batch_correlation")
    cm = data.get("crested_metrics")

    for i, topic in enumerate(topics):
        display = "block" if i == 0 else "none"
        html += f'<div id="drill-tab-{topic}" style="display:{display}; padding:12px 0;">\n'

        row = master[master["topic"] == topic].iloc[0]

        # Header
        cat = row.get("category", "---")
        label = row.get("label", "---")
        html += f'<h4>{topic}: {label} {category_badge(cat)}</h4>\n'

        # Metrics summary
        metrics_html = '<div style="display:flex; gap:24px; flex-wrap:wrap; margin:8px 0;">'
        for metric, val in [("auROC", row.get("auROC")), ("auPR", row.get("auPR")),
                            ("auPR FC", row.get("aupr_fc")), ("Regions", row.get("n_regions")),
                            ("Gini", row.get("Gini")), ("Batch", row.get("batch_flag"))]:
            if pd.notna(val):
                fmt = f"{val:.3f}" if isinstance(val, float) else str(int(val)) if isinstance(val, (int, np.integer, np.floating)) and metric == "Regions" else str(val)
                metrics_html += f'<div><b>{metric}:</b> {fmt}</div>'
        metrics_html += '</div>\n'
        html += metrics_html

        # Cell type z-score bar
        if ct_zscores is not None and topic in ct_zscores.index:
            zscores = ct_zscores.loc[topic].sort_values(ascending=True)
            fig, ax = plt.subplots(figsize=(5, 3))
            bar_colors = [ct_colors.get(ct, "#999") for ct in zscores.index]
            ax.barh(range(len(zscores)), zscores.values, color=bar_colors)
            ax.set_yticks(range(len(zscores)))
            ax.set_yticklabels(zscores.index, fontsize=8)
            ax.set_xlabel("Z-score")
            ax.set_title(f"{topic} Cell Type Specificity", fontsize=10)
            ax.axvline(0, color="black", linewidth=0.5)
            ax.grid(True, alpha=0.2, axis="x")
            fig.tight_layout()
            html += embed_img(fig, dpi=100, width="50%")

        # cCRE pie
        if ccre is not None:
            ccre_idx = ccre.set_index("topic")
            if topic in ccre_idx.index:
                ccre_row = ccre_idx.loc[topic]
                pie_cols = ["frac_PLS", "frac_pELS", "frac_dELS", "frac_CTCF-only", "frac_no_cCRE"]
                pie_labels = ["PLS", "pELS", "dELS", "CTCF-only", "no cCRE"]
                pie_colors_list = [CCRE_COLORS.get(l, "#ccc") for l in ["PLS", "pELS", "dELS", "CTCF-only", "no_cCRE"]]
                vals = [ccre_row.get(c, 0) for c in pie_cols]
                fig, ax = plt.subplots(figsize=(4, 3))
                ax.pie(vals, labels=pie_labels, colors=pie_colors_list, autopct="%1.0f%%",
                       textprops={"fontsize": 8}, startangle=90)
                ax.set_title(f"cCRE Composition (n={int(ccre_row.get('n_regions', 0))})", fontsize=10)
                fig.tight_layout()
                html += embed_img(fig, dpi=100, width="35%")

        # Top 5 HOMER motifs table
        if me is not None:
            topic_motifs = me[(me["topic"] == topic) & (me["rank"] <= 5)].copy()
            if len(topic_motifs) > 0:
                display_cols = ["rank", "motif_name", "motif_family", "log_p_value", "fold_enrichment", "pct_target", "pct_background"]
                avail_cols = [c for c in display_cols if c in topic_motifs.columns]
                html += '<h5>Top HOMER Motifs</h5>\n'
                html += table_html(topic_motifs[avail_cols])

        # Top 5 pathways
        if ge is not None:
            topic_genes = ge[(ge["topic"] == topic) & (ge["source"] == "GO:BP")].head(5).copy()
            if len(topic_genes) > 0:
                gcols = ["term_name", "p_value", "intersection_size", "query_size"]
                avail_gcols = [c for c in gcols if c in topic_genes.columns]
                html += '<h5>Top GO:BP Pathways</h5>\n'
                html += table_html(topic_genes[avail_gcols])

        html += '</div>\n'

    html += '</div>\n'  # close tab-container

    return html_section("7. Per-Topic Drill-Down", html)


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------
def build_html(sections, timestamp, title="Topic MultiQC Report"):
    toc_items = [
        ("1. Summary Table", "1.-summary-table"),
        ("2. Cell Type Identity", "2.-cell-type-identity"),
        ("3. Quality and Confidence", "3.-quality-and-confidence"),
        ("4. Regulatory Composition", "4.-regulatory-composition"),
        ("5. Motif and Gene Enrichment", "5.-motif-and-gene-enrichment"),
        ("6. Topic Categories", "6.-topic-categories"),
        ("7. Per-Topic Drill-Down", "7.-per-topic-drill-down"),
    ]
    toc_html = '<div id="toc"><h3>Contents</h3><ol>\n'
    for label, anchor in toc_items:
        toc_html += f'<li><a href="#{anchor}">{label}</a></li>\n'
    toc_html += '</ol></div>\n'

    content = toc_html + "\n".join(sections)

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
<script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
<script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
<style>
body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    max-width: 1500px; margin: 0 auto; padding: 24px; color: #333; line-height: 1.7; font-size: 16px;
}}
h1 {{ color: #1a1a2e; border-bottom: 3px solid #16213e; padding-bottom: 12px; font-size: 32px; }}
h2 {{ color: #16213e; border-bottom: 2px solid #ccc; padding-bottom: 8px; margin-top: 50px; font-size: 26px; }}
h3 {{ color: #0f3460; font-size: 22px; margin-top: 35px; }}
h4 {{ color: #1a237e; font-size: 18px; margin-top: 25px; }}
h5 {{ color: #333; font-size: 15px; margin-top: 16px; }}
p {{ font-size: 15px; line-height: 1.7; }}
img {{ display: block; margin: 18px auto; border: 1px solid #eee; border-radius: 6px; }}
table.data-table {{ border-collapse: collapse; width: 100%; font-size: 13px; margin: 18px 0; }}
table.data-table th {{ background: #16213e; color: white; padding: 10px 14px; text-align: left; font-size: 14px; }}
table.data-table td {{ padding: 8px 14px; border-bottom: 1px solid #eee; }}
table.data-table tr:nth-child(even) {{ background: #f8f9fa; }}
table.data-table tr:hover {{ background: #e8f0fe; }}
pre {{ background: #f5f5f5; padding: 14px; border-radius: 6px; overflow-x: auto; font-size: 14px; }}
code {{ background: #f0f0f0; padding: 2px 8px; border-radius: 4px; font-size: 13px; }}
#toc {{ background: #f8f9fa; padding: 24px; border-radius: 10px; margin: 24px 0; font-size: 16px; }}
#toc a {{ text-decoration: none; color: #0f3460; }}
#toc a:hover {{ text-decoration: underline; }}
#toc li {{ margin: 6px 0; }}
.timestamp {{ color: #888; font-size: 14px; }}
.tab-container {{ margin: 16px 0; }}
.tab-buttons {{ display: flex; flex-wrap: wrap; gap: 4px; margin-bottom: 12px; border-bottom: 2px solid #eee; padding-bottom: 4px; }}
.tab-btn {{ padding: 8px 16px; border: none; background: #f0f0f0; cursor: pointer; font-size: 13px; font-weight: 500; border-radius: 6px 6px 0 0; transition: background 0.2s; }}
.tab-btn:hover {{ background: #e0e0e0; }}
.tab-btn.active {{ background: #16213e; color: white; }}
</style>
<script>
function showDrillTab(topic) {{
    document.querySelectorAll('[id^="drill-tab-"]').forEach(el => el.style.display = 'none');
    document.getElementById('drill-tab-' + topic).style.display = 'block';
    // Update buttons
    var container = document.getElementById('drill-tab-' + topic).closest('.tab-container');
    if (container) {{
        container.querySelectorAll('.tab-btn').forEach(btn => btn.classList.remove('active'));
        container.querySelectorAll('.tab-btn').forEach(btn => {{
            if (btn.textContent === topic.replace('Topic', 'T')) btn.classList.add('active');
        }});
    }}
}}
</script>
</head>
<body>
<h1>{title}</h1>
<p class="timestamp">Generated: {timestamp} | Subset: endocrine_replicate | Topics: 50</p>

{content}

<hr>
<p class="timestamp">Generated by <code>generate_topic_report.py</code></p>

<script>
$(document).ready(function() {{
    $('table.data-table').each(function() {{
        if ($.fn.DataTable.isDataTable(this)) return;
        var headerCols = $(this).find('thead tr:first th').length;
        var firstRowCols = $(this).find('tbody tr:first td').length;
        if (headerCols === 0 || firstRowCols === 0 || headerCols !== firstRowCols) return;
        try {{
            $(this).DataTable({{
                pageLength: 25,
                lengthMenu: [10, 25, 50, 100],
                order: [],
                dom: 'lfrtip',
                scrollX: true,
            }});
        }} catch(e) {{
            console.warn('DataTables init failed for table:', e);
        }}
    }});
}});
</script>
</body>
</html>"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Generate Topic MultiQC HTML report")
    parser.add_argument("--base", default=DEFAULT_BASE, help="Base results directory for a topic model subset")
    parser.add_argument("--config", default=DEFAULT_CONFIG, help="Config directory with cell_type_metadata.tsv")
    parser.add_argument("--output-dir", default=None, help="Output directory (default: {base}/summary)")
    parser.add_argument("--manual-overrides", default=None, help="Path to topic_categories.tsv with manual overrides")
    args = parser.parse_args()

    output_dir = args.output_dir or os.path.join(args.base, "summary")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Topic MultiQC Report Generator")
    print(f"  Base: {args.base}")
    print(f"  Config: {args.config}")
    print(f"  Output: {output_dir}")
    print()

    # Phase 1: Load
    print("Phase 1: Loading data...")
    data = load_all_data(args.base, args.config)
    data["_base"] = args.base  # stash for later

    # Cell type colors
    ct_meta = data.get("ct_metadata")
    ct_colors = {}
    if ct_meta is not None:
        ct_colors = dict(zip(ct_meta["cell_type"], ct_meta["color"]))

    # Phase 2: Build master table + categorize
    print("\nPhase 2: Building master table...")
    master = build_master_table(data)
    master = assign_categories(master, data)

    # Apply manual overrides
    if args.manual_overrides:
        print(f"  Applying manual overrides from {args.manual_overrides}")
        master = apply_manual_overrides(master, args.manual_overrides)

    # Save TSVs
    master_path = os.path.join(output_dir, "topic_master_summary.tsv")
    master.to_csv(master_path, sep="\t", index=False)
    print(f"  Saved: {master_path}")

    cat_path = os.path.join(output_dir, "topic_categories.tsv")
    master[["topic", "category", "manual_override", "cell_types_list", "n_cell_types", "batch_flag"]].to_csv(
        cat_path, sep="\t", index=False)
    print(f"  Saved: {cat_path}")

    # Category summary stats
    agg_cols = {"topic": "count"}
    for c in ["auROC", "auPR", "aupr_fc", "n_regions", "Gini", "coherence", "top_ct_zscore"]:
        if c in master.columns:
            agg_cols[c] = "mean"
    cat_stats = master.groupby("category").agg(agg_cols).rename(columns={"topic": "n_topics"})
    cat_stats_path = os.path.join(output_dir, "category_summary_stats.tsv")
    cat_stats.to_csv(cat_stats_path, sep="\t")
    print(f"  Saved: {cat_stats_path}")

    # Category distribution
    print(f"\n  Category distribution:")
    for cat in CATEGORY_ORDER:
        n = (master["category"] == cat).sum()
        if n > 0:
            print(f"    {cat}: {n}")

    # Phase 3: Generate sections
    print("\nPhase 3: Generating report sections...")
    section_generators = [
        ("1. Summary Table", gen_section_summary),
        ("2. Cell Type Identity", gen_section_celltype),
        ("3. Quality & Confidence", gen_section_quality),
        ("4. Regulatory Composition", gen_section_regulatory),
        ("5. Motif & Gene Enrichment", gen_section_enrichment),
        ("6. Topic Categories", gen_section_categories),
        ("7. Per-Topic Drill-Down", gen_section_drilldown),
    ]

    sections = []
    for name, gen_func in section_generators:
        print(f"  {name}...")
        try:
            sections.append(gen_func(master, data, ct_colors))
        except Exception as e:
            tb = traceback.format_exc()
            print(f"    ERROR: {e}")
            sections.append(html_section(name, f'<p style="color:red;">Error generating section: {e}<pre>{tb}</pre></p>'))

    # Phase 4: Assemble HTML
    print("\nPhase 4: Assembling HTML...")
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    html = build_html(sections, timestamp)

    html_path = os.path.join(output_dir, "topic_report.html")
    with open(html_path, "w") as f:
        f.write(html)

    size_mb = os.path.getsize(html_path) / 1024 / 1024
    print(f"\nDone! Report saved: {html_path} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
