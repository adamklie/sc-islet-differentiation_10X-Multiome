#!/usr/bin/env python
"""Run selected model analysis on a cisTopic object.

Loads a cisTopic object and trained models, selects a model at the specified
n_topics, runs UMAP, computes marker topics, generates heatmaps and composition
plots. Saves the cisTopic object with selected model and a topic AnnData.

Adapted from Muto_mouse-kidney_10X-Multiome/bin/5_topic_models/2_selected_model_analysis.ipynb

Usage:
    python run_selected_model_analysis.py \
        --cistopic_obj <path>.pkl \
        --models <merged_dir> \
        --output_dir <analysis_dir> \
        --n_topics 50 \
        --groupby cell_type
"""

import argparse
import glob
import os
import pickle

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from pycisTopic.cistopic_class import CistopicObject
from pycisTopic.lda_models import evaluate_models
from pycisTopic.clust_vis import run_umap

import logging
logging.basicConfig(format="%(asctime)s %(levelname)-8s %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)


def load_models(models_path):
    """Load models from a single .pkl or a directory of .pkl files."""
    if os.path.isfile(models_path):
        logger.info("Loading models from file: %s", models_path)
        with open(models_path, "rb") as f:
            models = pickle.load(f)
        if not isinstance(models, list):
            models = [models]
    elif os.path.isdir(models_path):
        logger.info("Loading models from directory: %s", models_path)
        models = []
        for pkl_file in sorted(glob.glob(os.path.join(models_path, "*.pkl"))):
            m = pickle.load(open(pkl_file, "rb"))
            if isinstance(m, list):
                models.extend(m)
            else:
                models.append(m)
    else:
        raise FileNotFoundError(f"Models path not found: {models_path}")
    logger.info("Loaded %d model(s)", len(models))
    return models


def main():
    parser = argparse.ArgumentParser(description="Selected model analysis for cisTopic")
    parser.add_argument("--cistopic_obj", required=True, help="Path to cisTopic object .pkl")
    parser.add_argument("--models", required=True, help="Path to models dir or .pkl")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--n_topics", type=int, required=True, help="Number of topics to select")
    parser.add_argument("--groupby", default="cell_type", help="Cell metadata column for grouping")
    parser.add_argument("--colors_tsv", default=None, help="Optional TSV with cell_type, color columns")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load cisTopic object
    logger.info("Loading cisTopic object from %s", args.cistopic_obj)
    with open(args.cistopic_obj, "rb") as f:
        cistopic_obj = pickle.load(f)
    logger.info("  %d regions x %d cells", *cistopic_obj.fragment_matrix.shape)

    # Load models
    models = load_models(args.models)

    # Select and add model
    logger.info("Selecting model with n_topics=%d", args.n_topics)
    model = evaluate_models(
        models,
        select_model=args.n_topics,
        return_model=True,
        metrics=["Arun_2010", "Cao_Juan_2009", "Minmo_2011", "loglikelihood"],
        plot_metrics=False,
    )
    cistopic_obj.add_LDA_model(model)

    # Extract cell-topic matrix
    cell_topic_df = cistopic_obj.selected_model.cell_topic.copy()
    cell_topic_df.columns = (
        cell_topic_df.columns.str.split("___", expand=True).get_level_values(0)
    )

    # Run UMAP
    logger.info("Running UMAP...")
    run_umap(cistopic_obj, target="cell", scale=False)

    # Build AnnData
    logger.info("Building topic AnnData...")
    adata_topics = ad.AnnData(
        X=cell_topic_df.values.T,
        obs=pd.DataFrame(index=cell_topic_df.columns),
        var=pd.DataFrame(index=cell_topic_df.index),
    )

    # Add UMAP coordinates
    umap_df = cistopic_obj.projections["cell"]["UMAP"].copy()
    umap_df.index = umap_df.index.str.split("___", expand=True).get_level_values(0)
    adata_topics.obsm["X_cisTopic"] = adata_topics.X.copy()
    adata_topics.obsm["X_umap"] = umap_df.reindex(adata_topics.obs_names)[["UMAP_1", "UMAP_2"]].values

    # Add cell metadata
    cell_meta = cistopic_obj.cell_data.copy()
    cell_meta.index = cell_meta.index.str.split("___", expand=True).get_level_values(0)
    for col in cell_meta.columns:
        adata_topics.obs[col] = cell_meta.reindex(adata_topics.obs_names)[col].values

    # Load colors if provided
    group_colors = None
    if args.colors_tsv and os.path.exists(args.colors_tsv):
        colors_df = pd.read_csv(args.colors_tsv, sep="\t")
        if "cell_type" in colors_df.columns and "color" in colors_df.columns:
            group_colors = dict(zip(colors_df["cell_type"], colors_df["color"]))

    # UMAP plot
    logger.info("Plotting UMAP...")
    group_col = args.groupby
    if group_col in adata_topics.obs.columns:
        # Ensure categorical
        if not hasattr(adata_topics.obs[group_col], "cat"):
            adata_topics.obs[group_col] = adata_topics.obs[group_col].astype("category")

        if group_colors:
            adata_topics.uns[f"{group_col}_colors"] = [
                group_colors.get(ct, "#999999")
                for ct in adata_topics.obs[group_col].cat.categories
            ]

        sc.pl.umap(adata_topics, color=group_col, show=False, frameon=False)
        plt.savefig(os.path.join(args.output_dir, "umap_overview.pdf"), bbox_inches="tight", dpi=150)
        plt.close()

    # Cell-topic heatmap
    logger.info("Plotting cell-topic heatmap...")
    if group_col in adata_topics.obs.columns:
        ax = sc.pl.heatmap(
            adata_topics,
            var_names=adata_topics.var_names,
            groupby=group_col,
            standard_scale="var",
            show=False,
            figsize=(12, 8),
        )
        plt.savefig(os.path.join(args.output_dir, "cell_topic_heatmap.pdf"), bbox_inches="tight", dpi=150)
        plt.close()

    # Scaled layer + neighbors + Leiden
    logger.info("Running Leiden clustering...")
    adata_topics.layers["scaled"] = sc.pp.scale(adata_topics, copy=True).X
    sc.pp.neighbors(adata_topics, use_rep="X_cisTopic", n_neighbors=30, metric="cosine")
    for res in [0.2, 0.5, 0.8, 1.0]:
        sc.tl.leiden(adata_topics, resolution=res, key_added=f"leiden_{res}")

    # Matrix plots
    logger.info("Plotting matrix plots...")
    if group_col in adata_topics.obs.columns:
        sc.pl.matrixplot(adata_topics, var_names=adata_topics.var_names, groupby=group_col,
                         standard_scale="var", show=False)
        plt.savefig(os.path.join(args.output_dir, "matrixplot_scaled.pdf"), bbox_inches="tight", dpi=150)
        plt.close()

        sc.pl.matrixplot(adata_topics, var_names=adata_topics.var_names, groupby=group_col,
                         layer="scaled", show=False)
        plt.savefig(os.path.join(args.output_dir, "matrixplot_zscore.pdf"), bbox_inches="tight", dpi=150)
        plt.close()

    # Marker topics
    logger.info("Computing marker topics...")
    if group_col in adata_topics.obs.columns:
        sc.tl.rank_genes_groups(adata_topics, groupby=group_col, method="wilcoxon", n_genes=10)
        sc.pl.rank_genes_groups(adata_topics, n_genes=10, sharey=False, show=False)
        plt.savefig(os.path.join(args.output_dir, "marker_topics.pdf"), bbox_inches="tight", dpi=150)
        plt.close()

        sc.pl.rank_genes_groups_heatmap(adata_topics, groupby=group_col, n_genes=2,
                                        standard_scale="var", show=False)
        plt.savefig(os.path.join(args.output_dir, "marker_topics_heatmap.pdf"), bbox_inches="tight", dpi=150)
        plt.close()

    # Cluster composition
    logger.info("Plotting cluster composition...")
    if group_col in adata_topics.obs.columns:
        for meta_col in [group_col, "diff_batch", "diff_stage"]:
            if meta_col not in adata_topics.obs.columns:
                continue
            ct = pd.crosstab(adata_topics.obs["leiden_0.5"], adata_topics.obs[meta_col], normalize="index")
            ct.plot.bar(stacked=True, figsize=(10, 5))
            plt.title(f"Cluster composition by {meta_col}")
            plt.tight_layout()
            plt.savefig(os.path.join(args.output_dir, f"cluster_composition_{meta_col}.pdf"),
                        bbox_inches="tight", dpi=150)
            plt.close()

    # Save objects
    logger.info("Saving outputs...")
    obj_path = os.path.join(args.output_dir, "cistopic_obj_with_model.pkl")
    with open(obj_path, "wb") as f:
        pickle.dump(cistopic_obj, f)
    logger.info("  Saved cisTopic object: %s", obj_path)

    # Save h5ad
    adata_save = adata_topics.copy()
    for col in adata_save.obs.columns:
        dtype = adata_save.obs[col].dtype
        if hasattr(dtype, "categories"):
            continue
        if adata_save.obs[col].dtype == object:
            adata_save.obs[col] = adata_save.obs[col].astype(str)
    h5ad_path = os.path.join(args.output_dir, f"topics_{args.n_topics}.h5ad")
    adata_save.write_h5ad(h5ad_path)
    logger.info("  Saved topic AnnData: %s", h5ad_path)

    logger.info("Done!")


if __name__ == "__main__":
    main()
