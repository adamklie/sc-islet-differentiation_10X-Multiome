#!/usr/bin/env python
"""Compute pseudobulk mean expression per cell type for TF genes.

Reads the combined RNA AnnData (from build_rna_anndata.py), subsets to
transcription factor genes relevant to the motif analysis, and computes
mean expression per cell type.

Usage:
    python build_pseudobulk_tf_expression.py
"""

import os

import anndata as ad
import numpy as np
import pandas as pd


def get_tf_gene_set(motif_annotations_path):
    """Extract TF gene names from motif annotations + known islet TFs."""
    motifs = pd.read_csv(motif_annotations_path)

    tf_genes = set()
    for ann in motifs["annotation"]:
        parts = ann.split("_")
        name = parts[0].replace(".mouse", "").replace("HUMAN", "").replace("MOUSE", "")
        tf_genes.add(name.upper())
        tf_genes.add(parts[0].upper())

    # Known islet differentiation TFs
    islet_tfs = [
        "PDX1", "NKX2-2", "NKX6-1", "PAX6", "NEUROD1", "RFX6", "ISL1", "MAFA", "MAFB",
        "FOXA2", "FOXA1", "HNF1B", "HNF4A", "HNF4G", "GATA4", "GATA6",
        "SOX9", "SOX17", "ONECUT1", "ONECUT2",
        "TEAD1", "TEAD4", "JUN", "FOS", "FOSB", "ATF1",
        "CTCF", "REST", "SP1", "SP2", "NRF1",
        "ETV6", "ELK1", "ETS1",
        "RFX2", "RFX4", "RFX5",
        "ARX", "INS", "GCG", "SST", "GHRL",
        "SIX1", "BARHL2", "CUX2", "GRHL2",
        "FOXM1", "FOXD2", "FOXC1", "FOXG1",
        "CEBPB", "SRF", "TBX4", "ZIC3", "ZIC4",
        "OTX2", "CRX", "NFIC", "NFYB",
        "HNF1A", "CDX2",
        "NEUROG3", "NKX2-5", "PAX4",
    ]
    for tf in islet_tfs:
        tf_genes.add(tf.upper())

    return tf_genes


def main():
    # Paths
    data_dir = "/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome"
    rna_path = os.path.join(data_dir, "results/2_process_data/rna_matrix.h5ad")
    motif_annotations_path = os.path.join(
        data_dir, "results/3_single_task_models/motifs/curation/motif_annotations_full.csv"
    )
    out_path = os.path.join(data_dir, "results/2_process_data/pseudobulk_tf_expression.csv")

    # Load RNA
    print("Loading RNA anndata...")
    adata = ad.read_h5ad(rna_path)
    adata.obs_names_make_unique()
    print(f"Shape: {adata.shape}")

    # Get TF genes
    tf_genes = get_tf_gene_set(motif_annotations_path)

    # Find matching genes in anndata
    var_upper = pd.Series(adata.var_names).str.upper()
    found_genes = [v for v, vu in zip(adata.var_names, var_upper) if vu in tf_genes]
    print(f"TF genes requested: {len(tf_genes)}, found in anndata: {len(found_genes)}")

    # Subset and compute pseudobulk
    adata_tf = adata[:, found_genes].copy()
    print("Computing pseudobulk mean expression per cell type...")

    pseudobulk = pd.DataFrame(
        index=sorted(adata_tf.obs["cell_type"].dropna().unique()),
        columns=found_genes,
        dtype=float,
    )

    for ct in pseudobulk.index:
        mask = adata_tf.obs["cell_type"] == ct
        pseudobulk.loc[ct] = np.asarray(adata_tf[mask].X.mean(axis=0)).flatten()

    pseudobulk.index.name = "cell_type"

    # Save
    pseudobulk.to_csv(out_path)
    print(f"Saved: {pseudobulk.shape} to {out_path}")

    # Summary
    print(f"\nTop expressed TFs (mean across cell types):")
    top = pseudobulk.mean(axis=0).sort_values(ascending=False).head(20)
    for gene, val in top.items():
        print(f"  {gene:15s} {val:.3f}")


if __name__ == "__main__":
    main()
