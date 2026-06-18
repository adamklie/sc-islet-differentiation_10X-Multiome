#!/usr/bin/env python
"""Build a combined RNA AnnData from 10X Multiome CellRanger outputs.

Reads each sample's filtered_feature_bc_matrix.h5, extracts Gene Expression
features, subsets to barcodes present in the cell metadata (already annotated
with cell types from Seurat/WNN), concatenates into a single AnnData.

Adapted from Muto_mouse-kidney_10X-Multiome/bin/2_process_data/5_build_rna_anndata.py

Usage:
    python build_rna_anndata.py
"""

import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def main():
    # ------- paths -------
    data_dir = "/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome"
    cell_meta_path = os.path.join(data_dir, "results/1_get_data/cell_metadata.tsv")
    sample_meta_path = os.path.join(data_dir, "results/1_get_data/sample_metadata.tsv")
    processed_dir = os.path.join(data_dir, "processed")
    out_path = os.path.join(data_dir, "results/2_process_data/rna_matrix.h5ad")

    # ------- sample metadata (maps sample_name -> cellranger folder) -------
    sample_meta = pd.read_csv(sample_meta_path, sep="\t")
    sample_to_folder = dict(zip(sample_meta["sample"], sample_meta["cellranger_folder"]))
    sample_to_batch = dict(zip(sample_meta["sample"], sample_meta["batch"]))
    print(f"Sample metadata: {len(sample_meta)} samples")

    # ------- cell metadata (from Seurat, already annotated) -------
    cell_meta = pd.read_csv(cell_meta_path, sep="\t", low_memory=False)
    # The first unnamed column is the Seurat cell ID — extract barcode and sample
    cell_meta["gex_barcode"] = cell_meta["gex_barcode_cellranger"]
    cell_meta["cell_type"] = cell_meta["rna_annotation"]
    print(f"Cell metadata: {len(cell_meta)} annotated cells, {cell_meta['cell_type'].nunique()} cell types")

    # Build per-sample valid barcode sets
    valid_barcodes_per_sample = {}
    for sample_name in cell_meta["sample_name"].unique():
        sub = cell_meta[cell_meta["sample_name"] == sample_name]
        valid_barcodes_per_sample[sample_name] = set(sub["gex_barcode"])
    print(f"Samples with annotated cells: {len(valid_barcodes_per_sample)}")

    # ------- read each sample -------
    adatas = []
    for sample_name, folder_name in sample_to_folder.items():
        batch = sample_to_batch[sample_name]
        h5_path = os.path.join(
            processed_dir, batch, folder_name, "outs", "filtered_feature_bc_matrix.h5"
        )

        if not os.path.exists(h5_path):
            print(f"  SKIP {sample_name} — {h5_path} not found")
            continue

        print(f"Reading {sample_name} ({folder_name})...")

        # Read 10X h5 — extract GEX only
        adata = sc.read_10x_h5(h5_path, gex_only=True)
        adata.var_names_make_unique()

        # Subset to annotated barcodes only
        valid = valid_barcodes_per_sample.get(sample_name, set())
        if not valid:
            print(f"  SKIP {sample_name} — no annotated barcodes")
            continue

        keep_mask = adata.obs_names.isin(valid)
        n_raw = adata.n_obs
        adata = adata[keep_mask].copy()
        print(f"  {n_raw} filtered -> {adata.n_obs} annotated cells, {adata.n_vars} genes")

        if adata.n_obs == 0:
            print(f"  SKIP {sample_name} — 0 cells after filtering")
            continue

        # Add metadata
        adata.obs["sample_name"] = sample_name
        adata.obs["batch"] = batch

        # Add cell-level metadata from Seurat
        sample_cell_meta = cell_meta[cell_meta["sample_name"] == sample_name].set_index("gex_barcode")
        for col in ["cell_type", "diff_batch", "diff_stage", "nCount_RNA", "nFeature_RNA"]:
            if col in sample_cell_meta.columns:
                adata.obs[col] = sample_cell_meta.reindex(adata.obs_names)[col].values

        adatas.append(adata)

    # ------- concatenate -------
    print(f"\nConcatenating {len(adatas)} samples...")
    combined = ad.concat(adatas, join="outer", merge="same")
    print(f"Combined: {combined.n_obs} cells, {combined.n_vars} genes")

    # Summary
    n_annotated = combined.obs["cell_type"].notna().sum()
    print(f"Cells with cell type annotation: {n_annotated}")
    print(f"Cell types: {sorted(combined.obs['cell_type'].dropna().unique())}")
    print(f"Samples: {combined.obs['sample_name'].nunique()}")

    # ------- save -------
    print(f"\nSaving to {out_path}")
    combined.write_h5ad(out_path)
    print("Done!")


if __name__ == "__main__":
    main()
