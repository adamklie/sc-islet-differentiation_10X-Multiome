#!/usr/bin/env python

import os
import glob
import pandas as pd
import snapatac2 as snap
from tqdm.auto import tqdm

# --- Paths ---
PATH_SAMPLE_META = "/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/1_get_data/sample_metadata.tsv"
PATH_CELL_META = "/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/1_get_data/cell_metadata.tsv"
PATH_OUT = "/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data"
PATH_ADATAS = os.path.join(PATH_OUT, "adatas")
PATH_SUBSET = os.path.join(PATH_OUT, "subset")

if __name__ == "__main__":

    # --- Load metadata ---
    sample_metadata = pd.read_csv(PATH_SAMPLE_META, sep="\t")
    cell_metadata = pd.read_csv(PATH_CELL_META, sep="\t")
    fragments_dict = sample_metadata.set_index("sample")["path_frag_file"].to_dict()
    cell_labels = dict(zip(
        cell_metadata["sample_name"] + "#" + cell_metadata["gex_barcode_cellranger"],
        cell_metadata["rna_annotation"],
    ))

    # --- 1. Import fragments ---
    os.makedirs(PATH_ADATAS, exist_ok=True)
    samples = list(fragments_dict.keys())
    frag_paths = list(fragments_dict.values())
    out_paths = [os.path.join(PATH_ADATAS, f"{s}.h5ad") for s in samples]

    adatas = snap.pp.import_fragments(
        frag_paths,
        chrom_sizes=snap.genome.hg38,
        file=out_paths,
        sorted_by_barcode=False,
        n_jobs=14
    )
    for adata in adatas:
        adata.close()

    # --- 2. Add sample prefix to barcodes ---
    for fpath in tqdm(glob.glob(os.path.join(PATH_ADATAS, "*.h5ad"))):
        sample = os.path.basename(fpath).replace(".h5ad", "")
        adata = snap.read(fpath, backed=None)
        adata.obs.index = [f"{sample}#{bc}" for bc in adata.obs_names]
        adata.write_h5ad(fpath)

    # --- 3. Create dataset and subset to annotated cells ---
    files = sorted(glob.glob(os.path.join(PATH_ADATAS, "*.h5ad")))
    dataset = snap.AnnDataSet(
        adatas=[(os.path.basename(f).replace(".h5ad", ""), f) for f in files],
        filename=os.path.join(PATH_OUT, "dataset.h5ads"),
        add_key="sample",
    )
    keep = [x for x in dataset.obs_names if x in cell_labels]
    dataset, _ = dataset.subset(obs_indices=keep, out=PATH_SUBSET)
    dataset.obs["cell_type"] = [cell_labels[x] for x in dataset.obs_names]
    dataset.close()

    print("DONE!")
    