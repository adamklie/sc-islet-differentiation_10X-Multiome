"""
Step 0: Compute topic QC metrics and annotations.

Runs pycisTopic's compute_topic_metrics and topic_annotation on the cistopic object,
then saves the results as TSVs for downstream figure scripts.

Inputs:
    - cistopic_obj_with_model.pkl (from analysis/)
Outputs:
    - topic_qc_metrics.tsv (coherence, gini, assignments, etc.)
    - topic_annotations_auto.tsv (automatic topic-to-cell-type classification)
"""

import argparse
import pickle
import pandas as pd
from pycisTopic.topic_qc import compute_topic_metrics, topic_annotation


def main():
    parser = argparse.ArgumentParser(description="Compute topic QC metrics and annotations")
    parser.add_argument("--cistopic_obj", required=True, help="Path to cistopic object pickle")
    parser.add_argument("--annot_var", default="cell_type", help="Cell metadata column for annotation")
    parser.add_argument("--general_topic_thr", type=float, default=0.2, help="Threshold for general topic classification")
    parser.add_argument("--output_dir", required=True, help="Output directory for TSVs")
    args = parser.parse_args()

    # Load cistopic object
    print("Loading cistopic object...")
    with open(args.cistopic_obj, "rb") as f:
        cistopic_obj = pickle.load(f)
    print(f"  Loaded: {cistopic_obj}")

    # Compute topic QC metrics
    print("Computing topic QC metrics...")
    topic_qc_metrics = compute_topic_metrics(cistopic_obj, return_metrics=True)
    print(f"  Metrics shape: {topic_qc_metrics.shape}")
    print(topic_qc_metrics.head())

    # Save QC metrics
    qc_path = f"{args.output_dir}/topic_qc_metrics.tsv"
    topic_qc_metrics.to_csv(qc_path, sep="\t")
    print(f"  Saved: {qc_path}")

    # Compute topic annotations
    print("Computing topic annotations...")
    topic_annot = topic_annotation(
        cistopic_obj,
        annot_var=args.annot_var,
        general_topic_thr=args.general_topic_thr,
    )
    print(f"  Annotations shape: {topic_annot.shape}")
    print(topic_annot.head())

    # Save annotations
    annot_path = f"{args.output_dir}/topic_annotations_auto.tsv"
    topic_annot.to_csv(annot_path, sep="\t")
    print(f"  Saved: {annot_path}")

    print("Done!")


if __name__ == "__main__":
    main()
