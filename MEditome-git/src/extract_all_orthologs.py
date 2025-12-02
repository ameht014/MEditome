#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
from collections import defaultdict

# Default config
METAEDIT_TSV_FOLDER = "./merged_with_offset"
DEFAULT_SAMPLES_LIST = "./data/samples.txt"
OUTPUT_FILE = "./all_orthologs_with_metadata.tsv"

# Columns to extract
RELEVANT_COLUMNS = [
    "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
    "RNA_reads_A", "RNA_reads_C", "RNA_reads_G", "RNA_reads_T"
]

def parse_args():
    parser = argparse.ArgumentParser(description="""
    Extract all orthologs with metadata, only if they pass RNA coverage filters
    AND are observed in at least N samples. Adds a column with sample IDs.
    """)
    parser.add_argument("--samples", default=DEFAULT_SAMPLES_LIST,
                        help="Path to a file listing sample names (one per line).")
    parser.add_argument("--out", default=OUTPUT_FILE,
                        help="Output TSV file for ortholog metadata.")
    parser.add_argument("--coverage-threshold", type=int, default=5,
                        help="Minimum coverage threshold for (A & G) or (C & T). Default=5.")
    parser.add_argument("--min-sample-count", type=int, default=5,
                        help="Minimum number of samples in which ortholog must appear. Default=5.")
    return parser.parse_args()

def main():
    args = parse_args()
    coverage_threshold = args.coverage_threshold
    min_sample_count    = args.min_sample_count
    sample_list_file    = args.samples
    output_file         = args.out

    with open(sample_list_file, "r") as f:
        SAMPLE_LIST = [x.strip() for x in f if x.strip()]

    # Track which samples contain each ortholog
    ortholog_to_samples = defaultdict(set)
    ortholog_metadata   = {}

    for sample in SAMPLE_LIST:
        tsv_path = os.path.join(METAEDIT_TSV_FOLDER, f"{sample}_merged_output_with_offset.tsv")
        if not os.path.isfile(tsv_path):
            print(f"[Warning] File not found: {tsv_path}")
            continue

        try:
            df = pd.read_csv(tsv_path, sep="\t", dtype=str)
            if "seed_ortholog" not in df.columns:
                print(f"[Warning] Missing seed_ortholog col in {tsv_path}")
                continue
            df = df.dropna(subset=["seed_ortholog"])

            if "Gene_biotype" in df.columns:
                df = df[df["Gene_biotype"] != "Intergenic region"]

            coverage_cols = ["RNA_reads_A", "RNA_reads_C", "RNA_reads_G", "RNA_reads_T"]
            missing_cov = [col for col in coverage_cols if col not in df.columns]
            if missing_cov:
                print(f"[Warning] Missing coverage columns in {sample}: {missing_cov}")
                continue

            for c in coverage_cols:
                df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

            def coverage_ok(row):
                rA, rC, rG, rT = row["RNA_reads_A"], row["RNA_reads_C"], row["RNA_reads_G"], row["RNA_reads_T"]
                return ((rA >= coverage_threshold and rG >= coverage_threshold) or
                        (rC >= coverage_threshold and rT >= coverage_threshold))

            df = df[df.apply(coverage_ok, axis=1)]

            # Track which samples each ortholog appears in
            for orth in df["seed_ortholog"].unique():
                ortholog_to_samples[orth].add(sample)

            # Fill missing metadata columns
            missing_cols = set(RELEVANT_COLUMNS) - set(df.columns)
            for mc in missing_cols:
                df[mc] = ""

            subset_cols = [c for c in RELEVANT_COLUMNS if c in df.columns]
            subset = df[subset_cols].drop_duplicates(subset=["seed_ortholog"])

            for _, row in subset.iterrows():
                orth = row["seed_ortholog"]
                if orth not in ortholog_metadata:
                    ortholog_metadata[orth] = row.to_dict()

        except Exception as e:
            print(f"[Error] Failed to process {sample}: {e}")

    # === Filter orthologs by sample count ===
    orthologs_passing = [
        orth for orth, samples in ortholog_to_samples.items()
        if len(samples) >= min_sample_count
    ]

    print(f"[Info] Found {len(orthologs_passing)} orthologs present in ≥{min_sample_count} samples.")

    # === Create final dataframe with sample_ids ===
    filtered_metadata = {}

    for orth in orthologs_passing:
        if orth in ortholog_metadata:
            row = ortholog_metadata[orth].copy()
            sample_list = sorted(ortholog_to_samples[orth])
            row["sample_ids"] = ",".join(sample_list)
            filtered_metadata[orth] = row

    output_columns = RELEVANT_COLUMNS + ["sample_ids"]

    if not filtered_metadata:
        print("[Info] No orthologs passed all filters.")
        metadata_df = pd.DataFrame(columns=output_columns)
    else:
        metadata_df = pd.DataFrame.from_dict(filtered_metadata, orient="index")

    if "seed_ortholog" in metadata_df.columns:
        metadata_df = metadata_df.sort_values(by="seed_ortholog")

    # Reorder to keep only available columns in the desired order
    metadata_df = metadata_df[[col for col in output_columns if col in metadata_df.columns]]

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Save to file
    metadata_df.to_csv(output_file, sep="\t", index=False)
    print(f"[✓] Saved metadata for {len(metadata_df)} orthologs to: {output_file}")

if __name__ == "__main__":
    main()
