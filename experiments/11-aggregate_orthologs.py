#!/usr/bin/env python3

import os
import glob
import pandas as pd
from tqdm import tqdm

from collections import defaultdict

# Set folder containing MSA_merged_*.tsv files
INPUT_FOLDER = "./MSA_orthologs/"
OUTPUT_FILE  = "./MSA_merged/MSA_summary_aggregated.tsv"
coverage_columns = [
    "dna_reads_A", "dna_reads_C", "dna_reads_G", "dna_reads_T",
    "rna_reads_A", "rna_reads_C", "rna_reads_G", "rna_reads_T"
]

metadata_columns = [
    "eggNOG_OGs","max_annot_lvl","COG_category","Description",
    "Preferred_name","GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module",
    "KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"
]

preserve_columns = [
    "seed_ortholog", "gene_id", "alignment_position", "old_base", "new_base"
] + coverage_columns + ["avg_edit"] + metadata_columns

group_key = ["seed_ortholog", "alignment_position", "old_base", "new_base"]

# --- Read all MSA_merged files ---
all_files = glob.glob(os.path.join(INPUT_FOLDER, "MSA_merged_*.tsv"))
print(f"[Info] Found {len(all_files)} MSA_merged_*.tsv files")

all_rows = []
for file in tqdm(all_files, desc="Reading MSA_merged files"):
    try:
        df = pd.read_csv(file, sep="\t", dtype=str)
        df = df[df["alignment_position"] != "-1"]

        for col in coverage_columns + ["avg_edit", "alignment_position"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

        all_rows.append(df)
    except Exception as e:
        print(f"[Warning] Failed to read {file}: {e}")

if not all_rows:
    print("[Error] No data to process.")
    exit(1)

# --- Combine all data ---
combined = pd.concat(all_rows, ignore_index=True)

# --- Aggregation Step 1: Numeric + metadata ---
agg_dict = {col: "sum" for col in coverage_columns}
agg_dict["avg_edit"] = "mean"
agg_dict["gene_id"] = "first"
for col in metadata_columns:
    agg_dict[col] = "first"

summary = (
    combined
    .groupby(group_key, as_index=False)
    .agg(agg_dict)
)

# --- Aggregation Step 2: num_samples + sample_ids ---
sample_stats = (
    combined
    .groupby(group_key)
    .agg(
        num_samples=("sample", "nunique"),
        sample_ids=("sample", lambda x: ",".join(sorted(set(x))))
    )
    .reset_index()
)

# --- Merge ---
summary = summary.merge(sample_stats, on=group_key, how="left")

# --- Save ---
summary.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"[] Aggregated summary written to: {OUTPUT_FILE}")