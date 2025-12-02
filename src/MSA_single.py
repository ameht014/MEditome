#!/usr/bin/env python3

"""
merge_rna_edits_across_samples_parallel.py

1) Parse the *updated* MetaEdit + eggNOG TSV files from multiple samples,
   BUT filter out edits below coverage threshold in the parse stage.
2) Parse the .emapper.genepred.fasta for each sample
3) Parallelize the multiple-sequence alignments across orthologs (in-process)
4) Merge/compare the RNA-edit coordinates across samples in MSA space
   and write the final lines to 'merged_rna_edits.tsv'.
"""

import os
import argparse
import pandas as pd
import math
from multiprocessing import Pool, cpu_count
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO

# -------------------------------
# 1) PARAMETERS & PATHS
# -------------------------------
METAEDIT_TSV_FOLDER = "./merged_with_offset/"
BASE_DIR = "/data/metaEdit/experiments"
EGGNOG_FASTA_FOLDER = f"{BASE_DIR}/eggnog/eggnog_out"
OUTPUT_FOLDER       = "./MSA_merged_results"

DEFAULT_SAMPLES_LIST_FILE = "./data/samples.txt"

# 2) Extra EggNOG metadata columns (if present)
EGGNOG_METADATA_COLS = [
    "eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name","GOs","EC","KEGG_ko",
    "KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"
]

# 3) The final output columns
FINAL_COLUMNS = [
    "seed_ortholog",
    "sample",
    "gene_id",
    "raw_position_in_gene",
    "alignment_position",
    "old_base",
    "new_base",
    "avg_edit",
    "dna_reads_A",
    "dna_reads_C",
    "dna_reads_G",
    "dna_reads_T",
    "rna_reads_A",
    "rna_reads_C",
    "rna_reads_G",
    "rna_reads_T"
] + EGGNOG_METADATA_COLS

# Global structures
edits_data = {}  # e.g. edits_data[ortholog][sample] = [list_of_edit_dicts]
seq_data   = {}  # e.g. seq_data[sample][gene_id] = SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(description="Parallel MSA merging of RNA-edit data (coverage filtered at parse stage).")
    parser.add_argument("--threshold", type=int, default=10,
                        help="Coverage threshold: keep an edit only if (A & G >= threshold) or (C & T >= threshold).")
    parser.add_argument("--samples", default=DEFAULT_SAMPLES_LIST_FILE,
                        help="File with list of samples (one per line).")
    parser.add_argument("--nproc", type=int, default=0,
                        help="Number of processes (0 => half of available cores).")
    return parser.parse_args()


# -------------------------------
# 2) PARSE TSV FILES & FASTA FILES
# -------------------------------
def parse_metaedit_tsv(tsv_path, coverage_threshold):
    """
    Reads one sample's *merged_output_with_offset.tsv, returns a list of dicts (one per RNA-edit),
    but ONLY if coverage meets the threshold:
      (rnaA >= threshold && rnaG >= threshold) || (rnaC >= threshold && rnaT >= threshold).
    """
    edits_list = []

    df = pd.read_csv(tsv_path, sep="\t", dtype=str)

    # Filter out rows missing essential columns
    required_cols = ["Gene_id","edit_position_in_gene","seed_ortholog"]
    for col in required_cols:
        if col not in df.columns:
            print(f"[Warning] Missing required column '{col}' in {tsv_path}")
            return edits_list

    # Numeric columns
    numeric_cols = ["edit_position_in_gene","average_edit",
                    "DNA_reads_A","DNA_reads_C","DNA_reads_G","DNA_reads_T",
                    "RNA_reads_A","RNA_reads_C","RNA_reads_G","RNA_reads_T"]
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)
        else:
            df[c] = 0  # if missing, default 0

    # Exclude intergenic if Gene_biotype column exists
    if "Gene_biotype" in df.columns:
        df = df[df["Gene_biotype"] != "Intergenic region"]

    # Drop any rows missing required fields
    df = df.dropna(subset=["Gene_id","edit_position_in_gene","seed_ortholog"])

    # Coverage filtering here
    # We'll keep only those rows that satisfy coverage rule:
    # (rA >= threshold && rG >= threshold) || (rC >= threshold && rT >= threshold)
    def coverage_ok(row):
        rA = row["RNA_reads_A"]
        rC = row["RNA_reads_C"]
        rG = row["RNA_reads_G"]
        rT = row["RNA_reads_T"]
        return ((rA >= coverage_threshold and rG >= coverage_threshold) or
                (rC >= coverage_threshold and rT >= coverage_threshold))

    # convert these columns to int so we can compare quickly
    df["RNA_reads_A"] = df["RNA_reads_A"].astype(int)
    df["RNA_reads_C"] = df["RNA_reads_C"].astype(int)
    df["RNA_reads_G"] = df["RNA_reads_G"].astype(int)
    df["RNA_reads_T"] = df["RNA_reads_T"].astype(int)

    df = df[df.apply(coverage_ok, axis=1)]

    for idx, row in df.iterrows():
        gene_id       = str(row["Gene_id"]).strip()
        seed_ortholog = str(row["seed_ortholog"]).strip()
        if not seed_ortholog:
            continue

        old_base = row["Old_base"] if "Old_base" in df.columns else ""
        new_base = row["New_base"] if "New_base" in df.columns else ""

        e = {
            "gene_id": gene_id,
            "seed_ortholog": seed_ortholog,
            "old_base": old_base,
            "new_base": new_base,
            "avg_edit": float(row["average_edit"]) if "average_edit" in df.columns else 0.0,
            "dna_reads_A": int(row["DNA_reads_A"]),
            "dna_reads_C": int(row["DNA_reads_C"]),
            "dna_reads_G": int(row["DNA_reads_G"]),
            "dna_reads_T": int(row["DNA_reads_T"]),
            "rna_reads_A": int(row["RNA_reads_A"]),
            "rna_reads_C": int(row["RNA_reads_C"]),
            "rna_reads_G": int(row["RNA_reads_G"]),
            "rna_reads_T": int(row["RNA_reads_T"]),
            "edit_position_in_gene": int(row["edit_position_in_gene"]),
        }
        # Add any EggNOG columns
        for mc in EGGNOG_METADATA_COLS:
            e[mc] = row[mc] if mc in df.columns and not pd.isna(row[mc]) else ""

        edits_list.append(e)

    return edits_list


def parse_emapper_genepred_fasta(fasta_path):
    """Parse .emapper.genepred.fasta => { gene_id: SeqRecord }."""
    gene_dict = {}
    if not os.path.isfile(fasta_path):
        return gene_dict
    for record in SeqIO.parse(fasta_path, "fasta"):
        gene_dict[record.id] = record
    return gene_dict


def run_mafft_alignment(seq_records, tmp_prefix="tmp_mafft"):
    """Runs MAFFT on a list of SeqRecord objects, returns a MultipleSeqAlignment."""
    input_fasta  = tmp_prefix + "_input.fa"
    output_fasta = tmp_prefix + "_aligned.fa"

    with open(input_fasta, "w") as f:
        SeqIO.write(seq_records, f, "fasta")

    mafft_cline = MafftCommandline(input=input_fasta)
    stdout, stderr = mafft_cline()

    with open(output_fasta, "w") as f:
        f.write(stdout)

    alignment = AlignIO.read(output_fasta, "fasta")

    os.remove(input_fasta)
    os.remove(output_fasta)
    return alignment


def load_all_samples(coverage_threshold, samples_file):
    """Loads all samples into edits_data + seq_data, skipping edits below coverage threshold."""
    with open(samples_file, "r") as f:
        sample_list = [line.strip() for line in f if line.strip()]

    for sample in sample_list:
        tsv_file = os.path.join(METAEDIT_TSV_FOLDER, f"{sample}_merged_output_with_offset.tsv")
        fa_file  = os.path.join(EGGNOG_FASTA_FOLDER, sample, ".emapper.genepred.fasta")

        if not os.path.isfile(tsv_file):
            print(f"[Warning] TSV not found for sample {sample}: {tsv_file}")
            continue
        if not os.path.isfile(fa_file):
            print(f"[Warning] FASTA not found for sample {sample}: {fa_file}")
            continue

        # parse & filter by coverage threshold
        sample_edits = parse_metaedit_tsv(tsv_file, coverage_threshold=coverage_threshold)

        sample_gene_seqs = parse_emapper_genepred_fasta(fa_file)
        seq_data[sample] = sample_gene_seqs

        for entry in sample_edits:
            orth = entry["seed_ortholog"]
            if orth not in edits_data:
                edits_data[orth] = {}
            if sample not in edits_data[orth]:
                edits_data[orth][sample] = []
            edits_data[orth][sample].append(entry)

    # Optionally remove orthologs that ended up with no edits after filtering
    # But if an ortholog has 0 edits, it won't matter. We'll skip it in alignment anyway.
    print(f"[Info] Loaded {len(edits_data)} orthologs (some may have zero edits).")
    return sample_list


def process_one_ortholog(args):
    """
    Worker function to run MSA for a single ortholog, build raw->aln map,
    and return lines (as list of strings) for that ortholog.
    (We no longer do coverage filtering here, as it's done at parse time.)
    """
    (seed_orth, orth_data, seq_data_local) = args

    # Build seq_records
    seq_records = []
    unique_pairs = set()
    for sample in orth_data:
        for e in orth_data[sample]:
            g_id = e["gene_id"]
            unique_pairs.add((sample,g_id))

    for (sample,g_id) in unique_pairs:
        if g_id in seq_data_local[sample]:
            rec = seq_data_local[sample][g_id]
            new_id = f"{sample}::{g_id}"
            seq_records.append(SeqRecord(rec.seq, id=new_id, description=""))

    if len(seq_records) < 2:
        # means there's only 1 or 0 genes left for this ortholog => skip
        return []

    try:
        alignment = run_mafft_alignment(seq_records, tmp_prefix=f"tmp_mafft_{seed_orth}")
    except Exception as e:
        print(f"[Error] MAFFT alignment failed for ortholog: {seed_orth}, error: {e}")
        return []

    # build raw->aln map
    rawpos_map = {}
    for aln_rec in alignment:
        sp = aln_rec.id.split("::",1)
        if len(sp)!=2:
            continue
        smp,gnm = sp
        map_dict = {}
        raw_i=0
        aln_str = str(aln_rec.seq)
        for j, base in enumerate(aln_str):
            if base!='-':
                map_dict[raw_i]= j
                raw_i+=1
        rawpos_map[(smp,gnm)] = map_dict

    # Build lines
    lines = []
    for sample in orth_data:
        for e in orth_data[sample]:
            raw_pos_in_gene = e["edit_position_in_gene"]
            if raw_pos_in_gene is None:
                continue
            align_map = rawpos_map.get((sample,e["gene_id"]), {})
            alignment_pos = align_map.get(raw_pos_in_gene, -1)

            row_dict = {
                "seed_ortholog": seed_orth,
                "sample": sample,
                "gene_id": e["gene_id"],
                "raw_position_in_gene": raw_pos_in_gene,
                "alignment_position": alignment_pos,
                "old_base": e["old_base"],
                "new_base": e["new_base"],
                "avg_edit": f"{e['avg_edit']:.6f}",
                "dna_reads_A": e["dna_reads_A"],
                "dna_reads_C": e["dna_reads_C"],
                "dna_reads_G": e["dna_reads_G"],
                "dna_reads_T": e["dna_reads_T"],
                "rna_reads_A": e["rna_reads_A"],
                "rna_reads_C": e["rna_reads_C"],
                "rna_reads_G": e["rna_reads_G"],
                "rna_reads_T": e["rna_reads_T"],
            }
            # Extra EggNOG fields
            for mc in EGGNOG_METADATA_COLS:
                row_dict[mc] = e.get(mc,"")

            out_line = [ str(row_dict.get(col,"")) for col in FINAL_COLUMNS ]
            lines.append("\t".join(out_line))

    return lines


def main():
    parser = parse_args()
    coverage_threshold = parser.threshold
    samples_file       = parser.samples
    user_nproc         = parser.nproc

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    merged_output_path = os.path.join(OUTPUT_FOLDER, "merged_rna_edits.tsv")

    # 1) Load + filter coverage at parse stage
    sample_list = load_all_samples(coverage_threshold, samples_file)
    print(f"[Info] Done loading samples, we have {len(edits_data)} orthologs total.")
    print(f"[Info] Coverage threshold: {coverage_threshold}")

    # 2) Build tasks
    tasks = []
    # We'll only keep orthologs that have at least 2 total sequences after filtering
    # but let's do that check in the worker anyway, so let's just pass them all
    for seed_orth, orth_data in edits_data.items():
        tasks.append((seed_orth, orth_data, seq_data))

    # 3) Parallelize with Pool
    if user_nproc <= 0:
        nproc = max(1, cpu_count()//2)
    else:
        nproc = user_nproc
    print(f"[Info] Starting Pool with {nproc} processes.")
    with Pool(processes=nproc) as p:
        results_list = p.map(process_one_ortholog, tasks)

    # 4) Flatten + write
    with open(merged_output_path, "w") as out:
        out.write("\t".join(FINAL_COLUMNS) + "\n")
        for lines_for_ortholog in results_list:
            for line in lines_for_ortholog:
                out.write(line + "\n")

    print(f"[] Wrote final merged output to: {merged_output_path}")


if __name__=="__main__":
    main()
