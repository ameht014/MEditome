#!/usr/bin/env python3

"""
msa_for_chunk.py

Given a chunk file with two columns: seed_ortholog, sample
(e.g. ortholog_chunk_X.tsv), do:
  1) Read + group by ortholog => {orth: set_of_samples}
  2) Parse only the relevant sample files in ./merged_with_offset + eggNOG fasta
     *BUT* store large data in global variables to avoid re-pickling
  3) Coverage filter
  4) Each ortholog -> run MSA in a worker
  5) Worker writes its partial result to MSA_merged_{orth}.tsv immediately
     (no big "results" list in memory).
"""

import os
import sys
import argparse
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from multiprocessing import Pool, cpu_count

# ---------------------------------------
# Global dictionaries to avoid pickling overhead:
SEQ_DATA   = {}  # { sample: { gene_id: SeqRecord } }
EDITS_DATA = {}  # { ortholog: { sample: [list_of_edits] } }
OUT_DIR    = "./"  # We'll set this from initargs.
# ---------------------------------------

EGGNOG_METADATA_COLS = [
    "eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name",
    "GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass",
    "BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs"
]

FINAL_COLUMNS = [
    "seed_ortholog","sample","gene_id",
    "raw_position_in_gene","alignment_position",
    "old_base","new_base","avg_edit",
    "dna_reads_A","dna_reads_C","dna_reads_G","dna_reads_T",
    "rna_reads_A","rna_reads_C","rna_reads_G","rna_reads_T"
] + EGGNOG_METADATA_COLS


def parse_args():
    parser = argparse.ArgumentParser(description="Optimized MSA for chunked ortholog data (Approach A+B).")
    parser.add_argument("--chunk-file", required=True,
                        help="TSV with columns: seed_ortholog, sample (header optional).")
    parser.add_argument("--threshold", type=int, default=5,
                        help="Coverage threshold for RNA reads.")
    parser.add_argument("--nproc", type=int, default=0,
                        help="Number of processes for parallel MSA (default=half cores).")
    parser.add_argument("--metaedit-tsv-folder", default="./merged_with_offset/",
                        help="Folder containing <sample>_merged_output_with_offset.tsv")
    parser.add_argument("--eggnog-fasta-folder", default="./eggnog_out",
                        help="Folder with sample subfolders for .emapper.genepred.fasta")
    parser.add_argument("--outdir", default="./MSA_merged_results_jobs",
                        help="Output directory.")
    return parser.parse_args()


# ------------------------------
# Approach B: Initialize global data in workers
# ------------------------------
def init_global_data(seq_data_arg, edits_data_arg, outdir_arg):
    """
    This function is called once per worker to set the global variables
    without pickling them every time we call a worker.
    """
    global SEQ_DATA, EDITS_DATA, OUT_DIR
    SEQ_DATA   = seq_data_arg
    EDITS_DATA = edits_data_arg
    OUT_DIR    = outdir_arg


# ------------------------------
# coverage filter function
# ------------------------------
def coverage_ok(row, threshold):
    rA = int(row["RNA_reads_A"])
    rC = int(row["RNA_reads_C"])
    rG = int(row["RNA_reads_G"])
    rT = int(row["RNA_reads_T"])
    return ((rA >= threshold and rG >= threshold) or
            (rC >= threshold and rT >= threshold))


# ------------------------------
# parse chunk file => {orth: set_of_samples}
# ------------------------------
def load_chunk_map(chunk_file):
    orth_to_samples = defaultdict(set)
    relevant_samples = set()
    with open(chunk_file, "r") as f:
        header = True
        for line in f:
            line=line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if header and ("seed_ortholog" in parts and "sample" in parts):
                header=False
                continue
            header=False

            if len(parts)<2:
                continue
            orth = parts[0].strip()
            samp = parts[1].strip()
            orth_to_samples[orth].add(samp)
            relevant_samples.add(samp)

    return orth_to_samples, relevant_samples


# ------------------------------
# parse coverage data from each sample
# ------------------------------
def parse_metaedit_tsv(tsv_path, coverage_threshold, orth_set):
    import math
    if not os.path.isfile(tsv_path):
        return []

    df = pd.read_csv(tsv_path, sep="\t", dtype=str, na_filter=False)

    needed_cols = ["Gene_id","edit_position_in_gene","seed_ortholog",
                   "RNA_reads_A","RNA_reads_C","RNA_reads_G","RNA_reads_T"]
    for c in needed_cols:
        if c not in df.columns:
            return []

    numeric_cols = ["edit_position_in_gene","RNA_reads_A","RNA_reads_C","RNA_reads_G","RNA_reads_T","average_edit"]
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    if "Gene_biotype" in df.columns:
        df = df[df["Gene_biotype"] != "Intergenic region"]

    # filter by ortholog
    df = df[df["seed_ortholog"].isin(orth_set)]

    # coverage filter
    def row_cov_ok(row):
        return coverage_ok(row, coverage_threshold)
    df = df[df.apply(row_cov_ok, axis=1)]

    edits_list=[]
    for _, row in df.iterrows():
        e = {
            "gene_id": row["Gene_id"].strip(),
            "seed_ortholog": row["seed_ortholog"].strip(),
            "old_base": row.get("Old_base",""),
            "new_base": row.get("New_base",""),
            "avg_edit": float(row.get("average_edit",0)),
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
        # eggNOG metadata
        for mc in EGGNOG_METADATA_COLS:
            e[mc] = row.get(mc, "")
        edits_list.append(e)
    return edits_list


def parse_emapper_genepred_fasta(fa_path):
    from Bio import SeqIO
    gdict={}
    if os.path.isfile(fa_path):
        for rec in SeqIO.parse(fa_path,"fasta"):
            gdict[rec.id] = rec
    return gdict


# ------------------------------
# run MAFFT
# ------------------------------
def run_mafft_alignment(seq_records, tmp_prefix="tmp_mafft"):
    from Bio.Align.Applications import MafftCommandline
    from Bio import AlignIO
    from Bio import SeqIO
    import os

    input_fasta  = tmp_prefix+"_input.fa"
    output_fasta = tmp_prefix+"_aligned.fa"

    with open(input_fasta, "w") as f:
        SeqIO.write(seq_records, f, "fasta")

    cline = MafftCommandline(input=input_fasta)
    stdout, stderr = cline()

    with open(output_fasta,"w") as f:
        f.write(stdout)

    alignment = AlignIO.read(output_fasta,"fasta")
    os.remove(input_fasta)
    os.remove(output_fasta)
    return alignment


# ------------------------------
# Worker function: MSA for a single ortholog
# => We do not return lines, but write them directly to disk (Approach A)
# ------------------------------
def process_one_ortholog(orth):
    """
    Because we are using approach B (global data),
    we can reference global SEQ_DATA, EDITS_DATA, OUT_DIR.
    We'll run MSA for this ortholog, then write a separate file:
        <OUT_DIR>/MSA_merged_{orth}.tsv
    """
    global SEQ_DATA, EDITS_DATA, OUT_DIR
    from Bio.SeqRecord import SeqRecord
    from collections import defaultdict

    if orth not in EDITS_DATA:
        return  # no coverage data
    orth_data = EDITS_DATA[orth]  # { sample: [ list_of_edits ] }

    # collect unique (sample, gene_id)
    unique_pairs=set()
    for sample, e_list in orth_data.items():
        for e in e_list:
            unique_pairs.add((sample,e["gene_id"]))

    # gather seqRecords
    seq_records=[]
    for (smp, gid) in unique_pairs:
        if smp in SEQ_DATA and gid in SEQ_DATA[smp]:
            rec = SEQ_DATA[smp][gid]
            new_id = f"{smp}::{gid}"
            seq_records.append(SeqRecord(rec.seq, id=new_id, description=""))

    if len(seq_records)<2:
        # skip (no MSA needed, or not enough seqs)
        return

    # MSA
    try:
        alignment = run_mafft_alignment(seq_records, tmp_prefix=f"tmp_mafft_{orth}")
    except:
        return

    # build rawpos->alignmentpos
    rawpos_map={}
    for aln_rec in alignment:
        sp=aln_rec.id.split("::",1)
        if len(sp)!=2:
            continue
        smp, gnm=sp
        raw_i=0
        map_dict={}
        for j, base in enumerate(str(aln_rec.seq)):
            if base!='-':
                map_dict[raw_i] = j
                raw_i+=1
        rawpos_map[(smp,gnm)] = map_dict

    # we'll write output immediately
    out_path = os.path.join(OUT_DIR, f"MSA_merged_{orth}.tsv")
    # write header + lines
    with open(out_path,"w") as out:
        out.write("\t".join(FINAL_COLUMNS)+"\n")
        for sample, e_list in orth_data.items():
            for e in e_list:
                rp = e["edit_position_in_gene"]
                align_map = rawpos_map.get((sample,e["gene_id"]), {})
                alignment_pos = align_map.get(rp, -1)
                row_dict = {
                    "seed_ortholog": orth,
                    "sample": sample,
                    "gene_id": e["gene_id"],
                    "raw_position_in_gene": rp,
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
                for mc in EGGNOG_METADATA_COLS:
                    row_dict[mc] = e.get(mc,"")

                out_line = [str(row_dict.get(c,"")) for c in FINAL_COLUMNS]
                out.write("\t".join(out_line)+"\n")


# ------------------------------
def main():
    args = parse_args()

    # load chunk => orth->samples, relevant_samples
    orth_to_samples, relevant_samples = load_chunk_map(args.chunk_file)
    print(f"[Info] chunk {args.chunk_file}: {len(orth_to_samples)} orthologs, {len(relevant_samples)} samples, threshold={args.threshold}")

    # parse coverage + load seq data
    global SEQ_DATA, EDITS_DATA, OUT_DIR
    SEQ_DATA={}
    EDITS_DATA={}
    OUT_DIR=args.outdir
    os.makedirs(OUT_DIR,exist_ok=True)

    # parse coverage for relevant samples
    for sample in relevant_samples:
        tsv_file = os.path.join(args.metaedit_tsv_folder, f"{sample}_merged_output_with_offset.tsv")
        fa_file  = os.path.join(args.eggnog_fasta_folder, sample, ".emapper.genepred.fasta")
        if not os.path.isfile(tsv_file) or not os.path.isfile(fa_file):
            continue
        sample_edits = parse_metaedit_tsv(tsv_file, args.threshold, orth_to_samples.keys())
        sample_gene_seqs = parse_emapper_genepred_fasta(fa_file)
        SEQ_DATA[sample] = sample_gene_seqs

        for e in sample_edits:
            so = e["seed_ortholog"]
            if so not in EDITS_DATA:
                EDITS_DATA[so]={}
            if sample not in EDITS_DATA[so]:
                EDITS_DATA[so][sample]=[]
            EDITS_DATA[so][sample].append(e)

    # build tasks => just the orthologs
    tasks=[]
    for orth in orth_to_samples:
        if orth in EDITS_DATA:
            tasks.append(orth)

    total_orth=len(tasks)
    if total_orth==0:
        print("[Info] No orthologs to process.")
        return

    # parallel init with global data approach
    if args.nproc <=0:
        nproc = max(1, cpu_count()//2)
    else:
        nproc = args.nproc

    print(f"[Info] Starting MSA with {nproc} processes, {total_orth} ortholog tasks.")
    pool=Pool(nproc, initializer=init_global_data, initargs=(SEQ_DATA, EDITS_DATA, OUT_DIR))

    # track progress
    done=0
    for _ in pool.imap(process_one_ortholog, tasks):
        done+=1
        if done%10==0 or done==total_orth:
            print(f"[Info] Completed {done}/{total_orth} orthologs.")
    pool.close()
    pool.join()

    print(f"[✓] Done. Each ortholog result is in {OUT_DIR}/MSA_merged_<orth>.tsv")


if __name__=="__main__":
    main()

