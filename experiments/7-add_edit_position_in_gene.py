#!/usr/bin/env python3

import os
import pandas as pd


#######################
# 1) Utility Functions
#######################

def parse_gff(gff_file):
    """
    Parse a GFF (e.g. .emapper.genepred.gff) to build a dictionary:
        gene_dict[gene_id] = {
            'contig': <str>,
            'start':  <int>,
            'end':    <int>,
            'strand': <'+' or '-'>
        }

    Only parses lines that appear to be CDS (or gene) with an 'ID=' attribute.
    Adjust if needed for your actual GFF format.
    """
    gene_dict = {}
    if not os.path.isfile(gff_file):
        print(f"  [Warning] GFF file not found: {gff_file}")
        return gene_dict

    with open(gff_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            contig = parts[0]
            feature_type = parts[2].lower()
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            # Only parse 'CDS' lines; adjust if your file uses 'gene' or something else
            if feature_type not in ['cds', 'gene']:
                continue

            attributes = parts[8]  # e.g. "ID=k105_36249_42;some=..."
            gene_id = None

            # Extract "ID=..." from attributes
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('ID='):
                    gene_id = attr.split('=', 1)[1]
                    break

            if gene_id:
                gene_dict[gene_id] = {
                    'contig': contig,
                    'start': start,
                    'end': end,
                    'strand': strand
                }

    return gene_dict


def split_position(pos_str):
    """
    Given a string like 'k105_36249_13353',
    split into (contig_name='k105_36249', contig_position=13353).

    If format differs, adapt this function.
    """
    parts = pos_str.rsplit('_', 1)
    if len(parts) == 2:
        contig_name = parts[0]
        try:
            contig_pos = int(parts[1])
            return contig_name, contig_pos
        except ValueError:
            pass
    # If we cannot parse, return None
    return None, None


def compute_position_in_gene(contig_pos, gene_info):
    """
    Given a contig position (int) and a gene_info dict:
        {
          'contig': <str>,
          'start':  <int>,
          'end':    <int>,
          'strand': <'+ or '-'>
        }
    Return the 1-based offset of contig_pos within that gene.

    If strand == '+':
        offset = contig_pos - start + 1
    If strand == '-':
        offset = end - contig_pos + 1
    """
    if gene_info['strand'] == '-':
        offset = gene_info['end'] - contig_pos + 1
    else:
        offset = contig_pos - gene_info['start'] + 1
    return offset if offset > 0 else None  # If it's out of range, might be None


###################
# 2) Main Script
###################

# 2.1) Read the sample list
SAMPLE_LIST_FILE = '../eggnog/data/samples.txt'
with open(SAMPLE_LIST_FILE, 'r') as f:
    SAMPLE_LIST = [line.strip() for line in f if line.strip()]

# 2.2) Define input/output paths
BASE_DIR = "/data/metaEdit/experiments/eggnog"
BASE_GFF_DIR = f"{BASE_DIR}/eggnog_out"
BASE_INPUT_TSV_DIR = f"{BASE_DIR}/RNA_edit_result/RNA_filtered"
BASE_OUTPUT_DIR = "./merged_with_offset/"

# 2.3) Process each sample
for sample in SAMPLE_LIST:
    print(f"=== Processing sample {sample} ===")

    # (A) Paths
    gff_file = os.path.join(BASE_GFF_DIR, sample, ".emapper.genepred.gff")
    input_tsv = os.path.join(BASE_INPUT_TSV_DIR, sample, f"{sample}_merged_output.tsv")
    out_dir = BASE_OUTPUT_DIR
    out_csv = os.path.join(out_dir, f"{sample}_merged_output_with_offset.tsv")

    # Ensure output directory
    os.makedirs(out_dir, exist_ok=True)

    # (B) Parse GFF
    gene_dict = parse_gff(gff_file)
    if not gene_dict:
        print(f"  [Warning] gene_dict is empty for sample {sample}. Check GFF or feature type.")

    # (C) Read the merged output TSV
    if not os.path.isfile(input_tsv):
        print(f"  [Error] Input TSV not found: {input_tsv}. Skipping.")
        continue

    df = pd.read_csv(input_tsv, sep='\t')

    # (D) Create new columns
    new_cols = {
        'contig_name': [],
        'contig_position': [],
        'gene_start': [],
        'gene_end': [],
        'strand': [],
        'edit_position_in_gene': []
    }

    # (E) For each row, parse 'Position' and match 'Gene_id' to gene_dict
    for idx, row in df.iterrows():
        pos_str = str(row.get('Position', ''))
        gene_id = str(row.get('Gene_id', ''))

        contig_name, contig_pos = split_position(pos_str)
        g_start, g_end, strand, offset_in_gene = None, None, None, None

        if contig_name and contig_pos and (gene_id in gene_dict):
            ginfo = gene_dict[gene_id]
            # Check if contig matches, if you want to be strict:
            if ginfo['contig'] == contig_name:
                g_start = ginfo['start']
                g_end = ginfo['end']
                strand = ginfo['strand']
                offset_in_gene = compute_position_in_gene(contig_pos, ginfo)

        new_cols['contig_name'].append(contig_name)
        new_cols['contig_position'].append(contig_pos)
        new_cols['gene_start'].append(g_start)
        new_cols['gene_end'].append(g_end)
        new_cols['strand'].append(strand)
        new_cols['edit_position_in_gene'].append(offset_in_gene)

    # (F) Attach columns to the DataFrame
    for col_name in new_cols:
        df[col_name] = new_cols[col_name]

    # (G) Save the updated DataFrame
    df.to_csv(out_csv, sep='\t', index=False)
    print(f"  [OK] Wrote {out_csv}\n")

print("All samples processed.")


