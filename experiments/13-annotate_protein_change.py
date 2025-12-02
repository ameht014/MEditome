#!/usr/bin/env python3
from pathlib import Path
import sys, os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq

# ---------- 1. paths ----------
FILTER_CSV  = Path("./data/extracted_orthologs_filter5.csv")
BASE_DIR = Path("/data/metaEdit")
MSA_DIR     = BASE_DIR / "experiments/MSA_merge/MSA_orthologs"
EGGNOG_ROOT = BASE_DIR / "experiments/eggnog/eggnog_out"

OUT_DIR = Path("./output")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------- 2. constants ----------
GENETIC_CODE = 11
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
VERBOSE = True
def log(msg):
    if VERBOSE:
        print(msg, file=sys.stderr)

# ---------- 3. read filter ----------
flt = pd.read_csv(FILTER_CSV, header=0)
flt.rename(columns={flt.columns[0]: "seed_ortholog"}, inplace=True)
wanted = {o: set(map(int, g["alignment_position"]))
          for o, g in flt.groupby("seed_ortholog")}
log(f"[INFO] loaded filter for {len(wanted)} orthologs")

# ---------- 4. sample caches ----------
_cds, _strand, _loaded = {}, {}, set()
def load_sample(s):
    if s in _loaded: return
    fasta = EGGNOG_ROOT / s / ".emapper.genepred.fasta"
    gff   = EGGNOG_ROOT / s / ".emapper.genepred.gff"
    for rec in SeqIO.parse(str(fasta), "fasta"):
        _cds[(s, rec.id)] = rec.seq
    if gff.is_file():
        with open(gff) as fh:
            for l in fh:
                if l.startswith("#"): continue
                p=l.rstrip().split("\t")
                if len(p)<9 or p[2].lower()!="cds": continue
                gid = dict(kv.split("=",1) for kv in p[8].split(";") if "=" in kv).get("ID")
                if gid: _strand[(s, gid)] = p[6]
    _loaded.add(s)
def seq_of(s,g):  load_sample(s); return _cds.get((s,g))
def strand_of(s,g): load_sample(s); return _strand.get((s,g), "+")

# ---------- 5. helpers ----------
def comp_base(b):  return b.translate(COMP)
def comp_counts(r,p): return {f"{p}_A": int(r[f"{p}_T"]),
                              f"{p}_C": int(r[f"{p}_G"]),
                              f"{p}_G": int(r[f"{p}_C"]),
                              f"{p}_T": int(r[f"{p}_A"])}
def annotate(seq,pos0,new_b):
    alt = MutableSeq(str(seq)); alt[pos0]=new_b
    cs=(pos0//3)*3
    aa_ref = str(seq[cs:cs+3].translate(table=GENETIC_CODE))
    aa_alt = str(Seq(str(alt[cs:cs+3])).translate(table=GENETIC_CODE))
    eff = "synonymous" if aa_ref==aa_alt else "nonsense" if aa_alt=="*" else "missense"
    return aa_ref, aa_alt, eff, alt

# ---------- 6. main ----------
for orth, pos_set in wanted.items():
    tsv = MSA_DIR / f"MSA_merged_{orth}.tsv"
    if not tsv.is_file():
        log(f"[SKIP] {orth} (no TSV)")
        continue

    df = pd.read_csv(tsv, sep="\t", dtype=str)
    df["alignment_position"]   = pd.to_numeric(df["alignment_position"], errors="coerce")
    df["raw_position_in_gene"] = pd.to_numeric(df["raw_position_in_gene"], errors="coerce")
    df = df[df["alignment_position"].isin(pos_set)]
    if df.empty: continue

    outs=[]
    for _, row in df.iterrows():
        samp = row["sample"]
        gid  = row["gene_id"]
        pos1 = int(row["raw_position_in_gene"])
        pos0 = pos1 - 1

        strand = strand_of(samp,gid)
        strand = strand if strand in "+-" else "+"

        old_b = row["old_base"].upper()
        new_b = row["new_base"].upper()
        if strand=="-":
            old_b, new_b = comp_base(old_b), comp_base(new_b)
            dna = comp_counts(row,"dna_reads")
            rna = comp_counts(row,"rna_reads")
        else:
            dna = {k:int(row[k]) for k in ["dna_reads_A","dna_reads_C","dna_reads_G","dna_reads_T"]}
            rna = {k:int(row[k]) for k in ["rna_reads_A","rna_reads_C","rna_reads_G","rna_reads_T"]}

        seq = seq_of(samp,gid)
        if seq is None or pd.isna(pos1): continue

        aa_ref, aa_alt, eff, seq_alt = annotate(seq,pos0,new_b)

        out = row.to_dict()
        out.update({
            "strand": strand,
            "old_base_strand_corrected": old_b,
            "new_base_strand_corrected": new_b,
            "prot_seq_before": str(seq.translate(table=GENETIC_CODE)),   # NEW
            "prot_seq_after":  str(seq_alt.translate(table=GENETIC_CODE)), # NEW
            "aa_ref": aa_ref,
            "aa_alt": aa_alt,
            "effect": eff,
            **{k+"_strand_corrected": dna[k] for k in dna},
            **{k+"_strand_corrected": rna[k] for k in rna},
        })
        outs.append(out)

    if outs:
        out_file = OUT_DIR / f"RNA_edit_protein_effects_{orth}.tsv"
        pd.DataFrame(outs).to_csv(out_file, sep="\t", index=False)
        log(f"[OK] wrote {out_file} ({len(outs):,} rows)")

log("[✓] Done.")
