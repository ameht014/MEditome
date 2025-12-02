# MetaEdit2 Experiments

This repository contains scripts for detecting and analyzing RNA editing events in microbiome datasets. The `experiments/` folder provides a stepwise workflow used for the iHMP IBD samples. Below is a high-level overview of the numbered experiments.

## Experiment overview

1. **Select samples** – Jupyter notebook to identify iHMP samples that have both metagenomic (DNA) and metatranscriptomic (RNA) data.
2. **Download DNA** – Bash script to download and extract WGS reads for the selected samples.
3. **Download RNA** – Bash script to download and extract RNA reads for the same set of samples.
4. **Submit eggNOG** – SLURM submission wrapper running `eggnog-mapper` on contigs for each sample to obtain gene annotations.
5. **Submit MetaEdit** – SLURM script using `metaEdit_one.sh` to run the Adenosine‑to‑Inosine detection pipeline for RNA and DNA, filter results and merge with eggNOG annotations.
6. **Merge with annotations** – Python script that combines filtered edit calls with the corresponding eggNOG annotation tables.
7. **Add edit position in gene** – Computes the edit position relative to gene boundaries using the `.emapper.genepred.gff` output.
8. **Extract orthologs** – Calls `extract_all_orthologs.py` to gather orthologs passing coverage filters and appearing in multiple samples.
9. **Split orthologs by chunks** – Splits the list of ortholog/sample pairs into smaller chunk files for parallel processing.
10. **Submit MSA for chunk** – SLURM job running `msa_for_chunk.py` on each chunk to align gene sequences and map edits to alignment coordinates.
11. **Aggregate orthologs** – Aggregates all per‑ortholog MSA results into a single summary table.
12. **Extract orthologs for protein change** – Creates a list of ortholog/alignment positions to examine potential protein‑level effects.
13. **Annotate protein change** – Reads MSA results and gene sequences to predict amino‑acid changes caused by RNA edits.

The `experiments/README.md` file contains basic setup notes (conda environment and Singularity usage).
