#!/bin/bash
#SBATCH --job-name=extract_orthologs
#SBATCH --account=your_slurm_acc
#SBATCH --qos=your_qos
#SBATCH -p investor
#SBATCH -N 1
#SBATCH -o "./logs/extract_orthologs.out"
#SBATCH -e "./logs/extract_orthologs.err"

python3 ../src/extract_all_orthologs.py \
    --coverage-threshold 5 \
    --samples ./data/samples.txt \
      --min-sample-count 5 \
    --out ./orthologs_thresh5_samples5_with_samples/all_orthologs_with_metadata.tsv
#
#[Info] Found 9795 orthologs present in ≥5 samples.
#[✓] Saved metadata for 9795 orthologs to: ./orthologs_thresh5_samples5_with_samples/all_orthologs_with_metadata.tsv
#
