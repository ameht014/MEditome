#!/bin/bash
#SBATCH --job-name=annotate_protein_change
#SBATCH --account=your_slurm_acc
#SBATCH --qos=your_qos
#SBATCH -p investor
#SBATCH -N 1
#SBATCH -o "./logs/annotate_protein_change.out"
#SBATCH -e "./logs/annotate_protein_change.err"

# run the script
python3 ../experiments/13-annotate_protein_change.py

