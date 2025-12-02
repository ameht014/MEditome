#!/bin/bash
#SBATCH --job-name=metaEdit_chunks
#SBATCH --account=your_slurm_acc
#SBATCH --qos=your_qos
#SBATCH -p investor
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-99
#SBATCH -o "./logs/MSA_chunk_0330/msa_chunks_%A_%a.out"
#SBATCH -e "./logs/MSA_chunk_0330/msa_chunks_%A_%a.err"

# Let’s assume you have 100 chunk files: ortholog_chunk_0.txt .. ortholog_chunk_99.txt
BASE_DIR="/data/metaEdit"
CHUNKFILE="${BASE_DIR}/experiments/MSA_merge/orthologs_thresh5_samples5_with_samples/chunks/ortholog_chunk_${SLURM_ARRAY_TASK_ID}.tsv"
THREADS=32

mkdir -p ./logs/MSA_chunk_0330/

module load mafft-7.221-gcc-8.2.0-y6cgezm

OUT_DIR="${BASE_DIR}/experiments/MSA_merge/MSA_orthologs/"
mkdir -p $OUT_DIR
EGGNOG_DIR="${BASE_DIR}/experiments/eggnog/eggnog_out"
mkdir -p "$EGGNOG_DIR"

echo "[Info] Processing chunk $CHUNKFILE"
python ../src/msa_for_chunk.py \
   --chunk-file "$CHUNKFILE" \
   --threshold 5 \
   --nproc $THREADS \
   --metaedit-tsv-folder "${BASE_DIR}/experiments/MSA_merge/merged_with_offset/" \
   --eggnog-fasta-folder "$EGGNOG_DIR" \
   --outdir $OUT_DIR

echo "[Done] chunk $CHUNKFILE"


