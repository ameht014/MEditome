#!/bin/bash

IFS=$'\n' read -d '' -r -a ALL_SAMPLES < ./data/samples.txt
SAMPLE=${ALL_SAMPLES[${SLURM_ARRAY_TASK_ID}]}
OUTPUT_DIR=$(pwd)/eggnog_out/${SAMPLE}/
mkdir -p $OUTPUT_DIR

REF_DIR="../../data/References_contigs_raw/"
THREADS=$1
TMPDIR=$(pwd)/tmp/; export TMPDIR
mkdir -p $TMPDIR

OUT_FILE="${OUTPUT_DIR}/.emapper.genepred.gff"

if [ ! -f "$OUT_FILE" ]; then
    echo "File $OUT_FILE does not exist. Running emapper.py..."
    emapper.py -m mmseqs --itype metagenome -i $REF_DIR/${SAMPLE}_contigs.fna -o $OUTPUT_DIR --cpu $THREADS
else
    echo "File $OUT_FILE already exists. Skipping emapper.py for $SAMPLE."
fi
