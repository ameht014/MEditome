#!/bin/bash

# Set variables
INPUT="./orthologs_thresh5_samples5_with_samples/all_orthologs_with_metadata.tsv"
OUTDIR="./orthologs_thresh5_samples5_with_samples/chunks"
N_CHUNKS=100

# Make output directory
mkdir -p "$OUTDIR"

# Extract ortholog and sample_ids columns (as TSV), then explode samples into separate rows
tail -n +2 "$INPUT" | cut -f1,25 | awk 'BEGIN{FS=OFS="\t"} {split($2, samples, ","); for (i in samples) print $1, samples[i]}' > "$OUTDIR/all_ortholog_sample_pairs.tsv"

# Count total lines (ortholog-sample pairs)
TOTAL=$(wc -l < "$OUTDIR/all_ortholog_sample_pairs.tsv")
PER_CHUNK=$(( (TOTAL + N_CHUNKS - 1) / N_CHUNKS ))

# Split the exploded TSV into temp chunks
split -l "$PER_CHUNK" --numeric-suffixes=0 --suffix-length=2 "$OUTDIR/all_ortholog_sample_pairs.tsv" "$OUTDIR/tmp_chunk_"

# Rename chunks to ortholog_chunk_<i>.tsv with no leading zeros and add header
i=0
for file in "$OUTDIR"/tmp_chunk_*; do
    out_file="$OUTDIR/ortholog_chunk_${i}.tsv"
    echo -e "seed_ortholog\tsample" > "$out_file"
    cat "$file" >> "$out_file"
    rm "$file"
    i=$((i + 1))
done

# Cleanup intermediate file
rm "$OUTDIR/all_ortholog_sample_pairs.tsv"

echo "[✓] Wrote $i chunk files to $OUTDIR/"

