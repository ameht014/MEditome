#!/bin/bash

THREADS=$1
IFS=$'\n' read -d '' -r -a ALL_SAMPLES < ./data/samples.txt
SAMPLE=${ALL_SAMPLES[${SLURM_ARRAY_TASK_ID}]}

## Submit EGGNOG
OUTPUT_DIR=$(pwd)/eggnog_out/${SAMPLE}/
mkdir -p $OUTPUT_DIR

REF_DIR="../../data/References_contigs_raw/"
TMPDIR=$(pwd)/tmp/; export TMPDIR
mkdir -p $TMPDIR

OUT_FILE="${OUTPUT_DIR}/.emapper.genepred.gff"

if [ ! -f "$OUT_FILE" ]; then
    echo "File $OUT_FILE does not exist. Running emapper.py..."
    emapper.py -m mmseqs --itype metagenome -i $REF_DIR/${SAMPLE}_contigs.fna -o $OUTPUT_DIR --cpu $THREADS
else
    echo "File $OUT_FILE already exists. Skipping emapper.py for $SAMPLE."
fi

# Function to define variables based on TYPE
define_variables() {
    local type=$1  # Accept TYPE as a parameter
    BASE_DIR="/data/metaEdit"
    SAMPLES_DIR="${BASE_DIR}/samples/iHMP/IBD/${type}/${SAMPLE}/"
    TRIMMED_DIR=$(pwd)"/trimmed_${type}/${SAMPLE}"
    OUT_DIR=$(pwd)"/RNA_edit_result/${type}/${SAMPLE}"

    # Create parent directories explicitly
    mkdir -p "$(pwd)/trimmed_${type}"   # Ensure parent exists
    mkdir -p "$TRIMMED_DIR"            # Create target subdirectory

    mkdir -p "$(pwd)/RNA_edit_result/"  # Ensure parent exists
    mkdir -p "$(pwd)/RNA_edit_result/${type}"  # Ensure parent exists
    mkdir -p "$OUT_DIR"                       # Create target subdirectory
}

# Global constants
RNA_TYPE="RNA"
DNA_TYPE="DNA_updated"
ALIGNER="bowtie2"
AFFIX1="_R1.fastq.gz"
AFFIX2="_R2.fastq.gz"
AIMAP_PATH=$(pwd)"/../../src/third_party/aimap/"
chmod +x ${AIMAP_PATH}/bin/Adenosine_to_inosine.py
RNA_COV_THRESHOLD="5"
DNA_COV_THRESHOLD="5"
EDIT_THRESHOLD="0.03"
GFF=$(pwd)"/eggnog_out/${SAMPLE}/.emapper.genepred.gff"


# Define input/output variables for RNA
TYPE=$RNA_TYPE
define_variables $TYPE
RNA_DIR=$OUT_DIR

# Submit RNA Edit
singularity exec --bind $SAMPLES_DIR ../../env/metaEdit.sif python3 ${AIMAP_PATH}/bin/Adenosine_to_inosine.py  \
    -g ${REF_DIR}/${SAMPLE}_contigs.fna \
    -tr ${TRIMMED_DIR} \
    -l 18 \
    -a $GFF  \
    -t $THREADS \
    -o $OUT_DIR  \
    --outfile_name ${SAMPLE} \
    -m paired \
    -f1 ${SAMPLES_DIR}/${SAMPLE}${AFFIX1} \
    -f2 ${SAMPLES_DIR}/${SAMPLE}${AFFIX2} \
    --coverage $RNA_COV_THRESHOLD \
    --aligner $ALIGNER \
    -paired_sep \
    -contig_mode

rm -f ${OUT_DIR}/*.sam
rm -f ${OUT_DIR}/*.bam
rm -r $TRIMMED_DIR


# Define input/output variables for DNA
TYPE=$DNA_TYPE
define_variables $TYPE
DNA_DIR=$OUT_DIR

# Submit DNA Edit
singularity exec --bind $SAMPLES_DIR ../../env/metaEdit.sif python3 ${AIMAP_PATH}/bin/Adenosine_to_inosine.py  \
    -g ${REF_DIR}/${SAMPLE}_contigs.fna \
    -tr ${TRIMMED_DIR} \
    -l 18 \
    -a $GFF  \
    -t $THREADS \
    -o $OUT_DIR  \
    --outfile_name ${SAMPLE} \
    -m paired \
    -f1 ${SAMPLES_DIR}/${SAMPLE}${AFFIX1} \
    -f2 ${SAMPLES_DIR}/${SAMPLE}${AFFIX2} \
    --coverage $DNA_COV_THRESHOLD \
    --aligner $ALIGNER \
    -paired_sep \
    -contig_mode

rm -f ${OUT_DIR}/*.sam
rm -f ${OUT_DIR}/*.bam
rm -r $TRIMMED_DIR

OUT_DNA_DIR=$(pwd)"/RNA_edit_result/DNA_updated_filtered/${SAMPLE}"
OUT_RNA_DIR=$(pwd)"/RNA_edit_result/RNA_filtered/${SAMPLE}"

SCRIPT_DIR="$(dirname "$0")"
python3 "$SCRIPT_DIR/filter_RNA.py" $DNA_DIR $RNA_DIR $REF_DIR $SAMPLE $SAMPLE $RNA_COV_THRESHOLD $DNA_COV_THRESHOLD $EDIT_THRESHOLD $OUT_RNA_DIR $OUT_DNA_DIR $THREADS --sample $SAMPLE

rm -r $DNA_DIR

OUT_RNA_DIR=$(pwd)"/RNA_edit_result/RNA_filtered/"
EGGNONG_DIR=$(pwd)"/eggnog_out/"
python3 merge_with_annotations.py $SAMPLE $OUT_RNA_DIR $EGGNONG_DIR

