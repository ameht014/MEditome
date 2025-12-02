#!/bin/bash

BASE_DIR="/data/metaEdit"
SAMPLES_DIR="${BASE_DIR}/samples/DNA/"
mkdir -p $SAMPLES_DIR

mkdir -p $SAMPLES_DIR
#
while IFS= read -r sample_id
do
 if [[ $string == *"_P"* ]]; then
   wget https://downloads.hmpdacc.org/ihmp/ibd/genome/microbiome/wgs/raw/${sample_id}.fastq.gz -P $SAMPLES_DIR
 else
   wget https://downloads.hmpdacc.org/ihmp/ibd/genome/microbiome/wgs/raw/${sample_id}.tar -P $SAMPLES_DIR
 fi
done < data/samples.txt

while IFS= read -r sample_id
do
    mkdir -p ${SAMPLES_DIR}/${sample_id}
    cp ${SAMPLES_DIR}/${sample_id}.tar ${SAMPLES_DIR}/${sample_id}/${sample_id}.tar
    cd ${SAMPLES_DIR}/${sample_id} && tar -xvf ${SAMPLES_DIR}/${sample_id}/${sample_id}.tar
    rm ${SAMPLES_DIR}/${sample_id}/${sample_id}.tar
    
    if [[ $sample_id == *"_P"* ]]; then
        tmp_id="${sample_id%_P}"
        mv ${SAMPLES_DIR}/${sample_id}/${tmp_id}/${tmp_id}_1.fastq.gz ${SAMPLES_DIR}/${sample_id}/${sample_id}_1.fastq.gz
        mv ${SAMPLES_DIR}/${sample_id}/${tmp_id}/${tmp_id}_2.fastq.gz ${SAMPLES_DIR}/${sample_id}/${sample_id}_2.fastq.gz
        rm -r ${SAMPLES_DIR}/${sample_id}/${tmp_id}
    elif [[ $sample_id != *"_TR"* ]]; then
        mv *_1_sequence.txt.bz2 ${SAMPLES_DIR}/${sample_id}/${sample_id}_1_sequence.txt.bz2
        mv *_2_sequence.txt.bz2 ${SAMPLES_DIR}/${sample_id}/${sample_id}_2_sequence.txt.bz2
    fi
done < data/samples.txt


