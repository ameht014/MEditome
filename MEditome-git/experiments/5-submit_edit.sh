#!/bin/bash

module load singularity-3.5.3
source activate eggnog

THREADS=16
mkdir -p ./logs

# Generate a timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")  # Formats the date as YYYYMMDD_HHMMSS

# Create the directory with the timestamp
LOGS_DIR="./logs/metaEdit_test_multiple/"
mkdir -p $LOGS_DIR

######################### SLURM PARAMETERS #########################
SLURM_ACC="your_slurm_acc"
SLURM_QOS="your_qos"
SLURM_NODE_TYPE="default-partition"
####################################################################

# Total number of jobs
TOTAL_JOBS=749
BATCH_SIZE=99

# Loop through and submit jobs
for ((START=350; START<TOTAL_JOBS; START+=BATCH_SIZE)); do
    END=$((START + BATCH_SIZE - 1))
    if [[ $END -ge $TOTAL_JOBS ]]; then
        END=$((TOTAL_JOBS - 1))
    fi

    sbatch -J metaEdit -a ${START}-${END} \
           --account=$SLURM_ACC \
           --qos=$SLURM_QOS \
           -p $SLURM_NODE_TYPE \
           -N 1 \
           --cpus-per-task=$THREADS \
           -o "${LOGS_DIR}/stdout-%a.txt" \
           -e "${LOGS_DIR}/stderr-%a.txt" \
           ../src/metaEdit_one.sh $THREADS
done


