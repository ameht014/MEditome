#!/bin/bash

# nohup ./2-submit_eggnog.sh &

source activate eggnog

THREADS=16
MAX_JOBS=51
MAX_SUBMIT_LIMIT=50
BATCH_SIZE=20
INITIAL_BATCH=50
TOTAL_JOBS=748
SUBMITTED=200

mkdir -p ./logs

# Generate a timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")  # Formats the date as YYYYMMDD_HHMMSS

# Create the directory with the timestamp
LOGS_DIR="./logs/${TIMESTAMP}"
mkdir -p $LOGS_DIR

######################### SLURM PARAMETERS #########################
SLURM_ACC="your_slurm_acc"
SLURM_QOS="your_qos"
SLURM_NODE_TYPE="investor"

####################################################################

# Function to check the number of jobs already submitted
check_jobs() {
  squeue -u $USER | wc -l
}

# Flag to handle the initial submission
FIRST_BATCH=true

# Submit jobs in batches
while [ $SUBMITTED -lt $TOTAL_JOBS ]; do
  CURRENT_JOBS=$(check_jobs)

  # Subtract one to account for the header line in `squeue`
  CURRENT_JOBS=$((CURRENT_JOBS - 1))

  # Check if adding the next batch would exceed the submit limit
  if [ $((CURRENT_JOBS + BATCH_SIZE)) -gt $MAX_SUBMIT_LIMIT ]; then
    echo "AssocMaxSubmitJobLimit: Too many jobs in the queue. Waiting for jobs to finish..."
    sleep 60  # Wait and retry
    continue
  fi

  # Determine how many jobs can be submitted
  AVAILABLE_SLOTS=$((MAX_JOBS - CURRENT_JOBS))
  if [ $AVAILABLE_SLOTS -gt 0 ]; then
    if $FIRST_BATCH; then
      # For the first batch, submit up to INITIAL_BATCH
      JOBS_TO_SUBMIT=$((INITIAL_BATCH < AVAILABLE_SLOTS ? INITIAL_BATCH : AVAILABLE_SLOTS))
      FIRST_BATCH=false
    else
      # Subsequent batches are based on BATCH_SIZE
      JOBS_TO_SUBMIT=$((BATCH_SIZE < AVAILABLE_SLOTS ? BATCH_SIZE : AVAILABLE_SLOTS))
    fi

    END_JOB=$((SUBMITTED + JOBS_TO_SUBMIT - 1))

    # Ensure we don't exceed the total number of jobs
    if [ $END_JOB -ge $TOTAL_JOBS ]; then
      END_JOB=$((TOTAL_JOBS - 1))
      JOBS_TO_SUBMIT=$((END_JOB - SUBMITTED + 1))
    fi

    echo "Submitting jobs from $SUBMITTED to $END_JOB"
    sbatch -J metaEdit -a ${SUBMITTED}-${END_JOB} \
          --account=$SLURM_ACC \
          --qos=$SLURM_QOS \
          -p $SLURM_NODE_TYPE \
          -N 1 \
          --cpus-per-task=$THREADS \
          -o "${LOGS_DIR}/stdout-%a.txt" \
         -e "${LOGS_DIR}/stderr-%a.txt" \
          ../src/eggnog_submit.sh $THREADS

    if [ $? -eq 0 ]; then
      SUBMITTED=$((SUBMITTED + JOBS_TO_SUBMIT))
    else
      echo "Submission failed for jobs from $SUBMITTED to $END_JOB"
      sleep 60  # Wait and retry
    fi
  else
    echo "No slots available. Waiting for jobs to finish..."
  fi

  # Wait for a short period before checking again
  sleep 60
done

echo "All jobs submitted!"

